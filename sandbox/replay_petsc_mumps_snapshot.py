# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Replay a captured Saint-Venant Jacobian with PETSc LU/MUMPS only.

Use the SAME .petscbin file in the old and new conda environments.
This script imports no SCRIMP, GetFEM or Gmsh. It deliberately
recreates a PETSc KSP/PC for each cycle, then calls KSP.reset() and
KSP.destroy() to exercise the MUMPS-factor teardown.

Capture matrices in the OLD environment with
  SCRIMP_JACOBIAN_SNAPSHOT_PREFIX=sandbox/outputs/mumps_ab/jac \
  mpiexec -n 2 python -u sandbox/viscous_swe_2d_mpi.py ...

Then replay exactly the same file with, for example,
  mpiexec -n 2 python -u sandbox/replay_petsc_mumps_snapshot.py \
    --matrix sandbox/outputs/mumps_ab/jac_004.petscbin --cycles 3

PETSc binary files store the assembled matrix, NOT the Newton
iterate or the lifetime of the original TS/SNES/KSP objects.
A passing replay does not rule out earlier heap corruption in
GetFEM/SCRIMP or a bug specific to the TS factor lifecycle.

For a stronger simulation of in-place Newton Jacobian changes, pass:
  --matrix path/to/jac_001.petscbin --next-matrix path/to/jac_002.petscbin
The same KSP/PC and same mutable PETSc Mat then survive alternating
001 -> 002 -> 001 updates; only after all solves do we reset KSP.
This still does NOT instantiate SNES/TS or call GetFEM.
"""

import argparse
from pathlib import Path

from petsc4py import PETSc


def emit(comm, message):
    print(
        f"[snapshot MPI {comm.getRank()}/{comm.getSize()}] {message}",
        flush=True,
    )


def load_matrix(path, comm):
    viewer = PETSc.Viewer().createBinary(str(path), mode="r", comm=comm)
    matrix = PETSc.Mat().create(comm=comm)
    try:
        matrix.setType("aij")
        matrix.load(viewer)
    except BaseException:
        matrix.destroy()
        raise
    finally:
        viewer.destroy()
    return matrix


def replay(path, cycles, factor_only):
    comm = PETSc.COMM_WORLD
    # PETSc binary I/O is collective; all ranks read the same path.
    if not Path(path).is_file():
        raise FileNotFoundError(path)
    matrix = load_matrix(path, comm)
    try:
        shape = matrix.getSize()
        if shape[0] != shape[1]:
            raise ValueError(f"Expected square Jacobian, got {shape}")
        emit(
            comm,
            f"loaded {path}: size={shape}, "
            f"Frobenius={matrix.norm(PETSc.NormType.FROBENIUS):.12g}, "
            f"factor_only={factor_only}, cycles={cycles}",
        )

        x_true = matrix.createVecRight()
        b = matrix.createVecLeft()
        try:
            x_true.set(1.0)
            matrix.mult(x_true, b)
            emit(comm, f"constructed consistent RHS: ||b||={b.norm():.9g}")

            for cycle in range(cycles):
                x = matrix.createVecRight()
                ksp = PETSc.KSP().create(comm=comm)
                try:
                    x.set(0.0)
                    ksp.setOperators(matrix)
                    ksp.setType("preonly")
                    pc = ksp.getPC()
                    pc.setType("lu")
                    pc.setFactorSolverType("mumps")
                    emit(comm, f"cycle={cycle} before KSP.setUp")
                    ksp.setUp()
                    emit(comm, f"cycle={cycle} after KSP.setUp")
                    if not factor_only:
                        emit(comm, f"cycle={cycle} before KSP.solve")
                        ksp.solve(b, x)
                        emit(
                            comm,
                            f"cycle={cycle} after KSP.solve, "
                            f"reason={int(ksp.getConvergedReason())}, "
                            f"iterations={ksp.getIterationNumber()}",
                        )
                    emit(comm, f"cycle={cycle} before KSP.reset")
                    ksp.reset()
                    emit(comm, f"cycle={cycle} after KSP.reset")
                finally:
                    emit(comm, f"cycle={cycle} before KSP.destroy")
                    ksp.destroy()
                    emit(comm, f"cycle={cycle} after KSP.destroy")
                    x.destroy()
        finally:
            b.destroy()
            x_true.destroy()
    finally:
        matrix.destroy()
    emit(comm, "all snapshot cycles finished")


def replay_changing_operator(first_path, second_path, transitions):
    """Reuse a single KSP/PC and a single Mat across evolving sparsity.

    In contrast with replay(), this recreates the important Newton
    lifecycle: the solver has already factored A_001 when SCRIMP
    reassembles the SAME PETSc Mat handle with A_002's new structure.
    """
    comm = PETSc.COMM_WORLD
    if not Path(first_path).is_file() or not Path(second_path).is_file():
        raise FileNotFoundError(
            f"Snapshots must both exist: {first_path}; {second_path}"
        )
    initial = load_matrix(first_path, comm)
    later = load_matrix(second_path, comm)
    # Explicitly take a new mutable Mat; the two snapshots stay constant.
    working = initial.copy()
    ksp = PETSc.KSP().create(comm=comm)
    x_true = working.createVecRight()
    rhs = working.createVecLeft()
    x = working.createVecRight()
    check = working.createVecLeft()
    try:
        if initial.getSize() != later.getSize():
            raise ValueError("Snapshot dimensions differ")
        x_true.set(1.0)
        ksp.setType("preonly")
        pc = ksp.getPC()
        pc.setType("lu")
        pc.setFactorSolverType("mumps")
        ksp.setOperators(working)
        structure = PETSc.Mat.Structure.DIFFERENT_NONZERO_PATTERN

        emit(
            comm,
            f"persistent KSP: {initial.getSize()}, "
            f"first norm={initial.norm(PETSc.NormType.FROBENIUS):.12g}, "
            f"next norm={later.norm(PETSc.NormType.FROBENIUS):.12g}",
        )
        for index in range(transitions):
            source = initial if index % 2 == 0 else later
            if index > 0:
                # Keep the exact Mat handle referenced by KSP while
                # its sparsity changes, as happens during SWE Newton.
                emit(
                    comm,
                    f"transition={index} before in-place Mat.copy "
                    f"(pattern={'001' if index % 2 == 0 else '002'})",
                )
                source.copy(working, structure=structure)
                working.assemble()
                ksp.setOperators(working)
                emit(comm, f"transition={index} after in-place Mat.copy")

            working.mult(x_true, rhs)
            x.set(0.0)
            emit(comm, f"transition={index} before KSP.solve")
            ksp.solve(rhs, x)
            reason = int(ksp.getConvergedReason())
            working.mult(x, check)
            check.axpy(-1.0, rhs)
            relative_residual = check.norm() / max(rhs.norm(), 1.0)
            emit(
                comm,
                f"transition={index} after KSP.solve "
                f"reason={reason}, iterations={ksp.getIterationNumber()}, "
                f"relative_linear_residual={relative_residual:.9g}",
            )
            if reason <= 0 or relative_residual > 1e-5:
                raise RuntimeError(
                    "Persistent KSP failed after a nonlinear matrix update"
                )

        emit(comm, "persistent KSP before KSP.reset")
        ksp.reset()
        emit(comm, "persistent KSP after KSP.reset")
    finally:
        # If a PETSc exception was raised, destruction is still required.
        emit(comm, "persistent KSP before KSP.destroy")
        ksp.destroy()
        emit(comm, "persistent KSP after KSP.destroy")
        for obj in (x_true, rhs, x, check, working, later, initial):
            obj.destroy()
    emit(comm, "all persistent KSP transitions finished")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix", required=True,
                        help="PETSc binary Jacobian captured on OLD env")
    parser.add_argument("--cycles", type=int, default=3)
    parser.add_argument(
        "--next-matrix", default=None,
        help="Second PETSc snapshot; when provided, reuse a single "
             "KSP and Mat while alternating old/new Jacobian patterns",
    )
    parser.add_argument(
        "--factor-only", action="store_true",
        help="Set up/factor and reset, without KSP.solve()",
    )
    args = parser.parse_args()
    if args.cycles <= 0:
        parser.error("--cycles must be positive")
    if args.next_matrix:
        if args.factor_only:
            parser.error("--factor-only cannot be used with --next-matrix")
        replay_changing_operator(
            args.matrix, args.next_matrix, args.cycles
        )
    else:
        replay(args.matrix, args.cycles, args.factor_only)


if __name__ == "__main__":
    main()
