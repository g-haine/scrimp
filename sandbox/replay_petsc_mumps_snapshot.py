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


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix", required=True,
                        help="PETSc binary Jacobian captured on OLD env")
    parser.add_argument("--cycles", type=int, default=3)
    parser.add_argument(
        "--factor-only", action="store_true",
        help="Set up/factor and reset, without KSP.solve()",
    )
    args = parser.parse_args()
    if args.cycles <= 0:
        parser.error("--cycles must be positive")
    replay(args.matrix, args.cycles, args.factor_only)


if __name__ == "__main__":
    main()
