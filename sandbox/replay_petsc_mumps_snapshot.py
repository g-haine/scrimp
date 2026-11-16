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
import json
from pathlib import Path

from mpi4py import MPI
import numpy as np
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



def inspect_two_matrices(first_path, second_path):
    """Compare actual stored CSR and numerical values without factorization.

    The matrices are loaded independently: PETSc MatCopy cannot hide a
    removed nonzero by retaining the old destination's allocated slots.
    Report exact rank-local added/removed locations and the Frobenius
    difference of the matrices, reduced across all MPI ranks.
    """
    comm = PETSc.COMM_WORLD
    mpi = MPI.COMM_WORLD
    if not Path(first_path).is_file() or not Path(second_path).is_file():
        raise FileNotFoundError(
            f"Both snapshot files are required: {first_path}; {second_path}"
        )
    first = load_matrix(first_path, comm)
    second = load_matrix(second_path, comm)
    try:
        if (first.getSize() != second.getSize()
                or first.getOwnershipRange() != second.getOwnershipRange()):
            raise ValueError("Snapshot sizes or distributed ownership differ")
        ra, ca, va = first.getValuesCSR()
        rb, cb, vb = second.getValuesCSR()
        if len(ra) != len(rb):
            raise ValueError("Snapshot local row counts differ")

        added = removed = 0
        rows_pattern_changed = 0
        local_delta_sq = 0.0
        local_nonfinite_a = int(np.count_nonzero(~np.isfinite(va)))
        local_nonfinite_b = int(np.count_nonzero(~np.isfinite(vb)))
        for row in range(len(ra) - 1):
            ia0, ia1 = int(ra[row]), int(ra[row + 1])
            ib0, ib1 = int(rb[row]), int(rb[row + 1])
            cols_a, cols_b = ca[ia0:ia1], cb[ib0:ib1]
            vals_a, vals_b = va[ia0:ia1], vb[ib0:ib1]
            _shared, inda, indb = np.intersect1d(
                cols_a, cols_b, assume_unique=True, return_indices=True
            )
            row_added = len(cols_b) - len(indb)
            row_removed = len(cols_a) - len(inda)
            added += row_added
            removed += row_removed
            if row_added or row_removed:
                rows_pattern_changed += 1
            # Compare aligned values, including entries present in only one
            # of the two matrices, without any in-place PETSc operations.
            keep_a = np.ones(len(cols_a), dtype=bool)
            keep_b = np.ones(len(cols_b), dtype=bool)
            keep_a[inda] = False
            keep_b[indb] = False
            local_delta_sq += float(np.sum(np.abs(vals_a[inda] - vals_b[indb]) ** 2))
            local_delta_sq += float(np.sum(np.abs(vals_a[keep_a]) ** 2))
            local_delta_sq += float(np.sum(np.abs(vals_b[keep_b]) ** 2))

        norm_a = float(first.norm(PETSc.NormType.FROBENIUS))
        norm_b = float(second.norm(PETSc.NormType.FROBENIUS))
        diff_frob = float(np.sqrt(mpi.allreduce(local_delta_sq, op=MPI.SUM)))
        result = {
            "mode": "inspect_only_no_factorization",
            "mpi_size": mpi.Get_size(),
            "first": str(first_path),
            "second": str(second_path),
            "global_size": first.getSize(),
            "stored_nnz_first": int(mpi.allreduce(len(ca), op=MPI.SUM)),
            "stored_nnz_second": int(mpi.allreduce(len(cb), op=MPI.SUM)),
            "stored_nnz_added": int(mpi.allreduce(added, op=MPI.SUM)),
            "stored_nnz_removed": int(mpi.allreduce(removed, op=MPI.SUM)),
            "rows_with_pattern_change": int(
                mpi.allreduce(rows_pattern_changed, op=MPI.SUM)
            ),
            "stored_zeros_first": int(
                mpi.allreduce(np.count_nonzero(va == 0), op=MPI.SUM)
            ),
            "stored_zeros_second": int(
                mpi.allreduce(np.count_nonzero(vb == 0), op=MPI.SUM)
            ),
            "nonfinite_first": int(
                mpi.allreduce(local_nonfinite_a, op=MPI.SUM)
            ),
            "nonfinite_second": int(
                mpi.allreduce(local_nonfinite_b, op=MPI.SUM)
            ),
            "frobenius_first": norm_a,
            "frobenius_second": norm_b,
            "frobenius_difference": diff_frob,
            "relative_frobenius_difference": diff_frob / max(1.0, norm_a),
        }
        if mpi.Get_rank() == 0:
            print(json.dumps(result, indent=2), flush=True)
        if result["nonfinite_first"] or result["nonfinite_second"]:
            raise RuntimeError("At least one snapshot has non-finite coefficients")
    finally:
        second.destroy()
        first.destroy()


def relative_state_error(solution, expected):
    """Check the manufactured exact state, not only the algebraic residual.

    A tiny ||A*x-b|| / ||b|| alone may hide a large state error when
    the captured, potentially divergent Jacobian is ill-conditioned.
    """
    difference = solution.duplicate()
    try:
        solution.copy(difference)
        difference.axpy(-1.0, expected)
        return difference.norm() / max(expected.norm(), 1.0)
    finally:
        difference.destroy()


def replay(path, cycles, factor_only, solver="mumps"):
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
            f"factor_only={factor_only}, cycles={cycles}, solver={solver}",
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
                    pc.setFactorSolverType(solver)
                    emit(comm, f"cycle={cycle} before KSP.setUp")
                    ksp.setUp()
                    emit(comm, f"cycle={cycle} after KSP.setUp")
                    if not factor_only:
                        emit(comm, f"cycle={cycle} before KSP.solve")
                        ksp.solve(b, x)
                        check = matrix.createVecLeft()
                        try:
                            matrix.mult(x, check)
                            check.axpy(-1.0, b)
                            relative_residual = (
                                check.norm() / max(b.norm(), 1.0)
                            )
                        finally:
                            check.destroy()
                        state_error = relative_state_error(x, x_true)
                        emit(
                            comm,
                            f"cycle={cycle} after KSP.solve, "
                            f"reason={int(ksp.getConvergedReason())}, "
                            f"iterations={ksp.getIterationNumber()}, "
                            f"relative_linear_residual="
                            f"{relative_residual:.9g}, "
                            f"relative_state_error={state_error:.9g}",
                        )
                        if (int(ksp.getConvergedReason()) <= 0
                                or not (relative_residual <= 1e-5)
                                or not (state_error <= 1e-5)):
                            raise RuntimeError(
                                "LU solver failed to reproduce the "
                                "known consistent Jacobian solution"
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


def replay_changing_operator(first_path, second_path, transitions,
                             solver="mumps", transition_reset="none"):
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
        pc.setFactorSolverType(solver)
        ksp.setOperators(working)
        structure = PETSc.Mat.Structure.DIFFERENT_NONZERO_PATTERN

        emit(
            comm,
            f"persistent KSP: {initial.getSize()}, "
            f"first norm={initial.norm(PETSc.NormType.FROBENIUS):.12g}, "
            f"next norm={later.norm(PETSc.NormType.FROBENIUS):.12g}, "
            f"solver={solver}, transition_reset={transition_reset}",
        )
        for index in range(transitions):
            source = initial if index % 2 == 0 else later
            if index > 0:
                # Reset before changing the operator: after a symbolic
                # factorization the PC may still own LU symbolic
                # structures for the OLD sparse graph. This is an
                # opt-in diagnostic; the default reproduces the bug.
                if transition_reset != "none":
                    emit(
                        comm,
                        f"transition={index} BEFORE "
                        f"{transition_reset}.reset",
                    )
                    if transition_reset == "pc":
                        ksp.getPC().reset()
                    elif transition_reset == "ksp":
                        ksp.reset()
                    else:
                        raise ValueError(transition_reset)
                    emit(
                        comm,
                        f"transition={index} AFTER "
                        f"{transition_reset}.reset",
                    )
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
            state_error = relative_state_error(x, x_true)
            emit(
                comm,
                f"transition={index} after KSP.solve "
                f"reason={reason}, iterations={ksp.getIterationNumber()}, "
                f"relative_linear_residual={relative_residual:.9g}, "
                f"relative_state_error={state_error:.9g}",
            )
            if (reason <= 0 or not (relative_residual <= 1e-5)
                    or not (state_error <= 1e-5)):
                pc = ksp.getPC()
                details = {
                    "relative_linear_residual": relative_residual,
                    "relative_state_error": state_error,
                    "KSP_reason": reason,
                    "PC_type": pc.getType(),
                    "PC_failed_reason": str(pc.getFailedReason()),
                }
                if pc.getType() == "lu" and solver == "mumps":
                    try:
                        factor = pc.getFactorMatrix()
                        details["MUMPS_INFOG(1)"] = factor.getMumpsInfog(1)
                        details["MUMPS_INFOG(2)"] = factor.getMumpsInfog(2)
                    except Exception as exc:
                        details["factor_query_error"] = repr(exc)
                emit(
                    comm,
                    f"transition={index} factor_failure={details}"
                )
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
        "--inspect-only", action="store_true",
        help="Compare exact stored CSR patterns and numeric coefficients "
             "of two binary matrices without invoking KSP or MUMPS",
    )
    parser.add_argument(
        "--solver", choices=("mumps", "superlu"), default="mumps",
        help="Compare the exact same stored Jacobian with MUMPS versus "
             "SuperLU (one rank); defaults to MUMPS",
    )
    parser.add_argument(
        "--next-matrix", default=None,
        help="Second PETSc snapshot; when provided, reuse a single "
             "KSP and Mat while alternating old/new Jacobian patterns",
    )
    parser.add_argument(
        "--transition-reset", choices=("none", "pc", "ksp"),
        default="none",
        help="For a changing matrix, reset the factor PC or full KSP "
             "BEFORE in-place Mat.copy; isolate stale sparse symbolic "
             "factorization. Default none reproduces the failure.",
    )
    parser.add_argument(
        "--factor-only", action="store_true",
        help="Set up/factor and reset, without KSP.solve()",
    )
    args = parser.parse_args()
    if args.cycles <= 0:
        parser.error("--cycles must be positive")
    if args.inspect_only:
        if not args.next_matrix:
            parser.error("--inspect-only requires --next-matrix")
        if args.factor_only or args.transition_reset != "none":
            parser.error("--inspect-only is incompatible with factor-only/reset")
        inspect_two_matrices(args.matrix, args.next_matrix)
        return
    if args.next_matrix:
        if args.factor_only:
            parser.error("--factor-only cannot be used with --next-matrix")
        replay_changing_operator(
            args.matrix, args.next_matrix, args.cycles,
            solver=args.solver, transition_reset=args.transition_reset,
        )
    else:
        if args.transition_reset != "none":
            parser.error("--transition-reset requires --next-matrix")
        replay(args.matrix, args.cycles, args.factor_only,
               solver=args.solver)


if __name__ == "__main__":
    main()
