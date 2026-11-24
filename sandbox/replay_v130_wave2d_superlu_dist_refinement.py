# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# See LICENSE for license information.

"""Refine a SuperLU_DIST tiny-pivot solution on the ORIGINAL 2D CN matrix.

The prior independent replay confirmed the nonregularized SuperLU_DIST
factorization FAILS on MPI1/2 (PC_FACTOR_NUMERIC_ZEROPIVOT), while
ReplaceTinyPivot makes PREONLY converge with inadequate accuracy
(relative backward error ~5.08e-10, known-solution error ~7.3e-7).

On exactly the same SHA256-pinned factor matrix, try:
(1) native SuperLU_DIST IterRefine (preonly, option "iterrefine"),
(2) right-preconditioned GMRES with the original LU, checking the TRUE
    unpreconditioned residual and the known solution.

Both methods require the opt-in ReplaceTinyPivot; both are diagnostic
only. No GetFEM reassembly, changes to the Jacobian, shift, algebraic
port, physical energy or any SCRIMP implementation.
"""

import argparse
import hashlib
import json
from pathlib import Path
import traceback

import numpy as np
from mpi4py import MPI
from petsc4py import PETSc


N = 2230
EXPECTED_HASH = "24346856f358e4f81cf332ebc974e329dc084c60f5d6b98bfbcd788fdff37f54"


def load_original(path):
    with np.load(path, allow_pickle=False) as archived:
        values = {key: archived[key].copy() for key in archived.files}
    expected = {f"{m}_{k}" for m in ("P", "M", "K")
                for k in ("indptr", "indices", "data")}
    if set(values) != expected:
        raise RuntimeError("Not the original nine native CSR arrays")
    digest = hashlib.sha256()
    for key in sorted(values):
        digest.update(key.encode())
        digest.update(np.ascontiguousarray(values[key]).tobytes())
    if digest.hexdigest() != EXPECTED_HASH:
        raise RuntimeError("Wrong original M/K/P CSR SHA256 — abort")
    indptr = values["P_indptr"]
    indices = values["P_indices"]
    data = values["P_data"]
    if (indptr.shape != (N+1,)
            or int(indptr[0]) != 0
            or int(indptr[-1]) != len(indices)
            or len(indices) != len(data)
            or not np.all(indices >= 0)
            or not np.all(indices < N)
            or not np.all(np.isfinite(data))
            or np.max(np.abs(np.imag(data))) != 0):
        raise RuntimeError("Invalid original first-CN Jacobian CSR")
    return indptr, indices, data, digest.hexdigest()


def build_exact_owned_aij(comm, indptr, indices, data):
    # Each rank constructs only its PETSc-owned rows; no fake sequential
    # matrix cloned across ranks. Preserve all original numeric entries,
    # including the ten explicit zero diagonals of the Y_L DAE block.
    matrix = PETSc.Mat().createAIJ(
        size=(N, N), nnz=96, comm=PETSc.COMM_WORLD,
    )
    matrix.setOption(PETSc.Mat.Option.IGNORE_ZERO_ENTRIES, False)
    first, last = matrix.getOwnershipRange()
    for row in range(first, last):
        start = int(indptr[row])
        stop = int(indptr[row+1])
        if stop > start:
            matrix.setValues(
                row, np.asarray(indices[start:stop], dtype=np.int32),
                np.asarray(data[start:stop], dtype=PETSc.ScalarType),
            )
    matrix.assemble()
    # Verify the STORED graph AND each original coefficient. Exactly
    # ten explicit-zero diagonal entries belong to Y_L. Losing or
    # introducing structural zeros changes the SuperLU_DIST test!
    rp, ci, vv = matrix.getValuesCSR()
    if len(rp) != (last-first)+1:
        raise RuntimeError("Unexpected rank-local PETSc CSR row layout")
    row_max = 0.0
    graph_mismatches = 0
    for row in range(first, last):
        local_index = row-first
        a, b = int(indptr[row]), int(indptr[row+1])
        c, d = int(rp[local_index]), int(rp[local_index+1])
        original_order = np.argsort(indices[a:b])
        local_order = np.argsort(ci[c:d])
        orig_cols = indices[a:b][original_order]
        local_cols = ci[c:d][local_order]
        if not np.array_equal(orig_cols, local_cols):
            graph_mismatches += 1
            continue
        if len(orig_cols):
            delta = np.abs(
                np.asarray(vv[c:d])[local_order] -
                np.asarray(data[a:b])[original_order]
            )
            row_max = max(row_max, float(np.max(delta)))
    graph_errors = comm.allreduce(graph_mismatches, op=MPI.SUM)
    global_max = comm.allreduce(row_max, op=MPI.MAX)
    if graph_errors or global_max != 0.:
        raise RuntimeError(
            "Actual PETSc CSR differs from original factorized operator: "
            f"graph_mismatched_rows={graph_errors}, max_delta={global_max}"
        )
    return matrix, first, last


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix", type=Path, required=True)
    parser.add_argument("--refinement",
                        choices=("internal_refine", "right_gmres"),
                        required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    indptr, indices, data, checksum = load_original(args.matrix)
    mat, row_lo, row_hi = build_exact_owned_aij(
        comm, indptr, indices, data
    )
    known = mat.createVecRight()
    offset = np.arange(row_lo, row_hi, dtype=np.float64)
    known.setValues(
        np.arange(row_lo, row_hi, dtype=np.int32),
        np.sin(.013 * offset) + .25,
    )
    known.assemble()
    rhs = mat.createVecLeft()
    mat.mult(known, rhs)
    solution = mat.createVecRight()
    solution.set(0.)

    # One explicitly opt-in pivot replacement, identical to the prior
    # replay. Investigate ACCURACY recovery, not a different physical
    # operator or a hidden fallback to sequential LU/MUMPS.
    ksp = PETSc.KSP().create(comm=PETSc.COMM_WORLD)
    ksp.setOperators(mat)
    if args.refinement == "internal_refine":
        # SuperLU_DIST provides built-in iterative refinement (PETSc
        # option uses a bool, selects DOUBLE refinement internally).
        ksp.setType("preonly")
    else:
        # GMRES with the TRUE, unpreconditioned residual. LEFT
        # preconditioned residual could otherwise declare convergence
        # while original Px-b is still inaccurate.
        ksp.setType("gmres")
        ksp.setPCSide(PETSc.PC.Side.RIGHT)
        ksp.setNormType(PETSc.KSP.NormType.UNPRECONDITIONED)
        ksp.setTolerances(rtol=1e-12, atol=1e-14, max_it=30)
    pc = ksp.getPC()
    pc.setType("lu")
    pc.setFactorSolverType("superlu_dist")
    options = PETSc.Options()
    options.setValue("mat_superlu_dist_replacetinypivot", "true")
    if args.refinement == "internal_refine":
        options.setValue("mat_superlu_dist_iterrefine", "true")
    ksp.setFromOptions()
    # Reassert explicit true-residual convergence and preserved PETSc
    # PC_RIGHT after PETSc has processed the options database.
    if args.refinement == "right_gmres":
        ksp.setPCSide(PETSc.PC.Side.RIGHT)
        ksp.setNormType(PETSc.KSP.NormType.UNPRECONDITIONED)
        ksp.setTolerances(rtol=1e-12, atol=1e-14, max_it=30)

    error = None
    try:
        ksp.solve(rhs, solution)
    except Exception as exc:
        error = f"{type(exc).__name__}: {exc}"
    reason = int(ksp.getConvergedReason())
    failed = str(pc.getFailedReason())
    try:
        factor_solver = str(pc.getFactorSolverType()).lower()
    except Exception:
        factor_solver = "unavailable"

    rel_backward = None
    rel_solution = None
    if error is None and reason > 0:
        check = mat.createVecLeft()
        mat.mult(solution, check)
        check.axpy(-1., rhs)
        residual = float(check.norm(PETSc.NormType.INFINITY))
        rhs_norm = float(rhs.norm(PETSc.NormType.INFINITY))
        x_norm = float(solution.norm(PETSc.NormType.INFINITY))
        a_row_sum = np.bincount(
            np.repeat(np.arange(N), np.diff(indptr)),
            weights=np.abs(data), minlength=N
        )
        a_infinity = float(np.max(a_row_sum))
        rel_backward = residual / max(
            rhs_norm + a_infinity*x_norm, 1e-300
        )
        diff = solution.copy()
        diff.axpy(-1., known)
        rel_solution = float(diff.norm(PETSc.NormType.INFINITY)) / max(
            float(known.norm(PETSc.NormType.INFINITY)), 1e-300
        )
        check.destroy()
        diff.destroy()
    # Do not mistake a solver that ignored PC_RIGHT or the true
    # unpreconditioned norm for a valid GMRES correction. PREONLY
    # deliberately has no Krylov iterations; GMRES must have at least 1.
    right_true_norm = bool(
        ksp.getPCSide() == PETSc.PC.Side.RIGHT
        and ksp.getNormType() == PETSc.KSP.NormType.UNPRECONDITIONED
    )
    correct_method = bool(
        (args.refinement == "internal_refine"
         and str(ksp.getType()).lower() == "preonly")
        or (args.refinement == "right_gmres"
            and str(ksp.getType()).lower() == "gmres"
            and right_true_norm
            and int(ksp.getIterationNumber()) > 0)
    )
    candidate_ok = (
        error is None
        and reason > 0
        and factor_solver == "superlu_dist"
        and correct_method
        and rel_backward is not None
        and np.isfinite(rel_backward)
        and rel_backward < 1e-11
        and np.isfinite(rel_solution)
        and rel_solution < 1e-7
    )
    all_ok = bool(comm.allreduce(bool(candidate_ok), op=MPI.LAND))
    result = {
        "status": (
            "PASS_INDEPENDENT_ORIGINAL_PETSC_PC_REFINED_REPLAY_ONLY"
            if all_ok else "FAIL_OR_UNQUALIFIED_PETSC_PC_REFINEMENT"
        ),
        "diagnostic_only": True,
        "is_not_original_wave_CN_physics_qualification": True,
        "mode": args.refinement,
        "MPI": comm.Get_size(),
        "rank": rank,
        "owned_rows": [row_lo, row_hi],
        "matrix_n": N,
        "original_full_P_M_K_SHA256": checksum,
        "PETSc_owned_CSR_graph_and_values_identical_to_archive": True,
        "original_J_not_modified": True,
        "ksp_type": str(ksp.getType()).lower(),
        "ksp_iterations": int(ksp.getIterationNumber()),
        "ksp_norm_type": str(ksp.getNormType()),
        "ksp_pc_side": str(ksp.getPCSide()),
        "right_gmres_true_residual_norm_confirmed": (
            right_true_norm if args.refinement == "right_gmres" else None
        ),
        "correct_refinement_method_confirmed": correct_method,
        "ksp_convergence_tolerances": tuple(
            float(x) if i < 3 else int(x)
            for i,x in enumerate(ksp.getTolerances())
        ),
        "replace_tiny_pivot_requested": True,
        "internal_iterative_refinement_requested": (
            args.refinement == "internal_refine"
        ),
        "PC_type": "lu",
        "factor_solver": factor_solver,
        "KSP_reason": reason,
        "PC_failed": failed,
        "exception": error,
        "relative_backward_residual": rel_backward,
        "relative_known_solution_error": rel_solution,
        "factorization_passed_all_ranks": all_ok,
    }
    (args.output/f"rank{rank}.json").write_text(
        json.dumps(result, indent=2)+"\n", encoding="utf-8",
    )
    if rank == 0:
        print(json.dumps(result, indent=2), flush=True)
    ksp.destroy()
    rhs.destroy()
    known.destroy()
    solution.destroy()
    mat.destroy()
    if not all_ok:
        raise SystemExit(1)


if __name__ == "__main__":
    try:
        main()
    except Exception:
        traceback.print_exc()
        raise
