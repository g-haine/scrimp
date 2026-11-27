# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp


"""Analyze actual GetFEM/PETSc CN factor matrices independently with SciPy.

This program does not alter matrices to make SuperLU_DIST run. It only
reads the exact main-step PC operator P, mass M, stiffness K saved
before factorization for the MUMPS and SuperLU_DIST runs.

The principal scientific question is whether the TEN zero diagonal
entries of the mixed port-Hamiltonian saddle matrix imply structural
rank deficiency, or whether a nonregularized pivoted sparse LU can
factor the actual matrix. A successful SciPy LU is evidence for the
latter, NOT proof that PETSc SuperLU_DIST has been repaired.
"""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import structural_rank
from scipy.sparse.linalg import splu


N_EXPECTED = 2230


def load_snapshot(directory):
    directory = Path(directory)
    meta = json.loads(
        (directory/"first_main_CN_factor.json").read_text(encoding="utf-8")
    )
    if (meta.get("schema")
            != "scrimp-v130-first-main-CN-exact-PC-factor-matrix-v1"):
        raise ValueError("Not a qualified actual-factor-matrix snapshot")
    matrices = {}
    with np.load(directory/"first_main_CN_factor.npz",
                 allow_pickle=False) as archive:
        for name in ("P", "K", "M"):
            rowptr = archive[f"{name}_indptr"]
            indices = archive[f"{name}_indices"]
            values = archive[f"{name}_data"]
            dims = tuple(meta["matrix_shapes"][name])
            if (dims != (N_EXPECTED, N_EXPECTED)
                    or rowptr.shape != (N_EXPECTED+1,)
                    or not np.all(np.isfinite(values))):
                raise ValueError("Wrong native CSR dimensions/data for "+name)
            matrix = csr_matrix((values, indices, rowptr), shape=dims)
            matrix.sum_duplicates()
            matrices[name] = matrix
    if (abs(float(meta["TS_shift"])-250.) >= 1e-7
            or abs(float(meta["time"])-.008) >= 1e-9):
        raise ValueError("Not the genuine initial CN shift=250 matrix")
    return meta, matrices


def numeric_max(matrix):
    return float(np.max(np.abs(matrix.data))) if matrix.nnz else 0.


def evaluate_snapshot(meta, matrices):
    P, K, M = (matrices[x] for x in ("P", "K", "M"))
    sigma = float(meta["TS_shift"])
    missing = P - (K+sigma*M)
    missing.eliminate_zeros()
    assembly_error = numeric_max(missing)
    tol = 1e-10*max(1., numeric_max(P))
    # Actual stored numerical entries, not merely allocated zero edges.
    graph = P.copy()
    graph.eliminate_zeros()
    structural = int(structural_rank(graph))
    zero_diag = np.flatnonzero(P.diagonal() == 0).astype(int).tolist()
    zero_rows = np.flatnonzero(np.diff(graph.indptr) == 0).astype(int).tolist()
    intervals = meta["getfem_variable_intervals"]
    owners = {
        str(row):[
            v["name"] for v in intervals
            if v["lo"] <= row < v["hi"]
        ] for row in zero_diag
    }
    report = {
        "selected_solver": meta["selected_solver"],
        "shape": list(P.shape),
        "stored_nnz": int(P.nnz),
        "numerical_nnz": int(graph.nnz),
        "structural_rank": structural,
        "structural_full_rank": structural == N_EXPECTED,
        "numerical_zero_diagonal_ids": zero_diag,
        "zero_diagonal_variable_owners": owners,
        "zero_numerical_rows": zero_rows,
        "max_abs_P_minus_K_plus_sigma_M": assembly_error,
        "assembly_identity_tol": tol,
        "P_is_original_K_plus_shift_M": assembly_error <= tol,
        "exact_csr_sha256": meta["csr_SHA256"],
        "scipy_unregularized_lu": {},
    }
    # One common RHS, deliberately constructed from an exact
    # known solution. This probes a real pivoted factorization
    # independently of any time-integration nonlinear solver.
    x_true = np.sin(np.arange(N_EXPECTED, dtype=float)*.013)+.25
    b = P @ x_true
    for ordering in ("COLAMD", "MMD_AT_PLUS_A"):
        try:
            lu = splu(P.tocsc(), permc_spec=ordering,
                      diag_pivot_thresh=1.0)
            x = lu.solve(b)
            residual = P @ x - b
            norm = float(np.linalg.norm(residual, ord=np.inf))
            denominator = (float(np.max(np.abs(b))) +
                           numeric_max(P)*float(np.max(np.abs(x))) + 1e-30)
            rel = norm / denominator
            diag = np.abs(lu.U.diagonal())
            report["scipy_unregularized_lu"][ordering] = {
                "factorized": True,
                "U_min_abs_pivot": float(np.min(diag)),
                "U_max_abs_pivot": float(np.max(diag)),
                "relative_residual": rel,
                "physical_model_unchanged": True,
            }
        except Exception as error:
            report["scipy_unregularized_lu"][ordering] = {
                "factorized": False,
                "error": f"{type(error).__name__}: {error}",
            }
    return report


def compare_snapshots(left, right):
    result = {}
    for name in ("P", "K", "M"):
        A, B = left[name], right[name]
        if A.shape != B.shape:
            raise ValueError("Different physical Jacobian shape: "+name)
        diff = (A - B).tocsr()
        diff.eliminate_zeros()
        norm = numeric_max(diff)
        result[name] = {
            "max_abs_difference": norm,
            "same_to_numerical_precision":
                norm <= 1e-12*max(1., numeric_max(A), numeric_max(B)),
            "nnz_of_actual_difference": int(diff.nnz),
        }
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mumps", type=Path, required=True)
    parser.add_argument("--superlu-dist", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    m_meta, m_mat = load_snapshot(args.mumps)
    d_meta, d_mat = load_snapshot(args.superlu_dist)
    if (m_meta["selected_solver"] != "mumps"
            or d_meta["selected_solver"] != "superlu_dist"):
        raise ValueError("Backend names in original PETSc snapshots disagree")
    result = {
        "schema": "scrimp-v130-wave2d-sparse-Jacobian-diagnostic-v1",
        "status": "UNQUALIFIED_MATRIX_DIAGNOSTIC",
        "diagnostic_only": True,
        "does_not_claim_superlu_dist_backend_success": True,
        "mumps": evaluate_snapshot(m_meta, m_mat),
        "superlu_dist": evaluate_snapshot(d_meta, d_mat),
        "backend_matrix_parity": compare_snapshots(m_mat, d_mat),
    }
    checks = (
        all(result[x]["P_is_original_K_plus_shift_M"]
            and not result[x]["zero_numerical_rows"]
            for x in ("mumps", "superlu_dist"))
        and all(e["same_to_numerical_precision"]
                for e in result["backend_matrix_parity"].values())
    )
    # Report whether mathematically identical operators have a
    # perfect structural matching; only mark that fact, NEVER
    # infer a passing PETSc SuperLU_DIST run.
    if checks:
        result["status"] = "PASS_IDENTICAL_ACTUAL_MAIN_CN_OPERATORS_DIAGNOSTIC_ONLY"
    else:
        result["status"] = "FAIL_MATRIX_PARITY_OR_ASSEMBLY"
    result["structural_rank_not_numerical_rank"] = True
    args.output.write_text(json.dumps(result, indent=2)+"\n", encoding="utf-8")
    print(json.dumps(result, indent=2), flush=True)
    if not checks:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
