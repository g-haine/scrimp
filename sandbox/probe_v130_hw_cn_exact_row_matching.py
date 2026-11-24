# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Exact left-row permutation of the ORIGINAL first Heat-Wave CN Jacobian.

No FE reformulation, diagonal shift, pivot replacement, or relaxed KSP
threshold. A complete bipartite matching on truly NONZERO PETSc
coefficients defines the row permutation P. Solve P A x = P b with the
same physical solution x, rebuilding P A via PETSc COO for this
standalone diagnostic only. KSP/PC factorization uses fresh MPI1
instances (one Docker process per case).

Original A's CSR checksum, sparse values, and numerical action are
verified independently. Full PETSc TS and pH energy are NOT qualified.
"""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import maximum_bipartite_matching
from petsc4py import PETSc

from replay_base import matrix_digest

ORIGINAL_SHA = "c94c69c4e4307d7a5aa61d9eba186782c0495b49431cadcb7281d7a7b5bffaff"
N = 7689


def row_match(original):
    """Return a complete, nonzero-diagonal left permutation.

    maximum_bipartite_matching(..., perm_type='row') supplies for
    every COLUMN j the old ROW p[j] matching its nonzero coefficient.
    Thus B[j, :] = A[p[j], :] has B[j, j] != 0.
    """
    indptr, indices, data = original.getValuesCSR()
    A = csr_matrix(
        (np.asarray(data).copy(), np.asarray(indices).copy(),
         np.asarray(indptr).copy()),
        shape=(N, N),
    )
    A.sort_indices()
    graph = A.copy()
    graph.eliminate_zeros()  # exclude PETSc's 126 stored zero diagonal entries
    rows = np.asarray(
        maximum_bipartite_matching(graph, perm_type="row"), dtype=np.int64
    )
    if (len(rows) != N or np.any(rows < 0)
            or not np.array_equal(np.sort(rows), np.arange(N))):
        raise RuntimeError("No complete matching of nonzero numerical entries")
    B = A[rows, :].tocsr()
    diagonal = np.abs(B.diagonal())
    if np.any(diagonal == 0) or not np.all(np.isfinite(diagonal)):
        raise RuntimeError("Bipartite matched diagonal contains a zero")
    if A.nnz != B.nnz or not np.array_equal(
        np.sort(A.data), np.sort(B.data)
    ):
        raise RuntimeError("Structural row permutation altered source values")
    return A, B, rows, {
        "row_permutation_sha256":hashlib.sha256(
            np.asarray(rows, dtype="<i8").tobytes()
        ).hexdigest(),
        "rows_moved":int(np.count_nonzero(rows != np.arange(N))),
        "original_stored_nonzeros":int(A.nnz),
        "matched_stored_nonzeros":int(B.nnz),
        "original_numeric_zero_diagonal":int(np.count_nonzero(
            A.diagonal() == 0
        )),
        "matched_numeric_zero_diagonal":int(np.count_nonzero(
            B.diagonal() == 0
        )),
        "matched_smallest_absolute_diagonal":float(np.min(diagonal)),
        "matched_median_absolute_diagonal":float(np.median(diagonal)),
        "matching_uses_only_numerically_nonzero_entries":True,
        "bipartite_matching_complete":True,
    }


def petsc_from_permuted_coo(csr):
    """Diagnostic PETSc AIJ from COO: same values, row-permuted only."""
    B = PETSc.Mat().createAIJ(size=(N, N), comm=PETSc.COMM_WORLD)
    i = np.repeat(
        np.arange(N, dtype=PETSc.IntType), np.diff(csr.indptr)
    )
    j = np.ascontiguousarray(csr.indices, dtype=PETSc.IntType)
    v = np.ascontiguousarray(csr.data, dtype=PETSc.ScalarType)
    try:
        B.setPreallocationCOO(i, j)
        B.setValuesCOO(v, addv=PETSc.InsertMode.INSERT_VALUES)
        B.assemble()
        ia, ja, va = B.getValuesCSR()
        reread = csr_matrix((va.copy(), ja.copy(), ia.copy()),shape=(N,N))
        diff = reread-csr
        if diff.nnz and float(np.max(np.abs(diff.data))) > 1e-13:
            raise RuntimeError("PETSc COO changed row-permuted matrix entries")
        return B
    except BaseException:
        B.destroy()
        raise


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--matrix",type=Path,required=True)
    ap.add_argument("--output",type=Path,required=True)
    ap.add_argument("--mode",choices=("original","row_matching"),required=True)
    ap.add_argument("--solver",choices=("mumps","superlu_dist"),required=True)
    ap.add_argument("--colperm",
                    choices=("default","NATURAL","MMD_AT_PLUS_A",
                             "METIS_AT_PLUS_A"),default="default")
    args = ap.parse_args()
    if PETSc.COMM_WORLD.getSize() != 1:
        raise RuntimeError("Native row-matching diagnostic is MPI1 only")
    if not args.matrix.is_file() or args.output.exists():
        raise RuntimeError("Require exact original matrix and fresh evidence")
    if args.solver=="mumps" and args.colperm!="default":
        raise RuntimeError("SuperLU_DIST column options not valid for MUMPS")
    PETSc.Options()["mat_superlu_dist_replacetinypivot"]="false"
    PETSc.Options()["mat_superlu_dist_rowperm"]="NOROWPERM"
    if args.colperm != "default":
        PETSc.Options()["mat_superlu_dist_colperm"]=args.colperm

    viewer=PETSc.Viewer().createBinary(
        str(args.matrix),mode="r",comm=PETSc.COMM_WORLD
    )
    A=PETSc.Mat().create(comm=PETSc.COMM_WORLD)
    try:
        A.setType("aij")
        A.load(viewer)
    finally:
        viewer.destroy()

    B = xtrue = b0 = rhs = candidate_rhs = result = check = ksp = None
    try:
        if A.getSize()!=(N,N) or matrix_digest(A)!=ORIGINAL_SHA:
            raise RuntimeError("Original first CN Jacobian does not match pinned hash")
        A_sp, B_sp, rows, match = row_match(A)
        if args.mode=="row_matching":
            B=petsc_from_permuted_coo(B_sp)
            mat=B
            reordered=rows
        else:
            mat=A
            reordered=np.arange(N,dtype=np.int64)

        xtrue=A.createVecRight()
        b0=A.createVecLeft()
        rhs=mat.createVecLeft()
        candidate_rhs=mat.createVecLeft()
        result=mat.createVecRight()
        check=A.createVecLeft()
        order=np.arange(N,dtype=PETSc.IntType)
        xtrue.setValues(
            order,
            1.+0.05*np.sin(np.arange(N,dtype=float)),
            addv=PETSc.InsertMode.INSERT_VALUES,
        )
        xtrue.assemble()
        A.mult(xtrue,b0)
        old_rhs=np.asarray(b0.getArray(readonly=True)).copy()
        rhs.setValues(
            order,old_rhs[reordered],
            addv=PETSc.InsertMode.INSERT_VALUES,
        )
        rhs.assemble()
        mat.mult(xtrue,candidate_rhs)
        candidate_rhs.axpy(-1.,rhs)
        action_error=float(
            candidate_rhs.norm()/max(float(rhs.norm()),1e-30)
        )
        if action_error>1e-12:
            raise RuntimeError("PAx != P(Ax), invalid row permutation of operator")
        if not np.array_equal(
            np.sort(A_sp.data),np.sort(B_sp.data)
        ):
            raise RuntimeError("Matching altered original matrix values")
        if matrix_digest(A)!=ORIGINAL_SHA:
            raise RuntimeError("Source matrix mutated before factorization")
        ksp=PETSc.KSP().create(comm=PETSc.COMM_WORLD)
        ksp.setOperators(mat)
        ksp.setType(PETSc.KSP.Type.PREONLY)
        pc=ksp.getPC()
        pc.setType(PETSc.PC.Type.LU)
        pc.setFactorSolverType(args.solver)
        ksp.setFromOptions()
        error=None
        try:
            ksp.setUp()
            ksp.solve(rhs,result)
        except PETSc.Error as exc:
            error=str(exc)
        reason=int(ksp.getConvergedReason())
        rel_resid=rel_state=None
        if reason>0 and error is None:
            A.mult(result,check)
            check.axpy(-1.,b0)
            rel_resid=float(check.norm()/max(float(b0.norm()),1e-30))
            delta=result.copy()
            try:
                delta.axpy(-1.,xtrue)
                rel_state=float(
                    delta.norm()/max(float(xtrue.norm()),1e-30)
                )
            finally:
                delta.destroy()
        report={
            "schema":"scrimp-v130-original-CN-exact-row-permutation-v1",
            "status":"DIAGNOSTIC_REPLAY_NOT_PHYSICAL_TS",
            "mode":args.mode, "solver":args.solver,
            "colperm_requested":args.colperm,
            "rowperm_requested":"NOROWPERM",
            "original_CSR_SHA256_before":ORIGINAL_SHA,
            "original_CSR_SHA256_after":matrix_digest(A),
            "row_matched_CSR_SHA256":matrix_digest(B) if B is not None else None,
            "matched_graph":match,
            "operator_permutes_only_equations_not_unknowns":True,
            "relative_PA_action_error":action_error,
            "KSP_convergence_reason":reason,
            "PC_failure_reason":str(pc.getFailedReason()),
            "KSP_iterations":int(ksp.getIterationNumber()),
            "relative_residual_in_ORIGINAL_equations":rel_resid,
            "relative_known_state_error":rel_state,
            "PETSc_error":error,
            "tiny_pivot_replacement":False,
            "diagonal_shift_or_regularization":False,
            "main_CN_TS_completed":False,
            "Hamiltonian_power_balance_qualified":False,
        }
        if report["original_CSR_SHA256_after"]!=ORIGINAL_SHA:
            raise RuntimeError("Source changed during factorization")
        good = (
            reason>0 and error is None
            and rel_resid is not None and rel_resid<1e-9
            and rel_state is not None and rel_state<1e-6
        )
        report["status"]=(
            "PASS_PERMUTED_JACOBIAN_KNOWN_STATE_LU_ONLY"
            if good else
            "FACTOR_OR_ACCURACY_FAILURE_NOT_PHYSICAL_PASS"
        )
        args.output.write_text(
            json.dumps(report,indent=2)+"\n", encoding="utf-8"
        )
        print(
            f"mode={args.mode} solver={args.solver} colperm={args.colperm} "
            f"status={report['status']} KSP={reason} "
            f"PC={report['PC_failure_reason']} "
            f"original_residual={rel_resid} state_error={rel_state} "
            f"moved_rows={match['rows_moved']} "
            f"min_matching_diag={match['matched_smallest_absolute_diagonal']:.6g}",
            flush=True
        )
    finally:
        if ksp is not None:
            ksp.destroy()
        for obj in (check,result,candidate_rhs,rhs,b0,xtrue,B,A):
            if obj is not None:
                obj.destroy()


if __name__=="__main__":
    main()
