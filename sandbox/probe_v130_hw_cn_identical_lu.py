# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Reproduce two PETSc LU factorizations from ONE immutable CN Jacobian.

The binary snapshot is recorded by production SCRIMP during the
actual original Heat-Wave first CN IJacobian callback. No GetFEM
import, no reassembly, no matrix modification, and no artificial
diagonal perturbation. A single PETSc Mat is loaded and supplied in
turn to entirely new MUMPS and SuperLU_DIST KSP/PC instances.

A manufactured consistent RHS is used only AFTER factorization for
numerical residual/state-error checks; the LU failure is RHS-independent.
The output is a diagnostic, never a physical CN integration PASS.
"""

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from petsc4py import PETSc
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import structural_rank


def matrix_digest(A):
    ptr, col, values = A.getValuesCSR()
    h = hashlib.sha256()
    for arr in (
        np.asarray(ptr, dtype="<i8"),
        np.asarray(col, dtype="<i8"),
        np.asarray(values, dtype="<c16"),
    ):
        h.update(np.ascontiguousarray(arr).tobytes())
    return h.hexdigest()


def matrix_audit(A):
    ptr, col, values = A.getValuesCSR()
    n,m=map(int,A.getSize())
    if n!=m or PETSc.COMM_WORLD.getSize()!=1:
        raise RuntimeError("Original Heat-Wave comparison is serial and square")
    if (ptr.shape!=(n+1,) or int(ptr[-1])!=len(col) or
            len(col)!=len(values)):
        raise RuntimeError("Original PETSc CSR structure invalid")
    if np.any(~np.isfinite(values)):
        raise RuntimeError("Original CN Jacobian has non-finite entries")
    M=csr_matrix((np.asarray(values,dtype=np.complex128),
                  col.astype(np.int64),
                  ptr.astype(np.int64)),shape=(n,n))
    graph=csr_matrix((np.ones(len(col),dtype=np.int8),
                      col.astype(np.int64),
                      ptr.astype(np.int64)),shape=(n,n))
    row_nnz=np.diff(ptr)
    row_abs_sum=np.asarray(np.abs(M).sum(axis=1)).reshape(-1)
    column_abs_sum=np.asarray(np.abs(M).sum(axis=0)).reshape(-1)
    diagonal=M.diagonal()
    rank=int(structural_rank(graph))
    return {
        "size":n,
        "stored_nonzeros":int(len(col)),
        "numeric_nonzeros":int(np.count_nonzero(values)),
        "stored_zero_entries":int(np.count_nonzero(values==0)),
        "structural_rank":rank,
        "structural_rank_defect":int(n-rank),
        "structural_zero_rows":int(np.count_nonzero(row_nnz==0)),
        "numeric_zero_rows":int(np.count_nonzero(row_abs_sum==0)),
        "numeric_zero_cols":int(np.count_nonzero(column_abs_sum==0)),
        "numeric_zero_diagonal_entries":int(np.count_nonzero(diagonal==0)),
        "smallest_nonzero_row_1_norm":float(
            np.min(row_abs_sum[row_abs_sum>0])
        ) if np.any(row_abs_sum>0) else 0.,
        "largest_row_1_norm":float(np.max(row_abs_sum)),
        "Frobenius_norm":float(A.norm(PETSc.NormType.FROBENIUS)),
        "stored_CSR_SHA256":matrix_digest(A),
        "not_a_numerical_rank_proof":True,
    }


def solve_same_matrix(A, solver):
    x_true=A.createVecRight()
    b=A.createVecLeft()
    x=A.createVecRight()
    residual=A.createVecLeft()
    ksp=PETSc.KSP().create(comm=PETSc.COMM_WORLD)
    try:
        lo,hi=x_true.getOwnershipRange()
        idx=np.arange(lo,hi,dtype=PETSc.IntType)
        x_true.setValues(
            idx, 1.+0.05*np.sin(np.arange(lo,hi,dtype=float)),
            addv=PETSc.InsertMode.INSERT_VALUES
        )
        x_true.assemble()
        A.mult(x_true,b)
        rhs_norm=float(b.norm())
        ksp.setType(PETSc.KSP.Type.PREONLY)
        ksp.setOperators(A)
        pc=ksp.getPC()
        pc.setType(PETSc.PC.Type.LU)
        pc.setFactorSolverType(solver)
        status={
            "factor_solver":solver,
            "operator_CSR_before":matrix_digest(A),
            "consistent_manufactured_RHS_norm":rhs_norm,
        }
        try:
            ksp.setUp()
            ksp.solve(b,x)
        except PETSc.Error as exc:
            status["PETSc_error"]=str(exc)
        status.update({
            "KSP_reason":int(ksp.getConvergedReason()),
            "KSP_iterations":int(ksp.getIterationNumber()),
            "PC_failed":str(pc.getFailedReason()),
        })
        if int(ksp.getConvergedReason())>0 and "PETSc_error" not in status:
            A.mult(x,residual)
            residual.axpy(-1.,b)
            diff=x.copy()
            try:
                diff.axpy(-1.,x_true)
                status["relative_linear_residual"]=float(
                    residual.norm()/max(rhs_norm,1e-30)
                )
                status["relative_state_error"]=float(
                    diff.norm()/max(float(x_true.norm()),1e-30)
                )
            finally:
                diff.destroy()
        else:
            status["relative_linear_residual"]=None
            status["relative_state_error"]=None
        status["operator_CSR_after"]=matrix_digest(A)
        if status["operator_CSR_after"]!=status["operator_CSR_before"]:
            raise RuntimeError("Factorizer modified immutable source Mat")
        status["status"]=(
            "PASS_FACTOR_AND_CONSISTENT_MANUFACTURED_SOLVE"
            if int(ksp.getConvergedReason())>0
            and status["relative_linear_residual"] is not None
            and status["relative_linear_residual"]<1e-9
            and status["relative_state_error"]<1e-6
            else "FAIL_FACTOR_OR_CONSISTENT_MANUFACTURED_SOLVE"
        )
        print(
            f"REPLAY {solver}: {status['status']} "
            f"reason={status['KSP_reason']} PC={status['PC_failed']} "
            f"linear_residual={status['relative_linear_residual']} "
            f"state_error={status['relative_state_error']}",
            flush=True,
        )
        return status
    finally:
        ksp.destroy()
        residual.destroy()
        x.destroy()
        b.destroy()
        x_true.destroy()


def main():
    arg=argparse.ArgumentParser(description=__doc__)
    arg.add_argument("--matrix",type=Path,required=True)
    arg.add_argument("--output",type=Path,required=True)
    args=arg.parse_args()
    if not args.matrix.is_file() or not args.output.parent.is_dir():
        raise FileNotFoundError("Original binary matrix or evidence mount missing")
    if args.output.exists():
        raise RuntimeError("Refuse to overwrite original LU replay evidence")
    if PETSc.COMM_WORLD.getSize()!=1:
        raise RuntimeError("The first CN Jacobian LU replay must run in MPI1")
    viewer=PETSc.Viewer().createBinary(str(args.matrix),mode="r",
                                        comm=PETSc.COMM_WORLD)
    A=PETSc.Mat().create(comm=PETSc.COMM_WORLD)
    try:
        A.setType("aij")
        A.load(viewer)
    finally:
        viewer.destroy()
    try:
        report={
            "schema":"scrimp-v130-heat-wave-identical-CN-Jacobian-LU-replay-v1",
            "status":"DIAGNOSTIC_ONLY_NO_PHYSICAL_TS_CLAIM",
            "source_matrix_file":args.matrix.name,
            "matrix":matrix_audit(A),
            "factorizations":{},
            "same_original_matrix_for_both_factorizers":True,
            "no_GetFEM_assembly_in_replay":True,
            "no_diagonal_shift_or_pivot_replacement":True,
            "full_Crank_Nicolson_TS_completed":False,
            "scientific_energy_power_qualified":False,
        }
        digest=report["matrix"]["stored_CSR_SHA256"]
        for solver in ("mumps","superlu_dist"):
            result=solve_same_matrix(A,solver)
            if result["operator_CSR_before"]!=digest or (
                result["operator_CSR_after"]!=digest
            ):
                raise RuntimeError("LU solver got a different source matrix")
            report["factorizations"][solver]=result
        args.output.write_text(json.dumps(report,indent=2)+"\n",
                              encoding="utf-8")
        print(json.dumps(report,indent=2),flush=True)
    finally:
        A.destroy()


if __name__=="__main__":
    main()
