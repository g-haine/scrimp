# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Original CN matrix: actual KSP PREONLY + Python row-matched PC."""

import argparse
import json
from pathlib import Path

import numpy as np
from petsc4py import PETSc

from probe_v130_hw_matching_pc import RowMatchedSuperLUPC
from replay_base import matrix_digest


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--matrix",type=Path,required=True)
    ap.add_argument("--out",type=Path,required=True)
    args=ap.parse_args()
    if PETSc.COMM_WORLD.getSize()!=1 or not args.matrix.is_file():
        raise RuntimeError("Original first CN binary MPI1 required")
    if args.out.exists() or not args.out.parent.is_dir():
        raise RuntimeError("Refuse to overwrite existing evidence")
    viewer=PETSc.Viewer().createBinary(
        str(args.matrix),mode="r",comm=PETSc.COMM_WORLD
    )
    A=PETSc.Mat().create(comm=PETSc.COMM_WORLD)
    try:
        A.setType("aij")
        A.load(viewer)
    finally:
        viewer.destroy()
    xtrue=rhs=x=check=ksp=None
    ctx=RowMatchedSuperLUPC()
    try:
        old=matrix_digest(A)
        assert old=="c94c69c4e4307d7a5aa61d9eba186782c0495b49431cadcb7281d7a7b5bffaff"
        xtrue=A.createVecRight()
        rhs=A.createVecLeft()
        x=A.createVecRight()
        check=A.createVecLeft()
        n=A.getSize()[0]
        xtrue.setValues(np.arange(n,dtype=PETSc.IntType),
                        1.+0.05*np.sin(np.arange(n,dtype=float)),
                        addv=PETSc.InsertMode.INSERT_VALUES)
        xtrue.assemble()
        A.mult(xtrue,rhs)
        ksp=PETSc.KSP().create(comm=PETSc.COMM_WORLD)
        ksp.setOperators(A)
        ksp.setType(PETSc.KSP.Type.PREONLY)
        outer_pc=ksp.getPC()
        outer_pc.setType(PETSc.PC.Type.PYTHON)
        outer_pc.setPythonContext(ctx)
        ksp.setUp()
        ksp.solve(rhs,x)
        A.mult(x,check)
        check.axpy(-1.,rhs)
        rel=float(check.norm()/max(float(rhs.norm()),1e-30))
        delta=x.copy()
        try:
            delta.axpy(-1.,xtrue)
            state=float(delta.norm()/max(float(xtrue.norm()),1e-30))
        finally:
            delta.destroy()
        success=(
            int(ksp.getConvergedReason())>0 and
            rel<1e-9 and state<1e-6 and
            ctx.apply_count>=1 and
            ctx.setup_history and
            old==matrix_digest(A)
        )
        report={
            "schema":"scrimp-v130-hw-original-CN-matched-python-PC-v1",
            "status":"PASS_SNAPSHOT_PYTHON_PC_ORIGINAL_EQUATIONS_ONLY"
                     if success else "FAIL_SNAPSHOT_ROW_MATCHING_PC",
            "KSP_convergence_reason":int(ksp.getConvergedReason()),
            "PC_type":outer_pc.getType(),
            "residual_original_equations":rel,
            "known_state_relative_error":state,
            "original_CSR_SHA256_before":old,
            "original_CSR_SHA256_after":matrix_digest(A),
            "matched_PC":ctx.audit(),
            "actual_TS_executed":False,
            "physical_power_balance_qualified":False,
        }
        args.out.write_text(json.dumps(report,indent=2)+"\n",
                            encoding="utf-8")
        print(json.dumps(report,indent=2),flush=True)
        if not success:
            raise RuntimeError("Original KSP with row-matched PCPYTHON failed")
    finally:
        if ksp is not None:ksp.destroy()
        for obj in (check,x,rhs,xtrue,A):
            if obj is not None:obj.destroy()


if __name__=="__main__":
    main()
