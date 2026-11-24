# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Full original Heat-Wave CN with an opt-in matched Python PC.

Use the unchanged published original example, strict algebraic DAE
initialization, default GetFEM L2 and distributed COO, and existing
scientific Hamiltonian/dissipation/power audit. Only substitute the
*linear preconditioner* in the first real main TS IJacobian; original
SNES nonlinear residual F(z) and Jacobian A remain untouched.
MPI1 sandbox-only feasibility test, no production default change.
"""

import argparse
import json
import sys
from pathlib import Path

from petsc4py import PETSc

from scrimp.dphs import DPHS
from probe_v130_hw_matching_pc import RowMatchedSuperLUPC
import probe_v130_heat_wave_compatible_cn as original_cn


def main():
    arg=argparse.ArgumentParser(description=__doc__)
    arg.add_argument("--candidate",
                     choices=("clamp","constrained_l2"),required=True)
    arg.add_argument("--out",type=Path,required=True)
    args=arg.parse_args()
    if args.out.exists() and any(args.out.iterdir()):
        raise RuntimeError("Native evidence folder must be empty")

    original_ijacobian=DPHS.IJacobian
    original_argv=sys.argv
    contexts=[]

    def attach_matched_PC(system, TS, t, z, zd, sig, A, P):
        result=original_ijacobian(
            system, TS, t, z, zd, sig, A, P
        )
        if TS is not None:
            ksp=TS.getSNES().getKSP()
            pc=ksp.getPC()
            if pc.getType()!=PETSc.PC.Type.PYTHON:
                ctx=RowMatchedSuperLUPC()
                pc.setType(PETSc.PC.Type.PYTHON)
                pc.setPythonContext(ctx)
                contexts.append(ctx)
                print("MATCHED-PC installed in REAL original PETSc TS, "
                      "outer KSP="+str(ksp.getType()),flush=True)
            elif not contexts:
                raise RuntimeError("Unexpected preexisting PCPYTHON context")
        return result

    DPHS.IJacobian=attach_matched_PC
    try:
        sys.argv=[
            str(Path(__file__).name),
            "--candidate",args.candidate,
            "--backend","superlu_dist",
            "--out",str(args.out),
        ]
        original_cn.main()
    finally:
        DPHS.IJacobian=original_ijacobian
        sys.argv=original_argv

    if len(contexts)!=1:
        raise RuntimeError("Original CN must install exactly one matched PC")
    audit=contexts[0].audit()
    if audit["setup_count"]<1 or audit["apply_count"]<4:
        raise RuntimeError("Actual CN did not use row-matched LU in all steps")
    original_report=json.loads(
        (args.out/"compatible_cn_report.json").read_text()
    )
    if (not original_report["full_original_CN_TS_executed"]
            or original_report["original_numerical_report"]["steps"]!=4
            or original_report["DAE_correction"]["DAE_initializer"]["status"]!="PASS"):
        raise RuntimeError("Original physical CN acceptance not confirmed")
    report={
        "schema":"scrimp-v130-hw-original-CN-rowmatched-python-PC-v1",
        "status":"PASS_ACTUAL_4STEP_CN_NUMERICAL_ROW_MATCHED_PC"
                 if not original_report["physical_balance_gate_passed"]
                 else "PASS_ACTUAL_4STEP_CN_ROW_MATCHED_PC_AND_SCIENCE_GATES",
        "candidate":args.candidate,
        "original_example_modified":False,
        "true_original_Ts_4steps":True,
        "strict_original_DAE_initialization":True,
        "native_rowmatched_PC":audit,
        "actual_original_physics_thresholds":
            original_report["scientific_thresholds"],
        "scientific_energy_power_gate_passed":
            original_report["physical_balance_gate_passed"],
        "MPI2_qualified":False,
        "automatically_enabled_in_production":False,
    }
    f=args.out/"rowmatched_cn_audit.json"
    if f.exists():
        raise RuntimeError("Refuse to overwrite native PC report")
    f.write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(report,indent=2),flush=True)


if __name__=="__main__":
    main()
