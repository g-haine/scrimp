# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Genuine corrected Heat-Wave PETSc TS/CN with oriented 2D interface.

No published source changes. Original S.Domain/Control_Port/Brick are
patched only for the lifetime of this sandbox process to direct
thermal forms to heat outward faces (30) and wave forms to wave
outward faces (31). Original strict DAE boundary-compatible
initialization and original main TS solver/ParaView remain active.
Corrected PHYSICAL port-power measurement is explicitly separated
from old doubled region 10. Native failure is a diagnostic,
not an automatic physical PASS.
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from mpi4py import MPI
from petsc4py import PETSc

import scrimp as S
from v130_hw_interface_geometry import new_one_sided_regions
import probe_v130_main_2d_examples as main2d
import probe_v130_heat_wave_compatible_cn as original_cn
from measure_v130_hw_dual_facet_physics import (
    collect_heat_wave_dual_facet_physics,
)


HEAT_BRICKS=frozenset(("M_Y_T","B_T","B_T^T"))
WAVE_BRICKS=frozenset(("M_Y_w","B_w","B_w^T"))
PORTS={"Interface Heat":30,"Interface Wave":31}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--candidate",choices=("clamp","constrained_l2"),required=True
    )
    parser.add_argument("--out",type=Path,required=True)
    args=parser.parse_args()
    if MPI.COMM_WORLD.size!=1 or PETSc.COMM_WORLD.getSize()!=1:
        raise RuntimeError("Corrected physical Heat-Wave CN is MPI1 only")
    if args.out.exists() and any(args.out.iterdir()):
        raise RuntimeError("Refuse to reuse original CN native report directory")

    original_Domain=S.Domain
    original_Control_Port=S.Control_Port
    original_Brick=S.Brick
    original_collector=main2d.collect_heat_wave_physics
    original_argv=sys.argv
    calls={"domain":0,"heat_bricks":{},"wave_bricks":{},"ports":{}}
    captured={}

    def new_domain(*pos,**kw):
        domain=original_Domain(*pos,**kw)
        if pos[:1]!=("Concentric",):
            raise RuntimeError("Expected original Concentric mesh")
        calls["domain"]+=1
        if calls["domain"]!=1:
            raise RuntimeError("Unexpected duplicate Concentric domain")
        mesh=domain.get_mesh()[0]
        geometry=new_one_sided_regions(mesh)
        disk=np.asarray(mesh.region(30),dtype=np.int64)
        wave=np.asarray(mesh.region(31),dtype=np.int64)
        doubled=np.asarray(mesh.region(10),dtype=np.int64)
        cv_faces=lambda mat:set(map(tuple,mat.T.tolist()))
        if (disk.shape!=(2,78) or wave.shape!=(2,78)
                or doubled.shape!=(2,156)
                or cv_faces(doubled)!=(cv_faces(disk)|cv_faces(wave))
                or cv_faces(disk)&cv_faces(wave)):
            raise RuntimeError("Incorrect original 78+78 Concentric interface")
        captured["native_original_facet_topology"]=geometry
        return domain

    def control(*pos,**kw):
        if pos and pos[0] in PORTS:
            if kw.get("region")!=10:
                raise RuntimeError("Original control region changed")
            name=pos[0]
            kw["region"]=PORTS[name]
            calls["ports"][name]=kw["region"]
        return original_Control_Port(*pos,**kw)

    def brick(*pos,**kw):
        if pos and pos[0] in HEAT_BRICKS|WAVE_BRICKS:
            if len(pos)<3 or list(pos[2])!=[10]:
                raise RuntimeError("Original physical brick region changed")
            name=pos[0]
            region=30 if name in HEAT_BRICKS else 31
            pos=list(pos)
            pos[2]=[region]
            pos=tuple(pos)
            bucket="heat_bricks" if name in HEAT_BRICKS else "wave_bricks"
            calls[bucket][name]=region
        return original_Brick(*pos,**kw)

    def real_collector(system,times,H):
        if captured.get("physical_collector_called",False):
            raise RuntimeError("Corrected CN must collect physical states exactly once")
        captured["physical_collector_called"]=True
        return collect_heat_wave_dual_facet_physics(system,times,H)

    S.Domain=new_domain
    S.Control_Port=control
    S.Brick=brick
    main2d.collect_heat_wave_physics=real_collector
    try:
        sys.argv=[
            str(Path(__file__).name),
            "--candidate",args.candidate,
            "--backend","mumps","--out",str(args.out),
        ]
        original_cn.main()
    finally:
        S.Domain=original_Domain
        S.Control_Port=original_Control_Port
        S.Brick=original_Brick
        main2d.collect_heat_wave_physics=original_collector
        sys.argv=original_argv

    if (calls["domain"]!=1 or set(calls["ports"])!=set(PORTS)
            or set(calls["heat_bricks"])!=HEAT_BRICKS
            or set(calls["wave_bricks"])!=WAVE_BRICKS
            or not captured.get("physical_collector_called")):
        raise RuntimeError("Actual original Heat-Wave dual-facet CN not verified")
    physical=json.loads(
        (args.out/"heat_wave_physics_diagnostic.json").read_text()
    )
    report=json.loads(
        (args.out/"compatible_cn_report.json").read_text()
    )
    original=json.loads(
        (args.out/"original_2d_example_report.json").read_text()
    )
    if (not report["true_strict_DAE_pre_Ts"]
            or report["DAE_correction"]["DAE_initializer"]["status"]!="PASS"
            or not report["DAE_correction"]["DAE_initializer"][
                "differential_initial_exact"
            ]
            or original["steps"]!=4
            or original["SNES_failures"]!=0
            or original["rejected_steps"]!=0
            or original["PETSc_reason"]<=0
            or original["full_state_DOFs"]!=7689
            or physical["interface_heat_region"]!=30
            or physical["interface_wave_region"]!=31):
        raise RuntimeError("Original strict DAE / corrected CN native contract failed")
    nontrivial=physical["interface_exchange_nontrivial"]
    balance=report["physical_balance_gate_passed"]
    result={
        "schema":"scrimp-v130-HW-dual-facet-real-CN-v1",
        "status":(
            "DIAGNOSTIC_CN_CONVERGED_BUT_PHYSICAL_BALANCE_OPEN"
            if not balance else
            "PASS_CN_BALANCE_NONVACUOUS_INTERFACE_OPEN" if not nontrivial else
            "PASS_CORRECTED_FACET_CN_BALANCE_AND_NONTRIVIAL_PORT_EXCHANGE"
        ),
        "candidate":args.candidate,
        "source":"unchanged published original heat_wave_eq",
        "example_src_modified":False,
        "runtime_only_constructor_region_override":calls,
        "interface_heat_region":30,
        "interface_wave_region":31,
        "strict_DAE_initialization_PASSED":True,
        "actual_real_4_step_PETSc_CN_PASSED":True,
        "CN_discrete_midpoint_energy_gate_PASSED":balance,
        "CN_interface_exchange_nontrivial":nontrivial,
        "actual_energy_power_diagnostics":physical,
        "original_full_state_included_in_npz":True,
        "MPI2_or_parallel_export_qualified":False,
        "production_solver_policy_changed":False,
    }
    output=args.out/"dual_facet_cn_audit.json"
    if output.exists():
        raise RuntimeError("Refuse to overwrite corrected CN native evidence")
    output.write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps({
        "status":result["status"],
        "candidate":args.candidate,
        "steps":original["steps"],
        "DAE_residual":report["DAE_correction"][
            "DAE_initializer"]["native_final_residual_l2"],
        "max_CN_total_energy_balance_defect":
            physical["max_CN_step_H_balance_defect"],
        "max_CN_thermal_energy_balance_defect":
            physical["max_CN_step_heat_balance_defect"],
        "max_CN_wave_energy_balance_defect":
            physical["max_CN_step_wave_balance_defect"],
        "max_CN_interface_power_defect":
            physical["max_CN_interface_power_defect"],
        "max_individual_CN_interface_power":
            physical["max_individual_CN_interface_power"],
        "max_heat_dissipation":
            physical["max_CN_heat_dissipation"],
        "interface_exchange_ratio_to_heat_dissipation":
            physical["interface_exchange_ratio_to_heat_dissipation"],
    },indent=2),flush=True)


if __name__=="__main__":
    main()
