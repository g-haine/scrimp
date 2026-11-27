#!/usr/bin/env python3
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp


"""Independent host acceptance of the genuine inviscid shallow-water 2D pilot.

Native physical data must originate from the SAME unmodified original
examples/shallow_water.py, with opt-in nonlinear_affine initialization.
Two BDF2 time steps are measured at a common nonzero final time, on
MPI1 and MPI2. We test full FE coefficient parity, mass conservation,
original port work and independently recompute the Hamiltonian defect.

BDF2 + lagged explicit gyroscopic G does NOT conserve H exactly.
No numerical temporal order is asserted from a two-dt experiment,
particularly with the discontinuous original initial surface.
"""
import argparse
import json
import math
from pathlib import Path
import re

import numpy as np

IMAGE="1d6d0c90115213298047b28237a7f38087563a35"
CASES=((20,0.0001),(40,0.00005))
FINAL_T=0.002


def demand(condition, message):
    if not condition:
        raise RuntimeError(message)


def load_case(root, steps, dt, mpi, rank):
    name=f"steps{steps}_mpi{mpi}"
    directory=root/name
    demand(int((root/f"{name}.exitcode").read_text())==0,
           f"{name}: native process failed")
    report=json.loads((directory/"shallow_short_physical.json").read_text())
    demand(
        report.get("schema")=="scrimp-v130-issue28-original-shallow-short-physical-v1"
        and report.get("status")=="PASS_SHORT_PHYSICAL_TIME_GATE"
        and report.get("model")=="original shipped examples/shallow_water.py"
        and report.get("mpi_ranks")==mpi
        and report.get("mesh_h")==0.1
        and report.get("accepted_steps")==steps
        and report.get("rejected_steps")==0
        and report.get("nonlinear_solver_failures")==0
        and report.get("ts_converged_reason",0)>0
        and abs(report.get("t_final_actual",-1)-FINAL_T)<1e-10
        and abs(report.get("time_step_requested",-1)-dt)<1e-15
        and report.get("original_prescribed_h0_p0_bitwise") is True
        and report.get("original_boundary_controls_unchanged_zero") is True
        and report.get("original_explicit_gyro_unchanged") is True
        and report.get("Hamiltonian_GetFEM_and_registered_agree") is True
        and report.get("BDF2_exact_power_conservation_claimed") is False
        and report.get("original_rotating_viscous_Y0_qualified") is False,
        f"{name}: wrong original physical model or native TS metadata"
    )
    init=report.get("production_initializer_report",{})
    demand(init.get("strategy")=="nonlinear_affine"
           and init.get("status")=="PASS"
           and init.get("differential_initial_exact") is True
           and init.get("accuracy_target_met") is True
           and float(init.get("native_final_algebraic_residual_l2",float("inf")))
           <= float(init.get("accuracy_target",0)),
           f"{name}: strict GetFEM nonlinear-affine DAE initialization invalid")
    with np.load(directory/f"shallow_fe_rank{rank}.npz",allow_pickle=False) as z:
        data={key:z[key].copy() for key in z.files}
    demand(set(data)=={"t","z","H","mass","boundary_power"},
           f"{name} rank {rank}: missing actual native FE arrays")
    t,z,H,M,P=(data[key] for key in ("t","z","H","mass","boundary_power"))
    demand(
        t.ndim==1 and len(t)>=3
        and z.ndim==2 and z.shape[0]==len(t) and z.shape[1]>200
        and all(field.shape==t.shape for field in (H,M,P))
        and all(np.all(np.isfinite(arr)) for arr in data.values())
        and np.all(np.diff(t)>0)
        and abs(t[0])<1e-12 and abs(t[-1]-FINAL_T)<1e-10
        and np.all(H>0) and np.all(M>0),
        f"{name} rank {rank}: invalid genuine FE states, clock or physical integrals"
    )
    for key,published in (
        ("saved_times","t"),("hamiltonian","H"),("mass","mass"),
        ("total_boundary_power","boundary_power")
    ):
        array=np.asarray(report[key],dtype=np.float64)
        demand(array.shape==data[published].shape
               and np.allclose(array,data[published],atol=2e-9,rtol=2e-10),
               f"{name} rank {rank}: original FE trace does not match GetFEM report {key}")
    work=np.zeros_like(H)
    work[1:]=np.cumsum((P[:-1]+P[1:])*0.5*np.diff(t))
    defect=H-H[0]-work
    reported=np.asarray(report["energy_balance_defect"],dtype=np.float64)
    demand(
        reported.shape==defect.shape
        and np.allclose(defect,reported,atol=2e-8,rtol=2e-10)
        and abs(float(defect[-1])-float(report["signed_final_energy_balance_defect"]))<2e-8,
        f"{name}: BDF2 H(T)-H(0)-boundary work not reproduced"
    )
    relative_mass=float(np.max(abs(M-M[0]))/max(abs(M[0]),1e-15))
    relative_H=float(np.max(abs(H-H[0]))/max(abs(H[0]),1e-15))
    max_P=float(np.max(abs(P)))
    demand(relative_mass<=1e-5
           and relative_H<=5e-3
           and max_P<=1e-7
           and not report["physical_sanity_gate_failures"],
           f"{name}: original independent mass, energy, boundary power sanity gates FAIL")
    return data,report,{
        "MPI":mpi,"steps":steps,"dt":dt,"saved_states":len(t),
        "DOFs":int(z.shape[1]),"mass_relative_drift":relative_mass,
        "H_relative_drift":relative_H,"maximum_boundary_power":max_P,
        "H_initial":float(H[0]),"H_final":float(H[-1]),
        "signed_energy_balance_defect":float(defect[-1]),
        "max_absolute_energy_defect":float(np.max(abs(defect))),
        "exact_BDF2_H_conservation_claimed":False,
    }


def compare_states(a,b,*,tag,atol,rtol):
    demand(set(a)==set(b),f"{tag}: different original field lists")
    maxima={}
    for key in a:
        x,y=a[key],b[key]
        demand(x.shape==y.shape
               and np.all(np.isfinite(x)) and np.all(np.isfinite(y)),
               f"{tag}: differing state/trace shapes {key}")
        worst=float(np.max(abs(x-y)))
        if not np.allclose(x,y,atol=atol,rtol=rtol):
            raise RuntimeError(f"{tag}: actual FE/energy MPI parity FAIL for {key}: {worst}")
        maxima[key]=worst
    return maxima


def audit(root):
    result={
        "schema":"v130-original-inviscid-SWE-bounded-physical-MPI-v1",
        "status":"FAIL_OR_UNQUALIFIED_SWE_EXTENDED_BDF2",
        "cases":{},"MPI_parity":{},"errors":[],
        "BDF2_exact_power_conservation_claimed":False,
        "observed_temporal_order_qualified":False,
        "rotating_viscous_SWE_qualified":False,
        "original_default_0p5s_horizon_qualified":False,
        "MPI_parallel_ParaView_qualified":False,
        "mesh_convergence_qualified":False,
        "published_examples_modified":False,
    }
    try:
        manifest=json.loads((root/"manifest.json").read_text())
        demand(manifest.get("schema")=="v130-shallow-original-extended-physical-v1"
               and re.fullmatch(r"[a-f0-9]{40}",manifest.get("source",""))
               and manifest.get("image_revision")==IMAGE
               and manifest.get("model")=="examples/shallow_water.py"
               and manifest.get("fixed_final_time")==FINAL_T
               and manifest.get("steps")==[20,40]
               and manifest.get("MPI")==[1,2]
               and manifest.get("source_published_immutable") is True,
               "Native source/image/protocol manifest invalid")
        result["source"]=manifest["source"]
        results={}
        for steps,dt in CASES:
            v1,r1,stats1=load_case(root,steps,dt,1,0)
            v20,r20,stats2=load_case(root,steps,dt,2,0)
            v21,r21,_=load_case(root,steps,dt,2,1)
            # Both real MPI ranks must evaluate the original algebraic
            # constraints and preserve the same full finite element state.
            result["MPI_parity"][str(steps)]={
                "rank0_rank1_max_abs":compare_states(
                    v20,v21,tag=f"MPI2 ranks, {steps} steps",
                    atol=2e-9,rtol=2e-8),
                "MPI1_MPI2_max_abs":compare_states(
                    v1,v20,tag=f"MPI1 vs MPI2, {steps} steps",
                    atol=5e-7,rtol=2e-7),
            }
            result["cases"][f"{steps}_MPI1"]=stats1
            result["cases"][f"{steps}_MPI2"]=stats2
            results[steps]=stats1
            demand(stats1["DOFs"]==stats2["DOFs"],
                   "MPI partition altered physical DOF count")
        coarse=results[20]["signed_energy_balance_defect"]
        fine=results[40]["signed_energy_balance_defect"]
        ratio=(abs(coarse)/abs(fine) if abs(fine)>1e-15 else None)
        result["comparison"]={
            "same_final_time":FINAL_T,
            "coarse_signed_energy_defect":coarse,
            "fine_signed_energy_defect":fine,
            "abs_defect_ratio_coarse_over_fine":ratio,
            "energy_defect_decreases":bool(abs(fine)<abs(coarse)),
            "note":"Measured behavior only; 2 dt values and discontinuous "
                   "initial height do not establish convergence order",
        }
        result["status"]="PASS_SWE_BOUNDED_NATIVE_BDF2_MPI1_MPI2"
    except Exception as exc:
        result["errors"].append(f"{type(exc).__name__}: {exc}")
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("evidence",type=Path)
    args=parser.parse_args()
    root=args.evidence.resolve()
    summary=audit(root)
    (root/"summary.json").write_text(json.dumps(summary,indent=2)+"\n")
    print(json.dumps(summary,indent=2),flush=True)
    raise SystemExit(0 if summary["status"].startswith("PASS_") else 1)


if __name__=="__main__":
    main()
