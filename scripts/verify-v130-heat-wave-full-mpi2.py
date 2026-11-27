#!/usr/bin/env python3
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Independent full-horizon MPI1/MPI2 Heat-Wave CN scientific audit.

Checks ALL 3,001 accepted 7,689-DOF states, separate heat/wave
GetFEM powers, full 15-second fixed CN clock, independent
energy identities and strict complete MPI1/MPI2 parity. It
never equates this to parallel VTK writing or a gallery mode.
"""
import argparse
import json
from pathlib import Path
import sys

import numpy as np


SHORT=("mpi1_long","mpi2_long")
LONG="mpi1_long"
NEEDED=(
    "H_heat","H_wave","H_components","CN_D_heat_flux",
    "CN_P_interface_physical","CN_P_outer",
    "CN_per_step_heat_balance_error","CN_per_step_wave_balance_error",
    "CN_per_step_H_balance_error","CN_cumulative_total_balance_error",
)


def load(root,label,rank,steps,dt,profile):
    directory=root/label
    rc=int((root/f"{label}.exitcode").read_text())
    if rc:
        raise RuntimeError(f"{label} exited natively with code={rc}")
    meta=json.loads((directory/f"rank{rank}.json").read_text())
    if not (
        meta["schema"]=="scrimp-v130-heat-wave-MPI-and-long-CN-physical-v1"
        and meta["status"]=="PASS_NATIVE_STRICT_DAE_CN_ORIENTED_POWER"
        and meta["rank"]==rank
        and meta["MPI_ranks"]==(2 if label=="mpi2_long" else 1)
        and meta["steps"]==steps and abs(meta["dt"]-dt)<1e-15
        and meta["profile"]==profile
        and meta["DOFs"]==7689
        and meta["saved_all_accepted_steps"]
        and all(meta["physical_gates"].values())
        and not meta["parallel_ParaView_export_qualified"]
        and not meta["whispering_gallery_mode_qualified"]
        and meta["all_rank_native_assembly_and_quadratures"]
        and meta["initial_wave_boundary"]["method"]=="clamp"
        and meta["initial_wave_boundary"]["boundary_compatibility_enforced"]
        and meta["initial_wave_boundary"]["after_boundary_l2_coefficients"]==0.
    ):
        raise RuntimeError(f"{label}/rank{rank}: native physical metadata FAIL")
    with np.load(directory/f"rank{rank}.npz",allow_pickle=False) as npz:
        a={k:npz[k].copy() for k in npz.files}
    if not set(NEEDED).issubset(
        {k[5:] for k in a if k.startswith("phys_")}
    ):
        raise RuntimeError(f"{label}/rank{rank}: missing real GetFEM quadratures")
    t,H=a["t"],a["H"]
    if (t.shape!=(steps+1,) or H.shape!=(steps+1,)
        or a["z0"].shape!=(7689,) or a["zf"].shape!=(7689,)
        or not np.all(np.isfinite(H)) or not np.all(np.isfinite(t))
        or np.min(H)<=0
        or not np.allclose(t,np.arange(steps+1)*dt,
                           atol=1e-10,rtol=0)):
        raise RuntimeError(f"{label}/rank{rank}: FE time clock/H invalid")
    if label in SHORT and a.get("z",None) is None:
        raise RuntimeError("Short MPI comparison requires complete FE vectors")
    if label in SHORT and a["z"].shape!=(steps+1,7689):
        raise RuntimeError("Incorrect short genuine full FE history shape")
    if not all(np.all(np.isfinite(a[k])) for k in a):
        raise RuntimeError("Nonfinite original native FE or physical data")
    for key in NEEDED:
        if a["phys_"+key].shape!=t.shape:
            raise RuntimeError(f"Bad time-aligned original GetFEM trace {key}")

    def p(key):
        return a["phys_"+key]
    scale=max(1.,float(np.max(np.abs(H))))
    dt_real=np.diff(t)
    dissip=p("CN_D_heat_flux")[1:]
    interf=p("CN_P_interface_physical")[1:]
    outer=p("CN_P_outer")[1:]
    total_independent=np.diff(H)+dt_real*(dissip-interf-outer)
    # Guard independently against a bogus native "PASS": all balances
    # are recomputed on host from unmodified original FE quadratures.
    gates={
        "all_consecutive_true_CN_states_saved":True,
        "H_equals_its_independently_assembled_parts":
            bool(np.max(abs(H-p("H_components")))<1e-8*scale),
        "nonnegative_heat_and_wave":
            bool(min(np.min(p("H_heat")),np.min(p("H_wave")))>-1e-9*scale),
        "nontrivial_positive_dissipation":
            bool(np.min(dissip)>=-1e-10*scale
                 and np.max(dissip)>1e-10*scale),
        "oriented_interface_power_cancels":
            bool(np.max(abs(interf))<1e-8*scale),
        "outer_homogeneous_power_zero":
            bool(np.max(abs(outer))<1e-8*scale),
        "heat_step_CN_energy_power":
            bool(np.max(abs(p("CN_per_step_heat_balance_error")[1:]))
                 <1e-8*scale),
        "wave_step_CN_energy_power":
            bool(np.max(abs(p("CN_per_step_wave_balance_error")[1:]))
                 <1e-8*scale),
        "total_step_independent_CN_energy_power":
            bool(np.max(abs(total_independent))<1e-8*scale
                 and np.max(abs(p("CN_per_step_H_balance_error")[1:]))
                     <1e-8*scale),
        "full_cumulative_CN_energy_power":
            bool(np.max(abs(p("CN_cumulative_total_balance_error")))
                 <1e-8*scale),
        "total_energy_nonincreasing_each_step":
            bool(np.max(np.diff(H))<1e-8*scale),
    }
    if not all(gates.values()):
        raise RuntimeError(f"{label}/rank{rank} physical gates FAIL: "+
                           ",".join(k for k,v in gates.items() if not v))
    return a,meta,{
        "gates":gates,
        "H0":float(H[0]),"Hf":float(H[-1]),
        "max_independent_total_balance":float(np.max(abs(total_independent))),
        "max_thermal_balance":float(np.max(abs(p("CN_per_step_heat_balance_error")))),
        "max_wave_balance":float(np.max(abs(p("CN_per_step_wave_balance_error")))),
        "max_cumulative_balance":float(np.max(abs(p("CN_cumulative_total_balance_error")))),
        "max_interface_power":float(np.max(abs(interf))),
    }


def compare(reference,parallel):
    if reference.keys()!=parallel.keys():
        raise RuntimeError("MPI FE trace field list differs")
    comparisons={}
    for name in reference:
        x,y=reference[name],parallel[name]
        if x.shape!=y.shape or not np.all(np.isfinite(x)) or not np.all(np.isfinite(y)):
            raise RuntimeError("MPI FE trace shapes or finite values differ: "+name)
        defect=float(np.max(abs(x-y)))
        if not np.allclose(x,y,atol=5e-10,rtol=2e-8):
            raise RuntimeError(f"MPI1/MPI2 physical FE trace differs: {name} {defect}")
        comparisons[name]=defect
    return comparisons


def audit(root):
    result={
        "schema":"scrimp-v130-heat-wave-full-15s-MPI1-MPI2-physical-v1",
        "status":"FAIL_OR_UNQUALIFIED_FULL_HEAT_WAVE_MPI2",
        "errors":[],
        "MPI2_rank0_rank1_bitwise_identical":False,
        "MPI1_vs_MPI2_all_3001_FE_states_and_physical_traces":{},
        "complete_15s_physical_power_balance_qualified":False,
        "MPI2_full_horizon_qualified":False,
        "MPI_parallel_VTK_writers_qualified":False,
        "SuperLU_DIST_qualified":False,
        "whispering_gallery_mode_qualified":False,
        "cases":{},
    }
    try:
        manifest=json.loads((root/"manifest.json").read_text())
        dt=float(manifest["dt"])
        steps=int(manifest["steps"])
        profile=manifest["profile"]
        if not (dt==.005 and steps==3000 and profile=="gaussian"
                and manifest["source"]
                and manifest["original_model_unchanged"] is True
                and manifest["full_FE_states_saved_per_rank"] is True):
            raise RuntimeError("This is not the previously qualified Gaussian 15s problem")
        trajectories={}
        for label,size in (("mpi2_long",2),("mpi1_long",1)):
            for rank in range(size):
                arr,meta,summary=load(root,label,rank,steps,dt,profile)
                if arr["z"].shape!=(steps+1,7689):
                    raise RuntimeError(f"Full original FE history missing: {label}/{rank}")
                pre=json.loads((root/label/f"preflight-rank{rank}.json").read_text())
                if (pre["accepted_PETSc_TS_steps"]!=steps
                        or pre["actual_saved_FE_states"]!=steps+1
                        or pre["expected_saved_FE_states"]!=steps+1
                        or pre["FE_state_shape"]!=[steps+1,7689]
                        or pre["missing_accepted_grid_indices"]
                        or pre["TS_SNES_failures"]!=0
                        or pre["TS_step_rejections"]!=0):
                    raise RuntimeError(f"Original CN monitor preflight failed {label}/{rank}")
                result["cases"][f"{label}/rank{rank}"]=summary
                trajectories[label,rank]=arr
        a=trajectories["mpi2_long",0]
        b=trajectories["mpi2_long",1]
        if set(a)!=set(b) or not all(np.array_equal(a[k],b[k]) for k in a):
            raise RuntimeError("MPI2 original 15-second GetFEM states/powers differ across ranks")
        result["MPI2_rank0_rank1_bitwise_identical"]=True
        reference=trajectories["mpi1_long",0]
        comparison=compare(reference,a)
        if "z" not in comparison or "H" not in comparison:
            raise RuntimeError("Missing full original FE coefficient/power parity")
        result["MPI1_vs_MPI2_all_3001_FE_states_and_physical_traces"]=comparison
        result["max_FEM_coefficient_difference"]=comparison["z"]
        result["long_accepted_CN_steps"]=steps
        result["original_FE_DOFs"]=7689
        result["native_saved_states_per_rank"]=steps+1
        result["full_duration"]=float(steps*dt)
        result["H0_MPI1"]=result["cases"]["mpi1_long/rank0"]["H0"]
        result["Hfinal_MPI1"]=result["cases"]["mpi1_long/rank0"]["Hf"]
        result["H0_MPI2"]=result["cases"]["mpi2_long/rank0"]["H0"]
        result["Hfinal_MPI2"]=result["cases"]["mpi2_long/rank0"]["Hf"]
        # Independently archived original 15-second Gaussian
        # GetFEM/PETSc H values; not a new manufactured trajectory.
        if (abs(result["H0_MPI1"]-0.7836421211876778)>1e-8
                or abs(result["Hfinal_MPI1"]-0.008268389052222169)>2e-6):
            raise RuntimeError("Full MPI1 baseline disagrees with accepted native reference")
        result["complete_15s_physical_power_balance_qualified"]=True
        result["MPI2_full_horizon_qualified"]=True
        result["status"]="PASS_HEAT_WAVE_FULL_15S_MPI1_MPI2_CN_PHYSICAL_PARITY"
    except Exception as exc:
        result["errors"].append(f"{type(exc).__name__}: {exc}")
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("evidence",type=Path)
    args=p.parse_args()
    summary=audit(args.evidence)
    (args.evidence/"summary.json").write_text(json.dumps(summary,indent=2)+"\n")
    print(json.dumps(summary,indent=2),flush=True)
    raise SystemExit(0 if summary["status"].startswith("PASS_") else 1)


if __name__=="__main__":
    main()
