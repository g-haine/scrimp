#!/usr/bin/env python3
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Fail-closed verification of original MPI1/2 + 15s Heat–Wave CN records.

Checks original per-rank real FE time states, GetFEM heat/wave powers,
the fixed accepted-step clock and strict cross-MPI parity independently
of the native probe's own reported scientific gates. Does not label
the long case MPI2-qualified, and never claims modal eigenfunctions,
SuperLU_DIST or MPI-parallel VTK writer acceptance.
"""
import argparse
import json
from pathlib import Path
import sys

import numpy as np


SHORT=("mpi1_short","mpi2_short")
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
        and meta["MPI_ranks"]==(2 if label=="mpi2_short" else 1)
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
    if label==LONG and "z" in a:
        raise RuntimeError("Long acceptance must preserve compact original FE proof")
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


def audit(root,phase):
    result={
        "schema":"scrimp-v130-HeatWave-MPI2-and-long-physical-v1",
        "status":"FAIL_OR_UNQUALIFIED_HEAT_WAVE_MPI2_LONG",
        "source":json.loads((root/"manifest.json").read_text())["source"],
        "phase":phase,"cases":{},"MPI1_vs_MPI2":{},
        "MPI2_rank0_rank1_bitwise_identical":False,
        "full_long_time_power_balance_qualified":False,
        "MPI2_ParaView_qualified":False,
        "SuperLU_DIST_qualified":False,
        "whispering_gallery_eigenmode_qualified":False,
        "errors":[],
    }
    try:
        for label in SHORT:
            ranks=2 if label=="mpi2_short" else 1
            arrays=[]
            for rank in range(ranks):
                arr,meta,summary=load(root,label,rank,4,.001,"gaussian")
                arrays.append(arr)
                result["cases"][f"{label}/rank{rank}"]=summary
            if label=="mpi1_short":
                baseline=arrays[0]
            else:
                if set(arrays[0])!=set(arrays[1]) or not all(
                    np.array_equal(arrays[0][k],arrays[1][k])
                    for k in arrays[0]
                ):
                    raise RuntimeError("MPI2 physical rank trajectories not exactly equal")
                result["MPI2_rank0_rank1_bitwise_identical"]=True
                result["MPI1_vs_MPI2"]=compare(baseline,arrays[0])
        if phase=="short":
            result["status"]="PASS_HEAT_WAVE_MPI1_MPI2_SHORT_PHYSICAL"
        else:
            manifest=json.loads((root/"manifest.json").read_text())
            dt=float(manifest["long_dt"])
            steps=int(manifest["long_steps"])
            profile=manifest["long_profile"]
            if not (0.<dt<=.01 and 4<=steps<=10000
                    and profile in ("gaussian","gallery_packet")):
                raise RuntimeError("Native long-run physical manifest invalid")
            arr,meta,summary=load(root,LONG,0,steps,dt,profile)
            result["cases"][LONG]=summary
            result["full_long_time_power_balance_qualified"]=True
            result["long_duration"]=float(steps*dt)
            result["long_accepted_CN_steps"]=int(steps)
            # Optional cross-validation against ORIGINAL source-b9d56fc
            # long FE Hamiltonian, not a fabricated expected answer.
            if (profile=="gaussian" and steps==3000 and dt==.005):
                result["historical_long_H_comparison"]={
                    "H_initial_abs_difference":
                        float(abs(summary["H0"]-0.7836421211876778)),
                    "H_final_abs_difference":
                        float(abs(summary["Hf"]-0.008268389052222169)),
                }
                if (result["historical_long_H_comparison"]["H_initial_abs_difference"]>1e-8
                        or result["historical_long_H_comparison"]["H_final_abs_difference"]>2e-6):
                    raise RuntimeError(
                        "Full-step-sampled long Hamiltonian disagrees with independently "
                        "archived 15-second original GetFEM trajectory"
                    )
            result["status"]="PASS_HEAT_WAVE_MPI2_SHORT_AND_MPI1_LONG_CN_PHYSICAL"
    except Exception as exc:
        result["errors"].append(f"{type(exc).__name__}: {exc}")
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("evidence",type=Path)
    p.add_argument("--phase",choices=("short","complete"),required=True)
    args=p.parse_args()
    record=audit(args.evidence,args.phase)
    filename="short-verification.json" if args.phase=="short" else "summary.json"
    (args.evidence/filename).write_text(json.dumps(record,indent=2)+"\n")
    print(json.dumps(record,indent=2),flush=True)
    raise SystemExit(0 if record["status"].startswith("PASS_") else 1)


if __name__=="__main__":
    main()
