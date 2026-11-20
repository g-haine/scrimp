# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Strict MPI parity and accepted-state tests for REAL nonlinear SWE TS.

This is a short physical smoke; it does NOT assert exact pH energy
conservation or second order for the first-order lagged gyroscopic G.
The independent nonlinear mass-Jacobian FD gate has already PASSED.
"""

import argparse
import json
from pathlib import Path

import numpy as np


def compare(folder):
    reports={}
    for case in ("published","manufactured"):
        for ranks in (1,2):
            p=folder/f"shallow_explicit_TS_{case}_{ranks}r.json"
            r=json.loads(p.read_text(encoding="utf-8"))
            if (
                bool(r["manufactured_nonzero_G_initial_data"])!=
                (case=="manufactured")
                or r["TS_converged_reason"]<=0
                or r["TS_accepted_steps"]<3
                or r["TS_rejected_steps"]!=0
                or r["TS_SNES_failures"]!=0
                or r["max_frozen_explicit_RHS_drift"]!=0
                or r["main_IFunction_checks"]<r["TS_accepted_steps"]
                or len(r["main_explicit_prestep"])<
                    r["TS_accepted_steps"]
                or r["exact_H_balance_claimed"]
                or r["nonlinear_DAE_initialization_qualified"]
            ):
                raise RuntimeError(f"Invalid physical TS record {p.name}")
            expected=r["TS_accepted_steps"]*r["dt"]
            if abs(r["time_end"]-expected)>1e-10:
                raise RuntimeError(f"Wrong physical TS end time {p.name}")
            if case=="manufactured" and (
                r["max_nonzero_explicit_G_norm"]<=1e-12
            ):
                raise RuntimeError(
                    "Physical nonlinear G was vacuous: "+p.name
                )
            h=np.asarray(r["hamiltonian"],dtype=float)
            t=np.asarray(r["times"],dtype=float)
            z=np.asarray(r["full_state_trajectory"],dtype=float)
            if (
                len(t)<2 or not np.all(np.diff(t)>0)
                or z.shape[0]!=t.size or h.shape!=(len(t),)
                or not np.all(np.isfinite(h))
                or not np.all(np.isfinite(z))
                or abs(t[-1]-expected)>1e-10
            ):
                raise RuntimeError(f"Incomplete physical TS history {p.name}")
            times=np.array(
                [x["time"] for x in r["main_explicit_prestep"]]
            )
            if not np.all(np.diff(times)>=0):
                raise RuntimeError(
                    f"Explicit PreStep reversed accepted time {p.name}"
                )
            if not np.allclose(
                times,
                [x["frozen_rhs_time"] for x in
                 r["main_explicit_prestep"]],
                atol=1e-13,rtol=0.
            ):
                raise RuntimeError(
                    f"Frozen RHS time differs from accepted time {p.name}"
                )
            reports[case,ranks]=r
    max_scaled=0.
    max_h_rel=0.
    for case in ("published","manufactured"):
        a=reports[case,1]
        b=reports[case,2]
        t1=np.asarray(a["times"])
        t2=np.asarray(b["times"])
        z1=np.asarray(a["full_state_trajectory"])
        z2=np.asarray(b["full_state_trajectory"])
        h1=np.asarray(a["hamiltonian"])
        h2=np.asarray(b["hamiltonian"])
        if (t1.shape!=t2.shape or not np.allclose(
            t1,t2,atol=1e-13,rtol=0.
        ) or z1.shape!=z2.shape):
            raise RuntimeError(
                f"Shallow-water TS MPI samples incompatible: {case}"
            )
        scaled=max(
            np.linalg.norm(z1[i]-z2[i])/
            max(1.,np.linalg.norm(z1[i]),np.linalg.norm(z2[i]))
            for i in range(len(t1))
        )
        hrel=max(
            abs(u-v)/max(1.,abs(u),abs(v))
            for u,v in zip(h1,h2)
        )
        max_scaled=max(max_scaled,scaled)
        max_h_rel=max(max_h_rel,hrel)
        if scaled>5e-8 or hrel>5e-8:
            raise RuntimeError(
                f"Physical SWE MPI disparity {case}: "
                f"state={scaled:.9g}, Hamiltonian={hrel:.9g}"
            )
    summary={
        "schema":"scrimp-v130-original-SWE-short-explicit-TS-v1",
        "status":"PASS",
        "original_physical_weak_forms":True,
        "published_and_nonvacuous_manufactured_initial_data":True,
        "actual_petsc_TS_accepted_steps":True,
        "explicit_brick_frozen_during_main_newton":True,
        "MPI_max_scaled_full_state_difference":float(max_scaled),
        "MPI_max_relative_Hamiltonian_difference":float(max_h_rel),
        "exact_energy_preservation_claimed":False,
        "general_nonlinear_DAE_initializer_claimed":False,
        "published_signed_H_change":float(
            reports["published",1]["signed_H_final_minus_H_initial"]
        ),
        "manufactured_signed_H_change":float(
            reports["manufactured",1]["signed_H_final_minus_H_initial"]
        ),
    }
    (folder/"shallow_explicit_TS_comparison.json").write_text(
        json.dumps(summary,indent=2)+"\n",encoding="utf-8"
    )
    print("SCRIMP shallow-water REAL BDF2 short TS, fully explicit G")
    print("| Initial p | MPI | Accepted | G RHS max | signed H drift |")
    print("|---|---:|---:|---:|---:|")
    for case in ("published","manufactured"):
        for ranks in (1,2):
            r=reports[case,ranks]
            print(f"| {case} | {ranks} | {r['TS_accepted_steps']} | "
                  f"{r['max_nonzero_explicit_G_norm']:.9g} | "
                  f"{r['signed_H_final_minus_H_initial']:.9g} |")
    print(
        f"PASS: true PETSc TS steps, accepted-state explicit G cache, "
        f"MPI state parity={max_scaled:.9g}, "
        f"H parity={max_h_rel:.9g}. "
        "Signed H drift is MEASURED, not presumed zero."
    )
    return summary


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("output_dir",type=Path)
    args=p.parse_args()
    compare(args.output_dir)


if __name__=="__main__":
    main()
