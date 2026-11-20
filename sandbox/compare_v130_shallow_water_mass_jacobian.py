# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Require genuine nonzero-time-derivative mass-Jacobian accuracy in MPI.

Original shipped shallow_water physics, no TS solve; the corrected
default is evaluated against independent F(z +/- epsilon*v,
zd +/- shift*epsilon*v) and historical legacy on exactly the same
states. We do not claim general nonlinear DAE initialization here.
"""

import argparse
import json
from pathlib import Path

import numpy as np


def compare(root):
    data={}
    for ranks in (1,2):
        for mode in ("default","legacy"):
            src = root/f"shallow_{mode}_{ranks}r.json"
            r=json.loads(src.read_text(encoding="utf-8"))
            if (r["mode"]!=mode or r["ranks"]!=ranks
                or not r["nonlinear_mass"] or r["main_TS_executed"]
                or not r["independent_time_derivative"]
                or r["dofs"]<20
                or len(r["outcomes"])!=4):
                raise RuntimeError(f"Wrong native Jacobian record: {src.name}")
            if mode=="default" and r["mass_mode_env"] is not None:
                raise RuntimeError("Default correction was explicitly configured")
            if mode=="legacy" and r["mass_mode_env"]!="legacy":
                raise RuntimeError("Legacy correction not explicitly selected")
            actual=[(float(x["shift"]),float(x["zd_norm"]))
                    for x in r["outcomes"]]
            if ([x[0] for x in actual]!=[0.,1.,0.,1.]
                or actual[0][1]!=0 or actual[1][1]!=0
                or actual[2][1]<=0 or actual[3][1]<=0):
                raise RuntimeError("Missing independent zd or shift tests")
            errors=np.array([float(x["relative_error"]) for x in r["outcomes"]])
            if not np.all(np.isfinite(errors)):
                raise RuntimeError("Nonfinite relative Jacobian error")
            data[ranks,mode]=r
            if mode=="default" and np.max(errors)>1e-5:
                raise RuntimeError(
                    f"Incorrect default nonlinear Jacobian MPI={ranks}: "
                    f"max error={np.max(errors):.9g} > 1e-5"
                )
        correction = np.array([
            x["relative_error"] for x in data[ranks,"default"]["outcomes"][2:]
        ])
        historical = np.array([
            x["relative_error"] for x in data[ranks,"legacy"]["outcomes"][2:]
        ])
        if not np.max(historical)>max(
            10*np.max(correction),1e-7
        ):
            raise RuntimeError(
                "State-dependent mass correction not demonstrably "
                f"needed on original shallow water MPI={ranks}: "
                f"default={correction.tolist()}, legacy={historical.tolist()}"
            )
    errors1=np.array([x["relative_error"] for x in data[1,"default"]["outcomes"]])
    errors2=np.array([x["relative_error"] for x in data[2,"default"]["outcomes"]])
    if np.max(abs(errors1-errors2))>1e-5:
        raise RuntimeError(
            "Default Jacobian directional errors differ across MPI"
        )
    summary={
        "schema":"scrimp-v130-original-shallow-mass-jacobian-gate-v1",
        "status":"PASS",
        "original_shallow_water_weak_forms":True,
        "default_finite_difference_used_without_override":True,
        "nonlinear_time_derivative_checked":True,
        "max_error_default_M1":float(np.max(errors1)),
        "max_error_default_M2":float(np.max(errors2)),
        "max_error_legacy_M1":float(max(
            x["relative_error"] for x in data[1,"legacy"]["outcomes"]
        )),
        "max_error_legacy_M2":float(max(
            x["relative_error"] for x in data[2,"legacy"]["outcomes"]
        )),
        "native_state_dependent_mass_Jacobian_qualified":True,
        "nonlinear_DAE_initializer_qualified":False,
    }
    (root/"shallow_mass_jacobian_qualification.json").write_text(
        json.dumps(summary,indent=2)+"\n",encoding="utf-8"
    )
    print("Original shallow water: production default mass Jacobian")
    print("| MPI | Corrected max FD rel error | Legacy max FD rel error |")
    print("|---:|---:|---:|")
    for ranks in (1,2):
        print(
            f"| {ranks} | "
            f"{data[ranks,'default']['maximum_directional_relative_error']:.9g} | "
            f"{data[ranks,'legacy']['maximum_directional_relative_error']:.9g} |"
        )
    print("PASS: true shipped shallow water mass-dependent Jacobian, "
          "MPI1/2, nonzero zd, independent directional F and J; "
          "historical legacy control is worse")
    return summary


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("output_dir",type=Path)
    args=ap.parse_args()
    compare(args.output_dir)


if __name__=="__main__":
    main()
