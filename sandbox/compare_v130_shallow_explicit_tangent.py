# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Isolate the gyroscopic G tangent omission in the shipped shallow water.

This is an explanatory A/B diagnostic, NOT qualification of a production
change. It does NOT loosen the full residual/Jacobian <1e-5 gate:
scripts/run-v130-shallow-mass-jacobian.sh subsequently runs the original
strict production comparator, which must still fail if uncorrected.
"""

import argparse
import json
from pathlib import Path

import numpy as np


def compare(folder):
    records={}
    for ranks in (1,2):
        for tangent in ("production","diagnostic_full"):
            for mass in ("default","legacy"):
                tag=mass if tangent=="production" else f"{mass}_diagnostic_full"
                name=f"shallow_{tag}_{ranks}r.json"
                r=json.loads((folder/name).read_text(encoding="utf-8"))
                if (r.get("mode")!=mass or r.get("ranks")!=ranks
                    or r.get("tangent_policy")!=tangent
                    or r.get("gyro_G_explicit_restored") is not True
                    or r.get("nonlinear_mass") is not True
                    or r.get("main_TS_executed") is not False
                    or len(r.get("outcomes",[]))!=4):
                    raise RuntimeError(f"Inconsistent explicit-tangent A/B report {name}")
                if (mass=="default" and r.get("mass_mode_env") is not None):
                    raise RuntimeError(f"Correction not tested as actual default in {name}")
                if (mass=="legacy" and r.get("mass_mode_env")!="legacy"):
                    raise RuntimeError(f"Historical mass control not explicitly selected in {name}")
                rows=r["outcomes"]
                for i,row in enumerate(rows):
                    if row.get("shift")!=(0. if i%2==0 else 1.):
                        raise RuntimeError(f"Wrong shift sequence in {name}")
                    zd_norm=float(row.get("zd_norm",float("nan")))
                    if (zd_norm!=0 if i<2 else zd_norm<=0):
                        raise RuntimeError(f"Wrong time-derivative control in {name}")
                    if not np.isfinite(float(row["relative_error"])):
                        raise RuntimeError(f"Nonfinite Jacobian discrepancy in {name}")
                records[ranks,tangent,mass]=r

    # The direction, physical residual and FD inputs are IDENTICAL
    # between production and diagnostic tangent tests. The sole
    # change is inclusion of the already-compiled explicit G term
    # in the tangent, restored before IFunction is evaluated.
    for ranks in (1,2):
        for mass in ("default","legacy"):
            prod=records[ranks,"production",mass]["outcomes"]
            full=records[ranks,"diagnostic_full",mass]["outcomes"]
            for i in range(4):
                ap=float(prod[i]["finite_difference_norm"])
                bp=float(full[i]["finite_difference_norm"])
                if abs(ap-bp)>1e-8*max(1.,abs(ap),abs(bp)):
                    raise RuntimeError(
                        "Tangent-only hypothesis contaminated the physical "
                        f"residual MPI={ranks}, mass={mass}, case={i}"
                    )

    print("Original shallow-water full residual FD versus nonlinear Jacobian")
    print("| MPI | Tangent | Mass correction | max all | zd=0 shift=0 | zd!=0 shift=0 |")
    print("|---:|---|---|---:|---:|---:|")
    for ranks in (1,2):
        for tangent in ("production","diagnostic_full"):
            for mass in ("default","legacy"):
                r=records[ranks,tangent,mass]
                errors=[float(x["relative_error"]) for x in r["outcomes"]]
                print(f"| {ranks} | {tangent} | {mass} | "
                      f"{max(errors):.9g} | {errors[0]:.9g} | {errors[2]:.9g} |")

    prod0=records[1,"production","default"]["outcomes"][0]["relative_error"]
    full0=records[1,"diagnostic_full","default"]["outcomes"][0]["relative_error"]
    full_nonzero=max(
        x["relative_error"] for x in
        records[1,"diagnostic_full","default"]["outcomes"][2:]
    )
    legacy_nonzero=max(
        x["relative_error"] for x in
        records[1,"diagnostic_full","legacy"]["outcomes"][2:]
    )
    # Keep exactly the existing full-Jacobian acceptance threshold.
    threshold=1e-5
    hypothesis=(
        prod0>threshold
        and full0<=threshold
        and full_nonzero<=threshold
        and legacy_nonzero>max(10*full_nonzero,1e-7)
    )
    summary={
        "schema":"scrimp-v130-shallow-explicit-G-tangent-hypothesis-v1",
        "status":"DIAGNOSTIC_HYPOTHESIS_SUPPORTED" if hypothesis else
                 "DIAGNOSTIC_HYPOTHESIS_UNRESOLVED",
        "original_shallow_water_formulation_unchanged":True,
        "production_implicit_Jacobian_qualified":False,
        "production_zero_zd_error":float(prod0),
        "diagnostic_full_tangent_zero_zd_error":float(full0),
        "diagnostic_full_tangent_nonzero_zd_default_error":float(full_nonzero),
        "diagnostic_full_tangent_nonzero_zd_legacy_error":float(legacy_nonzero),
        "full_Jacobian_original_acceptance_threshold":threshold,
    }
    (folder/"shallow_explicit_tangent_diagnosis.json").write_text(
        json.dumps(summary,indent=2)+"\n",encoding="utf-8"
    )
    if hypothesis:
        print("DIAGNOSTIC HYPOTHESIS SUPPORTED: explicit nonlinear G "
              "present in IFunction but missing from fully implicit IJacobian; "
              "mass Hessian default improves the now-consistent tangent.")
        print("This does NOT modify production or qualify the release.")
    else:
        print("DIAGNOSTIC HYPOTHESIS NOT FULLY SUPPORTED: inspect native "
              "reports before making any production change.")
    return summary


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("output_dir",type=Path)
    args=ap.parse_args()
    compare(args.output_dir)


if __name__=="__main__":
    main()
