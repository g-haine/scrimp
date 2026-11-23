# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""#30 OFFLINE audit: apply the new public plot quadrature to original
immutable physical linear-wave DAE GetFEM/PETSc CN port powers.

This requires no Docker, no new simulation, no modification of original
v1/v2 physical evidence and never imports PETSc/GetFEM.
The new helper's raw signed integrals are turned into *supplied work*
with the wave model's separately demonstrated W=-int(P) orientation.
A public plotting smoke is still a distinct later acceptance gate.
"""

import argparse
import importlib.util
import json
from pathlib import Path

import numpy as np


ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location(
    "_v130_native_CN_power_integration",
    ROOT/"scrimp"/"utils"/"power_integration.py",
)
integration=importlib.util.module_from_spec(spec)
spec.loader.exec_module(integration)


CASES=(
    ("seed",1,None),("seed",2,None),
    ("continuous",1,None),("continuous",2,None),
    ("resume",1,1),("resume",2,1),
    ("resume",1,2),("resume",2,2),
)


def name(phase,ranks,source):
    return (f"wave_restart_from{source}_to{ranks}" if phase=="resume"
            else f"wave_{phase}_mpi{ranks}")


def segment(directory):
    report=json.loads(
        (directory/"linear_report.json").read_text(encoding="utf-8")
    )
    with np.load(directory/"linear_states.npz",allow_pickle=False) as z:
        t=np.asarray(z["times"],dtype=float)
        H=np.asarray(z["H"],dtype=float)
        p=np.asarray(z["port_power"],dtype=float)
    if (report["PETSc_reason"]<=0 or report["PETSc_rejected"]
            or report["SNES_failures"]
            or not (t.shape==H.shape==p.shape)):
        raise ValueError("Missing or failed original physical PETSc CN segment")
    raw=integration.integrate_port_power(t,p,convention="cn_step")
    legacy=integration.integrate_port_power(t,p,convention="legacy_trapezoid")
    work=-raw  # Original linear DAE left control port orientation, not generic.
    step=np.diff(H)+np.diff(t)*p[1:]
    defect=float(np.max(np.abs(step)))
    if defect>=1e-11:
        raise ValueError(f"Physical original CN step power mismatch: {defect}")
    if float(np.max(abs(H-H[0]-work)))>=1e-11:
        raise ValueError("Wrong signed physical work from original CN port")
    return t,H,p,raw,work,legacy,defect


def audit(root):
    root=Path(root)
    traces={}
    wrong_trapezoid=0.
    max_step=0.
    for phase,ranks,source in CASES:
        case=name(phase,ranks,source)
        folder=root/case
        t,H,p,raw,work,legacy,defect=segment(folder)
        traces[case]=(t,H,p,raw,work)
        wrong_trapezoid=max(wrong_trapezoid,
                            float(np.max(abs(legacy-raw))))
        max_step=max(max_step,defect)
    max_restart=0.
    for source in (1,2):
        for target in (1,2):
            seed=traces[name("seed",source,None)]
            full=traces[name("continuous",target,None)]
            resumed=traces[name("resume",target,source)]
            ts=np.r_[seed[0],resumed[0][1:]]
            if not np.allclose(ts,full[0],atol=1e-12,rtol=0):
                raise ValueError("Physical absolute time mismatch")
            split=np.r_[seed[4],seed[4][-1]+resumed[4][1:]]
            diff=float(np.max(abs(split-full[4])))
            if diff>=1e-10:
                raise ValueError(f"Physical CN port work restart mismatch: {diff}")
            max_restart=max(max_restart,diff)
    if wrong_trapezoid<=1e-7:
        raise ValueError(
            "False negative control: original CN midpoint/trapezoid "
            "conventions were not meaningfully distinguished"
        )
    return {
        "schema":"scrimp-v130-issue30-offline-native-port-quadrature-audit-v1",
        "status":"PASS_PUBLIC_CN_STEP_INTEGRAL_ON_ORIGINAL_NATIVE_PORT_ARRAYS",
        "native_source_segments_reused":len(CASES),
        "original_files_modified":False,
        "physical_integrator_run":False,
        "physical_port_sign_convention":"wave_left_control_W=-sum(dt*P[k>=1])",
        "maximum_original_native_CN_per_step_H_plus_P_dt":max_step,
        "max_reconstructed_continuous_vs_split_work_difference":max_restart,
        "max_old_wrong_trapezoid_vs_correct_CN_integral":wrong_trapezoid,
        "generic_other_TS_not_qualified":True,
        "public_matplotlib_plot_smoke_not_claimed":True,
    }


def main():
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--evidence",type=Path,required=True)
    cli.add_argument("--report",type=Path)
    args=cli.parse_args()
    result=audit(args.evidence)
    if args.report is not None:
        if args.report.exists():
            raise FileExistsError("Refusing to overwrite original evidence")
        args.report.parent.mkdir(parents=True,exist_ok=True)
        args.report.write_text(json.dumps(result,indent=2)+"\n",
                               encoding="utf-8")
    print(json.dumps(result,indent=2))
    print("PASS #30 original native GetFEM CN powers re-integrated "
          "without Docker (NOT a public matplotlib smoke)")


if __name__=="__main__":
    main()
