# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Native original Heat-Wave CN trajectory with DAE-compatible Gaussian.

This diagnostic intercepts only the original DPHS.init_step to correct
the projected p boundary trace, then calls the original *strict* DAE
initializer and continues into the true public Heat-Wave CN example,
MUMPS or SuperLU_DIST, using its unchanged mesh, weak forms, plot,
and ParaView outputs. No published source edits. MPI=1 only until
public ParaView export can be independently qualified for MPI>1.

Reuse the existing original main example evaluator for 4-step
numerical acceptance, Hamiltonian terms, GetFEM physical dissipation,
CN midpoint port powers, and numerical backend evidence.
Scientific power/energy qualifications remain separate from LU parity.
"""

import argparse
import json
import os
import sys
from pathlib import Path

import numpy as np
from mpi4py import MPI
from petsc4py import PETSc

from scrimp.dphs import DPHS
from probe_v130_heat_wave_boundary_compatibility import (
    wave_boundary_partition, constrained_l2_field,
)
import probe_v130_main_2d_examples as original_main


def add_boundary_compatible_state(system, candidate, report):
    """Change only the projected p field; call the real strict DAE gate."""
    original_p = np.asarray(
        system.gf_model.variable("p"), dtype=np.float64
    ).reshape(-1).copy()
    original_T = np.asarray(
        system.gf_model.variable("T")
    ).reshape(-1).copy()
    original_q = np.asarray(
        system.gf_model.variable("q")
    ).reshape(-1).copy()
    b, f, _, _ = wave_boundary_partition(
        system.ports["p"].get_fem(), field_size=len(original_p)
    )
    if len(b) != 126 or len(original_p) != 1105:
        raise RuntimeError("Unexpected original Heat-Wave CG1 boundary layout")
    if candidate == "clamp":
        corrected = original_p.copy()
        corrected[b] = 0.
        metrics = {}
        if not np.array_equal(corrected[f], original_p[f]):
            raise RuntimeError("Clamp changed an interior DOF")
    elif candidate == "constrained_l2":
        corrected, metrics = constrained_l2_field(
            system, original_p, b, f
        )
    else:
        raise ValueError("Unknown explicit boundary initial-data candidate")
    if np.any(corrected[b]) or not np.all(np.isfinite(corrected)):
        raise RuntimeError("Outer CG1 Dirichlet coefficients not zero")
    system.gf_model.set_variable(
        "p", np.ascontiguousarray(corrected, dtype=np.float64)
    )
    original_init = report.pop("_original_init")
    original_init(system)  # actual strict PETSc/LSQR algebraic initialization
    native = system.initialization_report
    if (native is None or native["status"] != "PASS"
            or not native["differential_initial_exact"]
            or native["pseudo_ts_executed"]
            or native["native_final_residual_l2"] >
               native["feasibility_tolerance"]):
        raise RuntimeError("Boundary-corrected original physical DAE not feasible")
    for name, initial in (("T",original_T),("q",original_q),("p",corrected)):
        current=np.asarray(system.gf_model.variable(name)).reshape(-1)
        if not np.array_equal(current,initial):
            raise RuntimeError(
                f"Strict DAE initializer modified differential state {name}"
            )
    if os.environ.get("SCRIMP_HW_TS_DIAGNOSTIC_MONITORS") == "1":
        # PETSc monitor options AFTER strict DAE, BEFORE main TS setup.
        # Observability only: CN method, dt and tolerances are unchanged.
        for name in ("snes_monitor", "snes_converged_reason",
                     "ksp_converged_reason", "ts_adapt_monitor"):
            system.time_scheme[name] = ""
    report.update({
        "candidate":candidate,
        "boundary_filtered_DOF_count":int(len(b)),
        "wave_region_DOF_count":int(len(original_p)),
        "boundary_initial_coefficient_l2":float(np.linalg.norm(original_p[b])),
        "boundary_corrected_coefficient_l2":float(np.linalg.norm(corrected[b])),
        "wave_coeff_relative_change":float(
            np.linalg.norm(corrected-original_p)/
            max(np.linalg.norm(original_p),1e-30)
        ),
        "initial_field_adjustment":metrics,
        "T_q_unchanged":True,
        "DAE_initializer":native,
        "physical_CN_main_TS_executed":False,
    })


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate",
                        choices=("clamp","constrained_l2"),required=True)
    parser.add_argument("--backend",
                        choices=("mumps","superlu_dist"),required=True)
    parser.add_argument("--out",type=Path,required=True)
    parser.add_argument("--capture-ts-failure",action="store_true")
    args=parser.parse_args()
    if MPI.COMM_WORLD.size != 1 or PETSc.COMM_WORLD.getSize() != 1:
        raise RuntimeError("Original public Heat-Wave export is MPI1 only")
    if not PETSc.Sys.hasExternalPackage(args.backend):
        raise RuntimeError("PETSc LU solver not available: "+args.backend)
    if args.out.exists() and any(args.out.iterdir()):
        raise RuntimeError("Native evidence output folder must be new")
    saved_argv=sys.argv
    original_init=DPHS.init_step
    captured={"_original_init":original_init}
    call_count=0
    def before_physical_TS(system):
        nonlocal call_count
        call_count+=1
        if call_count != 1:
            raise RuntimeError("Unexpected multiple initialization calls")
        return add_boundary_compatible_state(system,args.candidate,captured)
    DPHS.init_step=before_physical_TS
    try:
        sys.argv=[
            str(Path(__file__).name),
            "--example","heat_wave",
            "--backend",args.backend,
            "--out",str(args.out),
        ]
        try:
            original_main.main()
        except PETSc.Error as exc:
            if not args.capture_ts_failure:
                raise
            native = captured.get("DAE_initializer")
            if (native is None or native["status"]!="PASS"
                    or "DIVERGED_STEP_REJECTED" not in str(exc)):
                raise
            failed={
                "schema":"scrimp-v130-heat-wave-superlu-ts-failure-v1",
                "status":"DIAGNOSTIC_TS_REJECTED_NOT_PHYSICAL_PASS",
                "candidate":args.candidate,"backend":args.backend,
                "PETSc_error":str(exc),
                "strict_DAE_pass":True,
                "initial_algebraic_residual":
                    native["native_final_residual_l2"],
                "initial_algebraic_tolerance":
                    native["feasibility_tolerance"],
                "full_CN_completed":False,
                "physics_qualified":False,
                "original_example_modified":False,
            }
            output=args.out/"compatible_cn_failure.json"
            if output.exists():
                raise RuntimeError("Refuse to overwrite TS failure evidence")
            output.write_text(json.dumps(failed,indent=2)+"\n",
                              encoding="utf-8")
            print("DIAGNOSTIC CN TS REJECTED "+args.candidate+" "+
                  args.backend+" (strict DAE PASS, PHYSICS OPEN)",flush=True)
            return
    finally:
        DPHS.init_step=original_init
        sys.argv=saved_argv
        captured.pop("_original_init",None)
    if call_count != 1 or not captured.get("DAE_initializer"):
        raise RuntimeError("Original physical Heat-Wave DAE not verified")

    original_report_path=args.out/"original_2d_example_report.json"
    physics_path=args.out/"heat_wave_physics_diagnostic.json"
    evidence=json.loads(original_report_path.read_text(encoding="utf-8"))
    physics=json.loads(physics_path.read_text(encoding="utf-8"))
    if evidence["status"]!="PASS_NATIVE_INTEGRATION_HEAT_WAVE_PHYSICAL_AUDIT_PENDING":
        raise RuntimeError("Original numerical CN solver report not PASS")
    if (evidence["steps"]!=4 or evidence["rejected_steps"]!=0
            or evidence["SNES_failures"]!=0 or evidence["PETSc_reason"]<=0
            or evidence["full_state_DOFs"]!=7689):
        raise RuntimeError("Original Heat-Wave CN solver contract failed")
    if (evidence["factor_solver"]!=args.backend
            or not physics["CN_step_powers_are_assembled_at_mean_FE_state"]):
        raise RuntimeError("Original physical stage power measure missing")

    # Independent, explicit, auditable candidate scientific thresholds.
    # A failed energy/power test is reported, NOT silently turned into
    # a successful physical-model acceptance.
    energy_scale=max(1.,abs(evidence["Hamiltonian_initial"]))
    power_scale=max(
        1.,float(physics["dissipation_max"]),
        abs(float(physics["max_interface_port_power_sum"]))
    )
    limits={
        "max_CN_step_H_balance_defect":1.e-7*energy_scale,
        "max_CN_step_heat_balance_defect":1.e-7*energy_scale,
        "max_CN_step_wave_balance_defect":1.e-7*energy_scale,
        "max_CN_interface_power_defect":1.e-7*power_scale,
        "max_CN_positive_H_increment":1.e-7*energy_scale,
    }
    gates={}
    for key,threshold in limits.items():
        value=float(physics[key])
        gates[key]={
            "observed":value,"threshold":threshold,
            "passed":bool(np.isfinite(value) and value <= threshold),
        }
    nonvacuous={
        "positive_H_at_start":evidence["Hamiltonian_initial"]>0.,
        "nonzero_heat_dissipation":physics["dissipation_max"]>1.e-12,
        "nonzero_physical_dynamics":abs(evidence["Hamiltonian_delta"])>1.e-12,
    }
    science_pass=all(x["passed"] for x in gates.values()) and all(
        nonvacuous.values()
    )
    status=(
        "PASS_CN_HEAT_WAVE_SCIENTIFIC_THRESHOLDS_MPI1"
        if science_pass else
        "CN_NUMERICAL_PASS_PHYSICAL_BALANCE_NOT_QUALIFIED"
    )
    report={
        "schema":"scrimp-v130-heat-wave-boundary-corrected-cn-v1",
        "status":status,
        "candidate":args.candidate,
        "backend":args.backend,
        "fixed_Crank_Nicolson":True,
        "time_steps":4,
        "original_example_source_modified":False,
        "original_L2_projection_preserved":True,
        "full_original_CN_TS_executed":True,
        "true_strict_DAE_pre_Ts":True,
        "DAE_correction":captured,
        "original_numerical_report":evidence,
        "original_physics_report":physics,
        "scientific_thresholds":gates,
        "nonvacuous_controls":nonvacuous,
        "physical_balance_gate_passed":science_pass,
        "MPI2_or_parallel_export_qualified":False,
    }
    output=args.out/"compatible_cn_report.json"
    if output.exists():
        raise RuntimeError("Refuse to overwrite native CN evidence")
    output.write_text(json.dumps(report,indent=2)+"\n",
                      encoding="utf-8")
    print(
        f"ORIGINAL HEAT-WAVE CN {args.candidate} {args.backend} "
        f"4 steps: {status}; "
        f"DAE={captured['DAE_initializer']['native_final_residual_l2']:.4g}; "
        f"max_step_H_defect="
        f"{physics['max_CN_step_H_balance_defect']:.4g}; "
        f"interface_power="
        f"{physics['max_CN_interface_power_defect']:.4g}",
        flush=True,
    )
    # The outer comparison retains full diagnostic evidence even if
    # the scientific power balance thresholds are not met.


if __name__=="__main__":
    main()
