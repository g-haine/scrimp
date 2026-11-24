# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Compare complete original Heat-Wave CN trajectories and pH power gates.

Four independently executed native full CN runs, 4 steps each,
each with original strict algebraic DAE:
  - clamp: MUMPS versus row-matched SuperLU_DIST;
  - constrained_L2: MUMPS versus row-matched SuperLU_DIST.

Check complete PETSc state, all existing GetFEM FE observables, and
independent nonvacuous *interface exchange* (individual port power).
Small interface POWER SUM alone is not enough: that can vanish when
both interface powers vanish. Exact physical-balance PASS on weakly
coupled examples is NOT a release-wide multiphysics certification.

No matrix, model, or TS mutation in this comparison-only utility.
"""

import argparse
import json
import math
from pathlib import Path
import zipfile

import numpy as np

CASES=("clamp","constrained_l2")
SOLVERS=("mumps","matched_superlu")
REQUIRED_FIELDS=(
    "times", "states", "H", "T", "p", "q",
    "phys_H_heat", "phys_H_wave_kinetic",
    "phys_H_wave_potential", "phys_CN_D_heat_flux",
    "phys_CN_P_interface_heat", "phys_CN_P_interface_wave",
    "phys_CN_P_interface_physical",
    "phys_CN_per_step_H_balance_error",
    "phys_CN_per_step_heat_balance_error",
    "phys_CN_per_step_wave_balance_error",
    "phys_CN_cumulative_total_balance_error",
)


def load_case(root, initial, solver):
    folder=root/f"{initial}_{solver}"
    report=json.loads(
        (folder/"compatible_cn_report.json").read_text(encoding="utf-8")
    )
    numerical=report["original_numerical_report"]
    physical=report["original_physics_report"]
    if (report["candidate"]!=initial
            or report["backend"]!=(
                "mumps" if solver=="mumps" else "superlu_dist"
            )
            or numerical["steps"]!=4
            or numerical["SNES_failures"]!=0
            or numerical["rejected_steps"]!=0
            or numerical["PETSc_reason"]<=0
            or numerical["full_state_DOFs"]!=7689
            or report["DAE_correction"]["DAE_initializer"]["status"]!="PASS"
            or report["DAE_correction"]["DAE_initializer"][
                "native_final_residual_l2"
            ]>1e-11
            or report["DAE_correction"]["DAE_initializer"]["pseudo_ts_executed"]
            or not report["physical_balance_gate_passed"]
            or not physical["CN_step_powers_are_assembled_at_mean_FE_state"]):
        raise RuntimeError(
            f"Missing real, independently qualified strict CN run: {initial} {solver}"
        )

    matched_report=None
    if solver=="matched_superlu":
        matched_report=json.loads(
            (folder/"rowmatched_cn_audit.json").read_text(encoding="utf-8")
        )
        pc=matched_report["native_rowmatched_PC"]
        if (not matched_report["true_original_Ts_4steps"]
                or matched_report["automatically_enabled_in_production"]
                or matched_report["MPI2_qualified"]
                or pc["setup_count"]<4
                or pc["apply_count"]<4
                or pc["max_relative_residual_original_preconditioner"]>1e-9):
            raise RuntimeError("Row-matched SuperLU_DIST PC not qualified")

    with np.load(
        folder/"original_2d_example_trace.npz",allow_pickle=False
    ) as npz:
        arrays={name:np.asarray(npz[name]).copy()
                for name in npz.files}
    if not set(REQUIRED_FIELDS).issubset(arrays):
        raise RuntimeError("Missing original physical GetFEM observables")
    if not np.array_equal(
        arrays["times"],np.linspace(0.,0.004,5)
    ):
        # floating-point TS timestamps can differ at roundoff
        if not np.allclose(
            arrays["times"],np.arange(5)*0.001,rtol=0.,atol=1e-12
        ):
            raise RuntimeError("Original CN saved time grid changed")
    if arrays["states"].shape!=(5,7689):
        raise RuntimeError("Original full-state shape changed")
    for key,arr in arrays.items():
        if not np.issubdtype(arr.dtype, np.number) or not np.all(
            np.isfinite(arr)
        ):
            raise RuntimeError(f"Non-finite/non-numeric original trace: {key}")

    interface=max(
        float(np.max(np.abs(arrays["phys_CN_P_interface_heat"]))),
        float(np.max(np.abs(arrays["phys_CN_P_interface_wave"]))),
    )
    dissip=float(np.max(np.abs(arrays["phys_CN_D_heat_flux"])))
    interface_ratio=interface/max(dissip,1e-30)
    interface_exchange_nontrivial=interface_ratio>1e-8
    physical_delta=float(numerical["Hamiltonian_delta"])
    if (not abs(physical_delta)>1e-10 or
            not arrays["phys_H_wave_potential"][-1]>1e-9 or
            not dissip>1e-12):
        raise RuntimeError("Heat-Wave nontrivial wave/dissipation gate failed")
    observed={
        "Hamiltonian_initial":float(numerical["Hamiltonian_initial"]),
        "Hamiltonian_final":float(numerical["Hamiltonian_final"]),
        "Hamiltonian_delta":physical_delta,
        "max_CN_total_energy_defect":
            float(physical["max_CN_step_H_balance_defect"]),
        "max_CN_interface_power_sum_defect":
            float(physical["max_CN_interface_power_defect"]),
        "max_individual_CN_interface_power":interface,
        "max_CN_heat_dissipation":dissip,
        "interface_exchange_ratio_to_heat_dissipation":interface_ratio,
        "interface_exchange_nontrivial":bool(interface_exchange_nontrivial),
        "nodal_trapezoid_dissipation_defect_DIAGNOSTIC_ONLY":float(
            physical["max_abs_total_balance_nodal_trapz_defect"]
        ),
        "CN_midpoint_energy_balance_pass":True,
        "full_coupling_excitation_qualified":bool(interface_exchange_nontrivial),
        "strict_DAE_residual":float(report["DAE_correction"][
            "DAE_initializer"]["native_final_residual_l2"]),
        "rowmatched_PC":matched_report is not None,
    }
    return arrays,observed


def compare_complete_traces(left,right):
    if set(left)!=set(right):
        raise RuntimeError("MUMPS and matched LU have different observable fields")
    errors={}
    failed={}
    for name in sorted(left):
        a,b=left[name],right[name]
        if a.shape!=b.shape:
            failed[name]="shape mismatch"
            continue
        diff=np.asarray(a)-np.asarray(b)
        if name=="states":
            absolute=float(np.max(np.linalg.norm(diff,axis=1)))
            reference=float(np.max(np.linalg.norm(a,axis=1)))
        else:
            absolute=float(np.max(np.abs(diff)))
            reference=float(np.max(np.abs(a)))
        # Both absolute and magnitude-aware relative checks.
        scaled=absolute/max(1.,reference)
        okay=absolute<=1e-8*max(1.,reference)
        errors[name]={
            "max_absolute":absolute,
            "max_scaled":scaled,
            "reference_scale":reference,
            "passes":bool(okay),
        }
        if not okay:
            failed[name]={"max_absolute":absolute,"max_scaled":scaled}
    return errors,failed


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root",type=Path,required=True)
    args=parser.parse_args()
    root=args.root
    if not root.is_dir() or (root/"comparison.json").exists():
        raise RuntimeError("Require original output folder, no previous summary")

    summaries={}
    all_errors={}
    all_pass=True
    all_interface_excited=True
    for initial in CASES:
        report={}
        traces={}
        for solver in SOLVERS:
            traces[solver],report[solver]=load_case(root,initial,solver)
            all_interface_excited &= report[solver][
                "interface_exchange_nontrivial"
            ]
        errors,failed=compare_complete_traces(
            traces["mumps"],traces["matched_superlu"]
        )
        summaries[initial]={
            "MUMPS":report["mumps"],
            "matched_SuperLU_DIST":report["matched_superlu"],
            "trace_backend_parity":errors,
            "failing_observables":failed,
        }
        if failed:
            all_pass=False
            all_errors[initial]=failed

    result={
        "schema":"scrimp-v130-original-HW-4-CN-matched-LU-parity-v1",
        "status": (
            "FAIL_NATIVE_FULL_STATE_BACKEND_PARITY"
            if not all_pass else
            "PASS_FOUR_CN_BACKEND_PARITY_AND_ENERGY_INTERFACE_EXCITATION_OPEN"
            if not all_interface_excited else
            "PASS_FOUR_CN_BACKEND_PARITY_AND_ENERGY_WITH_INTERFACE_EXCITATION"
        ),
        "source":"original Heat-Wave revised physical RT/CN example",
        "initial_candidates":list(CASES),
        "four_actual_full_CN_TS_runs":True,
        "matched_LU_in_production":False,
        "MPI2_or_parallel_export_qualified":False,
        "both_backends_strict_DAE_and_CN_energy_gates_passed":True,
        "interface_energy_exchange_excited_all_runs":bool(
            all_interface_excited
        ),
        "interface_excitation_criterion":
            "max(|individual CN interface port power|)"
            "/max(|CN heat dissipation|) > 1e-8",
        "cases":summaries,
        "failed_observable_backend_parity":all_errors,
        "nodal_trapezoid_not_used_for_CN_energy_acceptance":True,
        "source_examples_unchanged_by_this_campaign":True,
    }
    (root/"comparison.json").write_text(
        json.dumps(result,indent=2)+"\n",encoding="utf-8"
    )
    archive=Path(str(root)+".zip")
    with zipfile.ZipFile(
        archive,"x",compression=zipfile.ZIP_DEFLATED,
        compresslevel=4,allowZip64=True
    ) as z:
        for path in sorted(root.rglob("*")):
            if path.is_file():
                z.write(path,arcname=str(path.relative_to(root)))
    concise={
        "status":result["status"],
        "four_actual_full_CN_TS_runs":True,
        "interface_energy_exchange_excited_all_runs":bool(
            all_interface_excited
        ),
        "candidates":{},
    }
    for initial in CASES:
        entry=summaries[initial]
        concise["candidates"][initial]={
            "MUMPS":entry["MUMPS"],
            "matched_SuperLU_DIST":entry["matched_SuperLU_DIST"],
            "max_state_scaled_backend_difference":entry[
                "trace_backend_parity"]["states"]["max_scaled"],
            "max_H_scaled_backend_difference":entry[
                "trace_backend_parity"]["H"]["max_scaled"],
            "failing_observables":entry["failing_observables"],
        }
    print(json.dumps(concise,indent=2))
    print("SHAREABLE ORIGINAL HEAT-WAVE FOUR-CN PARITY EVIDENCE:",archive)
    if not all_pass:
        raise SystemExit("Full physical-state backend parity failed; archive kept")


if __name__=="__main__":
    main()
