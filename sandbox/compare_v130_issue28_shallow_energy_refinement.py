# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Host-only audit of native #28 SWE energy defect time refinement, MPI 1/2.

No GetFEM/PETSc/Docker required. Compare full archived physical-time
reports at the SAME final time and SAME original physical mesh/initial
data. Measure, do not presume, a temporal convergence rate for the
signed port-power balance defect H(T)-H(0)-integrated boundary power.

The original initial free surface is nonsmooth, and its BDF2/explicit
gyro splitting is not a discrete-gradient scheme. We do NOT assume
a theoretical order from three short time steps, or call this a
full power-preservation certificate without time/mesh evidence.
"""

import argparse
import json
import math
from pathlib import Path
import re

EXPECTED_IMAGE_REV = "1d6d0c90115213298047b28237a7f38087563a35"
T = 0.0004
CASES = ((4, 0.0001), (8, 0.00005), (16, 0.000025))


def require(valid, message):
    if not valid:
        raise ValueError(message)


def read(path):
    require(path.is_file(), f"Missing original native evidence: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def finite(number, label):
    require(isinstance(number, (float, int)) and not isinstance(number, bool)
            and math.isfinite(number), f"Invalid {label}: {number!r}")
    return float(number)


def nearly(a, b, *, rel=0., abs_=1e-12):
    return math.isclose(a, b, rel_tol=rel, abs_tol=abs_)


def check(folder):
    folder = Path(folder).resolve()
    manifest = read(folder / "manifest.json")
    require(manifest.get("schema") == "issue28-shallow-energy-time-refinement-v1",
            "Unexpected archived native campaign type")
    require(manifest.get("image_revision") == EXPECTED_IMAGE_REV,
            "Unqualified Docker image revision")
    require(re.fullmatch(r"[0-9a-f]{40}", manifest.get("source", "")) is not None,
            "Missing full Git source SHA")
    require(manifest.get("original_physical_model") is True
            and manifest.get("spatial_mesh_step") == .1
            and nearly(manifest.get("time_final", -1), T),
            "The source physics, spatial mesh or final time changed")
    require(manifest.get("true_physical_step_counts") == [4, 8, 16],
            "Unexpected time-refinement test cases")

    report_by_case = {}
    max_h_mpi = 0.
    max_mass_mpi = 0.
    all_records = []

    for steps, dt in CASES:
        records = []
        for ranks in (1, 2):
            f = folder / f"result_steps{steps}_{ranks}" / "shallow_short_physical.json"
            r = read(f)
            require(r.get("schema") == "scrimp-v130-issue28-original-shallow-short-physical-v1"
                    and r.get("status") == "PASS_SHORT_PHYSICAL_TIME_GATE",
                    f"Native run was not qualified for {steps} steps MPI {ranks}")
            require(r.get("model") == "original shipped examples/shallow_water.py"
                    and r.get("mpi_ranks") == ranks and r.get("mesh_h") == .1,
                    f"Wrong physical case {steps}/{ranks}")
            require(r.get("original_prescribed_h0_p0_bitwise") is True
                    and r.get("original_boundary_controls_unchanged_zero") is True
                    and r.get("original_explicit_gyro_unchanged") is True,
                    f"Original physical model changed at {steps}/{ranks}")
            init = r.get("production_initializer_report", {})
            require(init.get("strategy") == "nonlinear_affine"
                    and init.get("status") == "PASS"
                    and init.get("differential_initial_exact") is True
                    and init.get("accuracy_target_met") is True,
                    f"Production initialization did not pass {steps}/{ranks}")
            fa = finite(init.get("native_final_algebraic_residual_l2"), "initial Fa")
            require(fa <= 1e-11, f"Initial algebraic residual too high {steps}/{ranks}")
            require(init.get("algebraic_null_columns") == 0,
                    f"Algebraic singularity undetected {steps}/{ranks}")
            require(r.get("expected_accepted_steps") == steps
                    and r.get("accepted_steps") == steps
                    and r.get("rejected_steps") == 0
                    and r.get("nonlinear_solver_failures") == 0
                    and r.get("ts_converged_reason", 0) > 0,
                    f"Actual PETSc TS failed to accept {steps} steps MPI {ranks}")
            require(nearly(finite(r.get("time_step_requested"), "dt"),
                           dt, abs_=1e-15)
                    and nearly(finite(r.get("t_final_actual"), "final time"),
                               T, abs_=1e-11)
                    and nearly(finite(r.get("t_final_requested"), "requested T"),
                               T, abs_=1e-15),
                    f"Temporal refinement did not keep fixed T {steps}/{ranks}")

            times = r.get("saved_times", [])
            h = r.get("hamiltonian", [])
            m = r.get("mass", [])
            p = r.get("total_boundary_power", [])
            work = r.get("cumulative_boundary_power", [])
            defect = r.get("energy_balance_defect", [])
            require(len(times) == len(h) == len(m) == len(p)
                    == len(work) == len(defect) >= 3,
                    f"Incomplete physical saved trajectory {steps}/{ranks}")
            require(all(math.isfinite(finite(t, "time")) for t in times)
                    and nearly(times[0], 0.)
                    and nearly(times[-1], T, abs_=1e-11)
                    and all(x < y for x, y in zip(times, times[1:])),
                    f"Missing or invalid time endpoints {steps}/{ranks}")
            # SCRIMP may save fewer states than accepted TS steps. Honor
            # ACTUAL sampling intervals, never invent missing values.
            require(r.get("saved_sample_count") == len(times),
                    "Saved trajectory count inconsistent")
            for key, values in (
                ("Hamiltonian", h), ("mass", m),
                ("boundary power", p), ("work", work), ("energy defect", defect),
            ):
                require(all(math.isfinite(finite(v, key)) for v in values),
                        f"Non-finite {key} for {steps}/{ranks}")
            require(min(h) > 0. and min(m) > 0.,
                    f"Nonphysical energy/mass for {steps}/{ranks}")
            require(r.get("Hamiltonian_GetFEM_and_registered_agree") is True,
                    f"GetFEM energy/registered energy mismatch {steps}/{ranks}")
            require(r.get("BDF2_exact_power_conservation_claimed") is False
                    and r.get("original_rotating_viscous_Y0_qualified") is False,
                    "Wrongly advertised exact energy conservation or rotating Y0")
            require(finite(r["relative_mass_drift"], "mass drift") <= 1e-5,
                    f"Excessive mass drift {steps}/{ranks}")
            require(finite(r["relative_hamiltonian_drift"], "H drift") <= 5e-3,
                    f"Excessive short-time Hamiltonian drift {steps}/{ranks}")
            require(finite(r["maximum_abs_boundary_power"], "max power") <= 1e-7,
                    f"Unexpected external work {steps}/{ranks}")
            for j, ts in enumerate(times):
                eq = h[j] - h[0] - work[j]
                require(nearly(eq, defect[j], rel=1e-10, abs_=2e-8),
                        f"Wrong actual energy defect {steps}/{ranks} at time {ts}")
            final = h[-1] - h[0] - work[-1]
            require(nearly(final, r["signed_final_energy_balance_defect"],
                           rel=1e-10, abs_=2e-8),
                    f"Final energy defect not measured {steps}/{ranks}")
            require(nearly(work[-1], 0., abs_=1e-10),
                    f"Unexpected total work with original zero controls {steps}/{ranks}")
            require(not r.get("physical_sanity_gate_failures", []),
                    f"Native physical sanity gate reported failures {steps}/{ranks}")
            records.append(r)

        a, b = records
        require(len(a["saved_times"]) == len(b["saved_times"])
                and all(nearly(x, y, abs_=1e-11)
                        for x, y in zip(a["saved_times"], b["saved_times"])),
                f"Different saved time grids MPI 1/2 for steps={steps}")
        for key in ("hamiltonian", "mass", "total_boundary_power"):
            for u, v in zip(a[key], b[key]):
                drift = abs(u-v)/max(1., abs(u), abs(v))
                if key == "hamiltonian":
                    max_h_mpi = max(max_h_mpi, drift)
                if key == "mass":
                    max_mass_mpi = max(max_mass_mpi, drift)
                require(drift <= 5e-8,
                        f"MPI mismatch {key} for {steps} steps: {drift:.9g}")

        report_by_case[steps] = records
        all_records.append({
            "requested_steps": steps,
            "dt": dt,
            "fixed_final_time": T,
            "signed_balance_defect_MPI1": float(a["signed_final_energy_balance_defect"]),
            "signed_balance_defect_MPI2": float(b["signed_final_energy_balance_defect"]),
            "relative_defect_to_H0_MPI1": float(
                abs(a["signed_final_energy_balance_defect"])/abs(a["hamiltonian"][0])
            ),
            "relative_defect_to_H0_MPI2": float(
                abs(b["signed_final_energy_balance_defect"])/abs(b["hamiltonian"][0])
            ),
            "max_mass_relative_drift_MPI1": float(a["relative_mass_drift"]),
            "max_mass_relative_drift_MPI2": float(b["relative_mass_drift"]),
            "number_of_saved_samples_MPI1": len(a["saved_times"]),
            "number_of_saved_samples_MPI2": len(b["saved_times"]),
            "work_final_MPI1": float(a["cumulative_boundary_power"][-1]),
        })

    # Record observed ratios/orders, do NOT force an unsupported
    # second-order claim for a discontinuous original initial height.
    errs = [abs(v["signed_balance_defect_MPI1"]) for v in all_records]
    ratios, observed = [], []
    for previous, current in zip(errs, errs[1:]):
        if min(previous, current) <= 1e-12:
            ratios.append(None)
            observed.append(None)
        else:
            ratio = previous/current
            ratios.append(float(ratio))
            observed.append(float(math.log2(ratio)))
    decreasing = bool(errs[0] > errs[1] > errs[2])
    classification = (
        "DECREASING_OVER_THREE_STEPS"
        if decreasing else "NO_MONOTONE_DECREASE_ON_TESTED_STEPS"
    )
    summary = {
        "schema": "scrimp-v130-issue28-shallow-time-refinement-v1",
        "status": "PASS_NATIVE_TIME_REFINEMENT_MEASUREMENTS",
        "original_physical_model_and_prescribed_initial_data_preserved": True,
        "mesh_h": 0.1,
        "fixed_final_time": T,
        "native_mpi_ranks": [1, 2],
        "native_true_accepted_step_counts": [4, 8, 16],
        "runs": all_records,
        "observed_abs_energy_defect_ratios_MPI1": ratios,
        "observed_energy_defect_orders_MPI1": observed,
        "observed_three_point_abs_defect_decreasing": decreasing,
        "trend_classification": classification,
        "max_mpi_relative_hamiltonian_difference": max_h_mpi,
        "max_mpi_relative_mass_difference": max_mass_mpi,
        "exact_BDF2_power_balance_claimed": False,
        "energy_defect_convergence_proven": False,
        "spatial_mesh_convergence_qualified": False,
        "rotating_viscous_initial_Y0_qualified": False,
    }
    output = folder/"shallow_energy_time_refinement_comparison.json"
    output.write_text(json.dumps(summary, indent=2)+"\n", encoding="utf-8")
    print("PASS: native original SWE time-refinement CAMPAIGN measured "
          "at fixed T, MPI1/2, no hidden initial changes.")
    print("Steps | dt | signed defect MPI1 | relative defect/H0 | saved samples")
    for v in all_records:
        print(
            f"{v['requested_steps']:>5d} | {v['dt']:.9g} | "
            f"{v['signed_balance_defect_MPI1']:.12g} | "
            f"{v['relative_defect_to_H0_MPI1']:.9g} | "
            f"{v['number_of_saved_samples_MPI1']}"
        )
    print(f"Observed adjacent defect ratios: {ratios}; orders: {observed}")
    print(f"Trend: {classification} (not an assumed order of BDF2).")
    print("MPI energy parity max relative=", f"{max_h_mpi:.9g}")
    print("BDF2 exact power conservation and rotating Y0 remain UNQUALIFIED.")
    print("Saved evidence:", folder)
    return summary


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--evidence", required=True, type=Path)
    args = ap.parse_args()
    check(args.evidence)


if __name__ == "__main__":
    main()
