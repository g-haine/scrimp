# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Compare full CN physics and MPI parity of two DAE initialization paths.

The model, prescribed q0/p0, drive and FE spaces must be identical.
MPI 1 vs 2 should reproduce each MODE's entire trajectory, energy,
port work and midpoint-power series. Pseudo vs projected physical
trajectories are NOT required to agree: different accepted t0
states represent different Cauchy problems. The difference itself
is reported, with Hamiltonian/port-work checks on BOTH.
"""

import argparse
import json
from pathlib import Path

import numpy as np


MODES = ("pseudo", "algebraic")
RANKS = (1, 2)


def _norm_scaled(a, b):
    a, b = np.asarray(a), np.asarray(b)
    return float(np.linalg.norm(a-b)/max(1.0, np.linalg.norm(a)))


def _read(directory, mode, ranks):
    stem = f"wave_{mode}_{ranks}r"
    report = json.loads(
        (directory / f"{stem}.json").read_text(encoding="utf-8")
    )
    with np.load(directory / f"{stem}.npz", allow_pickle=False) as src:
        states = {
            k: np.asarray(src[k])
            for k in (
                "times", "trajectory", "prescribed_initial",
                "pseudo_accepted", "main_initializer",
            )
        }
    if (report["initializer"] != mode or report["mpi_ranks"] != ranks):
        raise RuntimeError(f"Incorrect initializer/ranks in {stem}")
    steps = report["steps"]
    dofs = report["dofs"]
    if (
        states["times"].shape != (steps+1,)
        or states["trajectory"].shape != (steps+1, dofs)
        or states["main_initializer"].shape != (dofs,)
    ):
        raise RuntimeError(f"Incomplete/invalid physical trajectory: {stem}")
    return report, states


def compare(directory, *, atol=1e-9, rtol=5e-8, strict=True):
    records = {
        mode: {
            ranks: _read(directory, mode, ranks)
            for ranks in RANKS
        } for mode in MODES
    }
    reference = records["pseudo"][1][0]
    common = (
        "schema", "image_source_revision", "mesh_h", "dt",
        "steps", "final_time", "drive_left",
        "initial_p0_expression", "dofs",
    )
    summary = {
        "schema": "scrimp-v130-initializer-scientific-comparison-v1",
        "image_source_revision": reference["image_source_revision"],
        "mpi_parity": {},
        "physical_initial_state_drift": {},
        "power_balance": {},
    }
    if not reference["image_source_revision"]:
        raise RuntimeError("Missing image source provenance")
    for mode in MODES:
        a, az = records[mode][1]
        b, bz = records[mode][2]
        for ranks, report in ((1, a), (2, b)):
            for key in common:
                if report[key] != reference[key]:
                    raise RuntimeError(
                        f"Changed physics between modes/MPI: "
                        f"mode={mode}, ranks={ranks}, {key}"
                    )
            if (report["main_ts_converged_reason"] <= 0
                    or report["main_ts_rejected_steps"] != 0
                    or report["main_ts_snes_failures"] != 0):
                raise RuntimeError(
                    f"Physical CN solve not converged for {mode}/{ranks}"
                )
            if not report["nonvacuous_port_power"]:
                raise RuntimeError(
                    f"Nontrivial boundary port power missing for {mode}/{ranks}"
                )
            power = report["power"]
            if not (np.isfinite(power["absolute_port_work"])
                    and power["absolute_port_work"] > 1e-9):
                raise RuntimeError(
                    f"Invalid or vacuous power measurement for {mode}/{ranks}"
                )
            summary["power_balance"][f"{mode}_{ranks}r"] = {
                "pass": bool(report["physical_power_balance_within_tolerance"]),
                "maximum_absolute_balance_defect": power[
                    "maximum_absolute_balance_defect"
                ],
                "tolerance": report["absolute_balance_tolerance"],
                "integrated_work": power["port_work_integral"][-1],
                "absolute_port_work": power["absolute_port_work"],
                "energy_change": power["energy_change"],
            }
        for name in (
            "times", "prescribed_initial", "pseudo_accepted",
            "main_initializer",
        ):
            if not np.allclose(az[name], bz[name], atol=atol, rtol=rtol):
                raise RuntimeError(
                    f"MPI 1/2 differs in {name} for mode {mode}"
                )
        per_step = [
            _norm_scaled(x, y)
            for x, y in zip(az["trajectory"], bz["trajectory"])
        ]
        worst = max(per_step)
        if worst > atol + rtol:
            raise RuntimeError(
                f"Full-CN trajectory MPI parity fails {mode}: {worst}"
            )
        for field in (
            "hamiltonian", "port_work_integral",
            "step_balance_defects",
        ):
            if not np.allclose(
                a["power"][field], b["power"][field],
                atol=atol, rtol=rtol
            ):
                raise RuntimeError(
                    f"Physical {field} MPI parity fails for {mode}"
                )
        for name, power_series in a["power"][
            "boundary_midpoint_powers"
        ].items():
            if name not in b["power"]["boundary_midpoint_powers"]:
                raise RuntimeError("Different boundary-port sets")
            if not np.allclose(
                power_series,
                b["power"]["boundary_midpoint_powers"][name],
                atol=atol, rtol=rtol,
            ):
                raise RuntimeError(
                    f"Boundary-port {name} MPI parity fails {mode}"
                )
        summary["mpi_parity"][mode] = {
            "full_trajectory_max_scaled_l2": worst,
            "all_saved_steps_pass": True,
            "energy_and_port_power_pass": True,
        }

    # Same prescribed dynamic state, independent of which initializer
    # is allowed to publish the physical t0 state.
    for ranks in RANKS:
        p, pz = records["pseudo"][ranks]
        a, az = records["algebraic"][ranks]
        if not np.array_equal(
            pz["prescribed_initial"], az["prescribed_initial"]
        ):
            raise RuntimeError("Initializer modes changed prescribed data")
        if not np.array_equal(pz["pseudo_accepted"], az["pseudo_accepted"]):
            raise RuntimeError(
                "Pseudo-TS accepted state changed when introducing "
                "diagnostic algebraic-only initialization"
            )
        if not np.array_equal(
            az["main_initializer"], az["prescribed_initial"]
        ):
            # Expected difference in algebraic variables; check dynamic
            # blocks via explicit drift metadata below instead.
            pass
        for name, drift in a["main_initial_dynamic_drift_l2"].items():
            if drift > 1e-12:
                raise RuntimeError(
                    f"Projected initializer failed to preserve {name}: {drift}"
                )
        q_drift = p["main_initial_dynamic_drift_l2"]["q"]
        p_drift = p["main_initial_dynamic_drift_l2"]["p"]
        if q_drift < 1e-8 or p_drift < 1e-9:
            raise RuntimeError(
                "Pseudo-TS drift unexpectedly vanished; "
                "the sensitivity A/B experiment was not effective"
            )
        if a["algebraic_projection"] is None:
            raise RuntimeError("Missing independent native algebraic projection")
        tol = a["algebraic_projection"]["native_algebraic_tolerance"]
        residual = a["algebraic_projection"][
            "algebraic_residual_after_project_l2"
        ]
        if residual > tol:
            raise RuntimeError(
                "Projected state is not algebraically consistent"
            )
        step_differences = [
            _norm_scaled(x,y) for x,y in zip(
                pz["trajectory"], az["trajectory"]
            )
        ]
        summary["physical_initial_state_drift"][f"{ranks}r"] = {
            "pseudo_q0_drift_l2": q_drift,
            "pseudo_p0_drift_l2": p_drift,
            "algebraic_q0_drift_l2": a[
                "main_initial_dynamic_drift_l2"
            ]["q"],
            "algebraic_p0_drift_l2": a[
                "main_initial_dynamic_drift_l2"
            ]["p"],
            "full_trajectory_mode_difference_scaled_l2_by_time": (
                step_differences
            ),
            "hamiltonian_t0_difference": (
                a["power"]["hamiltonian"][0]
                - p["power"]["hamiltonian"][0]
            ),
        }

    checks = [
        value["pass"]
        for value in summary["power_balance"].values()
    ]
    summary["both_initializers_power_balance_pass"] = bool(all(checks))
    summary["full_trajectories_mpi_1_2_parity"] = "PASS"
    summary["diagnostic_projected_dynamic_initial_fidelity"] = "PASS"
    summary["pseudotime_dynamic_initial_fidelity"] = "FAIL"
    (directory / "wave_initializer_physics_comparison.json").write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )

    print(
        "SCRIMP v1.3.0 wave CN initializers: same compatible p0, "
        "1/2 ranks, full-state trajectories and signed port powers"
    )
    print(
        "| Initializer | MPI | initial q drift | initial p drift | "
        "Max |H-H0+W| | tolerance | absolute boundary work |"
    )
    print("| --- | ---: | ---: | ---: | ---: | ---: | ---: |")
    for mode in MODES:
        for ranks in RANKS:
            report, _states = records[mode][ranks]
            power = summary["power_balance"][f"{mode}_{ranks}r"]
            d = report["main_initial_dynamic_drift_l2"]
            print(
                f"| {mode} | {ranks} | {d['q']:.8g} | "
                f"{d['p']:.8g} | "
                f"{power['maximum_absolute_balance_defect']:.8g} | "
                f"{power['tolerance']:.8g} | "
                f"{power['absolute_port_work']:.8g} |"
            )
    for mode in MODES:
        print(
            f"PASS: complete CN MPI parity in mode={mode}, "
            "max scaled state discrepancy="
            f"{summary['mpi_parity'][mode]['full_trajectory_max_scaled_l2']:.9g}"
        )
    print(
        "PHYSICAL POWER BALANCE, BOTH INITIALIZERS: "
        + ("PASS" if all(checks) else "FAIL")
    )
    print(
        "PHYSICAL INITIAL-DATA FIDELITY: "
        "projected=PASS; production pseudo=FAIL (for this wave)"
    )
    if strict and not all(checks):
        raise RuntimeError(
            "A non-vacuous, signed Hamiltonian/port-power balance "
            "failed in at least one actual wave simulation"
        )
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument(
        "--diagnose-only", action="store_true",
        help="Write comparison and explicit scientific balance FAIL "
             "without terminating before all evidence is inspected",
    )
    args = parser.parse_args()
    compare(args.directory, strict=not args.diagnose_only)


if __name__ == "__main__":
    main()
