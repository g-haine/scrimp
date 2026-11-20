# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""A/B contrast: original left velocity vs p0/rho-matched boundary velocity.

Interpret as support/refutation of an *initial boundary compatibility*
hypothesis for one discrete wave model. Different boundary controls
are different mathematical initial-value problems, not MPI parity
or a production solver change. The dynamic initial q/p must agree.
"""

import argparse
import json
from pathlib import Path

import numpy as np


def _load(directory, ranks, scenario):
    suffix = "" if scenario == "published" else f"_{scenario}"
    path = directory / f"wave_algebraic_only_{ranks}r{suffix}.json"
    return json.loads(path.read_text(encoding="utf-8"))


def compare(directory):
    scenarios = ("published", "matched_left_control")
    data = {
        scenario: {n: _load(directory, n, scenario) for n in (1, 2)}
        for scenario in scenarios
    }
    results = {}
    for ranks in (1, 2):
        a = data["published"][ranks]
        b = data["matched_left_control"][ranks]
        for key in ("image_source_revision", "global_dofs",
                    "algebraic_dofs", "algebraic_intervals",
                    "mesh_h", "cn_steps", "time_step"):
            if a[key] != b[key]:
                raise ValueError(f"Different underlying model: {key}")
        for field, original in a["prescribed_dynamic_norm_l2"].items():
            variant = b["prescribed_dynamic_norm_l2"][field]
            if not np.isclose(original, variant, atol=1e-12, rtol=1e-12):
                raise ValueError(
                    f"Changed prescribed differential initial field {field}"
                )
        for scenario, current in (
            ("published", a), ("matched_left_control", b)
        ):
            if any(v != 0.0 for v in current[
                "dynamic_data_change_l2"
            ].values()):
                raise ValueError(f"Diagnostic modified q,p in {scenario}")
            if not current["linearization_validated"]:
                raise ValueError(
                    f"Invalid native IFunction linearization in {scenario}"
                )
        results[str(ranks)] = {
            scenario: {
                "projected_native_residual_l2": float(
                    data[scenario][ranks][
                        "independent_native_residual_after_project_l2"
                    ]
                ),
                "fixed_dynamics_feasibility": data[scenario][ranks][
                    "fixed_dynamics_actual_algebraic_feasibility"
                ],
                "dual_lower_bound_for_bounded_corrections": data[scenario][
                    ranks
                ]["linear_least_squares"][
                    "dual_lower_bound_for_corrections_within_radius_l2"
                ],
            } for scenario in scenarios
        }
    original = data["published"][1]
    alternate = data["matched_left_control"][1]
    supports_boundary_hypothesis = (
        original["fixed_dynamics_actual_algebraic_feasibility"]
        != "PASS"
        and alternate["fixed_dynamics_actual_algebraic_feasibility"]
        == "PASS"
    )
    report = {
        "schema": "scrimp-v130-wave-left-control-compatibility-v1",
        "image_source_revision": original["image_source_revision"],
        "mesh_h": original["mesh_h"],
        "dofs": original["global_dofs"],
        "algebraic_dofs": original["algebraic_dofs"],
        "prescribed_dynamic_norms_unchanged": True,
        "two_rank_consistency": True,
        "original_and_matched_left_control": results,
        "matched_control_removes_algebraic_obstruction": (
            bool(supports_boundary_hypothesis)
        ),
        "interpretation": (
            "An algebraic residual difference under different "
            "U_L(t) isolates a boundary-compatibility hypothesis "
            "for this discrete model. It does not prove a universal "
            "DAE theorem or approve a production initializer."
        ),
    }
    destination = directory / "wave_left_control_compatibility_summary.json"
    destination.write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    print("Wave DAE control-compatibility contrast: same prescribed q0/p0")
    print("| MPI | Published projected residual | Matched-left projected residual | Published | Matched |")
    print("| ---: | ---: | ---: | --- | --- |")
    for ranks in (1, 2):
        p = results[str(ranks)]["published"]
        c = results[str(ranks)]["matched_left_control"]
        print(
            f"| {ranks} | {p['projected_native_residual_l2']:.9g} | "
            f"{c['projected_native_residual_l2']:.9g} | "
            f"{p['fixed_dynamics_feasibility']} | "
            f"{c['fixed_dynamics_feasibility']} |"
        )
    print(
        "BOUNDARY COMPATIBILITY HYPOTHESIS: "
        + ("SUPPORTED (numerically)" if supports_boundary_hypothesis
           else "NOT YET ESTABLISHED")
    )
    print(f"Full report: {destination}")
    return report


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("directory", type=Path)
    args = p.parse_args()
    compare(args.directory)


if __name__ == "__main__":
    main()
