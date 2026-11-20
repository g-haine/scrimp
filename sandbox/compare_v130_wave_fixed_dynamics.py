# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Compare native fixed-dynamics algebraic projection probes at MPI 1 and 2.

Separate MPI numeric parity, actual algebraic compatibility and
ill-conditioned/incomplete iterative solve outcomes. A numeric
least-squares failure is not a theorem of DAE incompatibility.
"""

import argparse
import json
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--atol", type=float, default=1e-9)
    parser.add_argument("--rtol", type=float, default=1e-7)
    parser.add_argument(
        "--scenario",
        choices=(
            "published", "matched_left_control", "polynomial_compatible",
        ),
        default="published",
    )
    args = parser.parse_args()
    reports = {}
    suffix = (
        "" if args.scenario == "published"
        else f"_{args.scenario}"
    )
    for ranks in (1, 2):
        path = args.directory / f"wave_algebraic_only_{ranks}r{suffix}.json"
        obj = json.loads(path.read_text(encoding="utf-8"))
        if obj["mpi_ranks"] != ranks:
            raise RuntimeError(f"Invalid rank count in {path}")
        reports[ranks] = obj
    a, b = reports[1], reports[2]
    for key in (
        "schema", "image_source_revision", "global_dofs",
        "algebraic_dofs", "algebraic_intervals",
        "mesh_h", "cn_steps", "time_step",
        "scenario", "prescribed_dynamic_norm_l2",
    ):
        if a[key] != b[key]:
            raise RuntimeError(
                f"Mismatching scientific problem at MPI 1/2: {key}"
            )
    if a["scenario"] != args.scenario:
        raise RuntimeError(
            f"Wrong scenario in reports: {a['scenario']} != {args.scenario}"
        )
    if not a["image_source_revision"]:
        raise RuntimeError("Source image provenance missing")

    for ranks, r in reports.items():
        if not r["linearization_validated"]:
            raise RuntimeError(
                f"MPI={ranks}: extracted algebraic Jacobian does "
                "not predict actual IFunction; feasibility inference invalid"
            )
        if any(v != 0.0 for v in r["dynamic_data_change_l2"].values()):
            raise RuntimeError("Fixed q/p coefficients were modified")
        if r["main_ts_converged_reason"] <= 0 or r["cn_steps"] != 5:
            raise RuntimeError("Original main CN TS failed")
    if (a["fixed_dynamics_actual_algebraic_feasibility"]
            != b["fixed_dynamics_actual_algebraic_feasibility"]):
        raise RuntimeError("1/2 rank algebraic feasibility verdict differs")

    compare_keys = (
        "independent_native_residual_after_project_l2",
        "native_residual_after_pseudo_l2",
        "jacobian_prediction_mismatch_l2",
    )
    for key in compare_keys:
        first, second = float(a[key]), float(b[key])
        if abs(second - first) > args.atol + args.rtol * max(
            abs(first), abs(second)
        ):
            raise RuntimeError(
                f"1/2-rank discrepancy in {key}: {first:.9g}, "
                f"{second:.9g}"
            )

    print(
        f"DAE algebraic-only feasibility probe ({args.scenario}): "
        f"{a['global_dofs']} DOFs, {a['algebraic_dofs']} algebraic, "
        f"source {a['image_source_revision']}"
    )
    print(
        "| MPI ranks | Initial algebraic residual | "
        "Only-algebraic projected residual | "
        "Pseudo-TS residual | LSMR condition estimate |"
    )
    print("|---:|---:|---:|---:|---:|")
    for ranks, r in sorted(reports.items()):
        solver = r["linear_least_squares"]
        print(
            f"| {ranks} | "
            f"{solver['raw_constraint_residual_l2']:.8e} | "
            f"{r['independent_native_residual_after_project_l2']:.8e} | "
            f"{r['native_residual_after_pseudo_l2']:.8e} | "
            f"{solver['lsmr_condition_estimate']:.4g} |"
        )
    print(
        "PASS: same algebraic-feasibility outcome at MPI 1 and 2, "
        "dynamic q,p preserved exactly, validated linearization."
    )
    print("| Algebraic field | Initial residual | Projected residual | Pseudo residual |")
    print("| --- | ---: | ---: | ---: |")
    for name, values in a["algebraic_residual_by_unknown"].items():
        b_values = b["algebraic_residual_by_unknown"][name]
        for stage in ("before", "projected", "pseudo"):
            x = values["residual_l2"][stage]
            y = b_values["residual_l2"][stage]
            if abs(x-y) > args.atol + args.rtol * max(abs(x), abs(y)):
                raise RuntimeError(
                    f"Algebraic field {name} stage {stage} mismatch in MPI"
                )
        norm = values["residual_l2"]
        print(
            f"| {name} | {norm['before']:.9g} | "
            f"{norm['projected']:.9g} | {norm['pseudo']:.9g} |"
        )
    for ranks, record in sorted(reports.items()):
        witness = record["linear_least_squares"]
        print(
            f"MPI={ranks} LSMR stop={witness['lsmr_istop']} "
            f"iter={witness['lsmr_iterations']} "
            f"normal_residual={witness['normal_equation_residual_l2_recomputed']:.8g}; "
            f"dual w^H r={witness['residual_left_witness_pairing_with_initial']:.8g}, "
            f"||K_aa^H w||={witness['residual_left_witness_adjoint_l2']:.8g}; "
            f"bounded lower bound={witness['dual_lower_bound_for_corrections_within_radius_l2']:.8g} "
            f"for ||delta_a||<={witness['bounded_correction_radius_l2']:.8g}"
        )
    for ranks, record in sorted(reports.items()):
        changes = record["pseudo_dynamic_initial_data_drift_l2"]
        p0 = record["prescribed_momentum_expression"]
        print(
            f"Physical pseudo-TS differential drift at MPI={ranks}: "
            f"q={changes['q']:.9g}, p={changes['p']:.9g}; p0={p0}"
        )
    print(
        "FIXED-DYNAMICS ALGEBRAIC FEASIBILITY: "
        f"{a['fixed_dynamics_actual_algebraic_feasibility']}"
    )
    if a["fixed_dynamics_actual_algebraic_feasibility"] != "PASS":
        print(
            "WARNING: do not call this a proven incompatibility. "
            "Investigate algebraic rank, conditioning, boundary trace "
            "and iterative accuracy before modifying production."
        )


if __name__ == "__main__":
    main()
