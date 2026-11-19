# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Strict parity gates for a SHORT nonlinear Saint-Venant output test.

Check identical physical configuration, successful fixed-CN integration,
all states, Hamiltonian, boundary supplied power, viscous dissipation,
volume conservation and complete power/volume balance trajectories.
These are cross-backend COMPARISONS, not assertions that the physical
energy/mass defects themselves are acceptably small.
"""
import argparse
import json
from pathlib import Path

import numpy as np

REL_TOL = 1e-9
ABS_TOL = 1e-9


def _error_metrics(actual, expected):
    delta = np.asarray(actual) - np.asarray(expected)
    peak = float(np.max(np.abs(delta))) if delta.size else 0.0
    norm = float(np.linalg.norm(delta))
    l2_scale = max(
        float(np.linalg.norm(actual)),
        float(np.linalg.norm(expected)),
        float(np.finfo(float).tiny),
    )
    linf_scale = max(
        float(np.max(np.abs(actual))),
        float(np.max(np.abs(expected))),
        float(np.finfo(float).tiny),
    ) if delta.size else 1.0
    return {
        "max_abs": peak,
        "relative_l2": norm / l2_scale,
        "relative_linf": peak / linf_scale,
    }


def compare(output):
    output = Path(output)
    reports = {
        mode: json.loads(
            (output / f"{mode}_swe.json").read_text(encoding="utf-8")
        )
        for mode in ("memory", "async_binary")
    }
    ref, streamed = reports["memory"], reports["async_binary"]
    for key in (
        "source_commit", "ranks", "dofs", "snapshots", "mesh", "dt",
        "final_time", "spinup_duration", "gamma", "ts_type",
        "ts_adapt_type", "mumps_reset_policy",
    ):
        if ref[key] != streamed[key]:
            raise ValueError(f"Mismatched physical configuration {key}")

    for mode, report in reports.items():
        conv = report["convergence"]
        if (conv["reason"] <= 0 or
                conv["steps"] != round(ref["final_time"] / ref["dt"]) or
                conv["rejections"] or conv["snes_failures"] or
                abs(conv["t_final"] - ref["final_time"]) > 1e-12 or
                report["snapshots"] != conv["steps"] + 1 or
                report["minimum_water_height"] <= 0):
            raise ValueError(f"{mode}: nonlinear convergence/sampling failed")
    if ref["dofs"] <= 0 or ref["ranks"] <= 0:
        raise ValueError("Invalid problem dimensions")
    manifest_path = output / "binary" / "trajectory_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if (manifest["global_size"] != ref["dofs"] or
            manifest["number_of_states"] != ref["snapshots"] or
            manifest["ranks"] != ref["ranks"] or
            len(manifest["partitions"]) != ref["ranks"]):
        raise ValueError("Invalid distributed binary manifest")

    metrics = {}
    with np.load(
        output / "memory_swe.npz", allow_pickle=False
    ) as left, np.load(
        output / "async_binary_swe.npz", allow_pickle=False
    ) as right:
        if set(left.files) != set(right.files):
            raise ValueError("Different physical balance diagnostics")
        required = {
            "states", "t", "kinetic", "potential", "H",
            "d_grad", "d_div", "p_normal", "p_viscous",
            "supplied", "dissipated", "balance",
            "relative_balance", "relative_balance_exchange",
            "volume", "normal_volume_outflow",
            "cumulative_outflow", "volume_balance",
            "relative_volume_balance", "lyapunov",
        }
        if set(left.files) != required:
            raise ValueError("Incomplete power/volume balance diagnostics")
        for key in sorted(required):
            old = np.asarray(left[key])
            new = np.asarray(right[key])
            if (old.shape != new.shape or
                    not np.all(np.isfinite(old)) or
                    not np.all(np.isfinite(new))):
                raise ValueError(f"{key}: nonfinite or shape mismatch")
            if key == "t":
                if not np.array_equal(old, new):
                    raise ValueError("Different saved time grid")
            else:
                if not np.allclose(
                    old, new, rtol=REL_TOL, atol=ABS_TOL
                ):
                    m = _error_metrics(new, old)
                    raise ValueError(
                        f"{key}: output backend changed numerical result "
                        + json.dumps(m)
                    )
            metrics[key] = _error_metrics(new, old)

    summary = {
        "source_commit": ref["source_commit"],
        "dofs": ref["dofs"],
        "mpi_ranks": ref["ranks"],
        "snapshots": ref["snapshots"],
        "steps": ref["convergence"]["steps"],
        "spinup_duration": ref["spinup_duration"],
        "numerical_tolerances": {
            "relative_componentwise": REL_TOL,
            "absolute_componentwise": ABS_TOL,
        },
        "max_absolute_difference_by_quantity": {
            key: item["max_abs"] for key, item in metrics.items()
        },
        "errors": metrics,
        "physical_balance": {
            "memory_energy_defect_final":
                ref["energy_balance_defect_final"],
            "streamed_energy_defect_final":
                streamed["energy_balance_defect_final"],
            "memory_volume_defect_final":
                ref["volume_balance_defect_final"],
            "streamed_volume_defect_final":
                streamed["volume_balance_defect_final"],
        },
        "notes": [
            "This checks backend parity for a short nonlinear spin-up.",
            "Passing does NOT imply published-model energy/volume "
            "balances are physically accurate, only that both storage "
            "backends preserve the computed balances.",
            "No long simulation or solver-history changes were made.",
        ],
    }
    (output / "comparison_swe.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(
        "PASS: nonlinear Saint-Venant stream parity | "
        f"{ref['dofs']} DDL, {ref['ranks']} MPI, "
        f"{summary['steps']} CN steps, {ref['snapshots']} states"
    )
    for key in ("states", "H", "supplied", "dissipated",
                "balance", "volume", "volume_balance"):
        err = metrics[key]
        print(
            f"{key}: max_abs={err['max_abs']:.4e}, "
            f"rel_L2={err['relative_l2']:.4e}"
        )
    print(
        "Both energy and volume defect histories validated for "
        "numerical parity (not absolute physics accuracy)."
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", type=Path, required=True)
    args = parser.parse_args()
    compare(args.directory)


if __name__ == "__main__":
    main()
