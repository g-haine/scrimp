# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Compare complete 1-/2-rank wave DAE initial-state and CN trajectories.

Reports MPI numeric parity without equating it with physical fidelity.
The historical pseudo-TS may modify q(0) even when MPI trajectories
are indistinguishable; keep that failure visible in the final summary.
No PETSc/GetFEM dependency beyond the image's Python/NumPy is needed.
"""

import argparse
import json
from pathlib import Path

import numpy as np


def _scaled_l2(a, b):
    return float(
        np.linalg.norm(np.asarray(a) - np.asarray(b))
        / max(1.0, float(np.linalg.norm(a)))
    )


def _load(path, ranks):
    stem = f"wave_initial_{ranks}r"
    with (path / f"{stem}.json").open(encoding="utf-8") as handle:
        report = json.load(handle)
    with np.load(path / f"{stem}.npz", allow_pickle=False) as archive:
        states = {
            name: np.array(archive[name])
            for name in ("before", "after", "times", "trajectory")
        }
    if report["mpi_ranks"] != ranks:
        raise ValueError(f"Unexpected MPI rank count for {stem}")
    if len(states["times"]) != report["cn_steps"] + 1:
        raise ValueError(f"Incomplete saved time series for {stem}")
    if states["trajectory"].shape != (
        len(states["times"]), report["dofs"]
    ):
        raise ValueError(f"Incomplete global state at {stem}")
    return report, states


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--rtol", type=float, default=5e-8)
    parser.add_argument("--atol", type=float, default=1e-9)
    parser.add_argument(
        "--require-initial-fidelity", action="store_true",
        help="Fail if pseudo-initialization changes prescribed q,p "
             "(diagnostic mode allows this known open physics issue)",
    )
    args = parser.parse_args()
    ref, a = _load(args.directory, 1)
    parallel, b = _load(args.directory, 2)
    config = (
        "source_image_revision", "mesh_h", "dt", "final_time",
        "cn_steps", "dofs", "variable_intervals", "dynamic_names",
        "algebraic_names", "fidelity_atol", "fidelity_rtol",
    )
    for key in config:
        if ref[key] != parallel[key]:
            raise RuntimeError(f"Different model or source: {key}")
    if not ref["source_image_revision"]:
        raise RuntimeError("Image provenance label is missing")
    if not np.array_equal(a["times"], b["times"]):
        raise RuntimeError("MPI time grids do not match")
    tolerances = args.atol + args.rtol
    discrepancies = {}
    for name in ("before", "after", "trajectory"):
        scaled = _scaled_l2(a[name], b[name])
        discrepancies[name] = scaled
        if scaled > tolerances:
            raise RuntimeError(
                f"{name}: full-state MPI scaled L2={scaled:.9g} "
                f"> tolerance={tolerances:.9g}"
            )
    per_step = [
        _scaled_l2(x, y)
        for x, y in zip(a["trajectory"], b["trajectory"])
    ]
    if max(per_step) > tolerances:
        raise RuntimeError(
            f"Full-state MPI trajectory step differs: {max(per_step):.9g}"
        )
    algebraic_deviation = {}
    for stage in ("before", "after"):
        key = f"algebraic_residual_{stage}"
        va = ref[key]["algebraic_residual_l2"]
        vb = parallel[key]["algebraic_residual_l2"]
        diff = abs(va - vb)
        algebraic_deviation[stage] = diff
        if diff > args.atol + args.rtol * max(1.0, abs(va)):
            raise RuntimeError(
                f"Algebraic residual at {stage} disagrees in MPI"
            )
    if any(
        parallel["dynamic_initial_changes"][name][
            "preserved_within_tolerance"
        ] != ref["dynamic_initial_changes"][name][
            "preserved_within_tolerance"
        ]
        for name in ref["dynamic_names"]
    ):
        raise RuntimeError("Physical fidelity verdict differs across MPI")
    if ref["physical_initial_value_fidelity"] != parallel[
        "physical_initial_value_fidelity"
    ]:
        raise RuntimeError("MPI initial-data fidelity verdict mismatch")

    summary = {
        "schema": "scrimp-v130-wave-dae-initial-comparison-v1",
        "source_image_revision": ref["source_image_revision"],
        "mesh_h": ref["mesh_h"],
        "dt": ref["dt"],
        "final_time": ref["final_time"],
        "dofs": ref["dofs"],
        "steps": ref["cn_steps"],
        "mpi_ranks_checked": [1, 2],
        "full_state_scaled_l2": discrepancies,
        "maximum_scaled_l2_over_individual_times": max(per_step),
        "mpi_algebraic_residual_absolute_differences": algebraic_deviation,
        "dynamic_initial_change_l2": {
            name: ref["dynamic_initial_changes"][name][
                "absolute_change_l2"
            ] for name in ref["dynamic_names"]
        },
        "algebraic_residual_initial": ref[
            "algebraic_residual_before"
        ]["algebraic_residual_l2"],
        "algebraic_residual_accepted_pseudo_state": ref[
            "algebraic_residual_after"
        ]["algebraic_residual_l2"],
        "physical_initial_value_fidelity": ref[
            "physical_initial_value_fidelity"
        ],
        "mpi_full_state_parity": "PASS",
        "interpretation": (
            "MPI parity and algebraic-constraint residuals cannot "
            "prove preservation of prescribed differential initial data"
        ),
    }
    (args.directory / "wave_initial_dae_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print("SCRIMP wave DAE initialization and full CN trajectory")
    print(
        f"Image {ref['source_image_revision']}, h={ref['mesh_h']}, "
        f"dt={ref['dt']}, steps={ref['cn_steps']}, DoFs={ref['dofs']}"
    )
    print("| Property | 1 rank | 2 ranks |")
    print("| --- | ---: | ---: |")
    for field in ref["dynamic_names"]:
        a_d = ref["dynamic_initial_changes"][field]
        b_d = parallel["dynamic_initial_changes"][field]
        print(
            f"| Initial \\|\\|{field} - {field}_prescribed\\|\\|_2 | "
            f"{a_d['absolute_change_l2']:.9g} | "
            f"{b_d['absolute_change_l2']:.9g} |"
        )
    for key, label in (("before", "Initial algebraic residual"),
                       ("after", "Post-pseudo algebraic residual")):
        a_r = ref[f"algebraic_residual_{key}"]["algebraic_residual_l2"]
        b_r = parallel[f"algebraic_residual_{key}"]["algebraic_residual_l2"]
        print(f"| {label} L2 | {a_r:.9g} | {b_r:.9g} |")
    print(
        f"PASS: full-state MPI parity at every saved CN step, "
        f"max scaled L2={max(per_step):.9g}"
    )
    print(
        "PHYSICAL INITIAL-DATA FIDELITY: "
        f"{summary['physical_initial_value_fidelity']}"
    )
    if (args.require_initial_fidelity
            and summary["physical_initial_value_fidelity"] != "PASS"):
        raise RuntimeError("Prescribed q/p initial data were changed by pseudo-TS")


if __name__ == "__main__":
    main()
