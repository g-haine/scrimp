# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
"""Aggregate repeated 2D acoustic strong-scaling runs.

Example:
    python sandbox/analyze_acoustic_scaling.py \
      sandbox/outputs/scaling/wave_*ranks_*.json

Each report must describe the SAME spatial mesh, timestep and t_f.
The analysis fails on significant physical-state discrepancies and
uses median wall time to avoid overinterpreting a single noisy sample.
"""
import argparse
from collections import defaultdict
import json
from pathlib import Path
import statistics

import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("files", nargs="+", type=Path)
    parser.add_argument("--rtol", type=float, default=1e-5,
                        help="Allowed relative difference in field norms")
    parser.add_argument("--atol", type=float, default=1e-9)
    args = parser.parse_args()

    grouped = defaultdict(list)
    reference = None
    for path in args.files:
        run = json.loads(path.read_text(encoding="utf-8"))
        rank_count = int(run["mpi_size"])
        if run["ts_converged_reason"] <= 0 or run["ts_snes_failures"] != 0:
            raise SystemExit(f"{path}: solver did not converge cleanly")
        signature = (
            int(run["degrees_of_freedom"]),
            float(run["rectangle_mesh_h"]),
            float(run["requested_final_time"]),
            float(run["last_saved_time"]),
            float(run["actual_final_time"]),
            int(run["ts_steps"]),
            str(run["ts_type"]),
            str(run["ts_adapt_type"]),
        )
        if reference is None:
            reference = signature
        elif signature != reference:
            raise SystemExit(f"{path}: strong-scaling configuration mismatch {signature} != {reference}")
        grouped[rank_count].append((path, run))

    if 1 not in grouped:
        raise SystemExit("A one-rank baseline is required")
    baseline = grouped[1][0][1]
    field_reference = baseline["final_variable_norms"]
    for rank_count, runs in grouped.items():
        for path, run in runs:
            for name, expected in field_reference.items():
                if name not in run["final_variable_norms"]:
                    raise SystemExit(f"{path}: missing field {name}")
                value = float(run["final_variable_norms"][name])
                if not np.isclose(value, expected, atol=args.atol, rtol=args.rtol):
                    raise SystemExit(
                        f"{path}: MPI physical mismatch in {name}: "
                        f"{value} vs one-rank {expected}"
                    )

    baseline_solve = statistics.median(
        run["wall_solve_seconds_max"] for _, run in grouped[1]
    )
    print(
        f"Validated fixed-mesh acoustic strong scaling: "
        f"DOFs={reference[0]}, mesh={reference[1]}, "
        f"t_f={reference[2]}, integrator={reference[6]}"
    )
    print("| MPI ranks | repetitions | model build, s | solve, s | speedup | efficiency | max rank RSS, MiB |")
    print("|---:|---:|---:|---:|---:|---:|---:|")
    for p, runs in sorted(grouped.items()):
        build = statistics.median(
            run["wall_build_seconds_max"] for _, run in runs
        )
        solve = statistics.median(
            run["wall_solve_seconds_max"] for _, run in runs
        )
        peak = max(run["peak_rank_rss_kib"] for _, run in runs) / 1024.0
        speedup = baseline_solve / solve
        efficiency = speedup / p
        print(
            f"| {p} | {len(runs)} | {build:.4g} | {solve:.4g} | "
            f"{speedup:.3f} | {efficiency:.3f} | {peak:.1f} |"
        )
    print(
        "Speedup compares complete DPHS.solve() calls (assembly + DAE "
        "initialization + PETSc TS), NOT just MUMPS kernel throughput. "
        "Times are per-run MPI-max wall durations, aggregated by median; "
        "RSS is peak PER RANK, not sum."
    )


if __name__ == "__main__":
    main()
