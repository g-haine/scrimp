# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

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


def plot_policy_comparison(paths, prefix, rtol, atol):
    """Plot existing benchmark reports without changing their profiling API.

    Each tagged JSON report contains a series='auto' or 'always' field.
    Physically inconsistent or incomplete MPI campaigns are rejected.
    Median solve wall times, not individual noisy repetitions, are used
    for the speedup and time comparisons.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    groups = defaultdict(list)
    signature_ref = None
    field_ref = None
    for path in paths:
        run = json.loads(path.read_text(encoding="utf-8"))
        series = run.get("series")
        if series not in ("auto", "always"):
            raise SystemExit(
                f"{path}: comparison requires --series auto/always "
                "in each benchmark report"
            )
        if (run["ts_converged_reason"] <= 0 or
                run["ts_snes_failures"] != 0 or
                run["ts_rejected_steps"] != 0):
            raise SystemExit(f"{path}: solver failed the convergence gate")
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
        if signature_ref is None:
            signature_ref = signature
            field_ref = run["final_variable_norms"]
        elif signature != signature_ref:
            raise SystemExit(
                f"{path}: mismatched fixed-mesh benchmark: "
                f"{signature} vs {signature_ref}"
            )
        for name, expected in field_ref.items():
            measured = run["final_variable_norms"].get(name)
            if measured is None or not np.isclose(
                    float(measured), float(expected), rtol=rtol, atol=atol):
                raise SystemExit(
                    f"{path}: physical mismatch in {name}: "
                    f"{measured} vs reference {expected}"
                )
        groups[series, int(run["mpi_size"])].append(
            float(run["wall_solve_seconds_max"])
        )

    ranks_by_series = {
        series: {rank for label, rank in groups if label == series}
        for series in ("auto", "always")
    }
    if ranks_by_series["auto"] != ranks_by_series["always"]:
        raise SystemExit(
            f"Compare requires the same MPI ranks in both series: "
            f"{ranks_by_series}"
        )
    ranks = sorted(ranks_by_series["auto"])
    if len(ranks) < 2 or 1 not in ranks:
        raise SystemExit("Comparison requires rank 1 and >=1 parallel rank")
    medians = {
        series: {
            rank: statistics.median(groups[series, rank])
            for rank in ranks
        }
        for series in ("auto", "always")
    }
    # Missing or unequal repetitions would bias a comparison of medians.
    repetitions = {len(group) for group in groups.values()}
    if len(repetitions) != 1:
        raise SystemExit(
            "Each rank/policy combination needs the same number "
            f"of repetitions, found counts {sorted(repetitions)}"
        )
    n_reps = repetitions.pop()
    speedups = {
        series: [
            medians[series][1] / medians[series][rank]
            for rank in ranks
        ]
        for series in ("auto", "always")
    }
    prefix.parent.mkdir(parents=True, exist_ok=True)
    print(
        f"Validated policy comparison: {signature_ref[0]} DOFs, "
        f"mesh={signature_ref[1]}, t_f={signature_ref[2]}, "
        f"{n_reps} repetitions per point"
    )
    print(
        "| MPI | always, s | auto, s | reduction, % "
        "| S(always) | S(auto) |"
    )
    print("|---:|---:|---:|---:|---:|---:|")
    for idx, rank in enumerate(ranks):
        baseline = medians["always"][rank]
        automatic = medians["auto"][rank]
        saved = 100.0 * (1.0 - automatic / baseline)
        print(
            f"| {rank} | {baseline:.5g} | {automatic:.5g} "
            f"| {saved:+.2f} | {speedups['always'][idx]:.3f} "
            f"| {speedups['auto'][idx]:.3f} |"
        )

    colors = {"auto": "tab:orange", "always": "tab:green"}
    labels = {
        "auto": "Politique automatique (never si linéaire)",
        "always": "Réinitialisation systématique (always)",
    }

    fig, ax = plt.subplots(figsize=(11.0, 6.7))
    ax.plot(ranks, ranks, linestyle="--", linewidth=2,
            color="tab:blue", label="Accélération idéale")
    for series in ("always", "auto"):
        ax.plot(ranks, speedups[series], "o-", linewidth=2.5,
                markersize=7, color=colors[series],
                label=labels[series])
        for rank, value in zip(ranks, speedups[series]):
            ax.annotate(f"{value:.2f}×", (rank, value),
                        xytext=((-20, -18) if series == "always" and rank == 1
                                else (18, 11) if series == "auto" and rank == 1
                                else (0, -16) if series == "always"
                                else (0, 9)),
                        textcoords="offset points", ha="center",
                        fontsize=10)
    ax.set(
        xlabel="Nombre de rangs MPI (P-cores distincts)",
        ylabel="Accélération S(p) = T(1) / T(p)",
        title=(
            "SCRIMP — Scalabilité forte des ondes acoustiques 2D\n"
            f"{signature_ref[0]:,} DDL · "
            f"{signature_ref[5]} pas de Crank–Nicolson · "
            f"{n_reps} répétitions".replace(",", " ")
        ),
    )
    ax.set_xticks(ranks)
    ax.grid(alpha=0.3)
    ax.legend(loc="upper left")
    fig.tight_layout()
    speed_file = Path(str(prefix) + "_speedup.png")
    fig.savefig(speed_file, dpi=160)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(11.0, 6.7))
    for series in ("always", "auto"):
        values = [medians[series][rank] for rank in ranks]
        ax.plot(ranks, values, "o-", linewidth=2.5,
                markersize=7, color=colors[series],
                label=labels[series])
        for rank, value in zip(ranks, values):
            offset = (
                (-20, 14) if series == "always" and rank == 1
                else (20, -16) if series == "auto" and rank == 1
                else (0, 12) if series == "always"
                else (0, -16)
            )
            ax.annotate(f"{value:.3f} s", (rank, value),
                        xytext=offset, textcoords="offset points",
                        ha="center", fontsize=10)
    ax.set(
        xlabel="Nombre de rangs MPI (P-cores distincts)",
        ylabel="Temps médian de DPHS.solve() (s)",
        title=(
            "SCRIMP — Coût réel des politiques MUMPS\n"
            f"{signature_ref[0]:,} DDL · {n_reps} répétitions".replace(",", " ")
        ),
    )
    ax.set_xticks(ranks)
    ax.grid(alpha=0.3)
    ax.legend(loc="upper right")
    fig.tight_layout()
    times_file = Path(str(prefix) + "_times.png")
    fig.savefig(times_file, dpi=160)
    plt.close(fig)
    print(f"Saved speedup comparison: {speed_file}")
    print(f"Saved solve-time comparison: {times_file}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("files", nargs="+", type=Path)
    parser.add_argument("--rtol", type=float, default=1e-5,
                        help="Allowed relative difference in field norms")
    parser.add_argument("--atol", type=float, default=1e-9)
    parser.add_argument(
        "--plot-comparison", type=Path, metavar="PREFIX",
        help="Compare JSON series auto/always, validate their physical "
             "states, and save PREFIX_speedup.png and PREFIX_times.png. "
             "Normal single-campaign analysis is unchanged.",
    )
    args = parser.parse_args()
    if args.plot_comparison is not None:
        plot_policy_comparison(
            args.files, args.plot_comparison, args.rtol, args.atol
        )
        return

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

    profiles = [
        run.get("profile")
        for runs in grouped.values()
        for _, run in runs
    ]
    if any(p is not None for p in profiles):
        if not all(p is not None and p.get("enabled") for p in profiles):
            raise SystemExit(
                "Mixed profiled and unprofiled runs; repeat with --profile"
            )
        labels = [
            ("Full DPHS.solve", "solve.total"),
            ("GetFEM setup mass", "setup.getfem.mass"),
            ("GetFEM setup stiffness", "setup.getfem.stiffness"),
            ("GetFEM setup RHS", "setup.getfem.rhs"),
            ("Bridge setup mass", "setup.bridge.mass"),
            ("Bridge setup stiffness", "setup.bridge.stiffness"),
            ("Bridge setup RHS", "setup.bridge.rhs"),
            ("Initial DAE consistency (all)", "initial.consistency_total"),
            ("Initial PETSc TS.solve", "initial.petsc.ts_solve"),
            ("Initial IFunction callback", "initial.petsc.callback.ifunction"),
            ("Initial IJacobian callback", "initial.petsc.callback.ijacobian"),
            ("Main PETSc TS.solve", "main.petsc.ts_solve"),
            ("Main IFunction callback", "main.petsc.callback.ifunction"),
            ("Main IJacobian callback", "main.petsc.callback.ijacobian"),
            ("Main PETSc remainder", "main.petsc.non_callback_remainder"),
            ("Main MUMPS PC.reset (opt-in)", "main.petsc.pc_reset"),
            ("Main GetFEM RHS", "main.getfem.rhs"),
            ("Main GetFEM nonlinear mass", "main.getfem.nl_mass"),
            ("Main GetFEM nonlinear stiffness", "main.getfem.nl_stiffness"),
            ("Main bridge RHS", "main.bridge.rhs"),
            ("Main bridge nonlinear mass", "main.bridge.nl_mass"),
            ("Main bridge nonlinear stiffness", "main.bridge.nl_stiffness"),
        ]
        nrows = [row for row in labels if any(
            row[1] in run["profile"]["events"]
            for runs in grouped.values() for _, run in runs
        )]
        print("\nSCRIMP internal profiler (s, medians of MPI rank maxima):")
        print("| Event | " + " | ".join(str(p) for p in sorted(grouped)) + " |")
        print("|---|" + "|".join("---:" for _ in grouped) + "|")
        for title, key in nrows:
            values = []
            for ranks, runs in sorted(grouped.items()):
                durations = [
                    run["profile"]["events"].get(
                        key, {}
                    ).get("seconds_rank_max", 0.)
                    for _, run in runs
                ]
                values.append(f"{statistics.median(durations):.5g}")
            print("| " + title + " | " + " | ".join(values) + " |")
        print(
            "Timings are INCLUSIVE and nested: do not sum callback, "
            "GetFEM and PETSc TS durations. PETSc remainder = TS.solve "
            "minus Python IFunction/IJacobian callback times on each rank; "
            "it includes native factorization/solve, waiting and other "
            "PETSc work, NOT isolated MUMPS time."
        )

        print("\nPETSc callbacks (calls per rank, median):")
        print("| Callback | " + " | ".join(str(p) for p in sorted(grouped)) + " |")
        print("|---|" + "|".join("---:" for _ in grouped) + "|")
        for title, key in [
            ("Main IFunction", "main.petsc.callback.ifunction"),
            ("Main IJacobian", "main.petsc.callback.ijacobian"),
        ]:
            counts = []
            for ranks, runs in sorted(grouped.items()):
                counts.append(str(int(statistics.median([
                    run["profile"]["events"].get(key, {}).get(
                        "calls_rank_max", 0
                    ) for _, run in runs
                ]))))
            print("| " + title + " | " + " | ".join(counts) + " |")



if __name__ == "__main__":
    main()
