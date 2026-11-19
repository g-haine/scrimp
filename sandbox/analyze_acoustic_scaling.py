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



def plot_bridge_comparison(paths, prefix, rtol, atol, selected_policies):
    """Plot honest CSR/COO x MUMPS comparisons on identical problem setups.

    All reports must explicitly identify bridge_method, matrix mode and
    source commit. A complete Cartesian grid of ranks, policies and
    bridge methods with identical repetition counts is mandatory.
    Legacy unlabelled reports are NOT silently assumed to mean CSR.
    """
    import csv
    import math
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    allowed_policies = {"never", "on_change", "always", "auto"}
    if selected_policies is not None:
        selected = selected_policies.split(",")
        if (not selected or len(set(selected)) != len(selected)
                or any(x not in allowed_policies for x in selected)):
            raise SystemExit("Invalid --bridge-policies selection")
    else:
        selected = None

    groups = defaultdict(list)
    reference_signature = None
    reference_fields = None
    commit_ref = None
    for path in paths:
        run = json.loads(path.read_text(encoding="utf-8"))
        policy = run.get("series")
        if selected is not None and policy not in selected:
            continue
        bridge = run.get("bridge_method")
        if policy not in allowed_policies or bridge not in ("csr", "coo"):
            raise SystemExit(
                f"{path}: require explicit series=never/on_change/always "
                "and bridge_method=csr/coo (old unlabelled JSON "
                "cannot be mixed into this comparison)"
            )
        if run.get("getfem_matrix_mode") != "distributed":
            raise SystemExit(f"{path}: expected distributed GetFEM matrix mode")
        source_commit = run.get("benchmark_source_commit")
        if not source_commit:
            raise SystemExit(f"{path}: missing Git source commit metadata")
        if commit_ref is None:
            commit_ref = source_commit
        elif source_commit != commit_ref:
            raise SystemExit(
                f"{path}: compared runs must use same Git commit "
                f"({source_commit} vs {commit_ref})"
            )
        effective = "never" if policy == "auto" else policy
        if run.get("effective_mumps_reset_policy") != effective:
            raise SystemExit(
                f"{path}: MUMPS reset policy is not {effective}"
            )
        if (run.get("ts_converged_reason", 0) <= 0
                or run.get("ts_rejected_steps", 1) != 0
                or run.get("ts_snes_failures", 1) != 0):
            raise SystemExit(f"{path}: failed solver/convergence gate")
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
        if reference_signature is None:
            reference_signature = signature
            reference_fields = run["final_variable_norms"]
        elif signature != reference_signature:
            raise SystemExit(
                f"{path}: different physical benchmark signature: "
                f"{signature} vs {reference_signature}"
            )
        norms = run["final_variable_norms"]
        if set(norms) != set(reference_fields):
            raise SystemExit(f"{path}: different named physical fields")
        for name, expected in reference_fields.items():
            if not np.isclose(
                float(norms[name]), float(expected), rtol=rtol, atol=atol
            ):
                raise SystemExit(
                    f"{path}: physical mismatch in {name}: "
                    f"{norms[name]} vs {expected}"
                )
        wall = float(run["wall_solve_seconds_max"])
        if not math.isfinite(wall) or wall <= 0:
            raise SystemExit(f"{path}: invalid solve time")
        groups[policy, bridge, int(run["mpi_size"])].append(run)

    if not groups:
        raise SystemExit("No matching benchmark reports")

    policies = sorted({p for p, b, n in groups})
    if selected is not None and set(policies) != set(selected):
        raise SystemExit(
            f"Requested policies {selected} not all present: {policies}"
        )
    ranks = sorted({n for p, b, n in groups})
    if 1 not in ranks or len(ranks) < 2:
        raise SystemExit("Require MPI rank 1 and at least another rank")
    expected_grid = {
        (policy, bridge, rank)
        for policy in policies for bridge in ("csr", "coo")
        for rank in ranks
    }
    if set(groups) != expected_grid:
        missing = sorted(expected_grid - set(groups))
        raise SystemExit(
            "Incomplete CSR/COO x MUMPS x MPI grid; "
            f"missing {missing[:20]}"
        )
    counts = {len(group) for group in groups.values()}
    if len(counts) != 1:
        raise SystemExit(
            "Unequal repetition counts across series: "
            f"{sorted(counts)}"
        )
    nreps = counts.pop()
    medians = {
        key: statistics.median(
            float(run["wall_solve_seconds_max"]) for run in records
        )
        for key, records in groups.items()
    }

    def color(policy):
        return {
            "never": "tab:blue",
            "auto": "tab:blue",
            "on_change": "tab:orange",
            "always": "tab:red",
        }[policy]

    def label(policy, bridge):
        policy_label = "auto→never" if policy == "auto" else policy
        return f"{bridge.upper()} · MUMPS {policy_label}"

    title_extra = (
        f"{reference_signature[0]:,} DDL · "
        f"{reference_signature[5]} pas CN · "
        f"{nreps} répétition(s)/point · commit {commit_ref}"
    ).replace(",", " ")

    print(
        f"Validated CSR/COO x MUMPS comparison: {title_extra}; "
        f"ranks={ranks}"
    )
    print(
        "| MUMPS | MPI | CSR (s) | COO (s) | CSR/COO (×) | "
        "S_CSR(p) | S_COO(p) |"
    )
    print("|---|---:|---:|---:|---:|---:|---:|")
    rows = []
    for policy in policies:
        for n in ranks:
            old = medians[policy, "csr", n]
            new = medians[policy, "coo", n]
            s_old = medians[policy, "csr", 1] / old
            s_new = medians[policy, "coo", 1] / new
            print(
                f"| {policy} | {n} | {old:.5f} | {new:.5f} | "
                f"{old/new:.3f} | {s_old:.3f} | {s_new:.3f} |"
            )
            rows.append({
                "source_commit": commit_ref,
                "mumps_policy": policy,
                "mpi_ranks": n,
                "repetitions": nreps,
                "csr_median_s": old,
                "coo_median_s": new,
                "coo_acceleration": old / new,
                "csr_parallel_speedup": s_old,
                "coo_parallel_speedup": s_new,
                "csr_peak_rank_rss_mib": max(
                    x["peak_rank_rss_kib"]
                    for x in groups[policy, "csr", n]
                ) / 1024,
                "coo_peak_rank_rss_mib": max(
                    x["peak_rank_rss_kib"]
                    for x in groups[policy, "coo", n]
                ) / 1024,
            })

    prefix.parent.mkdir(parents=True, exist_ok=True)
    csv_file = Path(str(prefix) + "_summary.csv")
    with csv_file.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    # S(p) is ALWAYS relative to the rank-1 baseline of THAT
    # policy/bridge. This is the actual parallel strong-scaling metric.
    fig, ax = plt.subplots(figsize=(11.5, 7))
    ax.plot(
        ranks, ranks, "k--", linewidth=1.7,
        label="Scalabilité idéale",
    )
    for policy in policies:
        for bridge in ("csr", "coo"):
            values = [
                medians[policy, bridge, 1] / medians[policy, bridge, p]
                for p in ranks
            ]
            ax.plot(
                ranks, values, marker="o" if bridge == "csr" else "s",
                linestyle="--" if bridge == "csr" else "-",
                linewidth=2, markersize=6, color=color(policy),
                label=label(policy, bridge),
            )
    ax.set(
        xlabel="Nombre de rangs MPI (cœurs physiques distincts)",
        ylabel="Accélération parallèle S(p) = T(1) / T(p)",
        title="SCRIMP — Scalabilité forte CSR / COO\n" + title_extra,
    )
    ax.set_xticks(ranks)
    ax.grid(alpha=0.3)
    ax.legend(loc="best")
    fig.tight_layout()
    speed_file = Path(str(prefix) + "_speedup.png")
    fig.savefig(speed_file, dpi=160)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(11.5, 7))
    for policy in policies:
        for bridge in ("csr", "coo"):
            ax.plot(
                ranks, [medians[policy, bridge, p] for p in ranks],
                marker="o" if bridge == "csr" else "s",
                linestyle="--" if bridge == "csr" else "-",
                linewidth=2, markersize=6, color=color(policy),
                label=label(policy, bridge),
            )
    ax.set(
        xlabel="Nombre de rangs MPI (cœurs physiques distincts)",
        ylabel="Temps médian DPHS.solve() (s)",
        title="SCRIMP — Temps réels CSR / COO\n" + title_extra,
    )
    ax.set_xticks(ranks)
    ax.grid(alpha=0.3)
    ax.legend(loc="best")
    fig.tight_layout()
    times_file = Path(str(prefix) + "_times.png")
    fig.savefig(times_file, dpi=160)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(11.5, 7))
    for policy in policies:
        ax.plot(
            ranks,
            [
                medians[policy, "csr", p] / medians[policy, "coo", p]
                for p in ranks
            ],
            "o-", linewidth=2, color=color(policy),
            label=f"MUMPS {policy}",
        )
    ax.axhline(1, linestyle="--", color="black", alpha=0.5)
    ax.set(
        xlabel="Nombre de rangs MPI (cœurs physiques distincts)",
        ylabel="Gain COO / CSR = T_CSR(p) / T_COO(p)",
        title="SCRIMP — Gain du transfert COO à rang égal\n"
        + title_extra,
    )
    ax.set_xticks(ranks)
    ax.grid(alpha=0.3)
    ax.legend(loc="best")
    fig.tight_layout()
    gain_file = Path(str(prefix) + "_coo_gain.png")
    fig.savefig(gain_file, dpi=160)
    plt.close(fig)

    print(f"Saved: {speed_file}, {times_file}, {gain_file}")
    print(f"Summary: {csv_file}")
    print(
        "S(p) uses an independent rank-1 baseline for each CSR/COO and "
        "MUMPS combination. COO gain compares the SAME rank and policy. "
        "All times are median full DPHS.solve() wall times; this does "
        "not isolate MUMPS kernels or claim memory savings."
    )


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
    parser.add_argument(
        "--plot-bridge-comparison", type=Path, metavar="PREFIX",
        help="Compare explicit CSR/COO x MUMPS reports on the same "
             "commit, plot speedup/times/COO gain, write summary CSV.",
    )
    parser.add_argument(
        "--bridge-policies", type=str,
        help="Comma-separated MUMPS series to include (e.g. "
             "never,always or never,on_change,always).",
    )
    args = parser.parse_args()
    if args.plot_comparison is not None and args.plot_bridge_comparison is not None:
        parser.error("Choose only one comparison plot mode")
    if args.plot_bridge_comparison is not None:
        plot_bridge_comparison(
            args.files, args.plot_bridge_comparison,
            args.rtol, args.atol, args.bridge_policies,
        )
        return
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
