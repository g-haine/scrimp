# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Three reproducible publication-ready figures for serial v1.2.0 vs COO/MPI.

Uses only actual native per-run JSON and FULL initial/end PETSc vectors.
Requires both current reset policies (never/on_change), at the same
mesh, CN dt, end time, rank grid and source revision, with a native
v1.2.0 1-rank control. Never invent missing bars/timing measurements.
"""
import argparse
import csv
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

RANKS = (1, 2, 4, 6)
POLICIES = ("never", "on_change")
COLOR = {"never": "#1266A7", "on_change": "#D26D32"}
ROOT_COLOR = "#575C67"


def _load_reports(directory, repetitions, *, current_only=False):
    all_runs = {}
    cfg = None
    commits = {}
    reference_initial = reference_final = None
    numerical_errors = []
    series = [("current", RANKS, POLICIES)] if current_only else [
        ("release_120", (1,), ("release_native",)),
        ("current", RANKS, POLICIES),
    ]
    for label, ranks, policies in series:
        for policy in policies:
            for rank in ranks:
                values = []
                for repetition in range(1, repetitions + 1):
                    filename = directory / (
                        f"{label}_{policy}_{rank}r_rep{repetition}.json"
                    )
                    data = json.loads(filename.read_text(encoding="utf-8"))
                    if (data["variant"] != label or
                            data["mumps_policy"] != policy or
                            data["mpi_ranks"] != rank or
                            data["repetition"] != repetition or
                            data["steps"] != round(
                                data["time_final_requested"] / data["dt"]
                            ) or data["rejections"] != 0 or
                            data["snes_failures"] != 0 or
                            abs(
                                data["time_final_monitor"] -
                                data["time_final_requested"]
                            ) > 1e-10):
                        raise ValueError(
                            f"{filename}: invalid source/config/convergence"
                        )
                    if (label == "release_120" and (
                        data["source_tag"] != "v1.2.0"
                        or data["matrix_backend"] !=
                        "release_v1.2.0_native")):
                        raise ValueError("Baseline is NOT the original v1.2.0")
                    if (label == "current" and (
                        data["matrix_backend"] != "distributed_coo"
                        or data["ts_convergence_reason"] is None
                        or data["ts_convergence_reason"] <= 0)):
                        raise ValueError("COO run is not the expected current source")
                    for key in ("wall_build_s", "wall_solve_s", "wall_total_s"):
                        if not math.isfinite(data[key]) or data[key] <= 0:
                            raise ValueError(f"{filename}: nonfinite timing")
                    signature = (
                        data["benchmark"], data["mesh_h"], data["dofs"],
                        data["scheme"], data["dt"],
                        data["time_final_requested"], data["save_dt"],
                        data["steps"], data["ts_adapt_type"],
                    )
                    if cfg is None:
                        cfg = signature
                    elif signature != cfg:
                        raise ValueError(f"{filename}: mixed mesh/time or model")
                    if label not in commits:
                        commits[label] = data["source_sha"]
                    elif commits[label] != data["source_sha"]:
                        raise ValueError("Different source revisions in one series")
                    path = directory / data["state_file"]
                    with np.load(path, allow_pickle=False) as archive:
                        initial = np.asarray(archive["initial"])
                        final = np.asarray(archive["final"])
                    if initial.size != data["dofs"] or final.size != data["dofs"]:
                        raise ValueError(f"{filename}: incomplete PETSc state")
                    if not np.all(np.isfinite(final)):
                        raise ValueError("Nonfinite final state")
                    if reference_initial is None:
                        reference_initial = initial
                        reference_final = final
                    init_diff = initial - reference_initial
                    end_diff = final - reference_final
                    def err(delta, a, b):
                        if not np.size(delta):
                            return 0., 0., 0.
                        peak = float(np.max(np.abs(delta)))
                        l2 = float(np.linalg.norm(delta)) / max(
                            float(np.linalg.norm(a)),
                            float(np.linalg.norm(b)),
                            np.finfo(float).tiny,
                        )
                        linf = peak / max(
                            float(np.max(np.abs(a))),
                            float(np.max(np.abs(b))),
                            np.finfo(float).tiny,
                        )
                        return peak, l2, linf
                    init_error = err(init_diff, initial, reference_initial)
                    end_error = err(end_diff, final, reference_final)
                    if (init_error[1] > 1e-8 or init_error[2] > 1e-8
                            or end_error[1] > 1e-7
                            or end_error[2] > 1e-7):
                        raise ValueError(
                            f"{filename}: PHYSICAL STATE PARITY FAILED "
                            f"initial={init_error} final={end_error}"
                        )
                    numerical_errors.append({
                        "case": f"{label}/{policy}/{rank}r/rep{repetition}",
                        "initial_max_abs": init_error[0],
                        "final_max_abs": end_error[0],
                        "final_rel_l2": end_error[1],
                        "final_scaled_linf": end_error[2],
                    })
                    values.append(data)
                all_runs[(label, policy, rank)] = values
    if not current_only and len(set(commits.values())) != 2:
        raise ValueError("Baseline and current source SHAs must differ")
    if current_only and set(commits) != {"current"}:
        raise ValueError("Current-only plot must not load release measurements")
    return all_runs, cfg, commits, numerical_errors


def _stats(reports, key="wall_solve_s"):
    x = np.array([float(d[key]) for d in reports])
    return float(np.median(x)), float(x.min()), float(x.max())


def _save(fig, base):
    """Save only figures under the run's png/ and svg/ directories.

    All JSON/NPZ/CSV measurement files stay in the run root. This
    layout also applies when plotting an existing campaign again.
    """
    png = base.parent / "png" / (base.name + ".png")
    svg = base.parent / "svg" / (base.name + ".svg")
    png.parent.mkdir(parents=True, exist_ok=True)
    svg.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(png, dpi=240)
    fig.savefig(svg)
    plt.close(fig)


def _caption(ax, source, cfg, repetitions, *, current_only=False):
    _, h, ndof, _, dt, final, _, steps, _ = cfg
    caption = (
        f"Acoustic 2D | {ndof:,} DoFs | h={h:g} | CN dt={dt:g}, "
        f"T={final:g} ({steps} steps) | {repetitions} trial(s)/case\n"
        f"Same Docker / CPU pinning; median of complete DPHS solve times. "
        + (
            f"Current branch {source[:8]} only; no v1.2.0 reference "
            f"(MPI-only figures)."
            if current_only else
            f"Tag v1.2.0 vs commit {source[:8]}."
        )
    )
    ax.figure.text(
        .09, .035, caption, fontsize=8.4, va="bottom",
        color="#50555C",
    )


def plot(directory, repetitions, *, current_only=False):
    directory = Path(directory)
    runs, cfg, commits, errors = _load_reports(
        directory, repetitions, current_only=current_only
    )
    current_sha = commits["current"]
    baseline = (
        None if current_only else
        _stats(runs[("release_120", "release_native", 1)])[0]
    )
    current_1 = {
        policy: _stats(runs[("current", policy, 1)])[0]
        for policy in POLICIES
    }
    # The release gain figure is ONLY legitimate when complete
    # initial/final state parity with v1.2.0 has been demonstrated.
    if not current_only:
        fig, ax = plt.subplots(figsize=(10.8, 6.0))
        labels = ["SCRIMP v1.2.0\n1 processus", "Actuel + COO\n1 rang, never",
                  "Actuel + COO\n1 rang, on_change"]
        values = [baseline, current_1["never"], current_1["on_change"]]
        cols = [ROOT_COLOR, COLOR["never"], COLOR["on_change"]]
        bars = ax.bar(labels, values, color=cols, width=.62)
        for i, bar in enumerate(bars):
            h = bar.get_height()
            text = (f"{h:.2f} s" if i == 0 else
                    f"{h:.2f} s\n×{baseline/h:.1f} vs v1.2.0")
            ax.annotate(
                text, (bar.get_x() + bar.get_width()/2, h),
                xytext=(0, 7), textcoords="offset points",
                ha="center", va="bottom", fontsize=10, fontweight="bold",
            )
        ax.set_ylim(0, max(values) * 1.26)
        ax.set_ylabel("Temps de résolution (s) — médiane")
        ax.set_title(
            "Gain algorithmique : v1.2.0 vs SCRIMP actuel (1 seul rang)",
            fontsize=14, weight="bold", loc="left", pad=18,
        )
        ax.spines[["top","right"]].set_visible(False)
        ax.grid(axis="y", alpha=.22)
        ax.set_axisbelow(True)
        _caption(ax, current_sha, cfg, repetitions)
        fig.subplots_adjust(bottom=.21, left=.11, top=.88)
        _save(fig, directory / "01_version120_vs_coo_serial")

    fig, ax = plt.subplots(figsize=(10.8, 6.0))
    if not current_only:
        ax.axhline(
            baseline, color=ROOT_COLOR, lw=1.7, ls="--",
            label=f"v1.2.0 séquentielle ({baseline:.2f} s)"
        )
    for policy in POLICIES:
        summary = [_stats(runs[("current", policy, n)]) for n in RANKS]
        meds = np.array([x[0] for x in summary])
        low = np.array([x[1] for x in summary])
        hi = np.array([x[2] for x in summary])
        ax.plot(
            RANKS, meds, marker="o", markersize=8, lw=2.6,
            color=COLOR[policy], label=f"COO + MUMPS {policy}",
        )
        if repetitions > 1:
            ax.fill_between(
                RANKS, low, hi, color=COLOR[policy], alpha=.13
            )
        for rank, value in zip(RANKS, meds):
            ax.annotate(
                f"{value:.2f}", (rank, value), xytext=(0, 9),
                textcoords="offset points", ha="center", fontsize=9,
            )
    ax.set_xticks(RANKS)
    ax.set_xlabel("Nombre de rangs MPI (maillage fixe)")
    ax.set_ylabel("Temps de résolution (s), médiane")
    ax.set_title("Scalabilité forte — SCRIMP actuel, assemblage COO",
                 fontsize=14, weight="bold", loc="left", pad=18)
    ax.spines[["top","right"]].set_visible(False)
    ax.grid(alpha=.20)
    ax.legend(frameon=False, loc="upper right", fontsize=9.5)
    _caption(ax, current_sha, cfg, repetitions, current_only=current_only)
    fig.subplots_adjust(bottom=.19, left=.11, top=.89)
    _save(fig, directory / "02_strong_scaling_coo")

    fig, ax = plt.subplots(figsize=(10.8, 6.0))
    ax.plot(
        RANKS, RANKS, color="#90959E", ls=":", lw=2,
        label="Scalabilité idéale"
    )
    for policy in POLICIES:
        medians = np.array([
            _stats(runs[("current", policy, n)])[0]
            for n in RANKS
        ])
        speedup = medians[0] / medians
        ax.plot(
            RANKS, speedup, marker="o", markersize=8, lw=2.6,
            color=COLOR[policy], label=f"COO + {policy}",
        )
        for n, val in zip(RANKS, speedup):
            ax.annotate(
                f"{val:.2f}×", (n, val), xytext=(0, 9),
                textcoords="offset points", ha="center", fontsize=9,
            )
    ax.set_xticks(RANKS)
    ax.set_xlabel("Nombre de rangs MPI (maillage fixe)")
    ax.set_ylabel("Accélération forte T₁/Tₚ")
    ax.set_title("Accélération MPI : référence interne à chaque politique",
                 fontsize=14, weight="bold", loc="left", pad=18)
    ax.grid(alpha=.20)
    ax.spines[["top","right"]].set_visible(False)
    ax.legend(frameon=False, loc="upper left")
    _caption(ax, current_sha, cfg, repetitions, current_only=current_only)
    fig.subplots_adjust(bottom=.19, left=.11, top=.89)
    _save(fig, directory / "03_speedup_coo")

    fig, ax = plt.subplots(figsize=(10.8, 6.0))
    ax.axhline(100, ls=":", lw=2, color="#90959E",
               label="Efficacité idéale")
    for policy in POLICIES:
        medians = np.array([
            _stats(runs[("current", policy, n)])[0]
            for n in RANKS
        ])
        efficiency = 100. * medians[0] / (medians * np.asarray(RANKS))
        ax.plot(
            RANKS, efficiency, marker="o", markersize=8, lw=2.6,
            color=COLOR[policy], label=f"COO + {policy}",
        )
        for n, val in zip(RANKS, efficiency):
            ax.annotate(
                f"{val:.0f} %", (n, val), xytext=(0, 9),
                textcoords="offset points", ha="center", fontsize=9,
            )
    ax.set_xticks(RANKS)
    ax.set_xlabel("Nombre de rangs MPI (maillage fixe)")
    ax.set_ylabel("Efficacité parallèle : 100 × T₁/(p Tₚ) (%)")
    ax.set_ylim(0, max(115, ax.get_ylim()[1]))
    ax.set_title("Efficacité parallèle — COO + MUMPS never / on_change",
                 fontsize=14, weight="bold", loc="left", pad=18)
    ax.grid(alpha=.20)
    ax.spines[["top","right"]].set_visible(False)
    ax.legend(frameon=False)
    _caption(ax, current_sha, cfg, repetitions, current_only=current_only)
    fig.subplots_adjust(bottom=.19, left=.11, top=.89)
    _save(fig, directory / "04_efficiency_coo")

    csv_file = directory / "performance_medians.csv"
    with csv_file.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow([
            "variant", "policy", "ranks", "source_sha",
            "median_solve_s", "min_solve_s", "max_solve_s",
            "speedup_vs_v1.2.0", "strong_speedup", "mpi_efficiency",
        ])
        for (variant, policy, ranks), records in sorted(runs.items()):
            med, lo, hi = _stats(records)
            if variant == "current":
                strong = current_1[policy] / med
                eff = strong / ranks
            else:
                strong = 1.
                eff = 1.
            writer.writerow([
                variant, policy, ranks, commits[variant], med, lo, hi,
                baseline/med if baseline is not None else "", strong, eff,
            ])
    report = {
        "provenance": {
            "baseline_tag": (
                None if current_only else "v1.2.0"
            ),
            "baseline_commit": commits.get("release_120"),
            "current_commit": current_sha,
            "benchmark_configuration": {
                "example": cfg[0], "mesh_h": cfg[1], "dofs": cfg[2],
                "ts_type": cfg[3], "dt": cfg[4],
                "final_time": cfg[5], "save_dt": cfg[6],
                "steps": cfg[7], "adapt": cfg[8],
            },
            "repetitions": repetitions,
        },
        "baseline_median_s": baseline,
        "physical_parity_scope": (
            "current_COO_never_on_change_MPI_only"
            if current_only else "current_and_release_v1.2.0"
        ),
        "release_comparison_accepted": not current_only,
        "current_1rank_median_s": current_1,
        "initial_and_final_full_state_parity": {
            "max_final_relative_l2": max(e["final_rel_l2"] for e in errors),
            "max_final_scaled_linf":
                max(e["final_scaled_linf"] for e in errors),
            "all_runs": errors,
        },
        "notes": [
            "No COO or MUMPS never/on_change options retrofitted to v1.2.0.",
            "No v1.2.0 performance gain may be asserted if current_only.",
            "Current-only rank/policy comparisons are still checked for "
            "complete initial and final state parity.",
            "Same wave.py source, same Docker image and same fixed-step CN.",
            "One rank for old release: this is a serial algorithm comparison.",
            "Strong scaling uses current 1-rank T1 for each MUMPS policy.",
            "Speedup vs release is separate and may combine algorithm "
            "improvement with MPI parallelization.",
            "No statistical error bars if only one repetition was run.",
        ],
    }
    (directory / "performance_summary.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    print(
        f"PASS: validated full initial/final states across all "
        f"{len(errors)} runs; "
        f"{3 if current_only else 4} PNG+SVG charts created."
    )
    if current_only:
        print(
            "INFO: release v1.2.0 not part of this analysis. "
            "No cross-version speedup figure or claim."
        )
        print("; ".join(
            f"current/COO/{p} 1r {v:.3f}s"
            for p, v in current_1.items()
        ))
    else:
        print(
            f"v1.2.0 median {baseline:.3f}s; "
            + "; ".join(
                f"current/COO/{p} 1r {v:.3f}s "
                f"(×{baseline/v:.2f} vs release)"
                for p, v in current_1.items()
            )
        )
    for policy in POLICIES:
        print(
            "Current strong scaling " + policy + ": "
            + ", ".join(
                f"{r}r {_stats(runs[('current',policy,r)])[0]:.3f}s"
                for r in RANKS
            )
        )
    print("Saved:", directory / "performance_summary.json")
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", type=Path, required=True)
    parser.add_argument("--repetitions", type=int, default=2)
    parser.add_argument(
        "--current-only", action="store_true",
        help="Build ONLY COO/MPI scaling figures using current-release "
             "measurements. Release-v1.2.0 is not required or loaded; "
             "no cross-version speedup or baseline overlay.",
    )
    args = parser.parse_args()
    if args.repetitions < 1 or args.repetitions > 10:
        parser.error("repetitions must be 1..10")
    plot(args.directory, args.repetitions, current_only=args.current_only)


if __name__ == "__main__":
    main()
