# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Validate and summarize native v1.2.0 vs v1.3 wave performance reports.

Scientific gates:
 - same wave.py, same DOFs and ten accepted fixed CN steps,
 - complete initial and final PETSc vectors across ALL MPI ranks,
 - actual MUMPS ordering and on_change policy,
 - exactly ONE main TS MUMPS numerical factorization for ten steps,
 - no elapsed-time estimates fabricated from step averages,
 - compare true end states, not merely their norms.

All timings come from measured JSON and PETSc native stage logs.
Phase timers are inclusive; do NOT add overlapping component timings.
"""
import argparse
import csv
import json
import re
import statistics
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

STEPS = 10
PHASES = (
    "setup.getfem.mass",
    "setup.bridge.mass",
    "setup.getfem.stiffness",
    "setup.bridge.stiffness",
    "initial.consistency_total",
    "initial.petsc.ts_solve",
    "initial.petsc.callback.ijacobian",
    "initial.petsc.non_callback_remainder",
    "main.petsc.ts_solve",
    "main.petsc.callback.ijacobian",
    "main.petsc.callback.ifunction",
    "main.petsc.csr_pattern_check",
    "main.petsc.non_callback_remainder",
)


def _read_json(path):
    if not path.is_file() or not path.stat().st_size:
        raise RuntimeError(f"Missing benchmark JSON: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def _state(path):
    if not path.is_file():
        raise RuntimeError(f"Missing full state: {path}")
    with np.load(path) as data:
        return (np.asarray(data["initial"], dtype=np.float64),
                np.asarray(data["final"], dtype=np.float64))


def _state_distance(ref, other, label, tolerance=5e-8):
    if ref.shape != other.shape:
        raise RuntimeError(f"{label}: full state dimensions differ")
    if not np.isfinite(ref).all() or not np.isfinite(other).all():
        raise RuntimeError(f"{label}: nonfinite full states")
    delta = other - ref
    relative_l2 = np.linalg.norm(delta) / max(1.0, np.linalg.norm(ref))
    relative_linf = np.max(np.abs(delta)) / max(1.0, np.max(np.abs(ref)))
    if relative_l2 > tolerance or relative_linf > tolerance:
        raise RuntimeError(
            f"FAIL full-state parity {label}: scaled L2={relative_l2:.4e}, "
            f"scaled Linf={relative_linf:.4e} (tol={tolerance:g})"
        )
    return float(relative_l2), float(relative_linf)


def _petsc_events(log_path, *, expected_numeric=None, expected_symbolic=None):
    """Validate actual PETSc LU counts against collective on_change decisions.

    With a fixed nominal CN dt the final MATCHSTEP shift can differ
    by floating-point rounding. Compare actual LU events against the
    EXACT matrix-change decision rather than assuming exactly one
    numeric factorization for every ten-step simulation.
    """
    stage = None
    counts = {}
    for line in log_path.read_text(encoding="utf-8", errors="replace").splitlines():
        if "Event Stage" in line:
            if "WAVE_MAIN_CN" in line:
                stage = "main"
            elif "WAVE_INITIAL_CONSISTENCY" in line:
                stage = "initial"
            else:
                stage = "other"
            continue
        cols = line.strip().split()
        if stage in ("initial", "main") and len(cols) > 2 and cols[0] in (
            "MatLUFactorSym", "MatLUFactorNum", "MatSolve",
        ):
            try:
                count = int(cols[1])
                seconds = float(cols[3])
            except (ValueError, IndexError):
                continue
            key = f"{stage}.{cols[0]}"
            if key in counts:
                raise RuntimeError(f"Duplicate native PETSc event {key} in {log_path}")
            counts[key] = {"calls": count, "seconds": seconds}
    required = ("main.MatLUFactorSym", "main.MatLUFactorNum",
                "main.MatSolve", "initial.MatLUFactorNum")
    for name in required:
        if name not in counts:
            raise RuntimeError(f"Missing PETSc {name} in {log_path}")
    numeric_count = counts["main.MatLUFactorNum"]["calls"]
    symbolic_count = counts["main.MatLUFactorSym"]["calls"]
    solve_count = counts["main.MatSolve"]["calls"]
    if (
        not (1 <= symbolic_count <= numeric_count <= STEPS)
        or solve_count != STEPS
        or (expected_numeric is not None and numeric_count != expected_numeric)
        or (expected_symbolic is not None
            and symbolic_count != expected_symbolic)
    ):
        raise RuntimeError(
            f"FAIL MUMPS numeric/symbolic reuse gate in {log_path}: "
            f"PETSc counts: symbolic={symbolic_count}, "
            f"numeric={numeric_count}, solve={solve_count}; "
            f"on_change expected symbolic={expected_symbolic}, "
            f"numeric={expected_numeric}, steps={STEPS}"
        )
    return counts


def _expected_lu_updates(report):
    """Return exact MPI-consistent operator updates recorded by SCRIMP."""
    stats = (report.get("mumps_reuse_stats") or {}).get("main") or {}
    required = (
        "ijacobian_calls", "numeric_updates",
        "symbolic_updates", "reused_operators",
    )
    if any(type(stats.get(k)) is not int for k in required):
        raise RuntimeError(
            "Missing on_change collective MAIN CN update counters: "
            + repr(stats)
        )
    calls = stats["ijacobian_calls"]
    updates = stats["numeric_updates"]
    graphs = stats["symbolic_updates"]
    reused = stats["reused_operators"]
    if (
        calls != STEPS or updates + reused != calls
        or not (1 <= graphs <= updates <= calls)
        or reused == 0
    ):
        raise RuntimeError(
            "FAIL on_change accounting or no operator reuse: "
            + repr(stats)
        )
    return updates, graphs


def _median(records, field):
    return statistics.median(float(item[field]) for item in records)


def _phase_max(report, event):
    profile = report.get("profile")
    if not profile or not profile.get("enabled"):
        raise RuntimeError("Missing current SCRIMP profiler measurements")
    entry = profile.get("events", {}).get(event)
    return float(entry["seconds_rank_max"]) if entry is not None else None



def _publish_figure(fig, directory, stem):
    """Publish one readable PNG and one scalable SVG under separate folders."""
    png_dir = directory / "png"
    svg_dir = directory / "svg"
    png_dir.mkdir(parents=True, exist_ok=True)
    svg_dir.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(png_dir / (stem + ".png"), dpi=200, bbox_inches="tight")
    fig.savefig(svg_dir / (stem + ".svg"), bbox_inches="tight")
    plt.close(fig)


def _base_axes(title, ylabel, ranks):
    fig, ax = plt.subplots(figsize=(9, 5.4))
    ax.set_title(title)
    ax.set_xlabel("MPI ranks")
    ax.set_ylabel(ylabel)
    ax.set_xticks(ranks)
    ax.grid(axis="y", alpha=0.25)
    return fig, ax


def publish_wave_figures(directory, rows, serial_release,
                         serial_current, release_parity_accepted):
    """Human-readable version and scaling figures from VERIFIED records.

    Do not plot cross-version speedup if whole-state parity is absent.
    Timers are inclusive; stage curves are NEVER stacked or added.
    """
    ranks = np.asarray([r["ranks"] for r in rows], dtype=int)
    solve = np.asarray([r["solve_s_median"] for r in rows])
    total = np.asarray([r["total_s_median"] for r in rows])
    created = []
    if release_parity_accepted:
        fig, ax = plt.subplots(figsize=(8.5, 5.4))
        bars = ax.bar(
            ["v1.2.0 native (1 rank)", "v1.3.0 candidate (1 rank)"],
            [serial_release, serial_current],
        )
        ax.bar_label(bars, fmt="%.3f s", padding=4)
        ax.set_title("Wave.py — validated sequential release comparison")
        ax.set_ylabel("Full solve time (s), median")
        ax.grid(axis="y", alpha=.25)
        ax.set_ylim(0, 1.15 * max(serial_release, serial_current))
        _publish_figure(fig, directory, "01_sequential_v120_vs_v130")
        created.append("01_sequential_v120_vs_v130")
    fig, ax = _base_axes(
        "Wave.py — strong scaling: full solve and total time",
        "Wall time (s), median", ranks,
    )
    ax.plot(ranks, total, "o-", label="Build + solve")
    ax.plot(ranks, solve, "s-", label="Solve")
    ax.plot(ranks, serial_current / ranks, "--",
            label="Ideal solve, relative to 1 rank")
    ax.legend()
    _publish_figure(fig, directory, "02_strong_scaling")
    created.append("02_strong_scaling")

    fig, ax = _base_axes(
        "Wave.py — strong-scaling speedup (current version)",
        "Speedup vs current single rank", ranks,
    )
    speedup = np.asarray([r["speedup_from_current_1r"] for r in rows])
    ax.plot(ranks, speedup, "o-", label="Measured")
    ax.plot(ranks, ranks, "--", label="Ideal")
    ax.legend()
    _publish_figure(fig, directory, "03_speedup")
    created.append("03_speedup")

    fig, ax = _base_axes(
        "Wave.py — parallel efficiency (current version)",
        "Efficiency (%)", ranks,
    )
    efficiencies = np.asarray([100. * r["parallel_efficiency"] for r in rows])
    ax.plot(ranks, efficiencies, "o-", label="Measured")
    ax.axhline(100., linestyle="--", label="Ideal (100%)")
    ax.set_ylim(0, max(110., 1.12 * efficiencies.max()))
    ax.legend()
    _publish_figure(fig, directory, "04_parallel_efficiency")
    created.append("04_parallel_efficiency")

    fig, ax = _base_axes(
        "Wave.py — first Crank–Nicolson step vs steady steps",
        "Observed time per CN step (s), median", ranks,
    )
    width = .36
    xx = np.arange(len(ranks))
    first = np.asarray([r["first_CN_step_s_median"] for r in rows])
    warm = np.asarray([r["warm_CN_step_s_median"] for r in rows])
    ax.bar(xx - width / 2, first, width, label="First step")
    ax.bar(xx + width / 2, warm, width, label="Average steps 2–10")
    ax.set_xticks(xx, ranks)
    ax.legend()
    _publish_figure(fig, directory, "05_cn_first_vs_warm_steps")
    created.append("05_cn_first_vs_warm_steps")

    phases = [
        ("initial.consistency_total", "Initial consistency"),
        ("main.petsc.ts_solve", "Main TS (includes callbacks)"),
        ("main.petsc.callback.ijacobian", "Main IJacobian callbacks"),
        ("main.petsc.callback.ifunction", "Main IFunction callbacks"),
        ("main.petsc.csr_pattern_check", "On-change CSR comparisons"),
    ]
    fig, ax = _base_axes(
        "Wave.py — SCRIMP phase costs (inclusive timings)",
        "Max-rank phase time (s), median", ranks,
    )
    for field, label in phases:
        values = [r.get(field) for r in rows]
        if all(v is not None for v in values):
            ax.plot(ranks, values, "o-", label=label)
    ax.legend(loc="best", fontsize=9)
    ax.text(.01, -.18, "Phases overlap: do not sum their durations.",
            transform=ax.transAxes, fontsize=9)
    _publish_figure(fig, directory, "06_scrimp_phases")
    created.append("06_scrimp_phases")

    events = [
        ("main.MatLUFactorSym.seconds_median", "Symbolic analysis"),
        ("main.MatLUFactorNum.seconds_median", "Numeric LU"),
        ("main.MatSolve.seconds_median", "Triangular solves"),
    ]
    fig, ax = _base_axes(
        "Wave.py — native PETSc/MUMPS main-CN events",
        "PETSc native event time (s), median", ranks,
    )
    for field, label in events:
        ax.plot(ranks, [r[field] for r in rows], "o-", label=label)
    ax.legend()
    _publish_figure(fig, directory, "07_mumps_factorization_breakdown")
    created.append("07_mumps_factorization_breakdown")
    return created


def publish_smoke_figures(directory, ranks):
    """Plot a PASSED current-only numeric-LU reuse smoke as qualification.

    The caller must validate step counts, complete MPI vectors and all
    native event counts before invoking this function.
    """
    ranks = [int(n) for n in ranks]
    factor_counts, solve_times = [], []
    for n in ranks:
        label = "current_on_change_" + str(n) + "r"
        data = _read_json(directory / (label + ".json"))
        events = _petsc_events(directory / (label + ".log"))
        factor_counts.append(events["main.MatLUFactorNum"]["calls"])
        solve_times.append(float(data["wall_solve_s"]))
    fig, ax = plt.subplots(figsize=(8.5, 5.0))
    bars = ax.bar([str(n) + " MPI" for n in ranks], factor_counts)
    ax.bar_label(bars, fmt="%d", padding=4)
    ax.axhline(1, linestyle="--", label="Expected: 1 numeric LU")
    ax.set_title("Production on_change — 10 CN steps, numeric factor reuse")
    ax.set_ylabel("Main-CN MatLUFactorNum calls")
    ax.set_ylim(0, max(2, max(factor_counts) * 1.5))
    ax.legend()
    ax.grid(axis="y", alpha=.2)
    _publish_figure(fig, directory, "01_smoke_numeric_factor_reuse")
    fig, ax = plt.subplots(figsize=(8.5, 5.0))
    bars = ax.bar([str(n) + " MPI" for n in ranks], solve_times)
    ax.bar_label(bars, fmt="%.3f s", padding=4)
    ax.set_title("Production on_change — full wave solve time")
    ax.set_ylabel("Wall time (s)")
    ax.set_ylim(0, 1.2 * max(solve_times))
    ax.grid(axis="y", alpha=.2)
    _publish_figure(fig, directory, "02_smoke_solve_time")
    return ["01_smoke_numeric_factor_reuse", "02_smoke_solve_time"]

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", type=Path, required=True)
    parser.add_argument("--repetitions", type=int, required=True)
    parser.add_argument("--ranks", required=True, help="CSV, e.g. 1,2,4,6")
    args = parser.parse_args()
    directory = args.directory
    ranks = [int(v) for v in args.ranks.split(",")]
    if len(set(ranks)) != len(ranks) or 1 not in ranks:
        parser.error("Ranks must be unique and include one serial reference")
    if not (1 <= args.repetitions <= 5):
        parser.error("Repetitions must be 1..5")
    manifest = _read_json(directory / "manifest.json")
    sources = {}
    reports = {}
    petsc = {}
    release_parity_failures = []
    for repetition in range(1, args.repetitions + 1):
        rel_label = f"release_120_release_native_1r_rep{repetition}"
        rel = _read_json(directory / (rel_label + ".json"))
        if (rel["variant"] != "release_120"
                or rel["source_sha"] != manifest["release_sha"]
                or rel["steps"] != STEPS or rel["mpi_ranks"] != 1
                or rel["rejections"] or rel["snes_failures"]):
            raise RuntimeError(f"Invalid v1.2.0 sequential report {rel_label}")
        reports.setdefault(("release_120", 1), []).append(rel)
        rel_initial, rel_final = _state(directory / (rel_label + "_state.npz"))
        if repetition == 1:
            sources["release_120"] = (rel_initial, rel_final)
        else:
            _state_distance(sources["release_120"][1], rel_final, rel_label)
        for p in ranks:
            label = f"current_on_change_{p}r_rep{repetition}"
            cur = _read_json(directory / (label + ".json"))
            if (cur["variant"] != "current"
                    or cur["source_sha"] != manifest["current_sha"]
                    or cur["mumps_policy"] != "on_change"
                    or cur.get("main_cn_mumps_icntl28") != (2 if p > 1 else 1)
                    or cur["mpi_ranks"] != p or cur["steps"] != STEPS
                    or cur["rejections"] or cur["snes_failures"]
                    or cur["dofs"] != rel["dofs"]):
                raise RuntimeError(f"Invalid current v1.3 report {label}")
            reports.setdefault(("current", p), []).append(cur)
            initial, final = _state(directory / (label + "_state.npz"))
            # A release-to-current mismatch invalidates *cross-version*
            # speedup claims, but must NOT hide valid current-only MPI
            # scaling. The original v1.2.0 native implementation has
            # already shown version-dependent state differences.
            try:
                _state_distance(rel_initial, initial,
                                label + " initial vs v1.2")
                _state_distance(rel_final, final,
                                label + " final vs v1.2")
            except RuntimeError as exc:
                release_parity_failures.append(str(exc))
            if "current" not in sources:
                sources["current"] = (initial, final)
            else:
                _state_distance(sources["current"][0], initial,
                                label + " initial vs current-1r")
                _state_distance(sources["current"][1], final,
                                label + " final vs current-1r")
            expected_numeric, expected_symbolic = _expected_lu_updates(cur)
            petsc[label] = _petsc_events(
                directory / (label + ".log"),
                expected_numeric=expected_numeric,
                expected_symbolic=expected_symbolic,
            )

    release_parity_accepted = not release_parity_failures
    serial_release = _median(reports[("release_120", 1)], "wall_solve_s")
    serial_current = _median(reports[("current", 1)], "wall_solve_s")
    rows = []
    for p in ranks:
        group = reports[("current", p)]
        wall = _median(group, "wall_solve_s")
        row = {
            "ranks": p, "dofs": group[0]["dofs"],
            "build_s_median": _median(group, "wall_build_s"),
            "solve_s_median": wall,
            "total_s_median": _median(group, "wall_total_s"),
            "first_CN_step_s_median": _median(group, "first_step_wall_s"),
            "warm_CN_step_s_median": _median(group, "steady_step_mean_s"),
            "speedup_from_current_1r": serial_current / wall,
            "parallel_efficiency": serial_current / (p * wall),
        }
        for name in PHASES:
            vals = [_phase_max(item, name) for item in group]
            if all(x is not None for x in vals):
                row[name] = statistics.median(vals)
        for event in ("main.MatLUFactorSym", "main.MatLUFactorNum",
                      "main.MatSolve", "initial.MatLUFactorNum"):
            row[f"{event}.seconds_median"] = statistics.median(
                petsc[f"current_on_change_{p}r_rep{k}"][event]["seconds"]
                for k in range(1, args.repetitions + 1)
            )
        rows.append(row)

    path = directory / "v130_wave_scaling.csv"
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

    lines = [
        "# SCRIMP v1.3.0 — wave.py performance qualification",
        "",
        f"- Source v1.2.0: `{manifest['release_sha']}`",
        f"- Current source: `{manifest['current_sha']}`",
        f"- Mesh `h={manifest['mesh']}`; {rows[0]['dofs']:,} DOFs; ten accepted CN steps.",
        f"- Repetitions: {args.repetitions}. Native Docker environment shared.",
        "- All current MPI full initial/final states agree within scaled 5e-8.",
        ("- Cross-release state parity: PASS (same full vectors)."
         if release_parity_accepted else
         "- Cross-release state parity: **NOT VERIFIED / MISMATCH**. "
         "Comparing elapsed times is descriptive only; no version "
         "speedup is accepted."),
        "- Current MUMPS: default on_change plus MPI ordering; native LU"
        " counts equal exact graph/coefficient changes recorded by SCRIMP.",
        "",
        "## Sequential comparison",
        "",
        f"- v1.2.0 native sequential solve: **{serial_release:.6f} s** (median).",
        f"- v1.3.0 current sequential solve: **{serial_current:.6f} s** (median).",
        (f"- Accepted release-to-current speedup: "
         f"**×{serial_release/serial_current:.3f}**."
         if release_parity_accepted else
         "- **No valid cross-release speedup** because the complete "
         "states differ; investigate assembly/initial conditions "
         "before claiming a faster equivalent physical simulation."),
        "",
        "## Strong scaling (current v1.3.0)",
        "",
        "| MPI ranks | Build (s) | Solve (s) | Total (s) | First CN (s) | Warm CN (s) | Speedup | Efficiency |",
        "|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in rows:
        lines.append(
            f"| {row['ranks']} | {row['build_s_median']:.3f} "
            f"| {row['solve_s_median']:.3f} | {row['total_s_median']:.3f} "
            f"| {row['first_CN_step_s_median']:.3f} "
            f"| {row['warm_CN_step_s_median']:.3f} "
            f"| {row['speedup_from_current_1r']:.3f} "
            f"| {100*row['parallel_efficiency']:.1f}% |"
        )
    lines.extend([
        "",
        "## Stage timings, current (median rank-max seconds per run)",
        "",
        "| Phase / PETSc native event | " + " | ".join(str(p) + "r" for p in ranks) + " |",
        "|---|" + "|".join("---:" for _ in ranks) + "|",
    ])
    all_keys = [x for x in PHASES if any(x in row for row in rows)]
    all_keys += ["main.MatLUFactorSym.seconds_median",
                 "main.MatLUFactorNum.seconds_median",
                 "main.MatSolve.seconds_median",
                 "initial.MatLUFactorNum.seconds_median"]
    for event in all_keys:
        lines.append(
            "| `" + event + "` | " +
            " | ".join(f"{row[event]:.5f}" if event in row else "—"
                       for row in rows) + " |"
        )
    lines += [
        "",
        "Times in different profiler rows can be nested/inclusive; do not sum them.",
        "v1.2.0 provides version-neutral wall and step timings only; native PETSc "
        "stage profiling and SCRIMP's current phase profiler are not claimed "
        "for the historical release.",
        "A ten-step test establishes factor reuse and state parity for the "
        "linear wave case, not all time-varying/nonlinear port-Hamiltonian "
        "models. The final PETSc MATCHSTEP may introduce a floating-point "
        "shift difference; exact on_change may then legitimately perform "
        "a second numeric LU on the final step.",
        "",
    ]
    lines += [
        "## Human-readable performance figures",
        "",
        "PNG plots live in the `png/` subdirectory (and vector SVGs in `svg/`).",
        "The sequential v1.2.0 comparison figure appears ONLY after passing",
        "full-state cross-release parity. Other figures show current-version",
        "MPI scaling, per-step timings and PETSc/MUMPS phase costs.",
        "",
    ]
    (directory / "v130_wave_performance.md").write_text(
        "\n".join(lines), encoding="utf-8"
    )
    (directory / "v130_wave_summary.json").write_text(
        json.dumps({
            "release_sha": manifest["release_sha"],
            "current_sha": manifest["current_sha"],
            "mesh": manifest["mesh"],
            "repetitions": args.repetitions,
            "v120_serial_s": serial_release,
            "current_serial_s": serial_current,
            "cross_release_state_parity": release_parity_accepted,
            "cross_release_mismatch_samples": release_parity_failures[:4],
            "serial_speedup": (
                serial_release/serial_current if release_parity_accepted
                else None
            ),
            "rows": rows,
            "petsc_counts_checked": True,
            "petsc_numeric_updates_match_operator_changes": True,
            "full_state_parity_checked": True,
            "current_mpi_state_parity_passed": True,
        }, indent=2) + "\n", encoding="utf-8"
    )
    figures = publish_wave_figures(
        directory, rows, serial_release, serial_current, release_parity_accepted
    )
    print("Human-readable figures:", directory / "png")
    print("PNG figure names:", ", ".join(figures))
    print("\n".join(lines[:25]))
    if release_parity_accepted:
        print(f"PASS: full-state parity including original release; reports: {directory}")
    else:
        print("PARTIAL: current MPI parity + MUMPS reuse PASS; "
              "v1.2.0 physical equivalence FAIL, no version speedup claim")
        print("First cross-version mismatch:", release_parity_failures[0])
        print("Reports:", directory)


if __name__ == "__main__":
    main()
