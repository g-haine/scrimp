# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
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


def _petsc_events(log_path):
    """Extract native PETSc MatLUFactorSym/Num/MatSolve counts by stage."""
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
    if (counts["main.MatLUFactorNum"]["calls"] != 1
            or counts["main.MatLUFactorSym"]["calls"] != 1
            or counts["main.MatSolve"]["calls"] != STEPS):
        raise RuntimeError(
            f"FAIL MUMPS numeric factor reuse gate in {log_path}: "
            f"{counts} (expected one symbolic, one numeric, {STEPS} solves)"
        )
    return counts


def _median(records, field):
    return statistics.median(float(item[field]) for item in records)


def _phase_max(report, event):
    profile = report.get("profile")
    if not profile or not profile.get("enabled"):
        raise RuntimeError("Missing current SCRIMP profiler measurements")
    entry = profile.get("events", {}).get(event)
    return float(entry["seconds_rank_max"]) if entry is not None else None


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
            petsc[label] = _petsc_events(directory / (label + ".log"))

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
        "- Current MUMPS: default on_change plus MPI ordering, exact one numeric LU per ten CN steps.",
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
        "A ten-step test establishes correctness and factor reuse for the "
        "linear wave case, not all time-varying/nonlinear port-Hamiltonian models.",
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
            "full_state_parity_checked": True,
            "current_mpi_state_parity_passed": True,
        }, indent=2) + "\n", encoding="utf-8"
    )
    fig, ax = plt.subplots(figsize=(8, 5))
    x = np.asarray([r["ranks"] for r in rows])
    y = np.asarray([r["solve_s_median"] for r in rows])
    ax.plot(x, y, "o-", label="v1.3.0 on_change")
    ax.plot(x, serial_current/x, "--", label="Ideal scaling from current 1r")
    if release_parity_accepted:
        ax.axhline(serial_release, linestyle=":",
                   label="Equivalent native v1.2.0 1r")
    ax.set_xlabel("MPI ranks")
    ax.set_ylabel("SCRIMP solve time (s), median")
    ax.set_xticks(x)
    ax.grid(alpha=.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(directory / "v130_wave_scaling.png", dpi=170)
    plt.close(fig)
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
