# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Current-only large wave strong scaling, without any v1.2.0 baseline.

Validates all native 1/2/4/6-MPI CN states, exact PETSc LU events
against collective on_change decisions, full initial/final MPI
vectors across all repetitions, and produces PNG/SVG figures.
Use a single physical mesh for strong scaling; do not invent timings
or count overlapping profiler phases as disjoint costs.
"""
import argparse
import csv
import json
import math
import statistics
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from analyze_v130_wave_benchmark import (
    PHASES, _expected_lu_updates, _petsc_events, _phase_max,
    _publish_figure, _read_json, _state, _state_distance,
    publish_wave_figures,
)

STEPS = 10
TOLERANCE = 5e-8
EVENTS = ("main.MatLUFactorSym", "main.MatLUFactorNum",
          "main.MatSolve", "initial.MatLUFactorNum")


def median(items):
    return float(statistics.median(items))


def validate_record(directory, name, record, manifest, rank, rep, dofs):
    """Reject non-equivalent simulations before showing speedup figures."""
    if (
        record["variant"] != "current"
        or record["source_sha"] != manifest["current_sha"]
        or record["repetition"] != rep
        or record["mpi_ranks"] != rank
        or record["mesh_h"] != manifest["mesh"]
        or record["scheme"] != "cn"
        or record["ts_adapt_type"] != "none"
        or record["steps"] != STEPS
        or record["rejections"] != 0
        or record["snes_failures"] != 0
        or record["ts_convergence_reason"] <= 0
        or abs(record["time_final_monitor"] - .02) > 1e-10
        or abs(record["dt"] - .002) > 1e-14
        or record["matrix_backend"] != "distributed_coo"
        or record["matrix_memory_mode"] != "compact"
        or record["mumps_policy"] != "on_change"
        or record["main_cn_mumps_icntl28"] != (2 if rank > 1 else 1)
        or (dofs is not None and record["dofs"] != dofs)
    ):
        raise RuntimeError(f"FAIL: inconsistent physical problem or solver in {name}")
    for field in ("wall_solve_s", "wall_build_s", "wall_total_s",
                  "first_step_wall_s", "steady_step_mean_s"):
        if not math.isfinite(record[field]) or record[field] <= 0:
            raise RuntimeError(f"FAIL: missing/nonpositive {field} in {name}")
    times = record["step_times_max_rank_s"]
    if len(times) != STEPS or not all(
        math.isfinite(x) and x > 0 for x in times
    ):
        raise RuntimeError(f"FAIL: invalid 10-step timing vector in {name}")

    rss = record.get("rank_hwm_kib")
    if (
        not isinstance(rss, list) or len(rss) != rank
        or any(type(x) is not int or x <= 0 for x in rss)
        or record.get("max_rank_hwm_kib") != max(rss)
        or record.get("sum_rank_hwm_kib") != sum(rss)
    ):
        raise RuntimeError(f"FAIL: incomplete per-rank RSS high-water marks in {name}")

    updates_num, updates_sym = _expected_lu_updates(record)
    events = _petsc_events(
        directory / (name + ".log"),
        expected_numeric=updates_num,
        expected_symbolic=updates_sym,
    )
    return events


def load_and_validate(directory, manifest, ranks, repetitions):
    if (
        manifest.get("schema") != "scrimp-v130-current-large-wave-v1"
        or manifest.get("ranks") != ",".join(map(str, ranks))
        or manifest.get("repetitions") != repetitions
        or manifest.get("variant") != "current"
        or manifest.get("mumps_policy") != "on_change"
        or manifest.get("source") != "examples/wave.py"
        or manifest.get("steps") != STEPS
    ):
        raise RuntimeError("FAIL: invalid current-only large-wave manifest")
    reports = {}
    events = {}
    initial_ref, final_ref = None, None
    dofs = None
    errors = {}
    for rep in range(1, repetitions+1):
        for rank in ranks:
            label = f"current_on_change_{rank}r_rep{rep}"
            r = _read_json(directory / (label + ".json"))
            e = validate_record(
                directory, label, r, manifest, rank, rep, dofs,
            )
            if dofs is None:
                dofs = r["dofs"]
                if dofs <= 0:
                    raise RuntimeError("FAIL: invalid global DOF count")
            initial, final = _state(directory / (label + "_state.npz"))
            if initial.size != dofs or final.size != dofs:
                raise RuntimeError(f"FAIL: incomplete global PETSc vector in {label}")
            if initial_ref is None:
                initial_ref = initial
                final_ref = final
                diffs = {"scaled_L2_initial": 0., "scaled_Linf_initial": 0.,
                         "scaled_L2_final": 0., "scaled_Linf_final": 0.}
            else:
                il2, iinf = _state_distance(
                    initial_ref, initial, label + " initial", TOLERANCE,
                )
                fl2, finf = _state_distance(
                    final_ref, final, label + " final", TOLERANCE,
                )
                diffs = {"scaled_L2_initial": il2, "scaled_Linf_initial": iinf,
                         "scaled_L2_final": fl2, "scaled_Linf_final": finf}
            reports[(rank, rep)] = r
            events[(rank, rep)] = e
            errors[label] = diffs
            print(
                f"PASS {label}: DOFs={dofs:,} wall_solve={r['wall_solve_s']:.4f}s "
                f"LU=sym{e['main.MatLUFactorSym']['calls']}/"
                f"num{e['main.MatLUFactorNum']['calls']}, "
                f"RSS max={r['max_rank_hwm_kib']/1024**2:.3f} GiB, "
                f"sum={r['sum_rank_hwm_kib']/1024**2:.3f} GiB, "
                f"scaled final state L2={diffs['scaled_L2_final']:.3e}",
                flush=True,
            )
    return reports, events, errors, dofs


def aggregate(reports, events, ranks, repetitions):
    serial = median([reports[(1, rep)]["wall_solve_s"]
                     for rep in range(1, repetitions+1)])
    rows = []
    for rank in ranks:
        group = [reports[(rank, rep)] for rep in range(1, repetitions+1)]
        native = [events[(rank, rep)] for rep in range(1, repetitions+1)]
        solve = [r["wall_solve_s"] for r in group]
        row = {
            "ranks": rank,
            "dofs": group[0]["dofs"],
            "repetitions": repetitions,
            "build_s_median": median([r["wall_build_s"] for r in group]),
            "solve_s_median": median(solve),
            "solve_s_min": min(solve),
            "solve_s_max": max(solve),
            "total_s_median": median([r["wall_total_s"] for r in group]),
            "first_CN_step_s_median": median([
                r["first_step_wall_s"] for r in group
            ]),
            "warm_CN_step_s_median": median([
                r["steady_step_mean_s"] for r in group
            ]),
            "speedup_from_current_1r": serial / median(solve),
            "parallel_efficiency": serial / (rank * median(solve)),
            "rss_max_rank_gib_median": median([
                r["max_rank_hwm_kib"]/1024**2 for r in group
            ]),
            "rss_sum_rank_hwm_gib_median": median([
                r["sum_rank_hwm_kib"]/1024**2 for r in group
            ]),
            "numerical_updates_median": median([
                r["mumps_reuse_stats"]["main"]["numeric_updates"]
                for r in group
            ]),
            "symbolic_updates_median": median([
                r["mumps_reuse_stats"]["main"]["symbolic_updates"]
                for r in group
            ]),
        }
        for key in PHASES:
            values = [_phase_max(record, key) for record in group]
            if all(value is not None for value in values):
                row[key] = median(values)
        for event in EVENTS:
            row[event + ".seconds_median"] = median([
                n[event]["seconds"] for n in native
            ])
            row[event + ".calls_median"] = median([
                n[event]["calls"] for n in native
            ])
        row["cn_per_step_median_s"] = [
            median([r["step_times_max_rank_s"][i] for r in group])
            for i in range(STEPS)
        ]
        rows.append(row)
    return rows


def extra_figures(directory, rows, reports, ranks, repetitions):
    result = []
    x = np.asarray(ranks, dtype=float)
    names = [str(p) for p in ranks]

    # Repeated run distributions. Min/max whiskers, no unsupported
    # confidence interval for three measurements.
    fig, ax = plt.subplots(figsize=(9.2, 5.2))
    med = np.array([r["solve_s_median"] for r in rows])
    ylo = med - np.array([r["solve_s_min"] for r in rows])
    yhi = np.array([r["solve_s_max"] for r in rows]) - med
    ax.errorbar(x, med, yerr=np.vstack([ylo, yhi]), fmt="o-",
                capsize=6, label="Median and observed min/max")
    ax.plot(x, med[0]/x, "--", label="Ideal strong scaling")
    ax.set_xticks(x)
    ax.set_xlabel("MPI ranks")
    ax.set_ylabel("Full solve wall time (s)")
    ax.set_title("Large wave: repeated-solve spread (no confidence interval)")
    ax.grid(alpha=.25)
    ax.legend()
    _publish_figure(fig, directory, "01_repetition_spread")
    result.append("01_repetition_spread")

    fig, ax = plt.subplots(figsize=(9.2, 5.2))
    ax.plot(x, [r["rss_max_rank_gib_median"] for r in rows],
            "o-", label="Largest per-process RSS high-water mark")
    ax.plot(x, [r["rss_sum_rank_hwm_gib_median"] for r in rows],
            "s-", label="Sum of process RSS high-water marks")
    ax.set_xticks(x)
    ax.set_xlabel("MPI ranks")
    ax.set_ylabel("Process memory (GiB, medians)")
    ax.set_title("Large wave: MPI resident-memory capacity")
    ax.text(.01, -.20,
            "Summed independent peaks are NOT simultaneous peak physical memory.",
            transform=ax.transAxes, fontsize=9)
    ax.grid(alpha=.25)
    ax.legend()
    _publish_figure(fig, directory, "08_memory_high_water")
    result.append("08_memory_high_water")

    fig, ax = plt.subplots(figsize=(9.2, 5.2))
    stages = (
        ("setup.getfem.mass", "GetFEM mass"),
        ("setup.bridge.mass", "Mass to PETSc"),
        ("setup.getfem.stiffness", "GetFEM stiffness"),
        ("setup.bridge.stiffness", "Stiffness to PETSc"),
    )
    for key, label in stages:
        if all(key in row for row in rows):
            ax.plot(x, [row[key] for row in rows], "o-", label=label)
    ax.set_xlabel("MPI ranks")
    ax.set_ylabel("Inclusive profiler phase time (s, median rank maximum)")
    ax.set_title("Large wave: mass/stiffness construction and GetFEM→PETSc bridge")
    ax.set_xticks(x)
    ax.grid(alpha=.25)
    ax.legend()
    _publish_figure(fig, directory, "09_getfem_bridge_stages")
    result.append("09_getfem_bridge_stages")

    fig, ax = plt.subplots(figsize=(9.2, 5.2))
    for row in rows:
        ax.plot(np.arange(1, STEPS+1), row["cn_per_step_median_s"],
                "o-", label=f"{row['ranks']} MPI")
    ax.set_xticks(np.arange(1, STEPS+1))
    ax.set_xlabel("Accepted Crank–Nicolson step")
    ax.set_ylabel("Per-step maximum-rank wall time (s), median of repeats")
    ax.set_title("Large wave: first, intermediate and final CN step")
    ax.grid(alpha=.25)
    ax.legend()
    _publish_figure(fig, directory, "10_cn_per_step")
    result.append("10_cn_per_step")

    fig, ax = plt.subplots(figsize=(9.2, 5.2))
    ind = np.arange(len(ranks))
    width = .26
    for offset, key, lab in (
        (-1, "main.MatLUFactorSym.calls_median", "MUMPS symbolic"),
        (0, "main.MatLUFactorNum.calls_median", "MUMPS numeric"),
        (1, "main.MatSolve.calls_median", "MUMPS solves"),
    ):
        bars = ax.bar(ind + offset*width, [row[key] for row in rows],
                      width, label=lab)
        ax.bar_label(bars, fmt="%g", padding=3)
    ax.set_xticks(ind, names)
    ax.set_xlabel("MPI ranks")
    ax.set_ylabel("PETSc event calls in main CN, median")
    ax.set_title("Large wave: native LU reuse verified against exact matrix changes")
    ax.grid(axis="y", alpha=.25)
    ax.legend()
    _publish_figure(fig, directory, "11_native_lu_event_counts")
    result.append("11_native_lu_event_counts")

    fig, ax = plt.subplots(figsize=(9.2, 5.2))
    first = rows[0]["first_CN_step_s_median"]
    warm = rows[0]["warm_CN_step_s_median"]
    ax.plot(x, [first/r["first_CN_step_s_median"] for r in rows],
            "o-", label="First CN step speedup")
    ax.plot(x, [warm/r["warm_CN_step_s_median"] for r in rows],
            "s-", label="Warm CN step speedup")
    ax.plot(x, x, "--", label="Ideal")
    ax.set_xticks(x)
    ax.set_xlabel("MPI ranks")
    ax.set_ylabel("Speedup from current 1-MPI baseline")
    ax.set_title("Large wave: startup factorization vs reusable steady steps")
    ax.grid(alpha=.25)
    ax.legend()
    _publish_figure(fig, directory, "12_cn_startup_vs_warm_speedup")
    result.append("12_cn_startup_vs_warm_speedup")

    return result


def analyze(directory, ranks, repetitions):
    manifest = _read_json(directory / "manifest.json")
    if len(ranks) != len(set(ranks)) or 1 not in ranks:
        raise RuntimeError("MPI rank list must be unique and include rank 1")
    reports, events, errors, dofs = load_and_validate(
        directory, manifest, ranks, repetitions
    )
    rows = aggregate(reports, events, ranks, repetitions)

    with (directory / "v130_large_wave_runs.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        keys = ("repetition", "ranks", "dofs", "wall_build_s",
                "wall_solve_s", "wall_total_s",
                "first_step_wall_s", "steady_step_mean_s",
                "max_rank_hwm_kib", "sum_rank_hwm_kib")
        writer = csv.DictWriter(handle, fieldnames=keys)
        writer.writeheader()
        for rep in range(1, repetitions+1):
            for rank in ranks:
                record = reports[(rank, rep)]
                writer.writerow({
                    **{key: record[key] for key in keys
                       if key not in ("ranks", "repetition")},
                    "repetition": rep, "ranks": rank,
                })

    with (directory / "v130_large_wave_scaling.csv").open(
        "w", newline="", encoding="utf-8",
    ) as handle:
        keys = list(dict.fromkeys(
            key for r in rows for key in r
            if key != "cn_per_step_median_s"
        ))
        writer = csv.DictWriter(handle, fieldnames=keys)
        writer.writeheader()
        writer.writerows([
            {key: value for key, value in r.items()
             if key != "cn_per_step_median_s"} for r in rows
        ])

    figures = publish_wave_figures(
        directory, rows, serial_release=float("nan"),
        serial_current=rows[0]["solve_s_median"],
        release_parity_accepted=False,
    )
    figures.extend(extra_figures(
        directory, rows, reports, ranks, repetitions
    ))
    lines = [
        "# SCRIMP v1.3.0 — ~3-million-DOF current-only wave scaling",
        "",
        f"- Commit: `{manifest['current_sha']}`.",
        f"- Model: `{manifest['source']}`; mesh `h={manifest['mesh']}`.",
        f"- **Actual global DOFs: {dofs:,}**; 10 fixed CN steps, dt=0.002.",
        f"- MPI ranks: {', '.join(map(str, ranks))}.",
        f"- Repetitions per MPI size: {repetitions}; medians are measured.",
        "- **No v1.2.0 execution or comparison** in this campaign.",
        "- Validation: complete initial and final PETSc vector parity "
        "(scaled L2 and Linf <= 5e-8) across every run.",
        "- MUMPS: native symbolic and numeric LU counts verified against "
        "collective exact on_change graph/coefficients decisions.",
        "- Fixed Crank–Nicolson without time adaptation; matched physical problem.",
        "",
        "## Current-only strong scaling",
        "",
        "| MPI | Solve median (s) | Min–max (s) | Total median (s) | "
        "Speedup | Efficiency | Largest-rank RSS (GiB) | Sum of per-rank HWM (GiB) |",
        "|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for r in rows:
        lines.append(
            f"| {r['ranks']} | {r['solve_s_median']:.3f} "
            f"| {r['solve_s_min']:.3f}–{r['solve_s_max']:.3f} "
            f"| {r['total_s_median']:.3f} "
            f"| {r['speedup_from_current_1r']:.3f} "
            f"| {100*r['parallel_efficiency']:.1f}% "
            f"| {r['rss_max_rank_gib_median']:.2f} "
            f"| {r['rss_sum_rank_hwm_gib_median']:.2f} |"
        )
    lines.extend([
        "",
        "## Native PETSc factorization counts",
        "",
        "| MPI | Main symbolic LU | Main numeric LU | Main MatSolve |",
        "|---:|---:|---:|---:|",
    ])
    for r in rows:
        lines.append(
            f"| {r['ranks']} | {r['main.MatLUFactorSym.calls_median']:g} "
            f"| {r['main.MatLUFactorNum.calls_median']:g} "
            f"| {r['main.MatSolve.calls_median']:g} |"
        )
    lines.extend([
        "",
        "## Scientific and measurement limits",
        "",
        "- Three independent runs estimate observed spread, not a statistically "
        "justified confidence interval.",
        "- First-step and warm-step timings exclude monitor callbacks; "
        "phase profiler entries may nest and must not be summed.",
        "- RSS values are Linux per-process lifetime **high-water marks**; "
        "their sum is a capacity indicator, **not a synchronized whole-job "
        "instantaneous peak**. PETSc/MUMPS memory totals may include "
        "transient allocations not visible in a single sample.",
        "- Full initial/final state equality checks are not a complete "
        "trajectory or Hamiltonian/power proof; separate physical tests "
        "remain authoritative.",
        "- MPI speedup/efficiency are strictly relative to **current v1.3.0 "
        "single rank**, not to the historical release.",
        "- PNG figures are in `png/`; high-resolution vector figures "
        "are in `svg/`.",
        "",
    ])
    report = {
        "schema": "scrimp-v130-current-large-wave-analysis-v1",
        "manifest": manifest, "actual_dofs": dofs,
        "repetitions": repetitions, "ranks": ranks,
        "all_full_state_parity": True,
        "all_exact_mumps_factorization_gates": True,
        "tolerance": TOLERANCE,
        "vector_parity_by_run": errors,
        "rows": rows,
        "figures_png": ["png/" + name + ".png" for name in figures],
        "figures_svg": ["svg/" + name + ".svg" for name in figures],
        "rss_sum_is_not_instantaneous_peak": True,
    }
    (directory / "v130_large_wave_summary.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8",
    )
    (directory / "v130_large_wave_performance.md").write_text(
        "\n".join(lines), encoding="utf-8",
    )
    print("\n".join(lines[:22]))
    print("PASS: all current-only 3-million-DOF scaling and parity gates")
    print("Human figures:", directory / "png")
    print("Generated:", ", ".join(figures))
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", type=Path, required=True)
    parser.add_argument("--ranks", default="1,2,4,6")
    parser.add_argument("--repetitions", type=int, default=3)
    args = parser.parse_args()
    ranks = [int(n) for n in args.ranks.split(",")]
    if not 1 <= args.repetitions <= 5:
        parser.error("Repetitions must be 1..5")
    analyze(args.directory, ranks, args.repetitions)


if __name__ == "__main__":
    main()
