# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Compare native rank-resident RSS *phases* in very short A/B solves.

Purpose: locate sizeable changes across model assembly, PETSc
allocation and time integration (including first MUMPS factorization).
RSS deltas show correlations, not ownership; they cannot attribute
every byte to GetFEM, PETSc or MUMPS. No long experiment needed.
"""
import argparse
import csv
import json
from pathlib import Path


PHASES = (
    "solve_entry",
    "before_petsc_allocate_memory",
    "after_petsc_allocate_memory",
    "before_getfem_assembly",
    "after_mass_assembly",
    "after_stiffness_assembly",
    "after_rhs_assembly",
    "after_jacobian_and_vector_duplicates",
    "after_initial_consistency_ts_cleanup",
    "before_main_ts_solve",
    "main_monitor_step_0",
    "main_monitor_step_1",
    "after_main_ts_solve",
    "after_main_ts_cleanup",
    "solve_exit",
)
REQUIRED = (
    "solve_entry", "before_petsc_allocate_memory",
    "after_petsc_allocate_memory", "after_mass_assembly",
    "after_stiffness_assembly", "after_rhs_assembly",
    "after_jacobian_and_vector_duplicates",
    "before_main_ts_solve", "after_main_ts_solve",
    "after_main_ts_cleanup", "solve_exit",
)


def analyze(directory):
    root = Path(directory)
    rows = []
    details = {}
    commit = None
    for case in ("acoustic", "swe"):
        for mode in ("legacy", "compact"):
            sub = root / case / mode
            files = sorted(sub.glob("memory_phase_rank*.json"))
            if len(files) != 2:
                raise ValueError(
                    f"{case}/{mode}: expected exactly two phase reports"
                )
            ranks = []
            for path in files:
                report = json.loads(path.read_text(encoding="utf-8"))
                if (report["schema"] != "scrimp-phase-memory-v1"
                        or report["ranks"] != 2
                        or report["matrix_memory_mode"] != mode):
                    raise ValueError(f"{path}: invalid phase audit metadata")
                if commit is None:
                    commit = report["source_commit"]
                elif commit != report["source_commit"]:
                    raise ValueError(f"{path}: source commit mismatch")
                phases = report["samples"]
                names = [item["phase"] for item in phases]
                if len(names) != len(set(names)):
                    raise ValueError(f"{path}: duplicate phase sample")
                if not set(REQUIRED).issubset(names):
                    raise ValueError(
                        f"{path}: missing checkpoints "
                        + repr(set(REQUIRED) - set(names))
                    )
                for item in phases:
                    if (item["VmRSS_kib"] <= 0
                            or item["VmHWM_kib"] < item["VmRSS_kib"]):
                        raise ValueError(f"{path}: invalid RSS/HWM")
                ranks.append(report)
            if sorted(x["rank"] for x in ranks) != [0, 1]:
                raise ValueError(f"{case}/{mode}: invalid rank participation")
            table = {
                r["rank"]: {s["phase"]: s for s in r["samples"]}
                for r in ranks
            }
            details[f"{case}/{mode}"] = {}
            for phase in PHASES:
                if any(phase not in snapshots for snapshots in table.values()):
                    continue
                rss = [
                    table[rank][phase]["VmRSS_kib"] / 1024.
                    for rank in (0, 1)
                ]
                hwm = [
                    table[rank][phase]["VmHWM_kib"] / 1024.
                    for rank in (0, 1)
                ]
                value = {
                    "maximum_rank_rss_mib": max(rss),
                    "maximum_rank_hwm_mib": max(hwm),
                    "rank_rss_mib": dict(zip(("0", "1"), rss)),
                }
                details[f"{case}/{mode}"][phase] = value
                rows.append({
                    "case": case, "mode": mode, "phase": phase,
                    "max_rank_rss_mib": round(max(rss), 4),
                    "max_rank_hwm_mib": round(max(hwm), 4),
                    "rank0_rss_mib": round(rss[0], 4),
                    "rank1_rss_mib": round(rss[1], 4),
                })

    comparisons = {}
    for case in ("acoustic", "swe"):
        legacy = details[f"{case}/legacy"]
        compact = details[f"{case}/compact"]
        comparisons[case] = {}
        print(f"=== {case} phase RSS (max of 2 ranks), MiB ===")
        print(f"{'Phase':39s} {'legacy':>9s} {'compact':>9s} {'delta':>9s}")
        for phase in PHASES:
            if phase not in legacy or phase not in compact:
                continue
            l = legacy[phase]["maximum_rank_rss_mib"]
            c = compact[phase]["maximum_rank_rss_mib"]
            diff = c - l
            comparisons[case][phase] = {
                "legacy_mib": l, "compact_mib": c,
                "delta_mib": diff,
            }
            print(f"{phase:39s} {l:9.2f} {c:9.2f} {diff:+9.2f}")

        def jump(table, start, end):
            if start not in table or end not in table:
                return None
            return (table[end]["maximum_rank_rss_mib"]
                    - table[start]["maximum_rank_rss_mib"])

        print(
            "  RSS steps (legacy/compact), MiB: "
            "linear+nonlinear assembly "
            f"{jump(legacy, 'before_getfem_assembly', 'after_rhs_assembly')}/"
            f"{jump(compact, 'before_getfem_assembly', 'after_rhs_assembly')}; "
            "start->first accepted TS step "
            f"{jump(legacy, 'before_main_ts_solve', 'main_monitor_step_1')}/"
            f"{jump(compact, 'before_main_ts_solve', 'main_monitor_step_1')}"
        )
    filename = root / "matrix_phase_memory_summary.json"
    filename.write_text(
        json.dumps({
            "source_commit": commit,
            "cases": details,
            "comparison_by_phase": comparisons,
            "method": (
                "Local Linux VmRSS/VmHWM snapshots; phase-to-phase "
                "differences do NOT precisely attribute native memory "
                "ownership, allocator fragmentation or MUMPS factor bytes."
            ),
        }, indent=2) + "\n",
        encoding="utf-8",
    )
    with (root / "matrix_phase_memory_timeseries.csv").open(
        "w", encoding="utf-8", newline=""
    ) as stream:
        fields = [
            "case", "mode", "phase", "max_rank_rss_mib",
            "max_rank_hwm_mib", "rank0_rss_mib", "rank1_rss_mib",
        ]
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    print("Saved:", filename)
    return comparisons


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", type=Path, required=True)
    args = parser.parse_args()
    analyze(args.directory)


if __name__ == "__main__":
    main()
