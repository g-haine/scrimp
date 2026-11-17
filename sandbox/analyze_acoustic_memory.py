# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Validate and plot per-rank acoustic memory traces.

Require complete dense/sparse rank sets with identical source, PDE,
time scheme and full initial/final NumPy state coefficients. Save
machine-readable JSON/CSV and two plots. Never infer PETSc/MUMPS-only
memory from Linux RSS or count stored states as a memory leak.

Run inside the scrimp-mpi:local Docker image, WITHOUT mpiexec:
    python sandbox/analyze_acoustic_memory.py --directory /artifacts
"""
import argparse
import csv
import json
from pathlib import Path

import numpy as np


def analyze(directory, prefix=None, *, l2_rtol=1e-10, linf_rtol=1e-10):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    directory = Path(directory)
    prefix = Path(prefix) if prefix is not None else directory / "acoustic_memory"
    reports = {}
    sample_series = {}
    physical_signature = None
    commit_ref = None
    reference_states = None
    max_state_abs_difference = 0.
    max_state_relative_l2 = 0.
    max_state_relative_linf = 0.

    for mode in ("dense", "sparse"):
        files = sorted(directory.glob(f"memory_{mode}_rank*.json"))
        if not files:
            raise ValueError(f"Missing {mode} rank reports")
        for path in files:
            report = json.loads(path.read_text(encoding="utf-8"))
            rank = int(report["mpi_rank"])
            size = int(report["mpi_size"])
            if rank in reports.get(mode, {}):
                raise ValueError(f"Duplicate MPI rank: {mode} {rank}")
            if mode != report["save_mode"]:
                raise ValueError(f"{path}: save mode mismatch")
            if report["bridge_method"] != "coo":
                raise ValueError(f"{path}: memory probe requires COO")
            if report["getfem_matrix_mode"] != "distributed":
                raise ValueError(f"{path}: requires distributed GetFEM")
            if report["mumps_reset_policy"] != "never":
                raise ValueError(f"{path}: MUMPS policy mismatch")
            commit = report.get("source_commit")
            if not commit:
                raise ValueError(f"{path}: missing source commit")
            if commit_ref is None:
                commit_ref = commit
            elif commit != commit_ref:
                raise ValueError("Memory traces must come from same commit")
            signature = (
                size, int(report["dofs"]), float(report["mesh"]),
                float(report["dt"]), float(report["final_time"]),
                int(report["expected_steps"]), int(report["sample_every"])
            )
            if physical_signature is None:
                physical_signature = signature
            elif signature != physical_signature:
                raise ValueError(f"{path}: mixed physical/time configurations")
            convergence = report["ts_convergence"]
            if (convergence["reason"] <= 0
                    or convergence["steps"] != report["expected_steps"]
                    or convergence["rejections"] != 0
                    or convergence["snes_failures"] != 0
                    or abs(convergence["last_ts_time"]
                           - report["final_time"]) > 1e-12
                    or abs(report["last_saved_time"]
                           - report["final_time"]) > 1e-12):
                raise ValueError(f"{path}: convergence/endpoint gate failed")
            state_path = directory / (path.stem + "_state.npz")
            with np.load(state_path, allow_pickle=False) as data:
                initial = np.array(data["initial"], copy=True)
                final = np.array(data["final"], copy=True)
                t_end = float(data["final_time"])
            if final.size != int(report["dofs"]):
                raise ValueError(f"{state_path}: wrong global state size")
            if abs(t_end - report["final_time"]) > 1e-12:
                raise ValueError(f"{state_path}: saved final time mismatch")
            if reference_states is None:
                reference_states = (initial, final)
            else:
                for label, actual, expected in (
                    ("initial", initial, reference_states[0]),
                    ("final", final, reference_states[1]),
                ):
                    if actual.shape != expected.shape:
                        raise ValueError(
                            f"{state_path}: {label} state shape differs"
                        )
                    if not np.all(np.isfinite(actual)):
                        raise ValueError(
                            f"{state_path}: {label} state is non-finite"
                        )
                    delta = actual - expected
                    abs_max = float(np.max(np.abs(delta)))
                    scale_l2 = max(
                        float(np.linalg.norm(actual)),
                        float(np.linalg.norm(expected)),
                        float(np.finfo(float).tiny),
                    )
                    scale_linf = max(
                        float(np.max(np.abs(actual))),
                        float(np.max(np.abs(expected))),
                        float(np.finfo(float).tiny),
                    )
                    relative_l2 = float(np.linalg.norm(delta)) / scale_l2
                    relative_linf = abs_max / scale_linf
                    if (relative_l2 > l2_rtol
                            or relative_linf > linf_rtol):
                        raise ValueError(
                            f"{state_path}: {label} state differs across "
                            "rank or saving schedule: "
                            f"max_abs={abs_max:.6e}, "
                            f"rel_L2={relative_l2:.6e}, "
                            f"rel_Linf={relative_linf:.6e} "
                            f"(limits {l2_rtol:g}/{linf_rtol:g})"
                        )
                    max_state_abs_difference = max(
                        max_state_abs_difference, abs_max
                    )
                    max_state_relative_l2 = max(
                        max_state_relative_l2, relative_l2
                    )
                    max_state_relative_linf = max(
                        max_state_relative_linf, relative_linf
                    )
            reports.setdefault(mode, {})[rank] = report

        required = set(range(physical_signature[0]))
        if set(reports[mode]) != required:
            raise ValueError(
                f"{mode}: expected ranks {sorted(required)}, "
                f"found {sorted(reports[mode])}"
            )

    if min(x["snapshots"] for x in reports["dense"].values()) <= max(
        x["snapshots"] for x in reports["sparse"].values()
    ):
        raise ValueError("Dense and sparse save cadences did not diverge")

    rank_rows = []
    series = {}
    for mode in ("dense", "sparse"):
        for rank, report in sorted(reports[mode].items()):
            samples = report["samples"]
            if not samples:
                raise ValueError(f"{mode} rank {rank}: no memory samples")
            step_samples = {
                int(x["step"]): x for x in samples
                if x["phase"] == "main_monitor" and int(x["step"]) >= 0
            }
            if not step_samples or 0 not in step_samples:
                raise ValueError(f"{mode} rank {rank}: missing step samples")
            sample_series[mode, rank] = step_samples

            peak_hwm_kib = max(
                int(x["VmHWM_kib"]) for x in samples
            )
            peak_sampled_rss_kib = max(
                int(x["VmRSS_kib"]) for x in samples
            )
            end = next(
                x for x in reversed(samples)
                if x["phase"] == "after_solve"
            )
            rank_rows.append({
                "mode": mode,
                "rank": rank,
                "snapshots": report["snapshots"],
                "saved_payload_mib": report["stored_state_bytes"] / 2**20,
                "peak_rss_hwm_mib": peak_hwm_kib / 1024,
                "max_sampled_rss_mib": peak_sampled_rss_kib / 1024,
                "post_solve_rss_mib": end["VmRSS_kib"] / 1024,
                "post_solve_rss_anon_mib":
                    end.get("RssAnon_kib", 0) / 1024,
                "solve_wall_seconds": report["solve_wall_seconds_rank"],
                "petsc_mat_memory_local_mib": sum(
                    float(info.get("memory", 0.))
                    for info in report["mat_info_after_solve_local"].values()
                ) / 2**20,
            })

        series[mode] = {}
        steps_per_rank = [
            set(sample_series[mode, rank])
            for rank in sorted(reports[mode])
        ]
        if any(x != steps_per_rank[0] for x in steps_per_rank):
            raise ValueError(f"{mode}: ranks sampled different steps")
        for step in sorted(steps_per_rank[0]):
            local = [
                sample_series[mode, rank][step]
                for rank in sorted(reports[mode])
            ]
            series[mode][step] = {
                "step": step,
                "simulation_time": local[0]["simulation_time"],
                "max_rank_rss_mib": max(
                    x["VmRSS_kib"] for x in local
                ) / 1024,
                "max_rank_anon_mib": max(
                    x.get("RssAnon_kib", 0) for x in local
                ) / 1024,
                "max_rank_hwm_mib": max(
                    x["VmHWM_kib"] for x in local
                ) / 1024,
                "max_rank_payload_mib": max(
                    x["stored_state_bytes"] for x in local
                ) / 2**20,
                "max_rank_snapshots": max(
                    x["stored_snapshots"] for x in local
                ),
            }
    prefix.parent.mkdir(parents=True, exist_ok=True)
    rank_csv = Path(str(prefix) + "_rank_summary.csv")
    with rank_csv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rank_rows[0]))
        writer.writeheader()
        writer.writerows(rank_rows)

    curve_rows = [
        {"mode": mode, **value}
        for mode in ("dense", "sparse")
        for value in series[mode].values()
    ]
    curve_csv = Path(str(prefix) + "_timeseries.csv")
    with curve_csv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(curve_rows[0]))
        writer.writeheader()
        writer.writerows(curve_rows)

    summary = {
        "source_commit": commit_ref,
        "mpi_size": physical_signature[0],
        "dofs": physical_signature[1],
        "mesh": physical_signature[2],
        "dt": physical_signature[3],
        "final_time": physical_signature[4],
        "accepted_steps": physical_signature[5],
        "max_initial_final_state_abs_difference":
            max_state_abs_difference,
        "max_initial_final_state_relative_l2":
            max_state_relative_l2,
        "max_initial_final_state_relative_linf":
            max_state_relative_linf,
        "full_state_relative_error_limits": {
            "l2": l2_rtol,
            "linf": linf_rtol,
        },
        "rank_summary": rank_rows,
        "notes": [
            "VmRSS/VmHWM are rank-local Linux measurements including "
            "Python, NumPy, GetFEM, PETSc and native MUMPS.",
            "Saved payload is the NumPy array size, not an isolated "
            "RSS contribution; allocator caches and shared mappings "
            "prevent additive attribution.",
            "Mat.getInfo memory is PETSc matrix-local, excludes MUMPS LU "
            "factorization memory, and may undercount native allocations.",
            "Comparisons do not yet establish memory-leak behavior; "
            "a plateau test needs longer simulation horizons.",
            "Both runs have the same accepted steps, but a final state "
            "saved within TS.monitor may not be bitwise equal to the "
            "same-time state exported from TS after solve.",
            "Full initial/final state comparisons use strict global "
            "relative L2 and L-infinity error thresholds (each 1e-10), "
            "not per-coefficient np.allclose near-zero tests.",
            "No additional collectives were made during callbacks.",
        ],
    }
    summary_file = Path(str(prefix) + "_summary.json")
    summary_file.write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )

    def plot(field, ylabel, suffix):
        fig, ax = plt.subplots(figsize=(10.5, 6.4))
        for mode in ("dense", "sparse"):
            values = list(series[mode].values())
            ax.plot(
                [x["step"] for x in values],
                [x[field] for x in values],
                marker="." if mode == "dense" else "o",
                markersize=4, linewidth=1.8,
                label="Sauvegarde à chaque pas" if mode == "dense"
                else "Sauvegarde initiale/finale",
            )
        ax.set(
            xlabel="Pas de temps Crank–Nicolson",
            ylabel=ylabel,
            title=(
                "SCRIMP — Audit mémoire, ondes acoustiques 2D\n"
                f"{physical_signature[1]:,} DDL · "
                f"{physical_signature[0]} rangs MPI · "
                f"commit {commit_ref}"
            ).replace(",", " "),
        )
        ax.grid(alpha=0.3)
        ax.legend()
        fig.tight_layout()
        target = Path(str(prefix) + suffix)
        fig.savefig(target, dpi=160)
        plt.close(fig)
        return target

    rss_plot = plot(
        "max_rank_rss_mib", "RSS du rang le plus chargé (Mio)",
        "_rss.png",
    )
    saved_plot = plot(
        "max_rank_payload_mib",
        "Taille des états NumPy sauvegardés par rang (Mio)",
        "_saved_payload.png",
    )

    print(
        f"Validated memory audit: {physical_signature[1]} DoFs, "
        f"{physical_signature[0]} ranks, "
        f"{physical_signature[5]} CN steps, commit {commit_ref}"
    )
    print(
        "Initial/final full-state comparison: PASS "
        "(relative L2 and scaled L-infinity each <= 1e-10)"
    )
    print(
        f"Max absolute coefficient difference: "
        f"{max_state_abs_difference:.9e}; "
        f"relative L2: {max_state_relative_l2:.9e}; "
        f"relative Linf: {max_state_relative_linf:.9e}"
    )
    print(
        "| Save policy | Saved states/rank | Payload MiB/rank "
        "| Highest per-rank HWM MiB |"
    )
    print("|---|---:|---:|---:|")
    for mode in ("dense", "sparse"):
        rows = [r for r in rank_rows if r["mode"] == mode]
        print(
            f"| {mode} | {rows[0]['snapshots']} | "
            f"{max(r['saved_payload_mib'] for r in rows):.2f} | "
            f"{max(r['peak_rss_hwm_mib'] for r in rows):.2f} |"
        )
    print(f"Reports: {summary_file}, {rank_csv}, {curve_csv}")
    print(f"Plots: {rss_plot}, {saved_plot}")
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", required=True, type=Path)
    parser.add_argument("--prefix", type=Path)
    args = parser.parse_args()
    analyze(args.directory, args.prefix)


if __name__ == "__main__":
    main()
