# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Full short-run legacy/compact matrix-policy parity + peak-RSS analysis.

Do not mistake a one-shot whole-process VmHWM difference for a
precise per-matrix allocation measurement or for sustained scalability.
"""
import argparse
import json
from pathlib import Path

import numpy as np


def compare(directory, *, rtol=1e-9, atol=1e-10):
    directory = Path(directory)
    result = {
        "test": "short native matrix-memory A/B",
        "rtol": rtol, "atol": atol,
        "cases": {},
        "notes": [
            "All cases use the SAME physical model and time scheme.",
            "Whole-process VmHWM includes PETSc/MUMPS, GetFEM and Python.",
            "Peak RSS is not equal to allocated PETSc AIJ storage.",
            "One run per configuration; do not infer timing speedups.",
            "SWE balance parity does not establish absolute physics accuracy.",
        ],
    }
    commit = None
    for case, suffix in (
        ("acoustic", "_physical"), ("swe", "_swe"),
    ):
        folders = {
            mode: directory / case / mode
            for mode in ("legacy", "compact")
        }
        names = {
            mode: f"memory{suffix}" for mode in ("legacy", "compact")
        }
        reports = {
            mode: json.loads(
                (folder / f"{names[mode]}.json").read_text(
                    encoding="utf-8"
                )
            )
            for mode, folder in folders.items()
        }
        a, b = reports["legacy"], reports["compact"]
        if a["matrix_memory_mode"] != "legacy" or (
            b["matrix_memory_mode"] != "compact"
        ):
            raise ValueError(f"{case}: matrix modes are not legacy/compact")
        if a["source_commit"] != b["source_commit"]:
            raise ValueError(f"{case}: mixed source commits")
        if commit is None:
            commit = a["source_commit"]
        elif commit != a["source_commit"]:
            raise ValueError("Tests use different source commits")

        attributes = (
            ("mpi_size", "dofs", "n_snapshots", "steps",
             "rejected", "snes_failures", "final_time")
            if case == "acoustic" else
            ("ranks", "dofs", "snapshots", "mesh", "dt",
             "final_time", "spinup_duration", "ts_type",
             "ts_adapt_type", "mumps_reset_policy")
        )
        for key in attributes:
            if a[key] != b[key]:
                raise ValueError(f"{case}: changed model config {key}")
        if case == "swe":
            for mode in ("legacy", "compact"):
                cv = reports[mode]["convergence"]
                if (cv["reason"] <= 0 or cv["rejections"]
                        or cv["snes_failures"]
                        or cv["steps"] != 3
                        or abs(cv["t_final"] - a["final_time"]) > 1e-12):
                    raise ValueError(f"{case}/{mode}: nonlinear TS failed")
        else:
            if (a["rejected"] or a["snes_failures"] or a["steps"] != 3):
                raise ValueError("Acoustic TS did not run exactly 3 steps")

        comparison = {}
        with np.load(
            folders["legacy"] / f"{names['legacy']}.npz",
            allow_pickle=False,
        ) as old, np.load(
            folders["compact"] / f"{names['compact']}.npz",
            allow_pickle=False,
        ) as new:
            if set(old.files) != set(new.files):
                raise ValueError(f"{case}: different computed quantities")
            if "states" not in old.files:
                raise ValueError(f"{case}: missing full snapshots")
            for key in old.files:
                x = np.asarray(old[key])
                y = np.asarray(new[key])
                if x.shape != y.shape:
                    raise ValueError(f"{case}/{key}: different shape")
                if not np.all(np.isfinite(x)) or not np.all(np.isfinite(y)):
                    raise ValueError(f"{case}/{key}: nonfinite values")
                deviation = y - x
                maximum = float(np.max(np.abs(deviation))) if x.size else 0.
                xnorm = float(np.linalg.norm(x))
                ynorm = float(np.linalg.norm(y))
                relative = float(np.linalg.norm(deviation)) / max(
                    xnorm, ynorm, np.finfo(float).tiny
                )
                if key in ("t", "times"):
                    if not np.array_equal(x, y):
                        raise ValueError(f"{case}: different time samples")
                elif not np.allclose(x, y, rtol=rtol, atol=atol):
                    raise ValueError(
                        f"{case}/{key}: legacy vs compact output DIFFER: "
                        f"max_abs={maximum:.4e}, relL2={relative:.4e}"
                    )
                comparison[key] = {
                    "maximum_abs_difference": maximum,
                    "relative_l2_error": relative,
                }

        mem = {}
        for mode in ("legacy", "compact"):
            samples = reports[mode]["matrix_memory_per_rank"]
            ranks = (
                reports[mode]["mpi_size"]
                if case == "acoustic" else reports[mode]["ranks"]
            )
            if len(samples) != ranks or sorted(
                x["mpi_rank"] for x in samples
            ) != list(range(ranks)):
                raise ValueError(f"{case}/{mode}: incomplete rank metrics")
            for sample in samples:
                if sample["matrix_memory_mode"] != mode:
                    raise ValueError("Mismatched MPI matrix policy")
                if (sample["tangent_mass_is_mass"]
                        != (mode == "compact" and sample["linear_mass"])):
                    raise ValueError(
                        f"{case}/{mode}: unexpected mass tangent ownership"
                    )
                if (sample["tangent_stiffness_is_stiffness"]
                        != (mode == "compact" and sample["linear_stiffness"])):
                    raise ValueError(
                        f"{case}/{mode}: stiffness tangent ownership"
                    )
                if case == "swe" and mode == "compact":
                    if not sample["linear_mass"]:
                        if (sample["hessian_fd_matrix_allocations"] != 1
                                or sample["hessian_fd_probe_allocations"] != 1
                                or sample["hessian_fd_matrix_reuses"] < 1):
                            raise ValueError(
                                "SWE Hessian workspace was not reused"
                            )
            mem[mode] = {
                "maximum_rank_hwm_mib": max(
                    x["VmHWM_kib"] for x in samples
                ) / 1024.,
                "rank_rss_end_mib": {
                    str(x["mpi_rank"]): x["VmRSS_kib"] / 1024.
                    for x in samples
                },
                "per_rank": samples,
            }
        comparison_dict = {
            "dofs": a["dofs"],
            "mpi_ranks": len(reports["compact"]["matrix_memory_per_rank"]),
            "snapshots": (
                a["n_snapshots"] if case == "acoustic" else a["snapshots"]
            ),
            "quantities": comparison,
            "memory": mem,
            "peak_hwm_delta_mib_compact_minus_legacy":
                (mem["compact"]["maximum_rank_hwm_mib"]
                 - mem["legacy"]["maximum_rank_hwm_mib"]),
        }
        result["cases"][case] = comparison_dict
        print(
            f"PASS {case}: {a['dofs']} DDL, "
            f"{comparison_dict['mpi_ranks']} MPI ranks, "
            f"{comparison_dict['snapshots']} states; "
            f"full-state maxabs={comparison['states']['maximum_abs_difference']:.4e}"
        )
        print(
            "  rank-high-water VmHWM: "
            f"legacy={mem['legacy']['maximum_rank_hwm_mib']:.2f} MiB, "
            f"compact={mem['compact']['maximum_rank_hwm_mib']:.2f} MiB "
            f"(delta {comparison_dict['peak_hwm_delta_mib_compact_minus_legacy']:+.2f} MiB)"
        )
        if case == "swe":
            print(
                "  power/volume checks: "
                f"balance maxabs={comparison['balance']['maximum_abs_difference']:.4e}, "
                f"volume_balance maxabs={comparison['volume_balance']['maximum_abs_difference']:.4e}"
            )
    result["source_commit"] = commit
    target = directory / "matrix_memory_comparison.json"
    target.write_text(
        json.dumps(result, indent=2) + "\n", encoding="utf-8"
    )
    print(f"PASSED all physical parity and MPI matrix-identity gates: {target}")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", type=Path, required=True)
    args = parser.parse_args()
    compare(args.directory)


if __name__ == "__main__":
    main()
