# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Diagnose existing full-vector mismatch, WITHOUT rerunning any PDE.

The benchmark's "initial" vector is the CN step-zero state *after*
SCRIMP's DAE consistency initialization. It is not necessarily the
raw prescribed initial condition. Sorted coefficient comparisons can
rule out a pure global permutation; matching sorted values are only
suggestive, not proof of FE-degree-of-freedom equivalence.
"""
import argparse
import json
from pathlib import Path

import numpy as np


def _norms(ref, cur):
    if ref.shape != cur.shape:
        raise ValueError("Cannot compare different global vector dimensions")
    ref = np.asarray(ref, dtype=np.float64)
    cur = np.asarray(cur, dtype=np.float64)
    if not np.all(np.isfinite(ref)) or not np.all(np.isfinite(cur)):
        raise ValueError("Nonfinite state coefficients")
    delta = cur - ref
    peak = float(np.max(np.abs(delta)))
    return {
        "max_abs": peak,
        "relative_l2": float(np.linalg.norm(delta)) / max(
            float(np.linalg.norm(ref)),
            float(np.linalg.norm(cur)),
            np.finfo(float).tiny,
        ),
        "scaled_linf": peak / max(
            float(np.max(np.abs(ref))),
            float(np.max(np.abs(cur))),
            np.finfo(float).tiny,
        ),
        "norm_release": float(np.linalg.norm(ref)),
        "norm_current": float(np.linalg.norm(cur)),
        "zero_count_release": int(np.count_nonzero(ref == 0)),
        "zero_count_current": int(np.count_nonzero(cur == 0)),
    }


def analyze(directory):
    root = Path(directory)
    old_metadata_path = root / "release_120_release_native_1r_rep1.json"
    old_metadata = json.loads(
        old_metadata_path.read_text(encoding="utf-8")
    )
    if (old_metadata["variant"] != "release_120"
            or old_metadata["source_tag"] != "v1.2.0"
            or old_metadata["mpi_ranks"] != 1):
        raise ValueError("Missing authentic v1.2.0 serial reference")
    old_path = root / old_metadata["state_file"]
    result = {
        "release_source": old_metadata["source_sha"],
        "release_tag": old_metadata["source_tag"],
        "benchmark_config": {
            key: old_metadata[key] for key in
            ("benchmark", "mesh_h", "dofs", "dt",
             "time_final_requested", "steps")
        },
        "comparisons": {},
        "interpretation": {
            "initial_label": (
                "The CN main TS monitor step zero is AFTER the "
                "possibly different DAE-consistency pseudo-TS."
            ),
            "sorted_coefficients": (
                "Sorting eliminates permutations. Small sorted-error "
                "with large index-error is compatible with a global "
                "DOF permutation, but does not prove one. Large "
                "sorted-error rules out a pure permutation."
            ),
            "release_claim": (
                "Never claim a like-for-like performance gain unless "
                "the release and current physical solutions have "
                "been independently reconciled."
            ),
        },
    }
    with np.load(old_path, allow_pickle=False) as data:
        old_states = {
            key: np.asarray(data[key]) for key in ("initial", "final")
        }
    for policy in ("never", "on_change"):
        for ranks in (1, 2, 4, 6):
            path = root / f"current_{policy}_{ranks}r_rep1.json"
            meta = json.loads(path.read_text(encoding="utf-8"))
            if (
                meta["variant"] != "current"
                or meta["mumps_policy"] != policy
                or meta["mpi_ranks"] != ranks
                or meta["dofs"] != old_metadata["dofs"]
                or meta["steps"] != old_metadata["steps"]
                or meta["mesh_h"] != old_metadata["mesh_h"]
            ):
                raise ValueError(f"Different physical test: {path}")
            with np.load(root / meta["state_file"], allow_pickle=False) as v:
                current = {
                    key: np.asarray(v[key])
                    for key in ("initial", "final")
                }
            report = {}
            for key in ("initial", "final"):
                observed = _norms(old_states[key], current[key])
                sorted_error = _norms(
                    np.sort(old_states[key]), np.sort(current[key])
                )
                report[key] = {
                    "raw": observed,
                    "after_global_sort": sorted_error,
                    "pure_permutation_compatible": (
                        sorted_error["relative_l2"] < 1e-10
                        and sorted_error["scaled_linf"] < 1e-10
                    ),
                }
            result["comparisons"][f"{policy}_{ranks}r"] = report
            print(
                f"{policy:9s} {ranks}r | "
                f"t0 relL2={report['initial']['raw']['relative_l2']:.3e} "
                f"sorted={report['initial']['after_global_sort']['relative_l2']:.3e} | "
                f"t_f relL2={report['final']['raw']['relative_l2']:.3e} "
                f"sorted={report['final']['after_global_sort']['relative_l2']:.3e}"
            )
    result_path = root / "release_state_diagnostic.json"
    result_path.write_text(
        json.dumps(result, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print("Saved:", result_path)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", type=Path, required=True)
    args = parser.parse_args()
    analyze(args.directory)


if __name__ == "__main__":
    main()
