# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Compare COMPLETE short CN trajectories, Hamiltonians and port powers."""
import argparse
import json
from pathlib import Path

import numpy as np


def _relative_norms(actual, expected):
    delta = np.asarray(actual) - np.asarray(expected)
    abs_error = float(np.max(np.abs(delta))) if delta.size else 0.
    scale_l2 = max(
        float(np.linalg.norm(actual)), float(np.linalg.norm(expected)),
        np.finfo(float).tiny,
    )
    scale_linf = max(
        float(np.max(np.abs(actual))),
        float(np.max(np.abs(expected))),
        np.finfo(float).tiny,
    ) if np.size(actual) else 1.
    return abs_error, float(np.linalg.norm(delta)) / scale_l2, abs_error / scale_linf


def compare(directory):
    directory = Path(directory)
    reports = {
        mode: json.loads(
            (directory / f"{mode}_physical.json").read_text()
        )
        for mode in ("memory", "async_binary")
    }
    memory = reports["memory"]
    disk = reports["async_binary"]
    expected = (
        "source_commit", "mpi_size", "dofs", "n_snapshots",
        "steps", "rejected", "snes_failures", "final_time",
        "power_names", "snapshot_dtype",
    )
    for key in expected:
        if memory[key] != disk[key]:
            raise ValueError(f"Physically inconsistent benchmark {key}: "
                             f"{memory[key]!r} versus {disk[key]!r}")
    if (memory["rejected"] or memory["snes_failures"]
            or memory["n_snapshots"] < 3):
        raise ValueError("Acoustic benchmark failed convergence/save checks")

    metrics = {}
    with np.load(directory / "memory_physical.npz",
                 allow_pickle=False) as m, np.load(
                     directory / "async_binary_physical.npz",
                     allow_pickle=False) as d:
        if set(m.files) != set(d.files):
            raise ValueError("Different physical variables in backends")
        for key in m.files:
            left, right = np.asarray(m[key]), np.asarray(d[key])
            if left.shape != right.shape:
                raise ValueError(f"{key}: shape mismatch")
            if not np.all(np.isfinite(left)) or not np.all(np.isfinite(right)):
                raise ValueError(f"{key}: nonfinite values")
            abs_error, rel_l2, rel_linf = _relative_norms(right, left)
            if key == "times":
                # The identical fixed-CN grid and save scheduler should
                # produce the same exact timestamp representation.
                if not np.array_equal(left, right):
                    raise ValueError("Saved timestamp mismatch")
            elif rel_l2 > 1e-9 or rel_linf > 1e-9:
                raise ValueError(
                    f"{key}: different numerical trajectories "
                    f"max_abs={abs_error:.4e}, relL2={rel_l2:.4e}, "
                    f"relLinf={rel_linf:.4e}"
                )
            metrics[key] = {
                "max_abs": abs_error, "relative_l2": rel_l2,
                "relative_linf": rel_linf,
            }
    manifest = directory / "binary" / "trajectory_manifest.json"
    archive = json.loads(manifest.read_text())
    if (archive["number_of_states"] != memory["n_snapshots"]
            or archive["global_size"] != memory["dofs"]
            or archive["ranks"] != memory["mpi_size"]):
        raise ValueError("Incomplete or incorrect binary manifest")

    summary = {
        "commit": memory["source_commit"],
        "dofs": memory["dofs"],
        "mpi_ranks": memory["mpi_size"],
        "steps": memory["steps"],
        "snapshots": memory["n_snapshots"],
        "backends": ["memory", "async_binary"],
        "physical_parity": metrics,
        "binary_partitions": archive["partitions"],
        "notes": [
            "Same fixed CN, MUMPS and GetFEM model; all saved state "
            "coefficients, Hamiltonian and power values compared.",
            "The short smoke materializes the trajectory ONLY in "
            "the testing driver; production streaming remains lazy.",
            "A writer must be running before mpiexec. The backend "
            "never starts a child process from an initialized MPI rank.",
        ],
    }
    (directory / "comparison.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(
        f"PASS: memory and async_binary parity on "
        f"{memory['dofs']} DDL / {memory['mpi_size']} MPI ranks / "
        f"{memory['steps']} CN steps / {memory['n_snapshots']} states"
    )
    for key, data in metrics.items():
        print(
            f"{key}: max_abs={data['max_abs']:.4e}, "
            f"relative_L2={data['relative_l2']:.4e}, "
            f"relative_Linf={data['relative_linf']:.4e}"
        )
    print(f"Saved: {directory / 'comparison.json'}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", required=True, type=Path)
    args = parser.parse_args()
    compare(args.directory)


if __name__ == "__main__":
    main()
