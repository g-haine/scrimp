# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Verify byte-exact native MPI partitioned trajectory reconstruction."""
import argparse
from pathlib import Path

import numpy as np

from async_binary_trajectory import LazyTrajectory


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--directory", type=Path, required=True)
    parser.add_argument("--ranks", type=int, required=True)
    parser.add_argument("--n-global", type=int, default=17)
    args = parser.parse_args()

    reader = LazyTrajectory(args.directory / "trajectory_manifest.json")
    if reader.global_size != args.n_global or len(reader) != 4:
        raise RuntimeError("Wrong reconstructed trajectory dimensions")
    expected_steps = np.array([0, 1, 4, 5], dtype=np.int64)
    np.testing.assert_array_equal(reader.steps, expected_steps)
    np.testing.assert_array_equal(reader.times, 0.002 * expected_steps)
    for index, step in enumerate(expected_steps):
        expected = 100. * step + np.arange(
            args.n_global, dtype=reader.dtype
        )
        np.testing.assert_array_equal(reader[index], expected)
    print(
        f"PASS: {len(reader)} exact binary global snapshots reconstructed "
        f"from {args.ranks} MPI partitions; no global state saved on any rank",
        flush=True,
    )


if __name__ == "__main__":
    main()
