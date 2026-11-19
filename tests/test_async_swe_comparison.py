# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Validate the strict numerical-comparison logic without a native MPI run."""
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "_scrimp_swe_storage_compare",
    ROOT / "sandbox" / "compare_async_swe.py",
)
ANALYZE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(ANALYZE)

FIELDS = (
    "t", "kinetic", "potential", "H", "d_grad", "d_div",
    "p_normal", "p_viscous", "supplied", "dissipated",
    "balance", "relative_balance", "relative_balance_exchange",
    "volume", "normal_volume_outflow", "cumulative_outflow",
    "volume_balance", "relative_volume_balance", "lyapunov",
)


class SWEOutputParityTests(unittest.TestCase):
    def setUp(self):
        directory = tempfile.TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        self.directory = Path(directory.name)
        (self.directory / "binary").mkdir()
        (self.directory / "binary" / "trajectory_manifest.json").write_text(
            json.dumps({
                "ranks": 2, "global_size": 3, "number_of_states": 3,
                "partitions": [
                    {"rank": 0, "lo": 0, "hi": 2},
                    {"rank": 1, "lo": 2, "hi": 3},
                ],
            })
        )
        data = {"states": np.array([
            [5., 0., 0.],
            [5., .001, .002],
            [5., .003, .004],
        ])}
        for i, name in enumerate(FIELDS):
            data[name] = (
                np.array([0., .001, .002]) if name == "t"
                else np.array([1. + i, 1. + i + 1e-6, 1. + i + 2e-6])
            )
        for backend in ("memory", "async_binary"):
            np.savez_compressed(
                self.directory / f"{backend}_swe.npz", **data
            )
            report = {
                "mode": backend, "source_commit": "cafefeed",
                "ranks": 2, "dofs": 3, "snapshots": 3,
                "mesh": .35, "dt": .001, "final_time": .002,
                "spinup_duration": .03, "gamma": 0., "ts_type": "cn",
                "ts_adapt_type": "none", "mumps_reset_policy": "on_change",
                "minimum_water_height": 5.,
                "energy_balance_defect_final": 1e-6,
                "volume_balance_defect_final": 2e-6,
                "convergence": {
                    "reason": 2, "steps": 2, "rejections": 0,
                    "snes_failures": 0, "t_final": .002,
                },
            }
            (self.directory / f"{backend}_swe.json").write_text(
                json.dumps(report), encoding="utf-8"
            )

    def test_matching_balances_and_full_states_pass(self):
        ANALYZE.compare(self.directory)
        self.assertTrue(
            (self.directory / "comparison_swe.json").is_file()
        )

    def test_modified_volume_balance_is_rejected(self):
        file = self.directory / "async_binary_swe.npz"
        with np.load(file, allow_pickle=False) as arrays:
            data = {key: arrays[key] for key in arrays.files}
        data["volume_balance"] = np.array([0.1, 0.2, 0.3])
        np.savez_compressed(file, **data)
        with self.assertRaisesRegex(ValueError, "volume_balance"):
            ANALYZE.compare(self.directory)

    def test_incomplete_MPI_manifest_is_rejected(self):
        file = self.directory / "binary" / "trajectory_manifest.json"
        info = json.loads(file.read_text())
        info["number_of_states"] = 2
        file.write_text(json.dumps(info))
        with self.assertRaisesRegex(ValueError, "manifest"):
            ANALYZE.compare(self.directory)

    def test_unconverged_nonlinear_solve_is_rejected(self):
        file = self.directory / "memory_swe.json"
        data = json.loads(file.read_text())
        data["convergence"]["reason"] = -4
        file.write_text(json.dumps(data))
        with self.assertRaisesRegex(ValueError, "convergence"):
            ANALYZE.compare(self.directory)


if __name__ == "__main__":
    unittest.main()
