# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Synthetic regressions of the native acoustic memory-report analyzer.

No PETSc/GetFEM model required; run in scrimp-mpi:local:
    python -m unittest -v tests.test_acoustic_memory_analysis
"""
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "_scrimp_acoustic_memory_analysis",
    ROOT / "sandbox" / "analyze_acoustic_memory.py",
)
ANALYZER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(ANALYZER)


class AcousticMemoryAnalysisTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.path = Path(self.temp.name)
        for mode in ("dense", "sparse"):
            for rank in (0, 1):
                stem = self.path / f"memory_{mode}_rank{rank}"
                samples = [
                    {
                        "phase": "main_monitor",
                        "step": step, "simulation_time": 0.001 * step,
                        "VmRSS_kib": 102400 + 128 * step
                            + (800 if mode == "dense" else 0),
                        "VmHWM_kib": 102400 + 128 * step
                            + (800 if mode == "dense" else 0),
                        "RssAnon_kib": 90000 + 120 * step,
                        "stored_state_bytes": (
                            (step + 1) * 8 * 3
                            if mode == "dense" else 8 * 3
                        ),
                        "stored_snapshots": step + 1
                            if mode == "dense" else 1,
                    }
                    for step in (0, 5, 10)
                ]
                samples.append({
                    "phase": "after_solve", "step": -1,
                    "simulation_time": 0.01,
                    "VmRSS_kib": 105000,
                    "VmHWM_kib": 105000,
                    "RssAnon_kib": 92000,
                    "stored_state_bytes": (
                        11 * 24 if mode == "dense" else 2 * 24
                    ),
                    "stored_snapshots": 11 if mode == "dense" else 2,
                })
                report = {
                    "save_mode": mode,
                    "bridge_method": "coo",
                    "getfem_matrix_mode": "distributed",
                    "mumps_reset_policy": "never",
                    "source_commit": "abcdef12",
                    "mpi_rank": rank,
                    "mpi_size": 2,
                    "dofs": 3,
                    "mesh": 0.04,
                    "dt": 0.001,
                    "final_time": 0.01,
                    "expected_steps": 10,
                    "sample_every": 5,
                    "ts_convergence": {
                        "reason": 2, "steps": 10,
                        "rejections": 0, "snes_failures": 0,
                        "last_ts_time": 0.01,
                    },
                    "last_saved_time": 0.01,
                    "snapshots": 11 if mode == "dense" else 2,
                    "stored_state_bytes":
                        11 * 24 if mode == "dense" else 2 * 24,
                    "solve_wall_seconds_rank": 0.123,
                    "mat_info_after_solve_local": {
                        "mass": {"memory": 1024.0}
                    },
                    "samples": samples,
                }
                stem.with_suffix(".json").write_text(
                    json.dumps(report), encoding="utf-8"
                )
                np.savez_compressed(
                    str(stem) + "_state.npz",
                    initial=np.array([0., 1., 2.]),
                    final=np.array([3., 4., 5.]),
                    first_saved_time=0.,
                    final_time=0.01,
                )

    def test_valid_reports_generate_plots_and_csv(self):
        summary = ANALYZER.analyze(self.path)
        self.assertEqual(summary["mpi_size"], 2)
        self.assertEqual(summary["accepted_steps"], 10)
        self.assertEqual(
            summary["max_initial_final_state_abs_difference"], 0.,
        )
        for suffix in (
            "_rss.png", "_saved_payload.png", "_timeseries.csv",
            "_rank_summary.csv", "_summary.json",
        ):
            self.assertTrue(
                (self.path / ("acoustic_memory" + suffix)).is_file(),
                suffix,
            )

    def test_roundoff_level_final_difference_is_reported(self):
        """Do not mistake tiny global endpoint roundoff for divergence.

        PETSc's end-of-TS state can differ from the last TS.monitor
        sample when the save cadence changes. The archived 250-step
        campaign has relative L2 3.18e-11 and scaled Linf 2.60e-11.
        Accept and REPORT this under the independent 1e-10 gates.
        """
        stem = self.path / "memory_sparse_rank1_state.npz"
        # Match across *both* MPI ranks: only save schedule may differ.
        for rank in (0, 1):
            dest = self.path / f"memory_sparse_rank{rank}_state.npz"
            np.savez_compressed(
                dest, initial=np.array([0., 1., 2.]),
                final=np.array([3. + 2e-10, 4., 5.]),
                first_saved_time=0.,
                final_time=0.01,
            )
        # For this small vector 2e-10 / ||state|| ~= 3e-11.
        summary = ANALYZER.analyze(self.path)
        self.assertAlmostEqual(
            summary["max_initial_final_state_abs_difference"],
            2e-10, delta=1e-14,
        )
        self.assertGreater(
            summary["max_initial_final_state_relative_l2"], 0.
        )
        self.assertLess(
            summary["max_initial_final_state_relative_l2"], 1e-10
        )

    def test_saved_state_mismatch_rejected(self):
        stem = self.path / "memory_dense_rank1_state.npz"
        np.savez_compressed(
            stem, initial=np.array([0., 1., 2.]),
            final=np.array([999., 4., 5.]),
            final_time=0.01,
        )
        with self.assertRaisesRegex(ValueError, "final state differs"):
            ANALYZER.analyze(self.path)

    def test_missing_rank_rejected(self):
        (self.path / "memory_sparse_rank1.json").unlink()
        with self.assertRaisesRegex(ValueError, "expected ranks"):
            ANALYZER.analyze(self.path)

    def test_mixed_commit_rejected(self):
        file = self.path / "memory_sparse_rank1.json"
        data = json.loads(file.read_text())
        data["source_commit"] = "different"
        file.write_text(json.dumps(data))
        with self.assertRaisesRegex(ValueError, "same commit"):
            ANALYZER.analyze(self.path)

    def test_solver_failure_rejected(self):
        file = self.path / "memory_dense_rank0.json"
        data = json.loads(file.read_text())
        data["ts_convergence"]["rejections"] = 1
        file.write_text(json.dumps(data))
        with self.assertRaisesRegex(ValueError, "convergence"):
            ANALYZER.analyze(self.path)


if __name__ == "__main__":
    unittest.main()
