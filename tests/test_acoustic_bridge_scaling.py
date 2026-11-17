# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Pure JSON regression tests for the cross-method strong-scaling charts.

Run in scrimp-mpi:local without MPI or a GetFEM model:
    python -m unittest -v tests.test_acoustic_bridge_scaling
"""

import csv
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
ANALYSIS = ROOT / "sandbox" / "analyze_acoustic_scaling.py"
SPEC = importlib.util.spec_from_file_location(
    "_scrimp_acoustic_scaling", ANALYSIS
)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def _report(policy, bridge, mpi, repetition, *, commit="abc12345"):
    seconds = {
        ("never", "csr", 1): 10.,
        ("never", "csr", 2): 6.,
        ("never", "coo", 1): 4.,
        ("never", "coo", 2): 3.,
        ("always", "csr", 1): 12.,
        ("always", "csr", 2): 8.,
        ("always", "coo", 1): 6.,
        ("always", "coo", 2): 5.,
    }[policy, bridge, mpi]
    return {
        "series": policy,
        "bridge_method": bridge,
        "bridge_backend": bridge,
        "benchmark_source_commit": commit,
        "getfem_matrix_mode": "distributed",
        "effective_mumps_reset_policy": policy,
        "mpi_size": mpi,
        "degrees_of_freedom": 100,
        "rectangle_mesh_h": 0.1,
        "requested_final_time": 0.04,
        "actual_final_time": 0.04,
        "last_saved_time": 0.04,
        "ts_steps": 20,
        "ts_type": "cn",
        "ts_adapt_type": "none",
        "ts_converged_reason": 2,
        "ts_snes_failures": 0,
        "ts_rejected_steps": 0,
        "final_variable_norms": {"p": 2., "q": 3.},
        "wall_solve_seconds_max": seconds + repetition * 0.1,
        "peak_rank_rss_kib": 409600,
    }


class AcousticBridgePlotTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.files = []
        for policy in ("never", "always"):
            for bridge in ("csr", "coo"):
                for mpi in (1, 2):
                    for repetition in (0, 1):
                        p = self.root / (
                            f"wave_{bridge}_{policy}_{mpi}ranks_"
                            f"{repetition}.json"
                        )
                        p.write_text(
                            json.dumps(
                                _report(policy, bridge, mpi, repetition)
                            ), encoding="utf-8",
                        )
                        self.files.append(p)

    def run_analysis(self):
        MODULE.plot_bridge_comparison(
            self.files,
            self.root / "comparison",
            rtol=1e-5,
            atol=1e-9,
            selected_policies="never,always",
        )

    def test_plot_and_csv_have_independent_rank_1_speedups(self):
        self.run_analysis()
        for suffix in ("_speedup.png", "_times.png", "_coo_gain.png"):
            self.assertTrue(
                (self.root / ("comparison" + suffix)).is_file(),
                suffix,
            )
        with (self.root / "comparison_summary.csv").open(
            newline="", encoding="utf-8"
        ) as f:
            rows = list(csv.DictReader(f))
        self.assertEqual(len(rows), 4)
        sample = next(
            row for row in rows
            if row["mumps_policy"] == "never" and row["mpi_ranks"] == "2"
        )
        self.assertAlmostEqual(
            float(sample["coo_acceleration"]), 6.05 / 3.05
        )
        self.assertAlmostEqual(
            float(sample["csr_parallel_speedup"]), 10.05 / 6.05
        )
        self.assertAlmostEqual(
            float(sample["coo_parallel_speedup"]), 4.05 / 3.05
        )

    def test_missing_bridge_rank_is_rejected(self):
        self.files = [
            p for p in self.files
            if not ("coo_always_2ranks" in p.name)
        ]
        with self.assertRaisesRegex(SystemExit, "Incomplete"):
            self.run_analysis()

    def test_mixed_commits_are_rejected(self):
        p = self.files[0]
        run = json.loads(p.read_text())
        run["benchmark_source_commit"] = "other"
        p.write_text(json.dumps(run))
        with self.assertRaisesRegex(SystemExit, "same Git commit"):
            self.run_analysis()

    def test_cross_method_physical_mismatch_is_rejected(self):
        p = self.files[0]
        run = json.loads(p.read_text())
        run["final_variable_norms"]["p"] = 100.
        p.write_text(json.dumps(run))
        with self.assertRaisesRegex(SystemExit, "physical mismatch"):
            self.run_analysis()

    def test_unbalanced_repetitions_are_rejected(self):
        self.files = self.files[:-1]
        with self.assertRaisesRegex(SystemExit, "Unequal repetition"):
            self.run_analysis()

    def test_unlabelled_legacy_report_is_rejected(self):
        p = self.files[0]
        run = json.loads(p.read_text())
        del run["bridge_method"]
        p.write_text(json.dumps(run))
        with self.assertRaisesRegex(SystemExit, "explicit"):
            self.run_analysis()


if __name__ == "__main__":
    unittest.main()
