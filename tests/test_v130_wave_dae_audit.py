# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Source-only synthetic checks for the scientific DAE MPI comparator.

The fixtures do not claim a physical PETSc/GetFEM result: they verify
that complete-trajectory MPI mismatches fail and initial q/p drift
remains visible even when both MPI trajectories match perfectly.
"""

import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
COMPARATOR = ROOT / "sandbox" / "compare_v130_wave_dae_initial.py"


class DAEInitialComparatorTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.before = np.array([0.0, 1.0, 0.0, 0.0])
        self.after = np.array([0.001, 0.999, 0.002, 0.004])
        self.times = np.array([0.0, 0.002, 0.004])
        self.trajectory = np.stack([
            self.after, self.after + 0.01, self.after + 0.02
        ])
        self._write(1)
        self._write(2)

    def _write(self, rank, trajectory=None):
        report = {
            "mpi_ranks": rank,
            "source_image_revision": "synthetic-test-not-an-image",
            "mesh_h": 0.2,
            "dt": 0.002,
            "final_time": 0.004,
            "cn_steps": 2,
            "dofs": 4,
            "variable_intervals": {
                "q": [0, 1], "p": [1, 2],
                "e_q": [2, 3], "e_p": [3, 4],
            },
            "dynamic_names": ["q", "p"],
            "algebraic_names": ["e_q", "e_p"],
            "fidelity_atol": 1e-8,
            "fidelity_rtol": 1e-8,
            "dynamic_initial_changes": {
                "q": {
                    "absolute_change_l2": 0.001,
                    "preserved_within_tolerance": False,
                },
                "p": {
                    "absolute_change_l2": 0.001,
                    "preserved_within_tolerance": False,
                },
            },
            "algebraic_residual_before": {
                "algebraic_residual_l2": 1e-3,
            },
            "algebraic_residual_after": {
                "algebraic_residual_l2": 1e-12,
            },
            "physical_initial_value_fidelity": "FAIL",
        }
        stem = f"wave_initial_{rank}r"
        (self.root / f"{stem}.json").write_text(
            json.dumps(report), encoding="utf-8"
        )
        np.savez_compressed(
            self.root / f"{stem}.npz",
            before=self.before, after=self.after,
            times=self.times,
            trajectory=(self.trajectory if trajectory is None
                        else trajectory),
        )

    def _run(self, *opts):
        return subprocess.run(
            [sys.executable, str(COMPARATOR), str(self.root), *opts],
            capture_output=True, text=True,
        )

    def test_matching_mpi_trajectories_do_not_hide_physical_drift(self):
        result = self._run()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("PASS: full-state MPI parity", result.stdout)
        self.assertIn("PHYSICAL INITIAL-DATA FIDELITY: FAIL", result.stdout)
        summary = json.loads(
            (self.root / "wave_initial_dae_summary.json").read_text()
        )
        self.assertEqual(summary["mpi_full_state_parity"], "PASS")
        self.assertEqual(summary["physical_initial_value_fidelity"], "FAIL")

    def test_strict_physical_fidelity_gate_fails(self):
        result = self._run("--require-initial-fidelity")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Prescribed q/p initial data", result.stderr)

    def test_a_single_later_time_state_mismatch_fails(self):
        wrong = self.trajectory.copy()
        wrong[-1, 2] += 0.1
        self._write(2, trajectory=wrong)
        result = self._run()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("MPI scaled L2", result.stderr)


if __name__ == "__main__":
    unittest.main()
