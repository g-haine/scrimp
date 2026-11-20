# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Synthetic source-only tests for complete DAE initializers' physics gate.

These do NOT replace native GetFEM/PETSc. Ensure comparisons never
confuse MPI equivalence, physical initial-data fidelity and a
non-vacuous midpoint Hamiltonian/port-work identity.
"""

import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "sandbox" / "compare_v130_wave_initializers.py"


class InitializerComparisonTest(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.folder = Path(self.temporary.name)
        self.prescribed = np.array([0.0, 2.0, 0.0, 0.0])
        self.pseudo = np.array([0.002, 2.00002, 0.0, 0.0])
        self.projected = np.array([0.0, 2.0, 0.3, -0.2])
        self._write_all()

    def _record(self, mode, ranks, *, corrupt_trajectory=False,
                power_fail=False, bad_projected_q=False):
        start = (self.pseudo if mode == "pseudo"
                 else self.projected.copy())
        if mode == "algebraic" and bad_projected_q:
            start[0] += .01
        times = np.arange(3)*.002
        z = np.array([start + k*np.array([.001, -.005, 0.0, 0.0])
                      for k in range(len(times))])
        if corrupt_trajectory:
            z[-1, 2] += .1
        drift = {
            "q": float(abs(start[0]-self.prescribed[0])),
            "p": float(abs(start[1]-self.prescribed[1])),
        }
        report = {
            "schema": "scrimp-v130-wave-initializer-physics-v1",
            "image_source_revision": "synthetic-example-only",
            "initializer": mode,
            "mpi_ranks": ranks,
            "mesh_h": .2, "dt": .002, "steps": 2,
            "final_time": .004, "drive_left": "0.08*sin(20*t)",
            "initial_p0_expression": "x*(3-x)", "dofs": 4,
            "main_ts_converged_reason": 2,
            "main_ts_rejected_steps": 0,
            "main_ts_snes_failures": 0,
            "main_initial_dynamic_drift_l2": drift,
            "algebraic_projection": (
                None if mode == "pseudo"
                else {
                    "native_algebraic_tolerance": 1e-9,
                    "algebraic_residual_after_project_l2": 1e-12,
                }
            ),
            "nonvacuous_port_power": True,
            "physical_power_balance_within_tolerance": not power_fail,
            "absolute_balance_tolerance": 1e-9,
            "power": {
                "times": times.tolist(),
                "hamiltonian": [1.0, 1.02, 1.04],
                "port_work_integral": [0.0, -.02, -.04],
                "step_balance_defects": [0.0, 0.0],
                "boundary_midpoint_powers": {
                    "left": [-10.0, -10.0],
                },
                "absolute_port_work": .04,
                "maximum_absolute_balance_defect": (
                    .01 if power_fail else 0.0
                ),
                "energy_change": .04,
            },
        }
        stem = f"wave_{mode}_{ranks}r"
        (self.folder / f"{stem}.json").write_text(
            json.dumps(report), encoding="utf-8"
        )
        np.savez_compressed(
            self.folder / f"{stem}.npz",
            times=times, trajectory=z,
            prescribed_initial=self.prescribed,
            pseudo_accepted=self.pseudo,
            main_initializer=start,
        )

    def _write_all(self, **overrides):
        for mode in ("pseudo", "algebraic"):
            for ranks in (1, 2):
                options = overrides if (
                    mode == overrides.get("mode")
                    and ranks == overrides.get("ranks")
                ) else {}
                self._record(mode, ranks, **{
                    name: value for name, value in options.items()
                    if name not in ("mode", "ranks")
                })

    def _run(self):
        return subprocess.run(
            [sys.executable, str(SCRIPT), str(self.folder)],
            text=True, capture_output=True,
        )

    def test_all_mpi_and_science_gates_pass_except_known_pseudo_drift(self):
        result = self._run()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("projected=PASS; production pseudo=FAIL", result.stdout)
        summary = json.loads(
            (self.folder / "wave_initializer_physics_comparison.json"
             ).read_text(encoding="utf-8")
        )
        self.assertTrue(summary["both_initializers_power_balance_pass"])
        self.assertEqual(
            summary["diagnostic_projected_dynamic_initial_fidelity"],
            "PASS",
        )

    def test_any_later_trajectory_mpi_difference_fails(self):
        self._write_all(
            mode="algebraic", ranks=2, corrupt_trajectory=True
        )
        result = self._run()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("trajectory MPI parity fails", result.stderr)

    def test_bad_projected_initial_physics_fails(self):
        self._write_all(
            mode="algebraic", ranks=2, bad_projected_q=True
        )
        result = self._run()
        self.assertNotEqual(result.returncode, 0)

    def test_nonvacuous_but_unbalanced_energy_fails(self):
        self._write_all(
            mode="pseudo", ranks=2, power_fail=True
        )
        result = self._run()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("power balance failed", result.stderr)


if __name__ == "__main__":
    unittest.main()
