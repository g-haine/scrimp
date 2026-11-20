# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Pure-numpy fixtures for distributed PETSc DAE scientific comparator.

This is *source-only* acceptance logic, not a replacement for
native MPI, GetFEM's weak residual or PETSc LSQR qualification.
"""

import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
COMPARE = ROOT / "sandbox" / "compare_v130_petsc_mpi_initial.py"


class PETScDAEScientificComparatorTests(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.root = Path(tmp.name)
        self._write()

    def _write(self, *, bad_negative=False, drift_mpi=False,
               fail_balance=False):
        for scenario in ("polynomial_compatible", "published"):
            for ranks in (1, 2):
                positive = scenario == "polynomial_compatible"
                status = "PASS" if positive else (
                    "PASS" if bad_negative else "EXPECTED_REJECTION"
                )
                power = {
                    "absolute_port_work": 5.39545892e-5,
                    "port_work_integral": [0., -5.39545892e-5],
                    "maximum_absolute_balance_defect": (
                        .01 if fail_balance else 5e-15
                    ),
                    "hamiltonian": [1., 1.00005],
                    "step_balance_defects": [5e-15],
                    "boundary_midpoint_powers": {"left":[-.0027]},
                }
                report = {
                    "scenario": scenario,
                    "mpi_ranks": ranks,
                    "image_source_revision": "synthetic-fixture-only",
                    "mesh_h": .2, "dt": .002, "steps": 10,
                    "pseudo_ts_executed": False,
                    "getfem_stiffness_gathered_to_rank_0": False,
                    "ksp_type": "lsqr",
                    "algebraic_mass_row_max_abs": 0.,
                    "differential_state_exactly_preserved": True,
                    "algebraic_state_correction_l2": 1.,
                    "native_algebraic_before_l2": .1625,
                    "native_algebraic_after_l2": (
                        1e-12 if positive else 4.8e-6
                    ),
                    "native_residual_tolerance": 2e-9,
                    "scientific_status": status,
                    "physical_cn_executed": positive,
                    "expected_incompatible_gaussian_rejected": (
                        not positive and not bad_negative
                    ),
                }
                if positive:
                    report.update(
                        power=power,
                        power_balance_tolerance=1e-8,
                        main_initial_vs_projected_l2=0.,
                    )
                    times = np.arange(11)*.002
                    states = np.zeros((11, 4), dtype=float)
                    states[:,1] = 2.
                    states[:,0] = np.arange(11)*.001
                    if drift_mpi and ranks == 2:
                        states[-1, 2] += .1
                    np.savez_compressed(
                        self.root / (
                            f"distributed_{scenario}_{ranks}r.npz"
                        ),
                        times=times, states=states,
                        initial=np.array([0., 2., 0., 0.]),
                        projected=np.array([0., 2., 1., 1.]),
                    )
                (self.root / (
                    f"distributed_{scenario}_{ranks}r.json"
                )).write_text(
                    json.dumps(report), encoding="utf-8"
                )

    def _run(self):
        return subprocess.run(
            [sys.executable, str(COMPARE), str(self.root)],
            text=True, capture_output=True,
        )

    def test_source_fixture_positive_negative_mpi_pass(self):
        r = self._run()
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertIn("EXPECTED_REJECTION", r.stdout)
        summary = json.loads(
            (self.root / "distributed_dae_comparison.json").read_text()
        )
        self.assertEqual(
            summary["full_CN_trajectory_MPI_1_2"], "PASS"
        )

    def test_incompatible_gaussian_must_be_rejected(self):
        self._write(bad_negative=True)
        r = self._run()
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("Gaussian incompatible-data guard", r.stderr)

    def test_later_full_state_mpi_disparity_must_fail(self):
        self._write(drift_mpi=True)
        r = self._run()
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("trajectory MPI disparity", r.stderr)

    def test_nonvacuous_but_unbalanced_power_must_fail(self):
        self._write(fail_balance=True)
        r = self._run()
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("port power or physical start failed", r.stderr)


if __name__ == "__main__":
    unittest.main()
