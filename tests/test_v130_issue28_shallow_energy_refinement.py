# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Pure-Python regression tests for #28 archived temporal refinement."""

import json
from pathlib import Path
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "sandbox"))

from compare_v130_issue28_shallow_energy_refinement import (
    CASES, EXPECTED_IMAGE_REV, T, check,
)


class Issue28EnergyRefinementComparatorTest(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.root = Path(temp.name)
        self.data = {}
        manifest = {
            "schema": "issue28-shallow-energy-time-refinement-v1",
            "source": "a"*40,
            "image_revision": EXPECTED_IMAGE_REV,
            "original_physical_model": True,
            "spatial_mesh_step": .1,
            "time_final": T,
            "true_physical_step_counts": [4, 8, 16],
            "exact_discrete_power_conservation_claimed": False,
        }
        (self.root/"manifest.json").write_text(
            json.dumps(manifest), encoding="utf-8"
        )
        for k, (steps, dt) in enumerate(CASES):
            for ranks in (1, 2):
                defect = -0.004/(2**k)
                # Tiny rank-2 differences challenge MPI comparator
                # but are far below its requested tolerance.
                if ranks == 2:
                    defect += 1e-11
                h0 = 31636.2594029492
                mass = [2.499, 2.499, 2.499]
                energies = [h0, h0+defect/2, h0+defect]
                record = {
                    "schema": "scrimp-v130-issue28-original-shallow-short-physical-v1",
                    "status": "PASS_SHORT_PHYSICAL_TIME_GATE",
                    "model": "original shipped examples/shallow_water.py",
                    "mpi_ranks": ranks,
                    "mesh_h": .1,
                    "original_prescribed_h0_p0_bitwise": True,
                    "original_boundary_controls_unchanged_zero": True,
                    "original_explicit_gyro_unchanged": True,
                    "production_initializer_report": {
                        "strategy": "nonlinear_affine",
                        "status": "PASS",
                        "differential_initial_exact": True,
                        "accuracy_target_met": True,
                        "algebraic_null_columns": 0,
                        "native_final_algebraic_residual_l2": 2e-12,
                    },
                    "expected_accepted_steps": steps,
                    "accepted_steps": steps,
                    "rejected_steps": 0,
                    "nonlinear_solver_failures": 0,
                    "ts_converged_reason": 1,
                    "time_step_requested": dt,
                    "t_final_actual": T,
                    "t_final_requested": T,
                    "saved_times": [0., T/2, T],
                    "saved_sample_count": 3,
                    "hamiltonian": energies,
                    "mass": mass,
                    "total_boundary_power": [0.]*3,
                    "cumulative_boundary_power": [0.]*3,
                    "energy_balance_defect": [
                        0., defect/2, defect,
                    ],
                    "signed_final_energy_balance_defect": defect,
                    "relative_mass_drift": 0.,
                    "relative_hamiltonian_drift": abs(defect)/h0,
                    "maximum_abs_boundary_power": 0.,
                    "Hamiltonian_GetFEM_and_registered_agree": True,
                    "BDF2_exact_power_conservation_claimed": False,
                    "original_rotating_viscous_Y0_qualified": False,
                    "physical_sanity_gate_failures": [],
                }
                d = self.root/f"result_steps{steps}_{ranks}"
                d.mkdir()
                self.data[steps, ranks] = record
                self.write(steps, ranks)

    def write(self, steps, ranks):
        path = self.root/f"result_steps{steps}_{ranks}"/"shallow_short_physical.json"
        path.write_text(
            json.dumps(self.data[steps, ranks]), encoding="utf-8"
        )

    def test_qualified_six_run_campaign(self):
        result = check(self.root)
        self.assertEqual(result["status"], "PASS_NATIVE_TIME_REFINEMENT_MEASUREMENTS")
        self.assertEqual(result["native_true_accepted_step_counts"], [4, 8, 16])
        self.assertEqual(len(result["runs"]), 3)
        self.assertTrue(result["observed_three_point_abs_defect_decreasing"])
        self.assertFalse(result["energy_defect_convergence_proven"])
        self.assertFalse(result["exact_BDF2_power_balance_claimed"])

    def test_rejects_wrong_accepted_step_count(self):
        self.data[8, 2]["accepted_steps"] = 7
        self.write(8, 2)
        with self.assertRaisesRegex(ValueError, "PETSc TS failed"):
            check(self.root)

    def test_rejects_nonzero_external_work(self):
        self.data[16, 1]["cumulative_boundary_power"][-1] = 0.1
        self.write(16, 1)
        with self.assertRaises(ValueError):
            check(self.root)

    def test_rejects_wrong_original_differential_data(self):
        self.data[4, 1]["original_prescribed_h0_p0_bitwise"] = False
        self.write(4, 1)
        with self.assertRaisesRegex(ValueError, "Original physical model changed"):
            check(self.root)

    def test_rejects_wrong_final_time(self):
        self.data[16, 2]["t_final_actual"] = T/2
        self.write(16, 2)
        with self.assertRaisesRegex(ValueError, "fixed T"):
            check(self.root)


if __name__ == "__main__":
    unittest.main()
