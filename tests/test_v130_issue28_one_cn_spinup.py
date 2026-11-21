# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See LICENSE for license information.

"""Synthetic host checks for REAL one-step CN native gate.

These fixtures are NOT native SCRIMP evidence. They verify that the
independent comparator rejects false TS success, wrong power bookkeeping
and MPI-mismatched FULL initial/final state coefficient vectors.
"""

import ast
import json
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"sandbox"))
from compare_v130_issue28_one_cn_spinup import compare, OCI_REV


class OneCNHostTests(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.root = Path(tmp.name)
        (self.root/"manifest.json").write_text(json.dumps({
            "schema": "v130-issue28-one-cn-spinup-v1",
            "image_revision": OCI_REV,
            "physical_spinup_only": True,
            "production_source_unchanged": True,
            "one_physical_step_expected": True,
        }), encoding="utf-8")
        for mesh in (.2, .1):
            for ranks in (1, 2):
                self.store(mesh, ranks, self.fixture(mesh, ranks))

    @staticmethod
    def fixture(mesh, ranks):
        change = 1.2e-8
        work = 2e-8
        diss = 1e-8
        return {
            "schema": "scrimp-v130-issue28-one-cn-spinup-v1",
            "status": "NATIVE_SINGLE_PHYSICAL_CN_STEP_DIAGNOSTIC_ONLY",
            "source_mesh": mesh, "mpi_ranks": ranks,
            "new_physical_rest_spinup_experiment": True,
            "original_incompatible_rotation_initial_state_unchanged": True,
            "geometry_order": 2, "height_order": 2,
            "physical_facet_tangent_wall": True,
            "rotation_after_spinup": .2, "spinup_duration": .02,
            "physical_scheme": "cn", "physical_time_step": .002,
            "initial": {
                "strategy": "sandbox_exact_stationary_spinup_only",
                "status": "SANDBOX_ONLY_PUBLISHED_NOT_GENERAL_INITIALIZER",
                "initial_native_IFunction": {
                    "algebraic_l2": 1e-17, "dynamic_l2": 1e-16
                },
                "differential_q0_preserved_bitwise": True,
                "getfem_globalvec_published_bitwise": True,
                "no_pseudo_or_BDF_initialization": True,
                "analytic_gprime_t0_exact_zero": True,
                "physical_time_step_not_yet_accepted": True,
            },
            "ts_converged_reason": 2,
            "ts_accepted_steps": 1,
            "ts_rejected_steps": 0,
            "ts_snes_failures": 0,
            "ts_actual_final_time": .002,
            "original_IFunction_midpoint_l2": {
                "full_l2": 2e-12, "algebraic_l2": 1e-12,
                "dynamic_l2": 1e-12,
            },
            "original_IFunction_endpoint_l2": {
                "full_l2": 1e-4, "algebraic_l2": 1e-4,
                "dynamic_l2": 1e-5, "U_l2": 1e-4,
            },
            "two_endpoint_states_archived": True,
            "native_original_midpoint_equations_verified": True,
            "first_physical_state_changed": True,
            "dynamic_step_change_l2": 1e-3,
            "first_step_minimum_height_coefficient": 4.9,
            "hamiltonian_initial": 4.,
            "hamiltonian_final": 4.+change,
            "hamiltonian_change": change,
            "physical_boundary_normal_work": 0.,
            "physical_boundary_viscous_work": work,
            "integrated_boundary_work": work,
            "integrated_viscous_dissipation": diss,
            "integrated_dissipation_grad": .8e-8,
            "integrated_dissipation_div": .2e-8,
            "physical_power_balance_defect": change+diss-work,
            "physical_power_balance_defect_over_actual_exchange":
                (change+diss-work)/max(abs(change), work, diss, 1e-14),
            "physical_volume_balance_defect": 0.,
            "energy_balance_convergence_qualified": False,
            "physical_multistep_trajectory_qualified": False,
            "production_high_index_initializer_qualified": False,
        }

    def store(self, mesh, ranks, doc, *, bad_final=False):
        folder = self.root/(
            f"result_h{str(mesh).replace('.', 'p')}_mpi{ranks}"
        )
        folder.mkdir(exist_ok=True)
        (folder/"one_cn_spinup.json").write_text(
            json.dumps(doc), encoding="utf-8"
        )
        original = np.array([5., 5., 0., 0., .5, -.5])
        final = original.copy()
        final[2] = 1e-3
        if bad_final:
            final[3] = 1e-3
        np.savez_compressed(
            folder/"one_cn_spinup_endpoints.npz",
            times=np.array([0., .002]), initial=original, final=final,
            H=np.array([4., 4.+1.2e-8]),
            p_normal=np.zeros(2), p_viscous=np.zeros(2),
            d_grad=np.zeros(2), d_div=np.zeros(2),
        )

    def test_probe_syntax(self):
        path = ROOT/"sandbox"/"probe_v130_issue28_one_cn_spinup.py"
        ast.parse(path.read_text(encoding="utf-8"))

    def test_positive_fixture_does_not_certify_energy(self):
        r = compare(self.root)
        self.assertEqual(
            r["status"],
            "ONE_NATIVE_CN_STEP_ACCEPTED_POWER_BALANCE_UNQUALIFIED",
        )
        self.assertFalse(r["physical_energy_balance_refinement_qualified"])

    def test_rejects_false_ts_step_count(self):
        record = self.fixture(.2, 1)
        record["ts_accepted_steps"] = 0
        self.store(.2, 1, record)
        with self.assertRaisesRegex(ValueError, "exactly one physical"):
            compare(self.root)

    def test_rejects_initial_q0_mutation(self):
        record = self.fixture(.2, 1)
        record["initial"]["differential_q0_preserved_bitwise"] = False
        self.store(.2, 1, record)
        with self.assertRaisesRegex(ValueError, "Unsafe"):
            compare(self.root)

    def test_rejects_bad_original_midpoint_IFunction(self):
        record = self.fixture(.1, 2)
        record["original_IFunction_midpoint_l2"]["full_l2"] = 1e-2
        self.store(.1, 2, record)
        with self.assertRaisesRegex(ValueError, "midpoint GetFEM"):
            compare(self.root)

    def test_rejects_corrupted_boundary_power_bookkeeping(self):
        record = self.fixture(.2, 1)
        record["physical_boundary_viscous_work"] += 1e-4
        self.store(.2, 1, record)
        with self.assertRaisesRegex(ValueError, "weighted port work"):
            compare(self.root)

    def test_rejects_corrupted_power_balance_bookkeeping(self):
        record = self.fixture(.1, 1)
        record["physical_power_balance_defect"] = 1e-4
        self.store(.1, 1, record)
        with self.assertRaisesRegex(ValueError, "cannot be reproduced"):
            compare(self.root)

    def test_rejects_full_vector_mpi_mismatch_at_equal_scalar_reports(self):
        # Both reports have identical scalar metrics but one native
        # coefficient moves to a different DOF on the second rank.
        self.store(.1, 2, self.fixture(.1, 2), bad_final=True)
        with self.assertRaisesRegex(ValueError, "full final coefficient"):
            compare(self.root)

    def test_rejects_missing_actual_native_endpoint_file(self):
        folder = self.root/"result_h0p2_mpi2"/"one_cn_spinup_endpoints.npz"
        folder.unlink()
        with self.assertRaisesRegex(ValueError, "Missing original native"):
            compare(self.root)


if __name__ == "__main__":
    unittest.main()
