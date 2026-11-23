# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Synthetic host negative tests for native spinup first-derivative gate.

No fixture is evidence of a native simulation or of a general DAE initializer.
"""

import ast
import json
from pathlib import Path
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"sandbox"))
from compare_v130_issue28_spinup_first_derivative import (
    compare, EXPECTED_IMAGE, EXPECTED_DIMENSIONS,
)


class SpinupFirstDifferentiatedHostTest(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.root = Path(tmp.name)
        (self.root/"manifest.json").write_text(json.dumps({
            "schema": "v130-issue28-spinup-first-differentiated-v1",
            "image_revision": EXPECTED_IMAGE,
            "physical_spinup_only": True,
            "no_physical_time_steps": True,
            "no_accepted_Y0": True,
        }), encoding="utf-8")
        for h in (.2, .1):
            for rank in (1, 2):
                self.store(h, rank, self.fixture(h, rank))

    @staticmethod
    def fixture(h, rank):
        n, d, regular, y = EXPECTED_DIMENSIONS[h]
        a = .0001/.02
        g = lambda q: 3*q*q-2*q*q*q
        ratio = g(2*a)/g(a)
        return {
            "schema": "scrimp-v130-issue28-spinup-first-differentiated-v1",
            "status": "DIAGNOSTIC_SPINUP_T0_TRUE_RHS_SOLVED_NOT_INITIALIZED",
            "case": "spinup", "nominal_mesh": h, "mpi_ranks": rank,
            "original_P2_geometry": True,
            "rotation_after_spinup": .2, "initial_physical_p0_zero": True,
            "physical_wall_facet_tangent": True, "spinup_duration": .02,
            "analytic_partial_t_Fa_at_t0_exact_zero": True,
            "native_wall_smoothstep_ratio": ratio,
            "analytic_wall_smoothstep_ratio": ratio,
            "native_wall_time_change_U_delta_l2": 1e-5,
            "native_wall_time_change_U_2delta_l2": ratio*1e-5,
            "native_time_change_other_rows_l2": 1e-17,
            "dofs": n, "dynamic_dofs": d,
            "regular_algebraic_dofs": regular, "traction_Y_dofs": y,
            "initial_projection": {
                "type": "known-exact-analytical-spinup-rest",
                "native_final_Fa_l2": 1e-17,
            },
            "native_Fd_at_zero_zdot_before_Y0_l2": 2e-16,
            "native_algebraic_partial_t_Fa_l2": 0.,
            "dynamic_physical_rhs_l2": 2e-16,
            "bordered_dynamic_algebraic_operator_shape": [n, n],
            "bordered_true_rhs_l2": 2e-16,
            "bordered_KSP_reason": 1,
            "bordered_relative_true_rhs_residual": 1e-17,
            "matnest_original_row_and_column_IS_used": True,
            "matnest_physical_rhs_roundtrip_l2": 0.,
            "matnest_solution_reconstructed_by_block_IS": True,
            "candidate_qdot_coefficients": [0.]*d,
            "candidate_regular_algebraic_derivative_coefficients": [0.]*regular,
            "candidate_traction_Y0_coefficients": [0.]*y,
            "candidate_qdot_l2": 0.,
            "candidate_regular_algebraic_derivative_l2": 0.,
            "candidate_traction_Y0_l2": 0.,
            "original_native_IFunction_Fd_with_candidate_derivative_Y0_l2": 1e-18,
            "original_native_IFunction_Fa_with_candidate_Y0_l2": 1e-17,
            "independent_spinup_t0_algebraic_compatibility_qualified": True,
            "physical_rotating_Y0_accepted": False,
            "true_initialization_TS_executed": False,
            "true_physical_TS_executed": False,
            "first_differentiated_rank_proven": False,
            "no_pseudo_or_physical_TS": True,
            "full_original_IFunction_checked": True,
            "fixed_initial_h0_p0_preserved_bitwise": True,
            "original_q0_preserved_bitwise": True,
        }

    def store(self, mesh, ranks, report):
        path = self.root/(
            f"result_h{str(mesh).replace('.', 'p')}_mpi{ranks}"
        )
        path.mkdir(exist_ok=True)
        (path/"spinup_first_derivative.json").write_text(
            json.dumps(report), encoding="utf-8"
        )

    def test_probe_python_syntax(self):
        p = ROOT/"sandbox"/"probe_v130_issue28_spinup_first_derivative.py"
        ast.parse(p.read_text(encoding="utf-8"))

    def test_original_fixture_valid_diagnostic_only(self):
        r = compare(self.root)
        self.assertEqual(
            r["status"], "COMPLETE_DIFFERENTIATED_SPINUP_T0_DIAGNOSTIC_NO_TS"
        )
        self.assertFalse(r["physical_time_integration_qualified"])

    def test_rejects_time_derivative_central_fd_contamination(self):
        r = self.fixture(.2, 1)
        r["native_algebraic_partial_t_Fa_l2"] = 1e-6
        self.store(.2, 1, r)
        with self.assertRaisesRegex(ValueError, "C1 wall starts"):
            compare(self.root)

    def test_rejects_wrong_native_wall_smoothstep(self):
        r = self.fixture(.2, 1)
        r["native_wall_smoothstep_ratio"] = 2.
        self.store(.2, 1, r)
        with self.assertRaisesRegex(ValueError, "C1 smoothstep"):
            compare(self.root)

    def test_rejects_MPI_vector_permutation_with_matching_norms(self):
        # Both ranks report the SAME Y0 norm, so scalar MPI parity
        # passes. Only the full native-vector comparison can detect
        # the inconsistent coefficient locations.
        r1 = self.fixture(.1, 1)
        r2 = self.fixture(.1, 2)
        r1["candidate_traction_Y0_coefficients"][0] = 2e-9
        r2["candidate_traction_Y0_coefficients"][1] = 2e-9
        r1["candidate_traction_Y0_l2"] = 2e-9
        r2["candidate_traction_Y0_l2"] = 2e-9
        self.store(.1, 1, r1)
        self.store(.1, 2, r2)
        with self.assertRaisesRegex(ValueError, "vector MPI parity"):
            compare(self.root)

    def test_rejects_MPI_traction_norm_mismatch(self):
        # Deliberately fail the earlier scalar gate instead. These
        # are two separate invariants and should have separate tests.
        r2 = self.fixture(.1, 2)
        r2["candidate_traction_Y0_coefficients"][0] = 2e-9
        r2["candidate_traction_Y0_l2"] = 2e-9
        self.store(.1, 2, r2)
        with self.assertRaisesRegex(ValueError, "mismatch in Y0"):
            compare(self.root)

    def test_rejects_incompatible_initial_Fa(self):
        r = self.fixture(.2, 1)
        r["initial_projection"]["native_final_Fa_l2"] = 2e-5
        self.store(.2, 1, r)
        with self.assertRaisesRegex(ValueError, "initial algebraic"):
            compare(self.root)

    def test_rejects_false_TS_acceptance(self):
        r = self.fixture(.1, 1)
        r["true_physical_TS_executed"] = True
        self.store(.1, 1, r)
        with self.assertRaisesRegex(ValueError, "False physical"):
            compare(self.root)

    def test_rejects_wrong_original_mesh_layout(self):
        r = self.fixture(.1, 1)
        r["dynamic_dofs"] = 1452
        self.store(.1, 1, r)
        with self.assertRaisesRegex(ValueError, "P2 finite element"):
            compare(self.root)


if __name__ == "__main__":
    unittest.main()
