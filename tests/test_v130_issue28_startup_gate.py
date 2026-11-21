# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See LICENSE for license information.

"""Synthetic host regression tests for #28/#6 startup-gate comparator.

No test here is evidence of native GetFEM/PETSc execution.
"""

import ast
import json
from pathlib import Path
import tempfile
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"sandbox"))

from compare_v130_issue28_startup_gate import (
    compare, check, OCI_REVISION,
)


class StartupGateHostTests(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.root = Path(tmp.name)
        (self.root/"manifest.json").write_text(
            json.dumps({
                "schema": "v130-issue28-startup-gate-v1",
                "image_revision": OCI_REVISION,
                "production_code_unchanged": True,
                "accepted_Y0": False,
                "no_physical_time_steps": True,
            }), encoding="utf-8",
        )
        for h in (.2, .1):
            for case in ("compatible_rest_spinup", "incompatible_rotation"):
                for ranks in (1, 2):
                    self.store(h, case, ranks, self.fixture(h, case, ranks))

    @staticmethod
    def fixture(h, case, ranks):
        base = {
            "schema": "scrimp-v130-issue28-startup-gate-v1",
            "nominal_mesh": h, "mpi_ranks": ranks,
            "geometry_order": 2, "height_order": 2,
            "original_differential_h0_p0_preserved_bitwise": True,
            "entire_GetFEM_and_GlobalVec_rollback_bitwise": True,
            "no_pseudo_or_physical_TS": True,
            "physical_time_evolution_qualified": False,
            "production_nonlinear_initializer_qualified": False,
            "physical_rotating_Y0_accepted": False,
        }
        if case == "compatible_rest_spinup":
            a = .0001/.02
            smooth1 = 3*a*a-2*a*a*a
            smooth2 = 3*(2*a)**2-2*(2*a)**3
            ratio = smooth2/smooth1
            base.update({
                "status": "COMPATIBLE_REST_SPINUP_T0_DIAGNOSTIC_ONLY",
                "rotation_after_spinup": .2, "spinup_duration": .02,
                "original_h0": 5., "original_p0_zero": True,
                "original_physical_wall_facet_tangent": True,
                "initial_wall_rotation_exact_zero": True,
                "initial_wall_rotation_derivative_zero": True,
                "initial_algebraic_coenergy_e_h": .5,
                "initial_algebraic_normal_traction_Y_n": -.5,
                "native_IFunction_t0_Fd_l2": 1e-17,
                "native_IFunction_t0_Fa_l2": 1e-17,
                "native_IFunction_t0_full_l2": 1.5e-17,
                "native_IFunction_t0_blocks_l2":
                    {s: 1e-18 for s in (
                        "h", "p", "e_h", "e_p", "U", "U_n", "Y", "Y_n"
                    )},
                "native_F_U_at_delta_l2": 1e-6,
                "native_F_U_at_two_delta_l2": 1e-6 * ratio,
                "native_wall_smoothstep_ratio": ratio,
                "analytic_wall_smoothstep_ratio": ratio,
            })
        else:
            base.update({
                "status": "HIGH_INDEX_OPT_IN_REJECTED_BEFORE_TS",
                "rotation": .2,
                "original_rotating_q0_fixed": True,
                "original_wall_already_rotating_at_t0": True,
                "actual_production_init_step_dispatch": "nonlinear_affine",
                "expected_exception_type": "UnsupportedNonlinearAlgebraicDAE",
                "native_algebraic_null_columns": 8,
                "rejection_reason": (
                    "Algebraic Jacobian has 8 numerically null columns"
                ),
                "production_initialization_report": {
                    "strategy": "nonlinear_affine",
                    "status": "REJECTED",
                    "algebraic_null_columns": 8,
                    "no_pseudo_ts_executed": True,
                    "no_physical_ts_executed": True,
                    "high_index_boundary_multipliers_supported": False,
                },
            })
        return base

    def store(self, h, case, ranks, record):
        folder = self.root/(
            f"result_h{str(h).replace('.', 'p')}_{case}_{ranks}"
        )
        folder.mkdir(exist_ok=True)
        (folder/"startup_gate.json").write_text(
            json.dumps(record), encoding="utf-8",
        )

    def test_python_probe_syntax(self):
        code = (ROOT/"sandbox"/"probe_v130_issue28_startup_gate.py").read_text(
            encoding="utf-8"
        )
        ast.parse(code)

    def test_fixtures_pass_without_claiming_ts(self):
        summary = compare(self.root)
        self.assertEqual(
            summary["status"],
            "NATIVE_STARTUP_SAFETY_AND_T0_FEASIBILITY_DIAGNOSTIC_ONLY",
        )
        self.assertFalse(summary["physical_time_integration_qualified"])

    def test_rejects_false_physical_ts(self):
        p = self.fixture(.2, "compatible_rest_spinup", 1)
        p["physical_time_evolution_qualified"] = True
        self.store(.2, "compatible_rest_spinup", 1, p)
        with self.assertRaisesRegex(ValueError, "Unqualified TS"):
            compare(self.root)

    def test_rejects_wrong_rest_traction(self):
        p = self.fixture(.1, "compatible_rest_spinup", 1)
        p["initial_algebraic_normal_traction_Y_n"] = .5
        self.store(.1, "compatible_rest_spinup", 1, p)
        with self.assertRaisesRegex(ValueError, "analytical rest"):
            compare(self.root)

    def test_rejects_nonzero_rest_algebraic_residual(self):
        p = self.fixture(.1, "compatible_rest_spinup", 2)
        p["native_IFunction_t0_Fa_l2"] = 1e-5
        self.store(.1, "compatible_rest_spinup", 2, p)
        with self.assertRaisesRegex(ValueError, "Incompatible exact rest"):
            compare(self.root)

    def test_rejects_first_order_wall_impulse(self):
        p = self.fixture(.2, "compatible_rest_spinup", 2)
        p["native_wall_smoothstep_ratio"] = 2.
        self.store(.2, "compatible_rest_spinup", 2, p)
        with self.assertRaisesRegex(ValueError, "smoothstep"):
            compare(self.root)

    def test_rejects_fake_supported_high_index(self):
        p = self.fixture(.1, "incompatible_rotation", 1)
        p["production_initialization_report"][
            "high_index_boundary_multipliers_supported"
        ] = True
        self.store(.1, "incompatible_rotation", 1, p)
        with self.assertRaisesRegex(ValueError, "Production rejection"):
            compare(self.root)

    def test_rejects_missing_MPI2_case(self):
        file = (
            self.root/"result_h0p1_incompatible_rotation_2"/"startup_gate.json"
        )
        file.unlink()
        with self.assertRaisesRegex(ValueError, "Missing native startup"):
            compare(self.root)


if __name__ == "__main__":
    unittest.main()
