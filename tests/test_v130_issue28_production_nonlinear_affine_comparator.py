# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Source-only tests for archived native #28 production MPI evidence."""

import json
from pathlib import Path
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "sandbox"))
from compare_v130_issue28_production_nonlinear_affine import (
    check_evidence, EXPECTED_REVISION, SCHEMA,
)


class ProductionNonlinearAffineComparatorTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.base = Path(self.tmp.name)
        self.manifest = {
            "schema": "issue28-production-nonlinear-affine-v1",
            "source": "a" * 40,
            "immutable_revision": EXPECTED_REVISION,
            "opt_in_only": True,
            "physical_steps": 0,
        }
        (self.base / "manifest.json").write_text(
            json.dumps(self.manifest), encoding="utf-8"
        )
        self.results = {}
        for n, final in ((1, 2.199213556173762e-12),
                         (2, 2.3277177796294783e-12)):
            obj = {
                "schema": SCHEMA,
                "status": "PASS_NATIVE_PRODUCTION_INITIALIZER",
                "mpi_ranks": n,
                "model": "original shipped examples/shallow_water.py",
                "mesh": 0.1,
                "dofs": 4270,
                "original_dynamic_h0_p0_bitwise": True,
                "published_algebraic_dofs_modified": True,
                "actual_published_IFunction_Fa_l2": final,
                "explicit_gyro_from_original_shallow_water": True,
                "explicit_gyro_cache_ready": True,
                "main_physical_TS_executed": False,
                "pseudo_initial_TS_executed": False,
                "generic_nonlinear_algebraic_Newton_qualified": False,
                "high_index_rotating_Y_initialization_qualified": False,
                "real_production_init_step_report": {
                    "strategy": "nonlinear_affine",
                    "status": "PASS",
                    "dofs": 4270,
                    "differential_dofs": 2035,
                    "algebraic_dofs": 2235,
                    "differential_initial_exact": True,
                    "no_pseudo_ts_executed": True,
                    "no_physical_ts_executed": True,
                    "accuracy_target_met": True,
                    "algebraic_null_columns": 0,
                    "linear_solver": "PETSc PREONLY + MUMPS",
                    "linear_solver_converged_reason": 4,
                    "max_abs_algebraic_mass_rows": 0.0,
                    "native_initial_algebraic_residual_l2": 1281.781256,
                    "native_final_algebraic_residual_l2": final,
                    "accuracy_target": 1e-11,
                    "affine_algebraic_identity_error_l2": 1.6e-13,
                    "affine_identity_tolerance": 1e-7,
                },
            }
            self.results[n] = obj
            directory = self.base / f"result_{n}"
            directory.mkdir()
            (directory / "production_nonlinear_affine.json").write_text(
                json.dumps(obj), encoding="utf-8"
            )

    def write_result(self, n):
        directory = self.base / f"result_{n}"
        (directory / "production_nonlinear_affine.json").write_text(
            json.dumps(self.results[n]), encoding="utf-8"
        )

    def test_saved_successful_mpi1_2_data_passes(self):
        self.assertEqual(len(check_evidence(self.base)), 3)

    def test_rejects_false_fixed_differential_claim(self):
        self.results[2]["original_dynamic_h0_p0_bitwise"] = False
        self.write_result(2)
        with self.assertRaisesRegex(ValueError, "differential data changed"):
            check_evidence(self.base)

    def test_rejects_mpi2_algebraic_residual_failure(self):
        self.results[2]["actual_published_IFunction_Fa_l2"] = 1e-5
        self.write_result(2)
        with self.assertRaisesRegex(ValueError, "published native residual"):
            check_evidence(self.base)

    def test_rejects_missing_native_rank(self):
        (self.base / "result_2" / "production_nonlinear_affine.json").unlink()
        with self.assertRaisesRegex(ValueError, "Missing native evidence"):
            check_evidence(self.base)

    def test_rejects_wrong_oci_revision(self):
        self.manifest["immutable_revision"] = "0" * 40
        (self.base / "manifest.json").write_text(
            json.dumps(self.manifest), encoding="utf-8"
        )
        with self.assertRaisesRegex(ValueError, "Wrong Docker OCI revision"):
            check_evidence(self.base)


if __name__ == "__main__":
    unittest.main()
