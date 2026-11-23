# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Host-only negative controls for native exact algebraic elimination #28."""

import ast
import json
from pathlib import Path
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"sandbox"))
from compare_v130_issue28_exact_weak_elimination import (
    compare, EXPECTED_IMAGE, EXPECTED_WEAK_MOMENT_NORMS,
)


class NativeExactWeakEliminationComparatorTest(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.root = Path(tmp.name)
        (self.root/"manifest.json").write_text(json.dumps({
            "schema": "v130-issue28-exact-weak-elimination-v1",
            "image_revision": EXPECTED_IMAGE,
            "original_fixed_h0_p0": True,
            "accepted_Y0": False,
            "no_physical_time_steps": True,
        }), encoding="utf-8")
        for h in (.2, .1):
            for case in ("rest", "rotating"):
                for ranks in (1, 2):
                    self.store(h, case, ranks, self.fixture(h, case, ranks))

    @staticmethod
    def fixture(h, case, ranks):
        moment = 0. if case == "rest" else EXPECTED_WEAK_MOMENT_NORMS[h]
        b = [moment, 0.]
        Fy = [-5*moment, 0.]
        return {
            "schema": "scrimp-v130-issue28-exact-weak-elimination-v1",
            "status": "EXACT_BLOCK_ELIMINATION_DIAGNOSTIC_ONLY",
            "nominal_mesh": h, "mpi_ranks": ranks,
            "rotation": 0. if case == "rest" else .2,
            "dofs": 100, "original_geometry_P2": True,
            "original_p_e_p_same_FE_dimension": True,
            "original_U_Y_same_FE_dimension": True,
            "native_Fep_after_exact_copy_l2": 1e-14,
            "native_Fep_after_wall_elimination_l2": 1e-14,
            "native_FU_before_wall_elimination_l2": moment,
            "native_FU_after_wall_elimination_l2": 1e-14,
            "native_FY_after_both_exact_eliminations_l2": 5*moment,
            "native_FY_coefficients": Fy,
            "independent_native_weak_U_moment_coefficients": b,
            "native_FY_equals_signed_h0_wall_moment_sign": -1,
            "native_FY_weak_moment_identity_error_l2": 0.,
            "expected_signed_h0_wall_moment_l2": 5*moment,
            "original_fixed_h0_p0_preserved_bitwise": True,
            "original_traction_Y_preserved_bitwise": True,
            "entire_state_rollback_bitwise": True,
            "no_physical_TS": True,
            "no_pseudo_or_BDF_initializer": True,
            "consistent_rotating_initial_Y0_accepted": False,
            "differential_index_proven": False,
            "generic_initialization_qualified": False,
        }

    def store(self, h, case, ranks, report):
        folder = self.root/(
            f"result_h{str(h).replace('.', 'p')}_{case}_{ranks}"
        )
        folder.mkdir(exist_ok=True)
        (folder/"exact_weak_elimination.json").write_text(
            json.dumps(report), encoding="utf-8"
        )

    def test_native_probe_python_syntax(self):
        probe = ROOT/"sandbox"/"probe_v130_issue28_exact_weak_elimination.py"
        ast.parse(probe.read_text(encoding="utf-8"))

    def test_complete_diagnostic_does_not_qualify_Y0(self):
        summary = compare(self.root)
        self.assertFalse(summary["original_rotating_full_DAE_Y0_initialized"])
        self.assertEqual(
            summary["status"],
            "COMPLETE_BLOCK_ELIMINATION_DIAGNOSTIC_NOT_INITIALIZER",
        )

    def test_rejects_changed_original_data(self):
        r = self.fixture(.2, "rotating", 2)
        r["original_fixed_h0_p0_preserved_bitwise"] = False
        self.store(.2, "rotating", 2, r)
        with self.assertRaisesRegex(ValueError, "Original physical data"):
            compare(self.root)

    def test_rejects_false_Y_identity(self):
        r = self.fixture(.2, "rotating", 1)
        r["native_FY_coefficients"][0] *= .5
        r["native_FY_after_both_exact_eliminations_l2"] *= .5
        self.store(.2, "rotating", 1, r)
        with self.assertRaisesRegex(ValueError, "falsified|does not equal"):
            compare(self.root)

    def test_fine_mesh_native_1e_minus_16_roundoff_is_not_rejected(self):
        # Regression: old native gate failed at h=.1 because 1.0994e-16
        # was just above an arbitrary 1e-16 coefficient threshold.
        r = self.fixture(.1, "rotating", 1)
        r["native_FY_coefficients"][0] += 1.09940161e-16
        r["native_FY_after_both_exact_eliminations_l2"] = abs(
            r["native_FY_coefficients"][0]
        )
        r["native_FY_weak_moment_identity_error_l2"] = abs(
            r["native_FY_coefficients"][0]
            + 5*r["independent_native_weak_U_moment_coefficients"][0]
        )
        self.store(.1, "rotating", 1, r)
        self.assertEqual(
            compare(self.root)["status"],
            "COMPLETE_BLOCK_ELIMINATION_DIAGNOSTIC_NOT_INITIALIZER"
        )

    def test_rejects_actual_wall_identity_corruption(self):
        r = self.fixture(.1, "rotating", 1)
        r["native_FY_coefficients"][0] += 2e-12
        r["native_FY_after_both_exact_eliminations_l2"] = abs(
            r["native_FY_coefficients"][0]
        )
        r["native_FY_weak_moment_identity_error_l2"] = abs(
            r["native_FY_coefficients"][0]
            + 5*r["independent_native_weak_U_moment_coefficients"][0]
        )
        self.store(.1, "rotating", 1, r)
        with self.assertRaisesRegex(ValueError, "does not equal"):
            compare(self.root)

    def test_rejects_missing_MPI2_evidence(self):
        f = self.root/"result_h0p1_rotating_2"/"exact_weak_elimination.json"
        f.unlink()
        with self.assertRaisesRegex(ValueError, "Missing native elimination report"):
            compare(self.root)

    def test_rejects_false_nonzero_rest(self):
        r = self.fixture(.2, "rest", 1)
        r["native_FY_coefficients"] = [1e-7, 0.]
        r["native_FY_after_both_exact_eliminations_l2"] = 1e-7
        self.store(.2, "rest", 1, r)
        with self.assertRaises(ValueError):
            compare(self.root)


if __name__ == "__main__":
    unittest.main()
