# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp


"""Source-only, dependency-light interface brick plan regressions."""

import importlib.util
from pathlib import Path
import sys
import types
import unittest

ROOT = Path(__file__).resolve().parents[1]
pkg = types.ModuleType("scrimp")
pkg.__path__ = [str(ROOT / "scrimp")]
prior = sys.modules.get("scrimp")
sys.modules["scrimp"] = pkg
try:
    from scrimp.interface_ownership import InterfaceOwnershipError
    from scrimp.interface_assembly_plan import plan_brick_sides
finally:
    if prior is None:
        del sys.modules["scrimp"]
    else:
        sys.modules["scrimp"] = prior


PORTS = {
    "Interface Heat": ("Y_T", "U_T"),
    "Interface Wave": ("Y_w", "U_w"),
}
BRICKS = [
    ("M_Y_T", [10], "Y_T*Test_Y_T"),
    ("B_T", [10], "-U_T*Test_e_Q.Normal"),
    ("B_T^T", [10], "e_Q.Normal*Test_Y_T"),
    ("M_Y_w", [10], "Y_w*Test_Y_w"),
    ("B_w", [10], "U_w*Test_p"),
    ("B_w^T", [10], "p*Test_Y_w"),
    ("M_U_T", [10], "-U_T*Test_U_T"),
    ("U_T_feedback", [10], "Y_w*Test_U_T"),
    ("M_U_w", [10], "-U_w*Test_U_w"),
    ("U_w_feedback", [10], "Y_T*Test_U_w"),
    ("bulk", [1], "e_Q.Test_e_Q"),
]


class InterfaceAssemblyPlanTests(unittest.TestCase):
    def plan(self, heat, wave, bricks=BRICKS):
        return plan_brick_sides(
            10, (1, 2), PORTS,
            {"e_Q": heat, "T": heat, "p": wave, "q": wave},
            {"Interface Heat": heat, "Interface Wave": wave},
            bricks,
        )

    def test_original_heat_wave_all_ten_interface_bricks(self):
        result = {item["brick"]: item["owner"]
                  for item in self.plan(1, 2)}
        self.assertEqual(result, {
            "M_Y_T": 1, "B_T": 1, "B_T^T": 1,
            "M_Y_w": 2, "B_w": 2, "B_w^T": 2,
            "M_U_T": 1, "U_T_feedback": 1,
            "M_U_w": 2, "U_w_feedback": 2,
        })

    def test_permuted_subdomains(self):
        result = {item["brick"]: item["owner"]
                  for item in self.plan(2, 1)}
        self.assertEqual(
            {k for k, v in result.items() if v == 2},
            {"M_Y_T", "B_T", "B_T^T", "M_U_T", "U_T_feedback"},
        )
        self.assertEqual(
            {k for k, v in result.items() if v == 1},
            {"M_Y_w", "B_w", "B_w^T", "M_U_w", "U_w_feedback"},
        )

    def test_feedback_uses_target_not_source(self):
        form = [("feedback", [10], "Y_w*Test_U_T")]
        result = self.plan(1, 2, form)
        self.assertEqual(result[0]["owner"], 1)

    def test_no_test_fails_closed(self):
        with self.assertRaisesRegex(InterfaceOwnershipError, "no resolvable"):
            self.plan(1, 2, [("ambiguous", [10], "0.")])

    def test_mixed_test_sides_fails_closed(self):
        with self.assertRaisesRegex(InterfaceOwnershipError, "multiple subdomains"):
            self.plan(1, 2, [("mixed", [10], "Test_U_T+Test_U_w")])

    def test_ignores_non_interface_bricks(self):
        self.assertEqual(
            self.plan(1, 2, [("bulk", [1], "T*Test_T")]), []
        )

    def test_nonadjacent_port_owner_rejected(self):
        with self.assertRaisesRegex(InterfaceOwnershipError, "non-adjacent"):
            plan_brick_sides(
                10, (1, 2), PORTS, {"e_Q": 1, "p": 2},
                {"Interface Heat": 99, "Interface Wave": 2},
                BRICKS,
            )


if __name__ == "__main__":
    unittest.main()
