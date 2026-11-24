# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Source-only negative and positive tests for interface port ownership."""

import importlib.util
from pathlib import Path
import sys
import unittest

FILE = Path(__file__).resolve().parents[1] / "scrimp" / "interface_ownership.py"
SPEC = importlib.util.spec_from_file_location("scrimp_interface_ownership_test", FILE)
MOD = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MOD
SPEC.loader.exec_module(MOD)


class OwnershipTests(unittest.TestCase):
    def test_heat_wave_anchors_independent_of_names_or_region_order(self):
        ports = {
            "port_A": ("Y_T", "U_T"),
            "port_B": ("Y_w", "U_w"),
        }
        bulk = {"e_Q": 7, "p": 3, "q": 3, "T": 7}
        bricks = [
            ([10], "Y_T*Test_Y_T"),
            ([10], "-U_T*Test_e_Q.Normal"),
            ([10], "e_Q.Normal*Test_Y_T"),
            ([10], "Y_w*Test_Y_w"),
            ([10], "U_w*Test_p"),
            ([10], "p*Test_Y_w"),
            # Cross-subdomain feedback: MUST NOT assign sides.
            ([10], "(Y_w)*Test_U_T"),
            ([10], "(Y_T)*Test_U_w"),
            ([7], "T*Test_T"),
            ([3], "p*Test_p"),
        ]
        self.assertEqual(
            MOD.infer_interface_port_owners(10, ports, bulk, bricks),
            {"port_A": 7, "port_B": 3},
        )

    def test_no_guess_from_feedback_only(self):
        with self.assertRaisesRegex(MOD.InterfaceOwnershipError, "no anchored owner"):
            MOD.infer_interface_port_owners(
                10, {"unanchored": ("Y", "U")}, {"p": 2},
                [([10], "Y*Test_Y"), ([10], "(Y)*Test_U")],
            )

    def test_contradictory_anchor_rejected(self):
        with self.assertRaisesRegex(MOD.InterfaceOwnershipError, "conflicting owners"):
            MOD.infer_interface_port_owners(
                10, {"boundary_port": ("Y", "U")}, {"a": 1, "b": 2},
                [([10], "a*Test_Y"), ([10], "b*Test_U")],
            )

    def test_mixed_bulk_owners_in_one_form_rejected(self):
        with self.assertRaisesRegex(MOD.InterfaceOwnershipError, "multiple subdomains"):
            MOD.infer_interface_port_owners(
                10, {"boundary_port": ("Y", "U")}, {"a": 1, "b": 2},
                [([10], "(a+b)*Test_Y")],
            )

    def test_identifiers_are_not_substrings(self):
        used = MOD.variables_in_form(
            "Test_e_Q*e_Q2 + Test_Y_T", {"e_Q", "e_Q2", "Y_T", "Q"}
        )
        self.assertEqual(used, {"e_Q", "e_Q2", "Y_T"})

    def test_other_region_is_not_used_as_anchor(self):
        with self.assertRaisesRegex(MOD.InterfaceOwnershipError, "no anchored owner"):
            MOD.infer_interface_port_owners(
                10, {"port": ("Y", "U")}, {"a": 1},
                [([20], "a*Test_Y")],
            )


if __name__ == "__main__":
    unittest.main()
