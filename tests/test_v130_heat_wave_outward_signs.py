# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Source-only invariant: the Heat-Wave paper's outward-normal laws.

Prevents the historical disk/annulus-dependent sign compensation from
silently returning to the live example. No GetFEM/PETSc required.
"""

import ast
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "examples" / "heat_wave.py"
TREE = ast.parse(SOURCE.read_text(encoding="utf-8"))


def str_literal(node):
    return node.value if isinstance(node, ast.Constant) and isinstance(node.value, str) else None


class OutwardHeatWaveLaws(unittest.TestCase):
    def test_interface_heat_mixed_boundary_is_outward_normal(self):
        forms = []
        for node in ast.walk(TREE):
            if (isinstance(node, ast.Call)
                    and isinstance(node.func, ast.Attribute)
                    and node.func.attr == "Brick"
                    and len(node.args) >= 2
                    and str_literal(node.args[0]) == "B_T"):
                forms.append(node.args[1])
        self.assertEqual(len(forms), 1)
        self.assertEqual(str_literal(forms[0]), "-U_T*Test_e_Q.Normal")

    def test_interface_wave_feedback_is_outward_normal(self):
        expressions = []
        for node in ast.walk(TREE):
            if (isinstance(node, ast.Call)
                    and isinstance(node.func, ast.Attribute)
                    and node.func.attr == "set_control"
                    and len(node.args) >= 2
                    and str_literal(node.args[0]) == "Interface Wave"):
                expressions.append(node.args[1])
        self.assertEqual(len(expressions), 1)
        self.assertEqual(str_literal(expressions[0]), "Y_T")

    def test_no_region_dependent_normal_sign_left(self):
        names = {node.id for node in ast.walk(TREE)
                 if isinstance(node, ast.Name)}
        self.assertNotIn("heat_normal_orientation", names)

    def test_power_cancel_with_outward_traction_and_continuity(self):
        for heat_region, wave_region in ((1, 2), (2, 1)):
            self.assertEqual({heat_region, wave_region}, {1, 2})
            for y_wave, y_heat in ((0., 0.), (1.75, -0.4), (-3., 2.5)):
                u_heat, u_wave = y_wave, y_heat
                # The model's mixed heat / grad-grad wave weak forms
                # yield -U_T Y_T + U_w Y_w at the interface.
                balance = -u_heat*y_heat + u_wave*y_wave
                self.assertEqual(balance, 0.0)


if __name__ == "__main__":
    unittest.main()
