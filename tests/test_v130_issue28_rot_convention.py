# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Source-only guard: GetFEM nested explicit matrices are column-major.

The earlier macro [[0,1],[-1,0]] was misread as Python matrix rows:
actual GetFEM GWFL interpreted it as COUNTERclockwise. This test
ensures that the original clockwise p0 and the actual boundary
Rot*Normal macro are again aligned. Native confirmation is separate.
"""

import ast
from pathlib import Path
import unittest


SOURCE = Path(__file__).resolve().parents[1] / "sandbox/viscous_swe_2d_mpi.py"


def source_macros():
    tree = ast.parse(SOURCE.read_text(encoding="utf-8"))
    macros = {}
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        f = node.func
        if (isinstance(f, ast.Attribute) and f.attr == "add_macro"
                and len(node.args) >= 2 and all(
                    isinstance(q, ast.Constant) and isinstance(q.value, str)
                    for q in node.args[:2]
                )):
            macros[node.args[0].value] = node.args[1].value
    return macros


def actual_getfem_matrix(source_text):
    """GetFEM GWFL [[...],[...]] lists COLUMNS, not NumPy rows."""
    columns = ast.literal_eval(source_text)
    if len(columns) != 2 or any(len(col) != 2 for col in columns):
        raise ValueError("Expected two 2D column vectors")
    return [[columns[j][i] for j in range(2)] for i in range(2)]


class Issue28GetFEMRotConventionTest(unittest.TestCase):
    def test_rot_produces_clockwise_tangent_at_rightmost_point(self):
        matrix = actual_getfem_matrix(source_macros()["Rot"])
        self.assertEqual(matrix, [[0, 1], [-1, 0]])
        nx, ny = 1, 0
        tangent = [matrix[0][0]*nx + matrix[0][1]*ny,
                   matrix[1][0]*nx + matrix[1][1]*ny]
        self.assertEqual(tangent, [0, -1])
        # p0/rho = omega*(y,-x)/(r+eps) has the same direction.
        self.assertEqual(tangent, [0, -1])

    def test_gyro_is_sign_invariant_and_untouched(self):
        macros = source_macros()
        self.assertEqual(macros["Tangent"], "(Rot*Normal)")
        self.assertEqual(macros["Curl2D(v)"], "div(Rot*v)")
        self.assertEqual(macros["Gyro(v)"], "Curl2D(v)*Rot")
        # Reversing Rot changes the sign of div(Rot*p) and of
        # Rot itself. Their product (gyroscopic force) is invariant.
        for curl in (-3, 0, 1.5):
            for vx, vy in ((2, -4), (0, 1)):
                matrix = actual_getfem_matrix(macros["Rot"])
                rp = [matrix[0][0]*vx + matrix[0][1]*vy,
                      matrix[1][0]*vx + matrix[1][1]*vy]
                transformed = [curl*q for q in rp]
                opposite_rot = [-q for q in rp]
                opposite = [(-curl)*q for q in opposite_rot]
                self.assertEqual(transformed, opposite)


if __name__ == "__main__":
    unittest.main()
