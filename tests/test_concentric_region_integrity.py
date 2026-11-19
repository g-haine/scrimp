# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Source-level guard: never merge Concentric interfaces into 2D volumes.

This test deliberately imports neither GetFEM nor PETSc nor Gmsh;
native matrix/physics parity is validated by MPI integration audits.

The 2026-11-18 regression had 5,790 rather than 5,724 DoFs,
with 66 artificial dissipative DoFs and a 7.93e15 norm for
the e_r diagonal block from e_r/nu integrated on nu=0 at
the 1D interface r=0.6. Assembled M,K disagreed between
1 and 2 ranks; the PETSc COO bridge was coefficientwise exact.
"""

import ast
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "scrimp" / "utils" / "mesh.py"


def _region_calls(function):
    """Return (method, integer region args) for region manipulation."""
    for node in ast.walk(function):
        if not isinstance(node, ast.Call):
            continue
        if not isinstance(node.func, ast.Attribute):
            continue
        if node.func.attr not in ("region_merge", "region_subtract"):
            continue
        try:
            values = tuple(ast.literal_eval(a) for a in node.args)
        except (ValueError, TypeError):
            raise AssertionError("Concentric region operations must be static")
        yield (node.func.attr, values)


class ConcentricRegionIntegrityTests(unittest.TestCase):
    def test_interface_faces_are_not_merged_into_volume_regions(self):
        tree = ast.parse(SOURCE.read_text(encoding="utf-8"))
        functions = [
            node for node in tree.body
            if isinstance(node, ast.FunctionDef)
            and node.name == "Concentric"
        ]
        self.assertEqual(len(functions), 1)
        calls = list(_region_calls(functions[0]))
        self.assertIn(("region_merge", (10, 11)), calls)
        self.assertIn(("region_merge", (10, 12)), calls)
        self.assertIn(("region_merge", (20, 21)), calls)
        self.assertIn(("region_merge", (20, 22)), calls)
        self.assertFalse(
            any(
                method == "region_merge" and values[0] in (1, 2)
                for method, values in calls
            ),
            "1D interface elements must stay outside 2D bulk regions",
        )

    def test_concentric_keeps_distinct_interface_region(self):
        tree = ast.parse(SOURCE.read_text(encoding="utf-8"))
        concentric = next(
            n for n in tree.body
            if isinstance(n, ast.FunctionDef) and n.name == "Concentric"
        )
        # Concentric must still advertise a separate interface region
        # after the internal arc labels are merged into physical region 10.
        text = ast.get_source_segment(
            SOURCE.read_text(encoding="utf-8"), concentric
        )
        self.assertIn('"Interface": 10', text)
        self.assertIn('"Omega_Disk": 1', text)
        self.assertIn('"Omega_Annulus": 2', text)


if __name__ == "__main__":
    unittest.main()
