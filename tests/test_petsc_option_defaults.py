# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""PETSc diagnostics must not be enabled implicitly by importing SCRIMP.

Import-only MPI user projects do not construct a TS, KSP, or solve.
An old import-time PETSc.Options().setValue("info", None) in linalg.py
left an unused -info option in the PETSc database on every rank.
This source contract is runnable without PETSc native dependencies.
"""

import ast
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]
LINALG = ROOT / "scrimp" / "utils" / "linalg.py"


class PETScOptionDefaultsTest(unittest.TestCase):
    def test_linalg_does_not_force_petsc_info(self):
        tree = ast.parse(LINALG.read_text(encoding="utf-8"))
        for call in (node for node in ast.walk(tree)
                     if isinstance(node, ast.Call)
                     and isinstance(node.func, ast.Attribute)
                     and node.func.attr == "setValue"):
            if call.args and isinstance(call.args[0], ast.Constant):
                self.assertNotEqual(
                    call.args[0].value, "info",
                    "Do not activate PETSc -info unconditionally: "
                    "user/CLI configuration should control diagnostics",
                )


if __name__ == "__main__":
    unittest.main()
