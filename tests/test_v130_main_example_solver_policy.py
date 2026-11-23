# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Source-only guarantees for the four non-published example solver APIs.

No native GetFEM/PETSc simulation is claimed by these tests.
"""

import ast
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]
CASES = {
    "mass_spring_damper": {
        "name": "mass_spring_damper",
        "default_time": {"final_time":8.0, "time_step":.01},
    },
    "nonlinear_pendulum": {
        "name": "nonlinear_pendulum",
        "default_time": {"final_time":12.0, "time_step":.01},
    },
    "heat_wave": {
        "name": "heat_wave_eq",
        "default_time": {"final_time":2.0,"time_step":.001,"save_step":.05},
    },
    "shallow_water": {
        "name": "shallow_water_eq",
        "default_time": {"final_time":.5,"time_step":.0001,"save_step":.01},
    },
}


class MainExampleSolverPolicyTests(unittest.TestCase):
    def test_approved_four_example_signatures_and_defaults(self):
        for filename, item in CASES.items():
            with self.subTest(example=filename):
                source=(ROOT/"examples"/f"{filename}.py").read_text()
                tree=ast.parse(source)
                func=next(node for node in tree.body
                          if isinstance(node,ast.FunctionDef)
                          and node.name==item["name"])
                kw_args=[arg.arg for arg in func.args.kwonlyargs]
                self.assertIn("factor_solver",kw_args)
                index=kw_args.index("factor_solver")
                self.assertIsNone(
                    ast.literal_eval(func.args.kw_defaults[index])
                )
                defaults={
                    node.arg: ast.literal_eval(value)
                    for node,value in zip(func.args.args[-len(func.args.defaults):],
                                          func.args.defaults)
                }
                defaults.update({
                    arg.arg: ast.literal_eval(value)
                    for arg,value in zip(func.args.kwonlyargs,func.args.kw_defaults)
                })
                for key,value in item["default_time"].items():
                    self.assertEqual(defaults[key],value)
                # The examples must NEVER impose serial SuperLU.
                self.assertNotIn('pc_factor_mat_solver_type="superlu"',source)
                self.assertNotIn('pc_factor_mat_solver_type="superlu_dist"',source)
                self.assertIn('if factor_solver is not None:',source)
                self.assertIn('scheme["pc_factor_mat_solver_type"] = factor_solver',
                              source)
                self.assertIn("**scheme",source)

    def test_production_default_is_mumps_and_user_selection_is_not_masked(self):
        source=(ROOT/"scrimp"/"dphs.py").read_text()
        self.assertIn(
            'if not self.time_scheme.hasName("pc_factor_mat_solver_type"):',
            source
        )
        self.assertIn(
            'self.time_scheme["pc_factor_mat_solver_type"] = "mumps"',
            source
        )

    def test_published_snapshot_guard_not_bypassed(self):
        source=(ROOT/"scripts"/"check-v130-published-immutable.py").read_text()
        self.assertIn("EXPECTED_GIT_BLOBS",source)
        self.assertIn('PREFIX = "examples/published/"',source)
        self.assertIn("ARCHIVED WORKTREE BYTES changed",source)
        self.assertNotIn("UNTRACKED new source",source)


if __name__=="__main__":
    unittest.main()
