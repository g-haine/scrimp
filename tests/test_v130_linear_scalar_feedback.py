# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Dependency-free regression for exact constant-gain GetFEM feedback.

This loads only the real helper AST from dphs.py, not PETSc/GetFEM;
the physical CN gate still requires the native 2C MPI pilot. An
unknown or nonlinear law MUST retain the nonlinear GetFEM path.
"""
import ast
from pathlib import Path
import re
import unittest

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
PATH = ROOT/"scrimp"/"dphs.py"
TEXT = PATH.read_text(encoding="utf-8")
TREE = ast.parse(TEXT, filename=str(PATH))
FUNCTION = "_control_expression_is_direct_linear_feedback"
NODES = [node for node in TREE.body
         if isinstance(node, ast.FunctionDef) and node.name == FUNCTION]
if len(NODES) != 1:
    raise AssertionError("Actual DPHS feedback classifier is missing or duplicated")
source_module = ast.Module(body=[NODES[0]], type_ignores=[])
scope = {"np": np, "re": re}
exec(compile(source_module, str(PATH), "exec"), scope)
is_linear = scope[FUNCTION]


class ConstantFeedbackClassificationTests(unittest.TestCase):
    def setUp(self):
        self.variables = {"Y", "Y_T", "Y_w", "U", "p"}

    def test_exact_finite_literal_gain_is_linear(self):
        for expression in (
            "Y", "-Y", "+Y", " Y_w ", "- Y_T",
            "-0.5*Y", "+0.5*Y", ".5*Y", "-.5*Y",
            "2*Y_w", "-2.*Y_w", "2.5E-3*Y_T",
            "-1e-10 * p", "0.0*Y", "1e-2*Y",
        ):
            with self.subTest(expression=expression):
                self.assertTrue(is_linear(expression, self.variables))

    def test_reject_nonlinear_state_laws_unknown_symbols_and_nonliterals(self):
        for expression in (
            "Y*Y", "Y + Y_w", "Y*U", "Y/2", "Y**2",
            "sin(Y)", "sqrt(Y)", "exp(Y)", "abs(Y)",
            "(Y)", "(-0.5)*Y", "-k*Y", "gain*Y", "k*Y",
            "1/(1+Y)", "-0.5*Unknown", "Unknown",
            "0.5*Y*Y", "1e309*Y", "1e10000*Y",
            "nan*Y", "inf*Y", "", "0.5",
            "Y0", "--Y", "1e-3 * (Y)",
        ):
            with self.subTest(expression=expression):
                self.assertFalse(is_linear(expression, self.variables))

    def test_coenergy_feedback_is_exactly_classified_linear(self):
        example = (ROOT/"examples"/"wave_coenergy.py").read_text(encoding="utf-8")
        self.assertIn('FEEDBACK_GAIN = 0.5', example)
        self.assertIn('expressions = [f"{FEEDBACK_GAIN}*Y"]', example)
        coefficient = .5
        expression = f"{coefficient}*Y"
        self.assertTrue(is_linear(expression, self.variables))
        self.assertIn('direct_linear_feedback = (', TEXT)
        self.assertIn('linear=direct_linear_feedback,', TEXT)

    def test_no_hidden_nonlinear_initializer_or_explicit_zeroing(self):
        example = (ROOT/"examples"/"wave_coenergy.py").read_text(encoding="utf-8")
        self.assertIn('init_step_ts_type="algebraic", init_step=True', example)
        self.assertNotIn('init_step_ts_type="pseudo"', example)
        self.assertNotIn('init_step_ts_type="nonlinear_affine"', example)
        self.assertIn('f"{FEEDBACK_GAIN}*Y*Y"', example)
        self.assertIn('step_boundary = dt*P[1:]', example)
        self.assertIn('"cumulative_CN_dissipative_energy_balance"', example)


if __name__ == "__main__":
    unittest.main()
