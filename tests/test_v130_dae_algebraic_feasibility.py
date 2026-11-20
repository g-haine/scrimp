# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Pure SciPy numerical fixtures for the fixed-dynamic DAE feasibility probe.

These tiny synthetic matrices are tests of the decision/reporting
semantics, not evidence of the actual finite-element wave physics.
"""

import ast
import math
import sys
from pathlib import Path
import unittest

import numpy as np
from scipy.sparse import csr_matrix


SANDBOX = Path(__file__).resolve().parents[1] / "sandbox"
sys.path.insert(0, str(SANDBOX))
from solve_v130_dae_algebraic_block import solve_algebraic_block


class DAEAlgebraicFeasibilityTests(unittest.TestCase):
    def test_compatible_constraints_solve_without_dynamic_correction(self):
        Kaa = csr_matrix(np.diag([2.0, 3.0]))
        r = np.array([2.0, -3.0])
        correction, status = solve_algebraic_block(Kaa, r)
        self.assertTrue(status["linearized_compatible_within_tolerance"])
        self.assertEqual(
            status["linear_feasibility_decision"], "LINEAR_COMPATIBLE"
        )
        np.testing.assert_allclose(
            correction, [-1.0, 1.0], atol=1e-11, rtol=0
        )
        self.assertLess(
            np.linalg.norm(Kaa @ correction + r), 1e-11
        )

    def test_inconsistent_constraints_not_marked_compatible(self):
        # The second constraint is 0*x=1, impossible with fixed
        # differential data. Do not report a false projection PASS.
        Kaa = csr_matrix(np.diag([1.0, 0.0]))
        r = np.array([1.0, 1.0])
        correction, status = solve_algebraic_block(Kaa, r)
        self.assertFalse(status["linearized_compatible_within_tolerance"])
        self.assertGreater(
            np.linalg.norm(Kaa @ correction + r), 0.5
        )
        self.assertNotEqual(
            status["linear_feasibility_decision"], "LINEAR_COMPATIBLE"
        )

    def test_singular_but_consistent_is_allowed(self):
        Kaa = csr_matrix(np.diag([1.0, 0.0]))
        r = np.array([2.0, 0.0])
        correction, status = solve_algebraic_block(Kaa, r)
        self.assertTrue(status["linearized_compatible_within_tolerance"])
        np.testing.assert_allclose(Kaa @ correction, -r, atol=1e-11)

    def test_left_nullspace_witness_positive_for_inconsistent_block(self):
        Kaa = csr_matrix(np.diag([1.0, 0.0]))
        residual = np.array([1.0, 1.0])
        correction, status = solve_algebraic_block(Kaa, residual)
        self.assertFalse(status["linearized_compatible_within_tolerance"])
        self.assertAlmostEqual(
            status["residual_left_witness_adjoint_l2"], 0.0,
            delta=1e-10,
        )
        self.assertGreater(
            status["dual_lower_bound_for_corrections_within_radius_l2"],
            0.9,
        )
        self.assertGreater(
            status["residual_left_witness_pairing_with_initial"],
            0.9,
        )
        self.assertAlmostEqual(
            status["normal_equation_residual_l2_recomputed"], 0.0,
            delta=1e-10,
        )

    def test_dual_bound_zero_when_exact_linear_solution_exists(self):
        Kaa = csr_matrix(np.eye(2))
        residual = np.array([3.0, -2.0])
        _correction, status = solve_algebraic_block(
            Kaa, residual
        )
        self.assertTrue(status["linearized_compatible_within_tolerance"])
        self.assertLess(
            status["dual_lower_bound_for_corrections_within_radius_l2"],
            1e-9,
        )

    def test_matched_left_source_uses_getfem_not_python_power_syntax(self):
        # Do not import the probe: that requires native PETSc/GetFEM.
        # Check its literal diagnostic control independently of the
        # production wave example and verify its mathematical trace.
        source = (
            SANDBOX / "probe_v130_wave_fixed_dynamics.py"
        ).read_text(encoding="utf-8")
        tree = ast.parse(source)
        assigned = [
            node.value for node in tree.body
            if isinstance(node, ast.Assign)
            and any(
                isinstance(target, ast.Name)
                and target.id == "MATCHED_LEFT_VELOCITY_EXPRESSION"
                for target in node.targets
            )
        ]
        self.assertEqual(len(assigned), 1)
        expression = ast.literal_eval(assigned[0])
        self.assertNotIn("**", expression)
        self.assertIn("exp(", expression)
        for y in (0.0, 0.125, 0.25, 0.5, 0.75, 0.875, 1.0):
            value = eval(
                expression,
                {"__builtins__": {}, "exp": math.exp},
                {"y": y},
            )
            expected = (
                3.0 ** (-20.0 * (0.25 + (y - 0.5) ** 2))
                / 3.0
            )
            self.assertAlmostEqual(
                value, expected,
                delta=2e-14 * max(1.0, abs(expected)),
            )

    def test_exact_fe_momentum_positive_control_is_polynomial(self):
        # Pure source-level regression: prove the diagnostic's p0/rho
        # equals the polynomial effort x with exactly zero left trace.
        # Real GetFEM FE assembly is qualified separately in native MPI.
        source = (
            SANDBOX / "probe_v130_wave_fixed_dynamics.py"
        ).read_text(encoding="utf-8")
        tree = ast.parse(source)
        assigned = [
            node.value for node in tree.body
            if isinstance(node, ast.Assign)
            and any(
                isinstance(target, ast.Name)
                and target.id == "EXACT_FE_COMPATIBLE_P0"
                for target in node.targets
            )
        ]
        self.assertEqual(len(assigned), 1)
        p0_expression = ast.literal_eval(assigned[0])
        self.assertEqual(p0_expression, "x*(3-x)")
        for x in (0., .125, .5, 1., 1.5, 2.):
            p0 = x*(3.-x)
            rho = 3.-x
            self.assertEqual(p0/rho, x)
        self.assertEqual(0.*(3.-0.), 0.)

    def test_nonfinite_inputs_refused(self):
        with self.assertRaises(ValueError):
            solve_algebraic_block(
                csr_matrix(np.eye(2)),
                np.array([np.nan, 0.0]),
            )


if __name__ == "__main__":
    unittest.main()
