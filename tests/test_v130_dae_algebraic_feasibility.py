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

    def test_nonfinite_inputs_refused(self):
        with self.assertRaises(ValueError):
            solve_algebraic_block(
                csr_matrix(np.eye(2)),
                np.array([np.nan, 0.0]),
            )


if __name__ == "__main__":
    unittest.main()
