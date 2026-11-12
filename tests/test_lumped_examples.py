# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Smoke tests for 0D port-Hamiltonian release examples.

These tests need the full numerical SCRIMP environment and intentionally
run very short simulations to keep CI cost low. They do not replace
time-refinement or physical power-balance validation.
"""

import unittest
import numpy as np

try:
    from examples.mass_spring_damper import mass_spring_damper
    from examples.nonlinear_pendulum import nonlinear_pendulum
except ModuleNotFoundError as exc:
    if exc.name in {"getfem", "gmsh", "petsc4py", "slepc4py"}:
        raise unittest.SkipTest("SCRIMP numerical dependencies are unavailable") from exc
    raise


class LumpedExamplesSmokeTest(unittest.TestCase):
    def test_mass_spring_damper(self):
        system = mass_spring_damper(final_time=0.04, time_step=0.01)
        self.assertEqual(system.gf_model.interval_of_variable("q")[1], 1)
        self.assertEqual(system.gf_model.interval_of_variable("p")[1], 1)
        t = np.asarray(system.solution["t"], dtype=float)
        self.assertGreater(len(t), 1)
        self.assertAlmostEqual(float(t[-1]), 0.04, places=8)
        q = np.asarray(system.get_solution("q"), dtype=float).ravel()
        self.assertTrue(np.isfinite(q).all())
        self.assertLess(np.max(np.abs(q)), 2.0)

    def test_nonlinear_pendulum(self):
        system = nonlinear_pendulum(
            final_time=0.04, time_step=0.01,
            initial_angle=0.8, torque_amplitude=0.0
        )
        self.assertEqual(system.gf_model.interval_of_variable("theta")[1], 1)
        self.assertEqual(system.gf_model.interval_of_variable("p")[1], 1)
        t = np.asarray(system.solution["t"], dtype=float)
        self.assertGreater(len(t), 1)
        self.assertAlmostEqual(float(t[-1]), 0.04, places=8)
        theta = np.asarray(system.get_solution("theta"), dtype=float).ravel()
        self.assertTrue(np.isfinite(theta).all())
        self.assertLess(np.max(np.abs(theta)), 2.0)
        H = np.asarray(system.get_Hamiltonian(), dtype=float).ravel()
        self.assertTrue(np.isfinite(H).all())
        self.assertGreaterEqual(float(H[0]), 0.0)
        # Damped, unforced mechanics should not generate appreciable energy.
        self.assertLessEqual(float(H[-1]-H[0]), 1e-5)


if __name__ == "__main__":
    unittest.main()
