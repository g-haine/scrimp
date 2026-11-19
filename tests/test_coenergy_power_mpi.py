# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""No-PETSc tests for exact CN dissipation and MPI physical comparisons."""
import copy
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "_scrimp_compare_coenergy_power_mpi",
    ROOT / "sandbox" / "compare_coenergy_power_mpi.py",
)
CHECK = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(CHECK)


class CoenergyPowerMpiTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.directory = Path(temporary.name)
        self.steps = 5
        self.times = np.arange(self.steps+1)*.002
        self.states = np.ones((self.steps+1, 4), dtype=float)
        self.energy = np.arange(self.steps+1)*.001 + .1
        self.powers = np.full(self.steps, .01)
        self.base = {
            "schema": "scrimp-coenergy-physics-v1",
            "source_commit": "testsha",
            "mesh": .12,
            "dofs": 4, "steps": self.steps,
            "t_final": .01, "dt": .002,
            "full_state_snapshots": self.steps+1,
            "TS_converged_reason": 2,
            "formulation": "substituted_coenergy",
            "effective_mumps_policy": "on_change",
            "power_sign_convention": "synthetic CN balance",
            "maxima": {
                "cn_relative_residual": 1e-14,
                "algebraic_relative_residual": 1e-14,
                "relative_discrete_power_defect": 1e-14,
                "independent_GetFEM_H_relative_error": 1e-14,
            },
        }
        self.reports = {
            1: {**copy.deepcopy(self.base), "mpi_ranks": 1},
            2: {**copy.deepcopy(self.base), "mpi_ranks": 2},
        }
        self.data = {
            r: dict(
                times=self.times.copy(), states=self.states.copy(),
                energy=self.energy.copy(),
                energy_getfem=self.energy.copy(),
                power_damping_mid=self.powers.copy(),
                power_boundary_mid=self.powers.copy(),
                cn_relative_residuals=np.full(self.steps, 1e-14),
                algebraic_relative_residuals=np.full(self.steps, 1e-14),
                relative_power_defects=np.full(self.steps, 1e-14),
            )
            for r in (1, 2)
        }
        self.flush()

    def flush(self):
        for r in (1, 2):
            (self.directory/f"coenergy_{r}r_physics.json").write_text(
                json.dumps(self.reports[r]), encoding="utf-8",
            )
            np.savez_compressed(
                self.directory/f"coenergy_{r}r_physics.npz",
                **self.data[r]
            )

    def test_exact_cn_discrete_power_identity_on_quadratic_model(self):
        # SCRIMP's flow bricks enter the residual with a MINUS:
        # M=-M_physical, K=-R_physical.
        # Hamiltonian H=-0.5 z.T M z; P=-z.T K z is positive.
        mass = -np.diag([2., 3., 1.])
        stiffness = -np.diag([.5, 1., 2.])
        dt = .07
        z0 = np.array([1., 2., -.7])
        z1 = np.linalg.solve(
            mass + .5*dt*stiffness,
            (mass - .5*dt*stiffness)@z0,
        )
        midpoint = .5*(z0+z1)
        h0 = -.5*z0@mass@z0
        h1 = -.5*z1@mass@z1
        power = -midpoint@stiffness@midpoint
        self.assertAlmostEqual(
            h1-h0+dt*power, 0., places=13
        )
        self.assertGreater(power, 0.)
        np.testing.assert_allclose(
            mass@((z1-z0)/dt)+stiffness@midpoint,
            np.zeros(3), atol=1e-13,
        )

    def test_mixed_dae_signed_damping_and_feedback_balance(self):
        # Synthetic pH DAE in order [U,Y,e_r,p], reproducing
        # the source example's signs: M_p<0, K_ee<0,
        # K_UU<0; U=.5Y, e_r=p, Y=p.
        M = np.diag([0., 0., 0., -2.])
        K = np.array([
            [-1., .5, 0., 0.],
            [0., -1., 0., 1.],
            [0., 0., -1., 1.],
            [-1., 0., -1., 0.],
        ])
        dt = .04
        p0 = 2.
        p1 = ((2.-.75*dt)/(2.+.75*dt))*p0
        x0 = np.array([.5*p0, p0, p0, p0])
        x1 = np.array([.5*p1, p1, p1, p1])
        mid = .5*(x0+x1)
        h0 = -.5*x0@M@x0
        h1 = -.5*x1@M@x1
        pr = -mid[2]*K[2, 2]*mid[2]
        pb = -mid[0]*K[0, 0]*mid[1]
        self.assertGreater(pr, 0.)
        self.assertGreater(pb, 0.)
        self.assertAlmostEqual(
            h1-h0+dt*(pr+pb), 0., places=13
        )
        np.testing.assert_allclose(
            M@((x1-x0)/dt)+K@mid,
            np.zeros(4), atol=1e-13,
        )

    def test_complete_state_and_power_parity(self):
        report = CHECK.compare(self.directory)
        self.assertEqual(report["mpi_ranks"], [1, 2])
        self.assertEqual(
            report["metrics"]["states"]["relative_l2_difference"],
            0.,
        )

    def test_different_mumps_policy_is_rejected(self):
        self.reports[2]["effective_mumps_policy"] = "never"
        self.flush()
        with self.assertRaisesRegex(ValueError, "effective_mumps_policy"):
            CHECK.compare(self.directory)

    def test_state_drift_fails(self):
        self.data[2]["states"][2, 2] += .05
        self.flush()
        with self.assertRaisesRegex(ValueError, "physical drift"):
            CHECK.compare(self.directory)

    def test_power_drift_fails(self):
        self.data[2]["power_boundary_mid"][2] += .15
        self.flush()
        with self.assertRaisesRegex(ValueError, "physical drift"):
            CHECK.compare(self.directory)

    def test_nonphysical_power_defect_fails(self):
        self.reports[2]["maxima"][
            "relative_discrete_power_defect"
        ] = 1e-3
        self.flush()
        with self.assertRaisesRegex(ValueError, "power_defect"):
            CHECK.compare(self.directory)

    def test_missing_saved_step_fails(self):
        self.data[2]["times"] = self.data[2]["times"][:-1]
        self.flush()
        with self.assertRaisesRegex(ValueError, "saved accepted CN"):
            CHECK.compare(self.directory)


if __name__ == "__main__":
    unittest.main()
