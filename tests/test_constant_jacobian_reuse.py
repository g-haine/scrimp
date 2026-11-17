# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Pure-Python admission and ownership tests for constant-CN LU reuse."""
import importlib.util
from pathlib import Path
import unittest

SOURCE = (
    Path(__file__).resolve().parents[1]
    / "scrimp/utils/constant_jacobian_reuse.py"
)
SPEC = importlib.util.spec_from_file_location("_constant_lu", SOURCE)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class FakePC:
    def getType(self):
        return self.type

    def getFactorSolverType(self):
        return self.factor_solver

    def __init__(self, kind="lu", factor="mumps"):
        self.type = kind
        self.factor_solver = factor


class FakeKSP:
    def __init__(self, kind="preonly", pc=None):
        self.type = kind
        self.pc = pc or FakePC()
        self.calls = []

    def getType(self):
        return self.type

    def getPC(self):
        return self.pc

    def setReusePreconditioner(self, enabled):
        self.calls.append(bool(enabled))


class FakeTS:
    def __init__(self):
        self.ts_type = "cn"
        self.adapt = "none"
        self.dt = 0.002
        self.ksp = FakeKSP()

    def getType(self):
        return self.ts_type

    def getAdapt(self):
        return self

    def getTimeStep(self):
        return self.dt

    def getSNES(self):
        return self

    def getKSP(self):
        return self.ksp


class ConstantJacobianLuReuseTests(unittest.TestCase):
    def setUp(self):
        self.ts = FakeTS()
        self.kwargs = {
            "all_bricks_linear": True,
            "linear_mass": True,
            "linear_stiffness": True,
            "lumped_variables": {},
            "mumps_policy": "never",
            "initial_time": 0.,
            "final_time": .01,
            "configured_dt": .002,
        }

    def configure(self, **overrides):
        return MODULE.enable_constant_jacobian_lu_reuse(
            self.ts, **{**self.kwargs, **overrides}
        )

    def test_fixed_cn_linear_enables_ksp_reuse_exactly_once(self):
        details = self.configure()
        self.assertEqual(self.ts.ksp.calls, [True])
        self.assertEqual(details["method"], "KSPSetReusePreconditioner")
        self.assertEqual(details["steps"], 5)
        self.assertEqual(details["factor_solver"], "mumps")

    def test_nonlinear_or_lumped_cannot_reuse_lu(self):
        for change in (
            {"all_bricks_linear": False},
            {"linear_mass": False},
            {"linear_stiffness": False},
            {"lumped_variables": {"x": object()}},
        ):
            with self.subTest(change=change):
                with self.assertRaisesRegex(ValueError, "fully linear"):
                    self.configure(**change)
                self.assertEqual(self.ts.ksp.calls, [])

    def test_adaptive_or_wrong_integrator_cannot_reuse(self):
        self.ts.adapt = "basic"
        with self.assertRaisesRegex(ValueError, "adapt_type"):
            self.configure()
        self.ts.adapt = "none"
        self.ts.ts_type = "bdf"
        with self.assertRaisesRegex(ValueError, "Crank-Nicolson"):
            self.configure()
        self.assertEqual(self.ts.ksp.calls, [])

    def test_last_partial_step_rejected_to_preserve_shift(self):
        with self.assertRaisesRegex(ValueError, "partial final CN step"):
            self.configure(final_time=.009)
        self.assertEqual(self.ts.ksp.calls, [])

    def test_changed_ts_step_size_rejected(self):
        self.ts.dt = .004
        with self.assertRaisesRegex(ValueError, "Effective TS"):
            self.configure()
        self.assertEqual(self.ts.ksp.calls, [])

    def test_on_change_policy_rejected(self):
        with self.assertRaisesRegex(ValueError, "never"):
            self.configure(mumps_policy="on_change")
        self.assertEqual(self.ts.ksp.calls, [])

    def test_solver_must_be_preonly_lu_mumps(self):
        self.ts.ksp.type = "gmres"
        with self.assertRaisesRegex(ValueError, "KSPPREONLY"):
            self.configure()
        self.ts.ksp.type = "preonly"
        self.ts.ksp.pc.factor_solver = "superlu"
        with self.assertRaisesRegex(ValueError, "MUMPS"):
            self.configure()
        self.assertEqual(self.ts.ksp.calls, [])


if __name__ == "__main__":
    unittest.main()
