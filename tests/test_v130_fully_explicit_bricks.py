# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Isolated tests for the actual accepted-step explicit brick contract.

Pure fake-vector tests exercise the core pre-step freeze and do not
claim native PETSc/GetFEM/TS or power-balance qualification.
"""

import ast
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]
DPHS = ROOT / "scrimp/dphs.py"
SWE = ROOT / "examples/shallow_water.py"
PROBE = ROOT / "sandbox/check_v130_shallow_water_mass_jacobian.py"


def _isolated_methods():
    src = DPHS.read_text(encoding="utf-8")
    tree = ast.parse(src)
    cls = next(
        node for node in tree.body
        if isinstance(node, ast.ClassDef) and node.name == "DPHS"
    )
    wanted = {"_freeze_explicit_rhs", "_prestep_explicit_rhs"}
    functions = [
        node for node in cls.body
        if isinstance(node, ast.FunctionDef) and node.name in wanted
    ]
    assert len(functions) == 2
    exec_namespace = {}
    exec(
        compile(
            ast.Module(body=functions, type_ignores=[]),
            str(DPHS), "exec",
        ),
        exec_namespace,
    )
    return exec_namespace


class _Vector:
    def __init__(self, value=None):
        self.value = value

    def duplicate(self):
        return _Vector()

    def copy(self, other):
        other.value = self.value


class _Brick:
    def __init__(self, explicit):
        self.explicit = explicit

    def get_explicit(self):
        return self.explicit


class _FakeDPHS:
    def __init__(self):
        self.bricks = {"G": _Brick(True)}
        self.rhs = _Vector()
        self._explicit_rhs = None
        self._explicit_rhs_ready = False
        self._explicit_rhs_time = None
        self.time = None
        self.state = None
        self.assemblies = []
        methods = _isolated_methods()
        self._freeze_explicit_rhs = methods[
            "_freeze_explicit_rhs"
        ].__get__(self)
        self._prestep_explicit_rhs = methods[
            "_prestep_explicit_rhs"
        ].__get__(self)

    def _sync_lumped_trial_state(self, t, z):
        self.time = float(t)
        self.state = z

    def assemble_rhs(self, *, explicit_only=False):
        if not explicit_only:
            raise RuntimeError("Explicit cache assembled wrong brick set")
        self.assemblies.append((self.time, self.state))
        self.rhs.value = (self.time, self.state)


class _FakeTS:
    def __init__(self, t, z):
        self.t, self.z = t, z

    def getTime(self):
        return self.t

    def getSolution(self):
        return self.z


class TestFullyExplicitBricks(unittest.TestCase):
    def test_frozen_previous_accepted_step_is_stable_during_newton(self):
        d = _FakeDPHS()
        ts = _FakeTS(0.0, "accepted_z0")
        d._prestep_explicit_rhs(ts)
        self.assertTrue(d._explicit_rhs_ready)
        self.assertEqual(d._explicit_rhs.value, (0.0, "accepted_z0"))
        self.assertEqual(d.assemblies, [(0.0, "accepted_z0")])
        self.assertEqual(d._explicit_rhs_time, 0.0)

        # Trial states must NOT trigger any refresh: their derivative
        # is zero, because the *entire* G forcing was already frozen.
        d._sync_lumped_trial_state(0.005, "newton_trial_z1")
        d._sync_lumped_trial_state(0.005, "line_search_trial_z2")
        self.assertEqual(d._explicit_rhs.value, (0.0, "accepted_z0"))
        self.assertEqual(len(d.assemblies), 1)

        # A genuinely accepted new step refreshes forcing ONCE.
        ts = _FakeTS(0.005, "accepted_z1")
        d._prestep_explicit_rhs(ts)
        self.assertEqual(d._explicit_rhs.value, (0.005, "accepted_z1"))
        self.assertEqual(d._explicit_rhs_time, 0.005)
        self.assertEqual(len(d.assemblies), 2)

    def test_no_explicit_bricks_do_not_allocate_or_freeze(self):
        d = _FakeDPHS()
        d.bricks = {"G": _Brick(False)}
        d._prestep_explicit_rhs(_FakeTS(0.0, "accepted_z0"))
        self.assertIsNone(d._explicit_rhs)
        self.assertFalse(d._explicit_rhs_ready)
        self.assertEqual(d.assemblies, [])

    def test_explicit_gyro_sign_and_weak_form_unchanged(self):
        src = SWE.read_text(encoding="utf-8")
        self.assertIn(
            'S.Brick("G", "(Gyro(p) * e_p) . Test_p", [1]',
            src,
        )
        self.assertIn(
            'linear=False, explicit=True, position="effort"',
            src,
        )
        core = DPHS.read_text(encoding="utf-8")
        self.assertIn(
            'if position == "flow" or brick.get_explicit():',
            core,
        )

    def test_actual_implicit_rhs_and_J_exclude_explicit_brick(self):
        core = DPHS.read_text(encoding="utf-8")
        self.assertIn(
            'elif explicit_only:', core,
        )
        self.assertIn(
            'enabled = brick.get_explicit()', core,
        )
        self.assertIn(
            'and not brick.get_explicit()', core,
        )
        self.assertIn(
            'TS.setPreStep(self._prestep_explicit_rhs)', core,
        )
        self.assertEqual(
            core.count('TS.setPreStep(self._prestep_explicit_rhs)'), 2,
        )
        self.assertIn(
            'F.axpy(-1, self._explicit_rhs)', core,
        )
        self.assertIn(
            '"explicit=True nonlinear bricks require a frozen RHS "', core,
        )

    def test_directional_fd_freezes_nonzero_lagged_G(self):
        src=PROBE.read_text(encoding="utf-8")
        self.assertIn("system._freeze_explicit_rhs(0., accepted)",src)
        self.assertIn("frozen_norm <= 1e-13",src)
        self.assertIn("explicit_drift != 0.",src)
        self.assertNotIn("diagnostic_full_tangent=",src)

    def test_strict_runner_checks_production_only(self):
        runner=(ROOT/"scripts/run-v130-shallow-mass-jacobian.sh"
               ).read_text(encoding="utf-8")
        self.assertIn("for mode in default legacy;",runner)
        self.assertIn("STRICT PRODUCTION GATE",runner)
        self.assertNotIn("--tangent-policy", runner)
        self.assertNotIn(
            "compare_v130_shallow_explicit_tangent", runner
        )


if __name__=="__main__":
    unittest.main()
