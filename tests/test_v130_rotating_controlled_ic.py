# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Host-only controls for the NEW controlled rotating SWE preparation.

These tests verify the symbolic control schedule, the exact t0 normal
feedback cancellation, safe time grids and pilot source contracts.
They are NOT native GetFEM/PETSc or stationarity evidence.
"""

import ast
import math
from pathlib import Path
from types import SimpleNamespace
import unittest


ROOT = Path(__file__).resolve().parents[1]
MODEL = ROOT / "sandbox" / "viscous_swe_2d_mpi.py"
PILOT = ROOT / "sandbox" / "probe_v130_rotating_controlled_ic.py"


def function_source(path, names):
    """Load ONLY pure source functions; never import PETSc/GetFEM."""
    tree = ast.parse(path.read_text(encoding="utf-8"))
    nodes = [node for node in tree.body
             if isinstance(node, ast.FunctionDef) and node.name in names]
    if len(nodes) != len(names):
        raise AssertionError("Requested pure source function absent")
    namespace = {
        "np": SimpleNamespace(isfinite=math.isfinite), "math": math,
        "MAX_SEGMENT_STEPS": 100, "RAMP_DECAY_START": 75.0,
    }
    exec(compile(ast.Module(body=nodes, type_ignores=[]),
                 str(path), "exec"), namespace)
    return namespace


MODEL_FUNCS = function_source(
    MODEL, {"ramp_expression", "spinup_control_envelopes"}
)
PILOT_FUNCS = function_source(PILOT, {"checked_grid"})


def numerical(expression, t):
    return float(eval(
        expression, {"__builtins__": {}},
        {"t": float(t), "min": min, "max": max}
    ))


class RotatingControlledPreparationHostTest(unittest.TestCase):
    def test_wall_starts_smoothly_and_reaches_target(self):
        wall, feedback = MODEL_FUNCS["spinup_control_envelopes"](75., 200., .02)
        for t, target in [(0., 0.), (.01, .5), (.02, 1.), (.024, 1.),
                          (.04, 1.), (74.99, 1.)]:
            self.assertAlmostEqual(numerical(wall, t), target, places=12)
            self.assertAlmostEqual(numerical(feedback, t), 1., places=12)
        self.assertAlmostEqual(numerical(wall, -.01), 0., places=12)

    def test_feedback_law_is_original_decay_not_spinup_law(self):
        wall, feedback = MODEL_FUNCS["spinup_control_envelopes"](75., 200., .02)
        self.assertNotEqual(wall, feedback)
        for t in (0., .01, .02, .024, 75., 90., 150., 200., 250.):
            decay = numerical(MODEL_FUNCS["ramp_expression"](75., 200.), t)
            self.assertAlmostEqual(numerical(feedback, t), decay, places=12)
        self.assertAlmostEqual(numerical(feedback, 200.), 0., places=12)

    def test_at_zero_gamma_quarter_normal_feedback_is_off(self):
        wall, feedback = MODEL_FUNCS["spinup_control_envelopes"](75., 200., .02)
        gamma, rho, g = .25, 1., .1
        h0, h_target = 5., 4.
        yn0 = -rho*g*h0
        b = -0. - gamma*(1.-numerical(feedback, 0.)) * (
            yn0 + 0. + rho*g*h_target
        )
        self.assertEqual(b, 0.)
        # Regression: wrongly using the tangential rising f at t0
        # activates the normal feedback and destroys the exact rest.
        wrong = -gamma*(1.-numerical(wall, 0.)) * (
            yn0 + rho*g*h_target
        )
        self.assertGreater(abs(wrong), .01)

    def test_after_spinup_both_controls_recover_original_law(self):
        wall, feedback = MODEL_FUNCS["spinup_control_envelopes"](75., 200., .02)
        for t in (.02, .024, .1, 40., 75., 100., 200., 250.):
            self.assertAlmostEqual(
                numerical(wall, t), numerical(feedback, t), places=12
            )

    def test_without_spinup_both_controls_are_identical(self):
        wall, feedback = MODEL_FUNCS["spinup_control_envelopes"](75., 200., 0.)
        self.assertEqual(wall, feedback)
        self.assertEqual(wall, MODEL_FUNCS["ramp_expression"](75., 200.))

    def test_invalid_ramp_rejected(self):
        f = MODEL_FUNCS["spinup_control_envelopes"]
        for duration in (-.1, float("nan"), float("inf")):
            with self.subTest(duration=duration):
                with self.assertRaises(ValueError):
                    f(75., 200., duration)
        with self.assertRaises(ValueError):
            f(200., 75., .02)

    def test_boundary_bricks_use_the_separate_correct_envelopes(self):
        tree = ast.parse(MODEL.read_text(encoding="utf-8"))
        builder = next(node for node in tree.body
                       if isinstance(node, ast.FunctionDef)
                       and node.name == "build_system")
        assignments = {
            t.id: ast.unparse(node.value)
            for node in ast.walk(builder)
            if isinstance(node, ast.Assign)
            for t in node.targets if isinstance(t, ast.Name)
        }
        self.assertIn("f_feedback", assignments["normal_law"])
        self.assertNotIn("f_feedback", assignments["vector_law"])
        self.assertIn("{f}", assignments["vector_law"])

    def test_timestep_validation_and_bounded_archiving(self):
        f = PILOT_FUNCS["checked_grid"]
        self.assertEqual(f(0., .02, .002), 10)
        self.assertEqual(f(.02, .024, .002), 2)
        self.assertEqual(f(.1, .12, .002), 10)
        for a,b,dt in [(0., .023, .002),(0., .1, 0.),
                       (0., .202, .002),(.02, 0., .002),
                       (0., 75.01, .002)]:
            with self.subTest(grid=(a,b,dt)):
                with self.assertRaises(ValueError):
                    f(a,b,dt)

    def test_native_probe_keeps_accepted_state_and_not_steady_label(self):
        source = PILOT.read_text(encoding="utf-8")
        for required in (
            "publish_sandbox_rest(", "system.load_checkpoint(",
            "system.save_checkpoint(", "native_residual(",
            "diagnostics(system, gamma=args.gamma)",
            "stationary_rotating_state_certified", "False",
            "original_power_defect_over_exchange",
            "original_Fa_max_l2", "original_TSCN_max_l2",
        ):
            self.assertIn(required, source)

    def test_energy_failure_preserves_native_traces_before_rejection(self):
        source = PILOT.read_text(encoding="utf-8")
        self.assertIn('"physical_precheck_report.json"', source)
        self.assertIn('"physical_precheck_trace.npz"', source)
        self.assertIn('"step_balance_defect"', source)
        self.assertIn('"checkpoint_created": False', source)
        self.assertLess(
            source.index('"physical_precheck_trace.npz"'),
            source.index('if power_defect > 1e-3 or mass_defect > 1e-8:')
        )
        # The original physical acceptance gate is deliberately unchanged.
        self.assertIn("if power_defect > 1e-3 or mass_defect > 1e-8:", source)
        self.assertLess(
            source.index('if power_defect > 1e-3 or mass_defect > 1e-8:'),
            source.index("    system.save_checkpoint(")
        )

    def test_continuous_split_comparator_never_accepts_power_defect(self):
        source = (ROOT / "sandbox" /
                  "compare_v130_rotating_ic_continuous_split.py").read_text(
                      encoding="utf-8"
                  )
        for required in (
            "np.array_equal(s[\"full_state_final\"], r[\"full_state_start\"])",
            "continuous_vs_restarted_final_full_FE_l2",
            "continuous_vs_split_energy_defect_difference",
            "NATIVE_RESTART_PARITY_DIAGNOSTIC_NOT_ENERGY_ACCEPTANCE",
            "physical_segment_relative_energy_gate_not_overridden",
        ):
            self.assertIn(required, source)

    def test_no_pseudo_bdf_on_resume(self):
        source = PILOT.read_text(encoding="utf-8")
        tree = ast.parse(source)
        calls = [ast.unparse(x) for x in ast.walk(tree)
                 if isinstance(x, ast.Call)]
        self.assertIn("system.load_checkpoint(args.checkpoint_in)", calls)
        self.assertIn("system.solve()", calls)
        checkpoint_saves = [
            call for call in calls
            if call.startswith("system.save_checkpoint(args.checkpoint_out,")
        ]
        self.assertEqual(len(checkpoint_saves), 1)
        self.assertIn("'physical_clock': 'absolute'", checkpoint_saves[0])
        self.assertIn("'gamma': args.gamma", checkpoint_saves[0])


if __name__ == "__main__":
    unittest.main()
