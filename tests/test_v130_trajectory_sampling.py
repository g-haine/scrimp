# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Source-only regressions for dropped genuine PETSc TS monitor samples.

The long 1D wave run on 24 November 2026 accepted all 150 TSCN
steps without a solver failure, but monitor i=149, t≈1.49
did NOT save the accepted state because (previous+dt_save)>t
by one representable floating-point spacing. The monitor then
saved i=150, creating a spurious TWO-step interval among CN
step powers; physics validation correctly refused that archive.

These tests check the ACTUAL imported utility used by DPHS.monitor
without requiring GetFEM or PETSc. They are not native physics proof.
"""

import ast
import importlib.util
import math
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]
FILE = ROOT/"scrimp"/"utils"/"trajectory_sampling.py"
spec = importlib.util.spec_from_file_location(
    "scrimp_v130_actual_ts_sampling", FILE
)
sampling = importlib.util.module_from_spec(spec)
spec.loader.exec_module(sampling)


class AcceptedTSSampling(unittest.TestCase):
    def test_reproduce_exact_long_wave_monitor_missed_step_149(self):
        accepted_time = math.nextafter(1.49, -math.inf)
        next_scheduled = 1.48 + .01
        # Previous SCRIMP monitor used this exact, wrong comparison.
        self.assertGreater(next_scheduled-accepted_time, 0.0)
        self.assertFalse(next_scheduled-accepted_time <= 0.0)
        # New helper MUST save that actual accepted FE state.
        self.assertTrue(sampling.saved_time_is_due(
            accepted_time, next_scheduled, .01
        ))

    def test_151_real_accepted_endpoints_all_preserved(self):
        times = [i*.01 for i in range(151)]
        times[149] = math.nextafter(1.49, -math.inf)
        stored = []
        for t in times:
            next_scheduled = stored[-1]+.01 if stored else 0.0
            if not stored or sampling.saved_time_is_due(
                t, next_scheduled, .01
            ):
                stored.append(t)
        self.assertEqual(len(stored), 151)
        self.assertEqual(stored[149], times[149])
        self.assertEqual(stored[150], 1.5)

    def test_genuine_later_step_is_saved_as_before(self):
        self.assertTrue(sampling.saved_time_is_due(
            1.5, 1.49, .01
        ))

    def test_small_roundoff_under_a_due_time(self):
        saved = math.nextafter(.04, -math.inf)
        self.assertTrue(sampling.saved_time_is_due(saved, .04, .004))

    def test_early_real_physical_step_is_not_saved(self):
        self.assertFalse(sampling.saved_time_is_due(
            .019, .020, .01
        ))

    def test_not_a_time_interpolator(self):
        t = 1.49 - 1.e-6
        self.assertFalse(sampling.saved_time_is_due(t, 1.49, .01))

    def test_deliberate_subsampling_keeps_one_in_four(self):
        times = [i*.01 for i in range(9)]
        saved = []
        for t in times:
            target = saved[-1]+.04 if saved else 0.0
            if not saved or sampling.saved_time_is_due(t, target, .04):
                saved.append(t)
        self.assertEqual(saved, [0.0, 0.04, 0.08])

    def test_cap_tolerance_for_tiny_dt_save(self):
        target = 1.0
        # 16 ULPs would otherwise be huge compared with dt_save.
        self.assertFalse(sampling.saved_time_is_due(
            target-1.e-14, target, 1.e-10
        ))

    def test_nonfinite_and_nonpositive_spacing_fail_closed(self):
        for args in (
            (float("nan"), .02, .01),
            (.01, float("inf"), .01),
            (.01, .02, 0.),
            (.01, .02, -.01),
        ):
            with self.subTest(args=args):
                with self.assertRaisesRegex(ValueError, "finite"):
                    sampling.saved_time_is_due(*args)

    def test_fixed_CN_each_accepted_step_is_stored_at_real_long_step_2999(self):
        # Reproduce the specific observed Heat-Wave failure: accepted
        # step 2999 is slightly earlier than last_saved_t + dt_save.
        # The genuine PETSc monitor did NOT save step 2999 although
        # 3000 CN steps converged. This test changes no TS physics.
        dt=.005
        t_2998=14.99
        t_2999=14.995-1.e-12
        target=t_2998+dt
        self.assertFalse(sampling.saved_time_is_due(
            t_2999,target,dt
        ))
        self.assertTrue(sampling.fixed_cn_accepted_step_is_due(
            2999,dt=dt,dt_save=dt,
            ts_type="cn",ts_adapt_type="none"
        ))
        self.assertEqual(
            sum(sampling.fixed_cn_accepted_step_is_due(
                    i,dt=dt,dt_save=dt,
                    ts_type="cn",ts_adapt_type="none"
                ) for i in range(1,3001)),
            3000,
        )

    def test_fixed_CN_does_not_change_other_sampling_policies(self):
        cases = (
            (1,.005,.05,"cn","none"),   # deliberately subsampled
            (1,.005,.005,"cn","basic"), # adapted TS
            (1,.005,.005,"bdf","none"), # different scheme
            (1,.005,.004,"cn","none"),  # unequal sampling
            (-1,.005,.005,"cn","none"), # synthetic final callback
            (0,.005,.005,"cn","none"),  # initial TS state
        )
        for step,dt,dt_save,scheme,adapt in cases:
            with self.subTest(case=(step,dt,dt_save,scheme,adapt)):
                self.assertFalse(sampling.fixed_cn_accepted_step_is_due(
                    step,dt=dt,dt_save=dt_save,
                    ts_type=scheme,ts_adapt_type=adapt
                ))

    def test_fixed_CN_invalid_step_sizes_fail_closed(self):
        for bad in (0.,float("nan"),float("inf"),-1.):
            with self.subTest(bad=bad):
                with self.assertRaisesRegex(ValueError,"positive and finite"):
                    sampling.fixed_cn_accepted_step_is_due(
                        12,dt=.005,dt_save=bad,
                        ts_type="cn",ts_adapt_type="none"
                    )

    def test_actual_dphs_monitor_invokes_roundoff_safe_helper(self):
        source = (ROOT/"scrimp"/"dphs.py").read_text()
        tree = ast.parse(source)
        cls = next(x for x in tree.body
                   if isinstance(x, ast.ClassDef) and x.name == "DPHS")
        monitor = next(x for x in cls.body
                       if isinstance(x, ast.FunctionDef)
                       and x.name == "monitor")
        calls = [node for node in ast.walk(monitor)
                 if isinstance(node, ast.Call)
                 and isinstance(node.func, ast.Name)]
        self.assertIn("saved_time_is_due",
                      [node.func.id for node in calls])
        self.assertIn("fixed_cn_accepted_step_is_due",
                      [node.func.id for node in calls])
        self.assertNotIn("next_saved_at_t - t <= 0", ast.get_source_segment(
            source, monitor
        ))


if __name__ == "__main__":
    unittest.main()
