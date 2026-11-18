# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Synthetic visual regression for v1.3.0 PNG/SVG benchmark outputs.

No timings below are presented as actual measurements.
"""
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

SOURCE = Path(__file__).resolve().parents[1] / "sandbox" / "analyze_v130_wave_benchmark.py"
SPEC = importlib.util.spec_from_file_location("v130_wave_plots", SOURCE)
PLOT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(PLOT)


def _row(p):
    return {
        "ranks": p,
        "solve_s_median": 3.0 / p + 0.2,
        "total_s_median": 3.4 / p + 0.3,
        "first_CN_step_s_median": .1 / p,
        "warm_CN_step_s_median": .04 / p,
        "speedup_from_current_1r": 3.2 / (3.0 / p + .2),
        "parallel_efficiency": 3.2 / (p * (3.0 / p + .2)),
        "initial.consistency_total": .1 / p,
        "main.petsc.ts_solve": .8 / p,
        "main.petsc.callback.ijacobian": .12 / p,
        "main.petsc.callback.ifunction": .06 / p,
        "main.petsc.csr_pattern_check": .04 / p,
        "main.MatLUFactorSym.seconds_median": .25 / p,
        "main.MatLUFactorNum.seconds_median": .14 / p,
        "main.MatSolve.seconds_median": .09 / p,
    }


class WaveV130PlotTests(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.directory = Path(tmp.name)

    def test_scaling_without_release_parity_omits_version_claim(self):
        names = PLOT.publish_wave_figures(
            self.directory, [_row(p) for p in (1, 2, 4, 6)],
            serial_release=80., serial_current=3.2,
            release_parity_accepted=False,
        )
        self.assertEqual(len(names), 6)
        self.assertFalse(
            (self.directory / "png" / "01_sequential_v120_vs_v130.png").exists()
        )
        for name in names:
            self.assertGreater(
                (self.directory / "png" / (name + ".png")).stat().st_size,
                1000,
            )
            self.assertGreater(
                (self.directory / "svg" / (name + ".svg")).stat().st_size,
                1000,
            )

    def test_validated_sequential_release_comparison_is_visualized(self):
        names = PLOT.publish_wave_figures(
            self.directory, [_row(p) for p in (1, 2)],
            serial_release=80., serial_current=3.2,
            release_parity_accepted=True,
        )
        self.assertEqual(len(names), 7)
        self.assertGreater(
            (self.directory / "png" / "01_sequential_v120_vs_v130.png").stat().st_size,
            1000,
        )

    def test_short_numeric_reuse_smoke_has_two_readable_figures(self):
        for n in (1, 2):
            name = f"current_on_change_{n}r"
            (self.directory / (name + ".json")).write_text(
                json.dumps({"wall_solve_s": 1.3 / n}), encoding="utf-8",
            )
            (self.directory / (name + ".log")).write_text(
                "Event Stage 1: WAVE_INITIAL_CONSISTENCY\n"
                "MatLUFactorNum 1 1.0 .01 extra\n"
                "Event Stage 2: WAVE_MAIN_CN\n"
                "MatLUFactorSym 1 1.0 .02 extra\n"
                "MatLUFactorNum 1 1.0 .03 extra\n"
                "MatSolve 10 1.0 .01 extra\n",
                encoding="utf-8",
            )
        names = PLOT.publish_smoke_figures(self.directory, [1, 2])
        self.assertEqual(len(names), 2)
        for name in names:
            self.assertGreater(
                (self.directory / "png" / (name + ".png")).stat().st_size,
                1000,
            )


if __name__ == "__main__":
    unittest.main()
