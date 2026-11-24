# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Source-only audit of opt-in SuperLU_DIST TS failure capture."""
from pathlib import Path
import unittest

ROOT=Path(__file__).resolve().parents[1]


class OriginalSuperLUStepDiagnosticsTests(unittest.TestCase):
    def test_actual_ts_error_diagnostic_is_observer_only(self):
        source=(ROOT/"scrimp/dphs.py").read_text(encoding="utf-8")
        for text in (
            'os.environ.get("SCRIMP_TS_FAILURE_DIAGNOSTICS") == "1"',
            "MAIN TS FAILURE DIAGNOSTIC:",
            "snes.getConvergedReason",
            "ksp.getConvergedReason",
            "pc.getFailedReason",
            "pc.getFactorSolverType",
        ):
            self.assertIn(text,source)

    def test_failure_capture_does_not_approve_physics(self):
        source=(ROOT/"sandbox/probe_v130_heat_wave_compatible_cn.py").read_text(encoding="utf-8")
        for text in (
            "--capture-ts-failure","except PETSc.Error as exc:",
            '"DIVERGED_STEP_REJECTED" not in str(exc)',
            'native["status"]!="PASS"',
            '"physics_qualified":False',
            '"snes_monitor"',
        ):
            self.assertIn(text,source)

    def test_real_original_image_diagnostic_scope(self):
        runner=(ROOT/"scripts/run-v130-heat-wave-superlu-step-diagnostic.sh").read_text(encoding="utf-8")
        for text in (
            "for candidate in clamp constrained_l2",
            "--backend superlu_dist","--capture-ts-failure",
            "SCRIMP_TS_FAILURE_DIAGNOSTICS=1",
            "SCRIMP_HW_TS_DIAGNOSTIC_MONITORS=1",
            "check-v130-published-immutable.py",
            "DIAGNOSTIC_ONLY_TS_FAILURE_REPRODUCED",
        ):
            self.assertIn(text,runner)


if __name__=="__main__":
    unittest.main()
