# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Source-only guards for reversible original Heat-Wave interface pilot."""

from pathlib import Path
import unittest


ROOT=Path(__file__).resolve().parents[1]
PILOT=ROOT/"sandbox/probe_v130_hw_dual_sided_assembly.py"
RUN=ROOT/"scripts/run-v130-hw-dual-sided-assembly.sh"


class OriginalHeatWaveDualFaceAssemblyTests(unittest.TestCase):
    def test_unmodified_source_constructor_interception_is_restored(self):
        s=PILOT.read_text(encoding="utf-8")
        for token in (
            "original_Domain=S.Domain",
            "original_Control_Port=S.Control_Port",
            "original_Brick=S.Brick",
            "original_init=DPHS.init_step",
            "S.Domain=domain_wrapper",
            "S.Control_Port=control_wrapper",
            "S.Brick=brick_wrapper",
            "DPHS.init_step=before_native_DAE",
            "S.Domain=original_Domain",
            "S.Control_Port=original_Control_Port",
            "S.Brick=original_Brick",
            "DPHS.init_step=original_init",
            "raise AssemblyOnlyComplete",
        ):
            self.assertIn(token,s)

    def test_only_original_named_interface_ports_and_bricks_retargeted(self):
        s=PILOT.read_text(encoding="utf-8")
        for token in (
            '("M_Y_T", "B_T", "B_T^T")',
            '("M_Y_w", "B_w", "B_w^T")',
            '"Interface Heat": 30',
            '"Interface Wave": 31',
            "cells(old)!=(cells(disk)|cells(wave))",
            "new_one_sided_regions(mesh)",
            "compare_block_norms(",
            '"physics_energy_power_qualified":False',
            '"DAE_initialized":False',
        ):
            self.assertIn(token,s)

    def test_runner_compares_both_original_assemblies_but_no_ts(self):
        s=RUN.read_text(encoding="utf-8")
        for token in (
            "for variant in baseline dual_sided",
            "--variant \"$variant\"",
            "SCRIMP_GETFEM_BRIDGE_DISTRIBUTED=coo",
            "check-v130-published-immutable.py",
            "nonzero_bidirectional_feedback_blocks",
            "DAE_initialization_executed",
            "energy_or_physical_coupling_qualified",
        ):
            self.assertIn(token,s)


if __name__=="__main__":
    unittest.main()
