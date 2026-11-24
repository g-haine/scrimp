# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Source-only guards for original Heat-Wave dual-facet strict DAE gate."""

import ast
from pathlib import Path
import unittest

ROOT=Path(__file__).resolve().parents[1]
PILOT=ROOT/"sandbox/probe_v130_hw_dual_sided_dae.py"
RUN=ROOT/"scripts/run-v130-hw-dual-sided-dae.sh"


class NativeDualFacetsDAE(unittest.TestCase):
    def test_import_closure_and_strict_native_initializer_intercept(self):
        text=PILOT.read_text(encoding="utf-8")
        tree=ast.parse(text)
        found={
            node.module for node in tree.body
            if isinstance(node,ast.ImportFrom)
        }
        self.assertIn("v130_hw_interface_geometry",found)
        self.assertIn("probe_v130_heat_wave_boundary_compatibility",found)
        self.assertNotIn("probe_v130_main_2d_examples",found)
        for token in (
            "original_init=DPHS.init_step",
            "DPHS.init_step=native_DAE_gate",
            "original_init(system)",
            "except InconsistentDAEInitialData as exc:",
            '"INCOMPATIBLE_INITIAL_DAE_UNALTERED_PHYSICS"',
            '"PASS_NATIVE_STRICT_DAE_FEASIBILITY_ONLY"',
            "raise DAEOnlyComplete",
            "DPHS.init_step=original_init",
        ):
            self.assertIn(token,text)

    def test_corrected_p_only_no_tol_or_TS_changes(self):
        text=PILOT.read_text(encoding="utf-8")
        for token in (
            "wave_boundary_partition(",
            "constrained_l2_field(",
            "p_corrected[b]=0.",
            'system.gf_model.set_variable(',
            "np.array_equal(actual,initial)",
            '"differential_initial_exact"',
            '"pseudo_ts_executed"',
            '"main_CN_TS_executed":False',
            '"physical_power_balance_qualified":False',
            '"M_Y_T", "B_T", "B_T^T"',
            '"M_Y_w", "B_w", "B_w^T"',
        ):
            self.assertIn(token,text)
        self.assertNotIn("init_algebraic_atol",text)
        self.assertNotIn("init_algebraic_rtol",text)

    def test_independent_native_runs_and_all_mounts(self):
        sh=RUN.read_text(encoding="utf-8")
        for token in (
            "baseline:clamp dual_sided:clamp dual_sided:constrained_l2",
            "v130_hw_interface_geometry.py:/workspace/scrimp/sandbox/",
            "probe_v130_heat_wave_boundary_compatibility.py:/workspace/scrimp/sandbox/",
            "probe_v130_hw_dual_sided_dae.py:/workspace/scrimp/sandbox/",
            "scripts/check-v130-published-immutable.py",
            "strict_DAE.json",
            "DUAL_FACET_DAE_INITIAL_INCOMPATIBILITY_DIAGNOSED",
            '"physical_energy_balance_qualified":False',
        ):
            self.assertIn(token,sh)
        self.assertNotIn("probe_v130_main_2d_examples",sh)


if __name__=="__main__":
    unittest.main()
