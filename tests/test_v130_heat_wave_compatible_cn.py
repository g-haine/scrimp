# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Source-only checks of real CN candidate wrapper and physics gates."""

import ast
from pathlib import Path
import unittest

ROOT=Path(__file__).resolve().parents[1]
PILOT=ROOT/"sandbox/probe_v130_heat_wave_compatible_cn.py"
RUNNER=ROOT/"scripts/run-v130-heat-wave-compatible-cn.sh"


class HeatWaveCompatibleCNSourceTests(unittest.TestCase):
    def test_original_main_physics_pilot_is_called_after_strict_DAE(self):
        s=PILOT.read_text(encoding="utf-8")
        tree=ast.parse(s)
        names={n.name for n in tree.body if isinstance(n,ast.FunctionDef)}
        self.assertTrue({"main","add_boundary_compatible_state"} <= names)
        for literal in (
            "original_init(system)",
            "original_main.main()",
            "DPHS.init_step=before_physical_TS",
            "DPHS.init_step=original_init",
            '"T",original_T',
            '"q",original_q',
            'original_main.main()',
            "physical_balance_gate_passed",
        ):
            self.assertIn(literal,s)
        self.assertNotIn('init_step_ts_type="pseudo"',s)
        self.assertIn('native["pseudo_ts_executed"]',s)
    def test_both_original_LU_and_boundary_candidates_present(self):
        r=RUNNER.read_text(encoding="utf-8")
        for text in (
            "for candidate in clamp constrained_l2",
            "for backend in mumps superlu_dist",
            'SCRIMP_GETFEM_MATRIX_MODE=distributed',
            'SCRIMP_GETFEM_BRIDGE_DISTRIBUTED=coo',
            'original_2d_example_trace.npz',
            'phys_CN_per_step_H_balance_error',
            'check-v130-published-immutable.py',
            'mpiexec --oversubscribe --bind-to none -n 1',
        ):
            self.assertIn(text,r)
    def test_no_automatic_initial_projection_policy_change(self):
        model=(ROOT/"examples/heat_wave.py").read_text(encoding="utf-8")
        self.assertIn(
            'hw.set_initial_value("p", "5.*np.exp(-25*((x-0.6)*(x-0.6)+y*y))")',
            model,
        )
        self.assertIn('init_step_ts_type="algebraic"',model)
        self.assertIn("original_init(system)",PILOT.read_text(encoding="utf-8"))


if __name__=="__main__":
    unittest.main()
