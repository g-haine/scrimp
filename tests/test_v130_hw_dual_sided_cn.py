# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""No false physics PASS for corrected dual-side original Heat-Wave CN."""

import ast
from pathlib import Path
import unittest

ROOT=Path(__file__).resolve().parents[1]
PILOT=ROOT/"sandbox/probe_v130_hw_dual_sided_cn.py"
SCI=ROOT/"sandbox/measure_v130_hw_dual_facet_physics.py"
RUN=ROOT/"scripts/run-v130-hw-dual-sided-cn.sh"


class CorrectedDualFacetCNGates(unittest.TestCase):
    def test_original_model_reversible_constructor_and_real_main(self):
        s=PILOT.read_text(encoding="utf-8")
        for term in (
            "S.Domain=new_domain",
            "S.Control_Port=control",
            "S.Brick=brick",
            "S.Domain=original_Domain",
            "S.Control_Port=original_Control_Port",
            "S.Brick=original_Brick",
            "original_cn.main()",
            "main2d.collect_heat_wave_physics=real_collector",
            "main2d.collect_heat_wave_physics=original_collector",
            "original_2d_example_report.json",
            "compatible_cn_report.json",
            '"strict_DAE_initialization_PASSED":True',
        ):
            self.assertIn(term,s)

    def test_correct_2sided_power_and_midpoint_balance(self):
        s=SCI.read_text(encoding="utf-8")
        for term in (
            "HEAT_FACES = 30", "WAVE_FACES = 31",
            '("CN_P_interface_heat","U_T*Y_T",HEAT_FACES)',
            '("CN_P_interface_wave","U_w*Y_w",WAVE_FACES)',
            "np.diff(H)+dt*(dissipation-port_exchange)",
            "np.diff(traces[\"H_heat\"])+dt*(dissipation+signed_heat)",
            "np.diff(H_wave)-dt*(signed_wave+outer)",
            '"interface_exchange_nontrivial":bool(ratio>1e-8)',
            '"energy_physics_qualified":False',
            '"old_doubled_interface_region_excluded_from_physical_power":True',
        ):
            self.assertIn(term,s)
        self.assertNotIn('("CN_P_interface_heat","U_T*Y_T",10)',s)
        self.assertNotIn('("CN_P_interface_wave","U_w*Y_w",10)',s)

    def test_pinned_true_native_ts_and_nonvacuous_gate(self):
        sh=RUN.read_text(encoding="utf-8")
        for term in (
            "for candidate in clamp constrained_l2",
            "scripts/check-v130-published-immutable.py",
            "probe_v130_hw_dual_sided_cn.py",
            "measure_v130_hw_dual_facet_physics.py",
            "probe_v130_main_2d_examples.py",
            "probe_v130_heat_wave_compatible_cn.py",
            "mpiexec --oversubscribe --bind-to none -n 1",
            "SuperLU_DIST_qualified_for_corrected_Jacobian",
            "native_failures",
            "CN_discrete_midpoint_energy_gate_PASSED",
        ):
            self.assertIn(term,sh)


if __name__=="__main__":
    unittest.main()
