# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Source-only contracts of full original Heat-Wave CN LU parity campaign."""

from pathlib import Path
import unittest


ROOT=Path(__file__).resolve().parents[1]
RUNNER=ROOT/"scripts/run-v130-hw-cn-matched-lu-parity.sh"
COMPARE=ROOT/"sandbox/compare_v130_hw_cn_matched_lu_parity.py"


class OriginalHeatWaveFullCNParityTests(unittest.TestCase):
    def test_four_unmodified_original_heat_wave_runs_are_independent(self):
        sh=RUNNER.read_text(encoding="utf-8")
        for phrase in (
            "for candidate in clamp constrained_l2",
            "for solver in mumps matched_superlu",
            "--candidate \"$candidate\" --backend mumps",
            "--candidate \"$candidate\" --out /artifacts",
            "probe_v130_heat_wave_compatible_cn.py",
            "probe_v130_hw_matched_pc_cn.py",
            "check-v130-published-immutable.py",
            "SCRIMP_GETFEM_BRIDGE_DISTRIBUTED=coo",
            "sandbox/compare_v130_hw_cn_matched_lu_parity.py",
            "--workdir /scrimp-work",
        ):
            self.assertIn(phrase,sh)

    def test_complete_original_pH_traces_and_nonvacuity_separated(self):
        src=COMPARE.read_text(encoding="utf-8")
        for phrase in (
            '"times", "states", "H", "T", "p", "q"',
            '"phys_CN_D_heat_flux"',
            '"phys_CN_P_interface_heat"',
            '"phys_CN_P_interface_wave"',
            '"phys_CN_per_step_H_balance_error"',
            "compare_complete_traces(",
            "interface_exchange_ratio_to_heat_dissipation",
            "interface_exchange_nontrivial",
            "interface_energy_exchange_excited_all_runs",
            "nodal_trapezoid_not_used_for_CN_energy_acceptance",
            "failed_observable_backend_parity",
            '"rowmatched_PC":matched_report is not None',
        ):
            self.assertIn(phrase,src)

    def test_no_unqualified_physical_coupling_pass(self):
        src=COMPARE.read_text(encoding="utf-8")
        self.assertIn(
            "PASS_FOUR_CN_BACKEND_PARITY_AND_ENERGY_INTERFACE_EXCITATION_OPEN",
            src,
        )
        self.assertIn("interface_ratio>1e-8",src)
        self.assertIn("not all_interface_excited",src)
        self.assertIn('"MPI2_or_parallel_export_qualified":False',src)


if __name__=="__main__":
    unittest.main()
