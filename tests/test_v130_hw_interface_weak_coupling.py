# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Source-only regressions for original 2D mixed interface audit."""

import ast
from pathlib import Path
import unittest

ROOT=Path(__file__).resolve().parents[1]
PROBE=ROOT/"sandbox/probe_v130_hw_interface_weak_coupling.py"
RUNNER=ROOT/"scripts/run-v130-hw-interface-weak-coupling.sh"


class InterfaceWeakCouplingDiagnosticTests(unittest.TestCase):
    def test_actual_getfem_weak_port_and_power_traces(self):
        src=PROBE.read_text(encoding="utf-8")
        for phrase in (
            '"U_T", "Y_w", "U_w", "Y_T", "T", "p"',
            '"heat_interface_local_power_density_squared"',
            '"wave_interface_local_power_density_squared"',
            "region=REGION_INTERFACE,CN=midpoint",
            "original_collector(system,times,H,*a,**kw)",
            "original_main" if False else "base.collect_heat_wave_physics",
            '"full_multiphysics_coupling_qualified":False',
            '"raw_trace_mismatch_is_not_proof_of_weak_FEM_constraint_failure":True',
        ):
            self.assertIn(phrase,src)

    def test_sparse_assembled_original_port_operator_blocks(self):
        src=PROBE.read_text(encoding="utf-8")
        for phrase in (
            'system.stiffness', 'system.gf_model.interval_of_variable(name)',
            '"p_to_Y_w_trace"', '"e_Q_to_Y_T_flux"',
            '"Y_w_to_U_T_feedback"', '"Y_T_to_U_w_feedback"',
            '"U_T_to_heat_flux_equation"',
            '"U_w_to_wave_velocity_equation"',
            '"algebraic_Kz_row_residuals"',
            '"port_coefficient_norms"',
        ):
            self.assertIn(phrase,src)

    def test_runner_single_original_mumps_cn_and_immutability(self):
        sh=RUNNER.read_text(encoding="utf-8")
        for phrase in (
            "--candidate clamp --out /artifacts",
            "--env SCRIMP_GETFEM_BRIDGE_DISTRIBUTED=coo",
            "probe_v130_heat_wave_compatible_cn.py",
            "check-v130-published-immutable.py",
            "PASS_DIAGNOSTIC_CAPTURED_NOT_COUPLING_QUALIFIED",
            "mpiexec --oversubscribe --bind-to none -n 1",
        ):
            self.assertIn(phrase,sh)
        self.assertNotIn("matched_superlu",sh)


if __name__=="__main__":
    unittest.main()
