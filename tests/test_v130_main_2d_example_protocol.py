# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Negative controls for source-only 2D main-example test protocol.

No PETSc/GetFEM runtime is available in this CI. Native physics is
not called a PASS by these source-only contract tests.
"""

from pathlib import Path
import unittest


ROOT=Path(__file__).resolve().parents[1]


class Main2DExampleProtocolTests(unittest.TestCase):
    def test_numerical_solver_compare_is_real_original_examples(self):
        pilot=(ROOT/"sandbox"/"probe_v130_main_2d_examples.py").read_text(
            encoding="utf-8"
        )
        runner=(ROOT/"scripts"/"run-v130-example-2d-lu.sh").read_text(
            encoding="utf-8"
        )
        self.assertIn("from examples.heat_wave import heat_wave_eq",pilot)
        self.assertIn("from examples.shallow_water import shallow_water_eq",pilot)
        self.assertIn("system=heat_wave_eq(",pilot)
        self.assertIn("system=shallow_water_eq(",pilot)
        self.assertIn("if comm.size!=1:",pilot)
        self.assertIn("system.get_Hamiltonian()",pilot)
        self.assertIn("system.get_solution(name)",pilot)
        self.assertIn('system.get_quantity("h",region=1)',pilot)
        self.assertIn("original_2d_example_trace.npz",pilot)
        self.assertIn("check_output_files(args.example)",pilot)
        self.assertIn('for backend in mumps superlu_dist;',runner)
        self.assertIn('for example in "${examples[@]}";',runner)
        self.assertIn('scope="${SCRIMP_V130_2D_SCOPE:-heat_wave}"',runner)
        self.assertIn('all) examples=(heat_wave shallow_water) ;;',runner)
        self.assertIn('mpiexec --oversubscribe --bind-to none -n 1',runner)
        self.assertNotIn('mpiexec --oversubscribe --bind-to none -n 2',runner)
        self.assertIn("python3 scripts/check-v130-published-immutable.py",runner)
        self.assertNotIn("--volume \"$root/examples/published",runner)

    def test_adaptive_BDF2_accepted_steps_and_diagnostics_are_explicit(self):
        probe=(ROOT/"sandbox"/"probe_v130_main_2d_examples.py").read_text(
            encoding="utf-8"
        )
        physics=(ROOT/"sandbox"/"measure_v130_heat_wave_physics.py").read_text(
            encoding="utf-8"
        )
        runner=(ROOT/"scripts"/"run-v130-example-2d-lu.sh").read_text(
            encoding="utf-8"
        )
        self.assertIn("np.all(np.diff(times)>0.)",probe)
        self.assertIn("system.last_ts_steps < times.size-1",probe)
        self.assertIn('"requested_initial_dt":DT',probe)
        self.assertIn("collect_heat_wave_physics(",probe)
        # Only the heat-wave CN branch now deliberately requires
        # fixed four-step integration; shallow-water remains adaptive.
        self.assertIn('if args.example=="heat_wave":',probe)
        self.assertIn('system.last_ts_steps!=4',probe)
        self.assertIn('system.last_ts_rejections!=0',probe)
        self.assertIn('"H_heat": (".5*T*T", 1)',physics)
        self.assertIn('"H_wave_kinetic": (".5*p*p", 2)',physics)
        self.assertIn('"H_wave_potential": (".5*q.q", 2)',physics)
        self.assertIn('"D_heat_flux": ("e_Q.e_Q", 1)',physics)
        self.assertIn('"P_interface_heat": ("U_T*Y_T", 10)',physics)
        self.assertIn('"P_interface_wave": ("U_w*Y_w", 10)',physics)
        self.assertIn('"interface_T_minus_p_L2_sq"',physics)
        self.assertIn('"interface_T_plus_p_L2_sq"',physics)
        self.assertIn('python3 - "$out" "${examples[@]}"',runner)
        self.assertIn("ENDPOINTS_ONLY_DIFFERENT_ADAPTIVE_CLOCKS",runner)
        self.assertIn("PHYSICAL_AUDIT_OPEN",runner)
        self.assertNotIn('report["steps"]==4',runner)
        self.assertNotIn('trace["times"].shape==(5,)',runner)

    def test_original_2d_defaults_are_not_overridden_without_call_site(self):
        hw=(ROOT/"examples"/"heat_wave.py").read_text(encoding="utf-8")
        swe=(ROOT/"examples"/"shallow_water.py").read_text(encoding="utf-8")
        self.assertIn("final_time=2.0, time_step=0.001, save_step=0.05",hw)
        self.assertIn("factor_solver=None, final_time=0.5, time_step=0.0001,",swe)
        self.assertIn("init_step_ts_type=\"algebraic\"",hw)
        self.assertIn("init_step_ts_type=\"pseudo\"",swe)
        self.assertIn('ts_type="cn"',hw)
        self.assertIn('ts_adapt_type="none"',hw)
        self.assertIn("S.FEM(\"Heat flux\", k, \"RT\")",hw)
        self.assertIn("S.FEM(\"q\", k, \"RT\")",hw)
        self.assertIn("S.FEM(\"T\", k, \"DG\")",hw)
        self.assertIn("S.FEM(\"p\", k + 1, \"CG\")",hw)
        self.assertIn("S.FEM(\"Interface Heat\", k + 1, \"CG\")",hw)
        self.assertIn("S.FEM(\"Interface Wave\", k + 1, \"CG\")",hw)
        self.assertIn("S.FEM(\"Boundary\", k + 1, \"CG\")",hw)
        self.assertIn('k = 0',hw)
        self.assertIn('"U_bnd*Test_Y_bnd"',hw)
        self.assertIn('"Y_bnd*Test_p"',hw)
        self.assertIn('"-p*Test_Y_bnd"',hw)
        self.assertIn('"Y_T" if heat_normal_orientation == 1 else "-Y_T"',hw)
        self.assertIn('if FEM.get_name() in ("Heat flux", "q")',hw)
        self.assertIn('factor_solver=args.backend',(
            ROOT/"sandbox"/"probe_v130_main_2d_examples.py"
        ).read_text(encoding="utf-8"))
        self.assertIn("original_FEM_and_weak_forms_unchanged",(
            ROOT/"sandbox"/"probe_v130_main_2d_examples.py"
        ).read_text(encoding="utf-8"))


if __name__=="__main__":
    unittest.main()
