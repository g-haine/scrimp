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
        self.assertIn('for example in heat_wave shallow_water;',runner)
        self.assertIn('mpiexec --oversubscribe --bind-to none -n 1',runner)
        self.assertNotIn('mpiexec --oversubscribe --bind-to none -n 2',runner)
        self.assertIn("python3 scripts/check-v130-published-immutable.py",runner)
        self.assertNotIn("--volume \"$root/examples/published",runner)

    def test_original_2d_defaults_are_not_overridden_without_call_site(self):
        hw=(ROOT/"examples"/"heat_wave.py").read_text(encoding="utf-8")
        swe=(ROOT/"examples"/"shallow_water.py").read_text(encoding="utf-8")
        self.assertIn("final_time=2.0, time_step=0.001, save_step=0.05",hw)
        self.assertIn("factor_solver=None, final_time=0.5, time_step=0.0001,",swe)
        self.assertIn("init_step_ts_type=\"pseudo\"",hw)
        self.assertIn("init_step_ts_type=\"pseudo\"",swe)
        self.assertIn('factor_solver=args.backend',(
            ROOT/"sandbox"/"probe_v130_main_2d_examples.py"
        ).read_text(encoding="utf-8"))
        self.assertIn("original_FEM_and_weak_forms_unchanged",(
            ROOT/"sandbox"/"probe_v130_main_2d_examples.py"
        ).read_text(encoding="utf-8"))


if __name__=="__main__":
    unittest.main()
