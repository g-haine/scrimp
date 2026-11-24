# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# See LICENSE for license information.

"""Pure source guards for FULL physical refined SuperLU_DIST 2D wave.

These are not native PETSc/GetFEM tests and must not be described as
scientific acceptance. They guard the experiment's exact executable
code, including the embedded Python archive comparator.
"""

import ast
from pathlib import Path
import re
import unittest


ROOT = Path(__file__).resolve().parents[1]
DRIVER = ROOT/"sandbox"/"probe_v130_wave2d_refined_superlu_dist.py"
RUNNER = ROOT/"scripts"/"run-v130-wave2d-refined-superlu-dist-physical.sh"
CANONICAL = ROOT/"examples"/"wave.py"


class RefinedWave2DPhysicalSourceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.driver = DRIVER.read_text(encoding="utf-8")
        cls.runner = RUNNER.read_text(encoding="utf-8")
        cls.wave = CANONICAL.read_text(encoding="utf-8")

    def test_driver_valid_python_and_exactly_two_optin_options(self):
        ast.parse(self.driver, filename=str(DRIVER))
        self.assertIn("SCRIMP_V130_WAVE2D_ENABLE_REFINED_DIST", self.driver)
        self.assertIn(
            'os.environ.get("SCRIMP_V130_WAVE2D_ENABLE_REFINED_DIST") != "1"',
            self.driver,
        )
        self.assertIn('mat_superlu_dist_replacetinypivot', self.driver)
        self.assertIn('mat_superlu_dist_iterrefine', self.driver)
        self.assertEqual(
            self.driver.count('    "mat_superlu_dist_replacetinypivot",'), 1
        )
        self.assertEqual(self.driver.count('    "mat_superlu_dist_iterrefine",'), 1)
        self.assertIn("system.time_scheme[option] = \"true\"", self.driver)
        self.assertIn("system.time_scheme.getBool(option, False)", self.driver)

    def test_injection_after_original_time_scheme_before_solve(self):
        self.assertIn("original = S.DPHS.set_time_scheme", self.driver)
        self.assertIn("result = original(system, **kwargs)", self.driver)
        self.assertIn('actual_solver != "superlu_dist"', self.driver)
        self.assertIn('system.time_scheme.getString(', self.driver)
        self.assertIn('S.DPHS.set_time_scheme = enable_refined_superlu_dist', self.driver)
        self.assertIn("from examples.wave import main as canonical_wave_main",
                      self.driver)
        self.assertIn("canonical_wave_main()", self.driver)
        self.assertIn("invoked[0] != 1", self.driver)
        self.assertIn('getattr(system, "_v130_refined_dist_enabled", False)',
                      self.driver)
        self.assertIn("system._v130_refined_dist_enabled = True", self.driver)
        self.assertIn("S.DPHS.set_time_scheme = original", self.driver)
        self.assertNotIn("S.Brick(",self.driver)
        self.assertNotIn("set_initial_value(",self.driver)

    def test_runner_has_valid_python_reporting_heredoc(self):
        # One embedded reporter after the two genuine native runs.
        matches = re.findall(
            r"<<'PY'\n(.*?)\nPY\n", self.runner, flags=re.DOTALL
        )
        self.assertEqual(len(matches),1)
        ast.parse(matches[0],filename=str(RUNNER)+" <PY report>")
        self.assertIn("import zipfile",matches[0])
        self.assertIn('status',matches[0])

    def test_runner_reuses_exact_native_mumps_baseline_and_only_dist(self):
        self.assertIn(
            "v130-wave2d-conservative-five-lu-99d99f2d-tEOjqgHy.zip",
            self.runner,
        )
        self.assertIn('for ranks in 1 2; do',self.runner)
        self.assertIn('name="refined_dist_mpi',self.runner)
        self.assertIn('SCRIMP_V130_WAVE2D_ENABLE_REFINED_DIST=1',self.runner)
        self.assertIn('--solver superlu_dist',self.runner)
        self.assertIn('--final-time 0.08 --time-step 0.008',self.runner)
        self.assertIn('--check --export --output /artifacts',self.runner)
        self.assertIn('mpiexec --oversubscribe --bind-to none -n "$ranks"',
                      self.runner)
        self.assertNotIn('--solver mumps',self.runner)
        self.assertNotIn('--solver superlu ',self.runner)
        self.assertNotIn('mat_superlu_dist_rowperm',self.runner)

    def test_true_FE_4ports_energy_and_independent_mumps_parity(self):
        for text in (
            '"power_bottom"', '"power_right"', '"power_top"',
            '"power_left"', '"signed_work"', '"energy_step_defects"',
            '"t"', '"z"', '"H"',
            'max_abs_step_defect', 'max_abs_cumulative_defect',
            'max_abs_left_power','physical_gates',
            'DAE_initialization','differential_initial_exact',
            'requested_LU', 'selected_LU',
            'np.arange(11)*.008', 'state_SHA256',
            'np.ascontiguousarray(z).tobytes()',
            'np.allclose(ref[key],baseline[key],rtol=2e-8,atol=5e-9)',
            'rank_trajectories_agree',
        ):
            with self.subTest(token=text):
                self.assertIn(text,self.runner)

    def test_mpi_VTK_all_original_fields_and_geometry(self):
        for expected in (
            'for field in ("q","p","e_q","e_p")',
            'native_case_pvd(case,field,ref["t"])',
            'baseline_pvd(arch,field,baseline["t"])',
            'if grids != base_geometry:',
            'if len(entries) != 11:',
            'piece.find(f"./PointData/DataArray[@Name=',
            'if len(set(signatures)) != 1:',
            'report["status"]="PASS_OPTIN_REFINED_SUPERLU_DIST_FULL_WAVE2D_MPI1_MPI2_CANDIDATE_ONLY"',
        ):
            with self.subTest(expected=expected):
                self.assertIn(expected,self.runner)

    def test_immutable_public_example_and_failure_reporting(self):
        ast.parse(self.wave,filename=str(CANONICAL))
        self.assertIn("def wave_eq(", self.wave)
        self.assertIn('"0.5*q.T.q"', self.wave)
        self.assertIn('init_step_ts_type="algebraic"', self.wave)
        self.assertIn('"0.2*sin(4*pi*t)*sin(pi*y)"', self.wave)
        self.assertIn('superlu_dist_default_backend_qualified": False',
                      self.runner)
        self.assertIn('original_five_default_backends_gate": "3/5"',
                      self.runner)
        self.assertIn('if report["status"]!="PASS_OPTIN_REFINED_SUPERLU_DIST_FULL_WAVE2D_MPI1_MPI2_CANDIDATE_ONLY":',
                      self.runner)
        self.assertIn('zipfile.ZipFile(archive,"x"',self.runner)


if __name__ == "__main__":
    unittest.main()
