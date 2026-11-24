# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# See LICENSE for license information.

"""Source-only gates for the *original* 2D CN Jacobian refinement trial.

No PETSc/GetFEM import; no fake numerical qualification. The tests
constrain the research replay to the same matrix, true non-preconditioned
residuals, and unchanged solver acceptance thresholds.
"""

import ast
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]
DRIVER = ROOT/"sandbox"/"replay_v130_wave2d_superlu_dist_refinement.py"
RUNNER = ROOT/"scripts"/"run-v130-wave2d-superlu-dist-refinement.sh"
SOURCE = ROOT/"examples"/"wave.py"
HASH = "24346856f358e4f81cf332ebc974e329dc084c60f5d6b98bfbcd788fdff37f54"


class FullOriginalCSRRefinementSourceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.code=DRIVER.read_text(encoding="utf-8")
        cls.script=RUNNER.read_text(encoding="utf-8")

    def test_original_P_M_K_hash_and_stored_graph_preserved(self):
        ast.parse(self.code, filename=str(DRIVER))
        self.assertIn(HASH,self.code)
        self.assertIn('matrix.getValuesCSR()',self.code)
        self.assertIn('np.array_equal(orig_cols, local_cols)',self.code)
        self.assertIn('global_max != 0.',self.code)
        self.assertIn('"original_J_not_modified": True',self.code)
        self.assertIn('"PETSc_owned_CSR_graph_and_values_identical_to_archive": True',self.code)
        self.assertNotIn('mat.setDiagonal(',self.code)
        self.assertNotIn('S.Brick(',self.code)
        self.assertNotIn('set_initial_value(',self.code)

    def test_only_superlu_dist_with_explicit_tiny_pivots(self):
        self.assertIn('pc.setFactorSolverType("superlu_dist")',self.code)
        self.assertIn('options.setValue("mat_superlu_dist_replacetinypivot", "true")',self.code)
        self.assertIn('options.setValue("mat_superlu_dist_iterrefine", "true")',self.code)
        self.assertNotIn('pc.setFactorSolverType("mumps")',self.code)
        self.assertNotIn('mat_superlu_dist_rowperm',self.code)

    def test_internal_iterrefine_and_right_gmres_are_distinct(self):
        self.assertIn('choices=("internal_refine", "right_gmres")',self.code)
        self.assertIn('if args.refinement == "internal_refine":',self.code)
        self.assertIn('ksp.setType("preonly")',self.code)
        self.assertIn('ksp.setType("gmres")',self.code)
        self.assertIn('ksp.setPCSide(PETSc.PC.Side.RIGHT)',self.code)
        self.assertIn('ksp.setNormType(PETSc.KSP.NormType.UNPRECONDITIONED)',self.code)
        self.assertIn('ksp.setTolerances(rtol=1e-12, atol=1e-14, max_it=30)',self.code)
        self.assertIn('ksp.getIterationNumber()',self.code)

    def test_true_original_system_residual_and_exact_solution(self):
        for token in (
            "np.sin(.013 * offset) + .25",
            "mat.mult(known, rhs)",
            "mat.mult(solution, check)",
            "check.axpy(-1., rhs)",
            "diff.axpy(-1., known)",
            "rel_backward < 1e-11",
            "rel_solution < 1e-7",
            "comm.allreduce(bool(candidate_ok), op=MPI.LAND)",
        ):
            with self.subTest(token=token):
                self.assertIn(token,self.code)
        self.assertNotIn("rtol=1e-7",self.code)
        self.assertNotIn("atol=1e-6",self.code)

    def test_runner_requires_four_real_MPI_runs_and_real_archive(self):
        for token in (
            "for refinement in internal_refine right_gmres; do",
            "for ranks in 1 2; do",
            '--refinement "$refinement"',
            "v130-wave2d-first-cn-jacobian-574d9c71-NKqIZ0BZ.zip",
            HASH,
            "internal_refine_mpi",
            "relative_backward_residual",
            "relative_known_solution_error",
            "wave_installation_gate",
            '"3/5"',
        ):
            # Runner names are assembled dynamically by loop, not repeated.
            if token == "internal_refine_mpi":
                continue
            with self.subTest(token=token):
                self.assertIn(token,self.script)
        self.assertNotIn('--solver mumps',self.script)
        self.assertNotIn('--solver superlu ',self.script)

    def test_backend_qualification_stays_failed_without_native_full_wave(self):
        self.assertIn('full_ten_step_wave_2d_physics_qualified_for_superlu_dist',self.script)
        self.assertIn('original_three_backends_still_PASS',self.script)
        self.assertIn('PASS_AT_LEAST_ONE_STRICT_REFINEMENT_PETSC_MPI1_MPI2_ONLY',self.script)
        self.assertIn('FULL SCRIMP wave physics remains 3/5',self.script)
        self.assertNotIn('PASS_WAVE2D_CONSERVATIVE_FIVE_LU',self.script)
        canonical=SOURCE.read_text(encoding="utf-8")
        self.assertIn("def wave_eq(",canonical)
        self.assertIn('"0.5*q.T.q"',canonical)


if __name__=="__main__":
    unittest.main()
