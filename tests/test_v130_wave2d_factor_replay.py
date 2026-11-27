# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp


"""Source-only checks for the exact 2D first-CN PETSc factor replay.

They must never infer a successful five-LU wave simulation from a
direct standalone factorization. Native PETSc is NOT imported here.
"""

import ast
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]
REPLAY = ROOT/"sandbox"/"replay_v130_wave2d_superlu_dist_factor.py"
RUNNER = ROOT/"scripts"/"run-v130-wave2d-superlu-dist-factor-replay.sh"
SOURCE = ROOT/"examples"/"wave.py"
HASH = "24346856f358e4f81cf332ebc974e329dc084c60f5d6b98bfbcd788fdff37f54"


class FirstCNFactorReplaySourceTests(unittest.TestCase):
    def test_unchanged_real_native_input_operator_sha256(self):
        src = REPLAY.read_text(encoding="utf-8")
        ast.parse(src)
        self.assertIn(HASH, src)
        self.assertIn('matrix.getOwnershipRange()', src)
        self.assertIn("PETSc.Mat.Option.IGNORE_ZERO_ENTRIES, False", src)
        self.assertIn("original_full_P_M_K_SHA256", src)
        self.assertIn(
            "PETSc_owned_CSR_graph_and_values_identical_to_archive", src
        )
        self.assertIn('"P_indptr"', src)
        self.assertIn('"P_indices"', src)
        self.assertIn('"P_data"', src)
        self.assertIn("matrix.getValuesCSR()", src)
        self.assertIn("graph_mismatches", src)
        self.assertIn("np.array_equal(orig_cols, local_cols)", src)
        self.assertIn('global_max != 0.', src)
        self.assertNotIn("S.Brick(", src)
        self.assertNotIn("set_initial_value(", src)

    def test_true_native_PC_without_fallback(self):
        src = REPLAY.read_text(encoding="utf-8")
        self.assertIn('ksp.setType("preonly")', src)
        self.assertIn('pc.setType("lu")', src)
        self.assertIn('pc.setFactorSolverType("superlu_dist")', src)
        self.assertIn('"mat_superlu_dist_replacetinypivot", "true"', src)
        self.assertIn("ksp.setFromOptions()", src)
        self.assertIn("ksp.solve(rhs, solution)", src)
        self.assertNotIn("mumps", src.lower())
        self.assertNotIn('pc_factor_shift', src)
        self.assertNotIn('mat.setDiagonal(', src)
        self.assertNotIn('A.diagonal(', src)

    def test_original_Ax_true_and_real_backward_error(self):
        src = REPLAY.read_text(encoding="utf-8")
        for needle in (
            "np.sin(.013 * offset) + .25",
            "mat.mult(known, rhs)",
            "mat.mult(solution, check)",
            "check.axpy(-1., rhs)",
            "solution.copy()",
            "diff.axpy(-1., known)",
            "relative_backward_residual",
            "relative_known_solution_error",
            "rel_backward < 1e-11",
            "rel_solution < 1e-7",
            "comm.allreduce(bool(candidate_ok), op=MPI.LAND)",
        ):
            with self.subTest(needle=needle):
                self.assertIn(needle,src)

    def test_replay_runs_both_MPI_and_both_options(self):
        script = RUNNER.read_text(encoding="utf-8")
        for needle in (
            "for mode in default replace_tiny_pivot; do",
            "for ranks in 1 2; do",
            '--mode "$mode"',
            '--solver',  # Negative check below, no full wave CLI.
        ):
            if needle == '--solver':
                self.assertNotIn(needle, script)
            else:
                self.assertIn(needle, script)
        self.assertIn(HASH, script)
        self.assertIn("v130-wave2d-first-cn-jacobian-574d9c71-NKqIZ0BZ.zip",
                      script)
        self.assertIn("rank_json_present", script)
        self.assertIn('"wave_installation_gate":"3/5"', script)
        self.assertNotIn("SCRIMP_NONLINEAR_MASS_JACOBIAN", script)
        self.assertNotIn("mat_superlu_dist_rowperm", script)

    def test_original_full_wave_gates_are_explicitly_not_accepted(self):
        script=RUNNER.read_text(encoding="utf-8")
        src=REPLAY.read_text(encoding="utf-8")
        self.assertIn(
            "PASS_OPTIONAL_PETSC_LU_REPLAY_MPI1_MPI2_NOT_WAVE_ACCEPTANCE",
            script,
        )
        self.assertIn(
            '"full_wave_superlu_dist_solver_qualified":False', script
        )
        self.assertIn(
            '"original_J_not_modified": True', src
        )
        self.assertIn(
            '"is_not_original_wave_CN_physics_qualification": True',
            src,
        )
        self.assertNotIn('PASS_WAVE2D_CONSERVATIVE_FIVE_LU', script)

    def test_can_keep_model_and_dissipative_pedagogical_sequence_untouched(self):
        src=SOURCE.read_text(encoding="utf-8")
        self.assertIn("def wave_eq(", src)
        self.assertIn('S.Term(', src)
        self.assertIn('"0.5*q.T.q"', src)
        self.assertIn('init_step_ts_type="algebraic"', src)


if __name__ == "__main__":
    unittest.main()
