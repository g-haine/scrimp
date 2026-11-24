# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Source-only contract for the exact first-CN SuperLU_DIST LU sweep."""

from pathlib import Path
import unittest

ROOT=Path(__file__).resolve().parents[1]
PROBE=ROOT/"sandbox/probe_v130_hw_superlu_permutation_sweep.py"
RUNNER=ROOT/"scripts/run-v130-hw-superlu-permutation-sweep.sh"


class IdenticalCNMatrixSuperLUPermutationTests(unittest.TestCase):
    def test_options_are_documented_and_actual_petsc_uses_them(self):
        src=PROBE.read_text()
        for expected in (
            "mat_superlu_dist_rowperm", "mat_superlu_dist_colperm",
            "mat_superlu_dist_replacetinypivot",
            "ksp.setFromOptions()",
            "F.view(ascii_viewer)",
            "selected_rowperm_verified",
            "selected_colperm_verified",
            "relative_linear_residual", "relative_state_error",
            "matrix_digest(A)",
        ):
            self.assertIn(expected,src)
        self.assertNotIn("A.shift(",src)
        self.assertNotIn("Mat.Option.NEW_NONZERO_ALLOCATION_ERR",src)

    def test_all_cases_are_isolated_native_processes(self):
        s=RUNNER.read_text()
        for name in (
            "default:default:baseline",
            "default:MMD_AT_PLUS_A:default_mmd",
            "default:MMD_ATA:default_ata",
            "NOROWPERM:MMD_AT_PLUS_A:norow_mmd",
            "NOROWPERM:MMD_ATA:norow_ata",
            "NOROWPERM:METIS_AT_PLUS_A:norow_metis",
            "NOROWPERM:NATURAL:norow_natural",
            "LargeDiag_AWPM:MMD_AT_PLUS_A:awpm_mmd",
            "LargeDiag_AWPM:METIS_AT_PLUS_A:awpm_metis",
            "first_cn_J_001.petscbin",
            "check-v130-published-immutable.py",
        ):
            self.assertIn(name,s)

    def test_nonredistributable_mc64_never_invoked_and_optional_awpm_isolated(self):
        s=RUNNER.read_text(encoding="utf-8")
        self.assertNotIn('"LargeDiag_MC64:MMD_AT_PLUS_A:mc64_mmd"',s)
        self.assertNotIn('"LargeDiag_MC64:MMD_ATA:mc64_mmd_ata"',s)
        self.assertIn("NOT_ATTEMPTED_HSL_NONREDISTRIBUTABLE",s)
        self.assertIn("CombBLAS",s)
        self.assertIn('printf \'%s\\n\' "$rc" > "$out/$label.exitcode"',s)
        self.assertIn("DIAGNOSTIC_PARTIAL_UNEXPECTED_PROCESS_FAILURE",s)
        self.assertIn("if unexpected:",s)
        self.assertIn('"heat_wave_CN_trajectory_qualified": False',s)


if __name__=="__main__":
    unittest.main()
