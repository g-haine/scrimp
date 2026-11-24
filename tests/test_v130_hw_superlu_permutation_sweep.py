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
            "LargeDiag_MC64:MMD_AT_PLUS_A:mc64_mmd",
            "LargeDiag_AWPM:MMD_AT_PLUS_A:awpm_mmd",
            "LargeDiag_AWPM:METIS_AT_PLUS_A:awpm_metis",
            "LargeDiag_MC64:MMD_ATA:mc64_mmd_ata",
            "NOROWPERM:MMD_AT_PLUS_A:norow_mmd",
            "first_cn_J_001.petscbin",
            "check-v130-published-immutable.py",
        ):
            self.assertIn(name,s)


if __name__=="__main__":
    unittest.main()
