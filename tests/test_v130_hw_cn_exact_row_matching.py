# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Source-only regression of exact original CN row-permutation experiment."""

import ast
from pathlib import Path
import unittest

ROOT=Path(__file__).resolve().parents[1]
PILOT=ROOT/"sandbox/probe_v130_hw_cn_exact_row_matching.py"
RUNNER=ROOT/"scripts/run-v130-hw-cn-exact-row-matching.sh"


class OriginalCNRowMatchingContracts(unittest.TestCase):
    def test_real_bipartite_matching_excludes_structural_zeros(self):
        src=PILOT.read_text(encoding="utf-8")
        for key in (
            'maximum_bipartite_matching(graph, perm_type="row")',
            'graph.eliminate_zeros()',
            "A[rows, :].tocsr()",
            "np.sort(A.data)",
            "np.sort(B.data)",
            "matched_smallest_absolute_diagonal",
            "rows_moved",
            "bipartite_matching_complete",
        ):
            self.assertIn(key,src)
    def test_equivalent_petsc_coo_and_original_residual(self):
        src=PILOT.read_text(encoding="utf-8")
        for key in (
            "B.setPreallocationCOO(i, j)",
            "B.setValuesCOO(v, addv=PETSc.InsertMode.INSERT_VALUES)",
            "old_rhs[reordered]",
            "candidate_rhs.axpy(-1.,rhs)",
            "A.mult(result,check)",
            "check.axpy(-1.,b0)",
            "relative_residual_in_ORIGINAL_equations",
            "row_matched_CSR_SHA256",
            "matrix_digest(A)",
            "pc.setFactorSolverType(args.solver)",
        ):
            self.assertIn(key,src)
        self.assertNotIn("A.shift(",src)
    def test_pinned_native_runner_no_disallowed_build_features(self):
        src=RUNNER.read_text(encoding="utf-8")
        for key in (
            "row_matching:superlu_dist:NATURAL:matching_natural",
            "row_matching:superlu_dist:MMD_AT_PLUS_A:matching_mmd",
            "row_matching:superlu_dist:METIS_AT_PLUS_A:matching_metis",
            "row_matching:mumps:default:matching_mumps",
            "original:superlu_dist:NATURAL:original_natural",
            "SCRIMP_V130_DAE_IMAGE",
            "check-v130-published-immutable.py",
            "ORIGINAL CN STRUCTURAL MATCHING LU EVIDENCE",
        ):
            self.assertIn(key,src)
        self.assertNotIn("LargeDiag_MC64",src)
        self.assertNotIn("LargeDiag_AWPM",src)


if __name__=="__main__":
    unittest.main()
