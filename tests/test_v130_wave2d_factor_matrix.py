# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# See LICENSE for license information.

"""Synthetic regression tests of the real sparse-Jacobian analyzer.

A saddle matrix can have TEN zero diagonal entries and NO zero rows
and yet be exactly nonsingular. These tests prevent us from mistaking
a diagonal-zero count for proof that a port-Hamiltonian DAE is singular.
Real native matrix capture and native backend acceptance remain separate.
"""

import ast
import importlib.util
from pathlib import Path
import unittest

import numpy as np
from scipy.sparse import csr_matrix, eye


ROOT = Path(__file__).resolve().parents[1]
ANALYZER = ROOT/"sandbox"/"analyze_v130_wave2d_factor_matrix.py"
PROBE = ROOT/"sandbox"/"probe_v130_wave2d_factor_matrix.py"
RUNNER = ROOT/"scripts"/"run-v130-wave2d-superlu-dist-matrix-audit.sh"
SPEC = importlib.util.spec_from_file_location(
    "v130_actual_wave2d_scipy_sparse_matrix_analyzer", ANALYZER
)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def make_matrix_pair(*, singular=False):
    matrix = eye(MODULE.N_EXPECTED, format="lil")
    for k in range(80, 90, 2):
        matrix[k, k] = 0.
        matrix[k+1, k+1] = 0.
        matrix[k, k+1] = 1.
        matrix[k+1, k] = 1.
    if singular:
        matrix[80, :] = 0.
    P = matrix.tocsr()
    K = P.copy()
    M = csr_matrix(P.shape, dtype=float)
    meta = {
        "schema": "scrimp-v130-first-main-CN-exact-PC-factor-matrix-v1",
        "time": .008,
        "TS_shift": 250.,
        "selected_solver": "superlu_dist",
        "matrix_shapes": {k:list(P.shape) for k in ("P","K","M")},
        "csr_SHA256": "synthetic-test-only-not-a-native-proof",
        "getfem_variable_intervals": [
            {"name":"test_algebraic_constraint", "lo":80, "hi":90,
             "dofs":10},
        ],
    }
    return meta, {"P":P,"K":K,"M":M}


class NativeJacobianAuditSourceTests(unittest.TestCase):
    def test_ten_zero_diagonals_can_be_structurally_full_rank(self):
        meta,mats=make_matrix_pair()
        report=MODULE.evaluate_snapshot(meta,mats)
        self.assertEqual(report["numerical_zero_diagonal_ids"],
                         list(range(80,90)))
        self.assertEqual(report["zero_numerical_rows"], [])
        self.assertEqual(report["structural_rank"], MODULE.N_EXPECTED)
        self.assertTrue(report["structural_full_rank"])
        self.assertTrue(report["P_is_original_K_plus_shift_M"])
        self.assertEqual(report["max_abs_P_minus_K_plus_sigma_M"],0.)
        self.assertEqual(report["zero_diagonal_variable_owners"]["80"],
                         ["test_algebraic_constraint"])
        for ordering in ("COLAMD", "MMD_AT_PLUS_A"):
            self.assertTrue(
                report["scipy_unregularized_lu"][ordering]["factorized"]
            )
            self.assertLess(
                report["scipy_unregularized_lu"][ordering][
                    "relative_residual"], 1e-12
            )

    def test_true_structural_rank_deficiency_is_detected(self):
        meta,mats=make_matrix_pair(singular=True)
        report=MODULE.evaluate_snapshot(meta,mats)
        self.assertFalse(report["structural_full_rank"])
        self.assertEqual(report["structural_rank"], MODULE.N_EXPECTED-1)
        self.assertEqual(report["zero_numerical_rows"],[80])
        self.assertFalse(report["scipy_unregularized_lu"]["COLAMD"][
            "factorized"
        ])

    def test_mumps_and_superlu_dist_operator_mismatch_is_detected(self):
        _meta,left=make_matrix_pair()
        _meta,right=make_matrix_pair()
        perturbed=right["P"].tolil()
        perturbed[400,400]=1.1
        right["P"]=perturbed.tocsr()
        checks=MODULE.compare_snapshots(left,right)
        self.assertFalse(checks["P"]["same_to_numerical_precision"])
        self.assertTrue(checks["K"]["same_to_numerical_precision"])
        self.assertTrue(checks["M"]["same_to_numerical_precision"])

    def test_capture_probes_original_pc_operator_without_model_edits(self):
        source=PROBE.read_text(encoding="utf-8")
        ast.parse(source)
        self.assertIn("original = S.DPHS.IJacobian",source)
        self.assertIn('"P": P',source)
        self.assertIn('"M": system.tangent_mass',source)
        self.assertIn('"K": system.tangent_stiffness',source)
        self.assertIn('system.gf_model.interval_of_variable(variable)',source)
        self.assertIn('"_profile_phase", ""',source)
        self.assertIn('abs(float(shift)-250.)',source)
        self.assertIn('S.DPHS.IJacobian = original',source)
        self.assertNotIn("replacetinypivot",source)
        self.assertNotIn("S.Brick(",source)
        self.assertNotIn("set_initial_value(",source)

    def test_runner_is_serial_same_model_stores_failure_and_native_csr(self):
        script=RUNNER.read_text(encoding="utf-8")
        self.assertIn("for solver in mumps superlu_dist; do",script)
        self.assertIn('mpiexec --oversubscribe --bind-to none -n 1',script)
        self.assertIn('SCRIMP_V130_WAVE2D_CAPTURE_DIR=/artifacts',script)
        self.assertIn('--final-time 0.008 --time-step 0.008',script)
        self.assertIn('first_main_CN_factor.npz',script)
        self.assertIn('sparse_matrix_analysis.json',script)
        self.assertIn('SUPERLU_DIST',script)
        self.assertIn('three_original_native_LU_passes_preserved',script)
        self.assertIn('"full_five_backend_gate_still": "FAIL"',script)
        self.assertNotIn('--solver superlu ',script)
        self.assertNotIn('ReplaceTinyPivot',script)


if __name__ == "__main__":
    unittest.main()
