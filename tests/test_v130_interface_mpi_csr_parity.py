# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Dependency-light regression of canonical MPI sparse K parity."""

import importlib.util
from pathlib import Path
import unittest

import numpy as np
from scipy.sparse import csr_matrix, vstack

ROOT = Path(__file__).resolve().parents[1]
FILE = ROOT / "sandbox" / "v130_interface_csr_parity.py"
spec = importlib.util.spec_from_file_location("v130_interface_csr_parity", FILE)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
compare_assembled_csr = module.compare_assembled_csr


class SparseMPIParityTests(unittest.TestCase):
    def setUp(self):
        self.reference = csr_matrix(np.array([
            [2., 0., 0., 0.],
            [0., 3., 0., 0.],
            [0., 0., 4., 0.],
            [0., 0., 0., 5.],
        ]))
        # Additional STORED zeros on both MPI2 ranks.
        first = csr_matrix((
            np.array([2., 0., 3., 0.]),
            np.array([0, 2, 1, 3]),
            np.array([0, 2, 4]),
        ), shape=(2, 4))
        second = csr_matrix((
            np.array([0., 4., 0., 5.]),
            np.array([0, 2, 1, 3]),
            np.array([0, 2, 4]),
        ), shape=(2, 4))
        self.parallel = vstack([first, second], format="csr")

    def test_explicit_stored_zeros_are_not_numerical_differences(self):
        result = compare_assembled_csr(self.reference, self.parallel, 4)
        self.assertTrue(result["passed"])
        self.assertEqual(result["MPI1_stored_nnz"], 4)
        self.assertEqual(result["MPI2_stored_nnz"], 8)
        self.assertEqual(result["MPI1_explicit_zeros"], 0)
        self.assertEqual(result["MPI2_explicit_zeros"], 4)
        self.assertEqual(result["MPI1_numerical_nnz"], 4)
        self.assertEqual(result["MPI2_numerical_nnz"], 4)
        self.assertFalse(result["raw_stored_nnz_identical"])
        self.assertTrue(result["canonical_CSR_bitwise_identical"])
        self.assertEqual(result["difference_nnz"], 0)
        self.assertEqual(result["max_abs_difference"], 0.)

    def test_original_matrices_not_modified(self):
        compare_assembled_csr(self.reference, self.parallel, 4)
        self.assertEqual(self.parallel.nnz, 8)
        self.assertEqual(self.reference.nnz, 4)

    def test_real_numerical_difference_is_rejected(self):
        B = self.parallel.copy()
        B.data[0] += 0.01
        result = compare_assembled_csr(self.reference, B, 4)
        self.assertFalse(result["passed"])
        self.assertGreater(result["max_abs_difference"], result["threshold"])

    def test_changed_sparsity_rejected_even_for_tiny_entry(self):
        B = self.parallel.tolil()
        B[0, 1] = 1e-14
        result = compare_assembled_csr(self.reference, B.tocsr(), 4)
        self.assertFalse(result["passed"])
        self.assertFalse(result["canonical_sparsity_identical"])

    def test_incorrect_dimension_rejected(self):
        self.assertFalse(
            compare_assembled_csr(self.reference, self.parallel, 5)["passed"]
        )

    def test_nonfinite_matrix_rejected(self):
        B = self.parallel.copy()
        B.data[0] = np.nan
        with self.assertRaisesRegex(ValueError, "Nonfinite"):
            compare_assembled_csr(self.reference, B, 4)

    def test_generator_is_not_used_as_vstack_sequence(self):
        rows = [self.parallel[:2], self.parallel[2:]]
        full = vstack([part for part in rows], format="csr")
        self.assertTrue(compare_assembled_csr(
            self.reference, full, 4,
        )["passed"])


if __name__ == "__main__":
    unittest.main()
