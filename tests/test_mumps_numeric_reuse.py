# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp
"""Regression coverage for v1.3.0 MUMPS numeric + symbolic on_change.

Actual LU reuse counts and state parity are verified separately by the
full CN wave benchmark under the native MPI/MUMPS container.
"""
import os
import unittest

import numpy as np

from scrimp.dphs import (
    _get_mumps_reset_policy,
    _local_stored_csr_snapshot,
    _local_csr_change_flags,
)


class MockCSR:
    def __init__(self, *, rows=(0, 2), indptr=(0, 2, 3),
                 indices=(0, 1, 1), values=(2., -1., 4.)):
        self.rows = rows
        self.indptr = np.asarray(indptr, dtype=np.int32)
        self.indices = np.asarray(indices, dtype=np.int32)
        self.values = np.asarray(values)

    def getOwnershipRange(self):
        return self.rows

    def getValuesCSR(self):
        return self.indptr, self.indices, self.values


class TestExactFactorReuseDecision(unittest.TestCase):
    def test_default_is_on_change_even_for_linear_bricks(self):
        names = ("SCRIMP_MUMPS_RESET_POLICY",
                 "SCRIMP_MUMPS_RESET_PC_ON_JACOBIAN")
        originals = {name: os.environ.pop(name, None) for name in names}
        try:
            self.assertEqual(
                _get_mumps_reset_policy(all_bricks_linear=True), "on_change"
            )
            self.assertEqual(
                _get_mumps_reset_policy(all_bricks_linear=False), "on_change"
            )
            os.environ["SCRIMP_MUMPS_RESET_POLICY"] = "never"
            self.assertEqual(_get_mumps_reset_policy(), "never")
            os.environ["SCRIMP_MUMPS_RESET_POLICY"] = "always"
            self.assertEqual(_get_mumps_reset_policy(), "always")
        finally:
            for name in names:
                os.environ.pop(name, None)
                if originals[name] is not None:
                    os.environ[name] = originals[name]

    def test_exact_numeric_coefficient_change(self):
        mat = MockCSR()
        first = _local_stored_csr_snapshot(mat)
        self.assertEqual(_local_csr_change_flags(None, first), 3)
        self.assertEqual(
            _local_csr_change_flags(first, _local_stored_csr_snapshot(mat)), 0
        )
        mat.values[0] += 1.e-15
        self.assertEqual(
            _local_csr_change_flags(first, _local_stored_csr_snapshot(mat)), 2
        )
        self.assertEqual(first[3][0], 2.0)

    def test_structural_change_requires_both_factors(self):
        mat = MockCSR()
        before = _local_stored_csr_snapshot(mat)
        mat.indices[1] = 2
        self.assertEqual(
            _local_csr_change_flags(before, _local_stored_csr_snapshot(mat)), 3
        )

    def test_ownership_change_requires_reanalysis(self):
        mat = MockCSR()
        before = _local_stored_csr_snapshot(mat)
        mat.rows = (2, 4)
        self.assertEqual(
            _local_csr_change_flags(before, _local_stored_csr_snapshot(mat)), 3
        )

    def test_nonfinite_never_considered_unchanged(self):
        mat = MockCSR(values=(2., float("nan"), 4.))
        before = _local_stored_csr_snapshot(mat)
        self.assertEqual(
            _local_csr_change_flags(before, _local_stored_csr_snapshot(mat)), 2
        )

    def test_complex_entries_preserved(self):
        mat = MockCSR(values=(2+1j, -1j, 4+0j))
        before = _local_stored_csr_snapshot(mat)
        self.assertEqual(
            _local_csr_change_flags(before, _local_stored_csr_snapshot(mat)), 0
        )
        mat.values[2] += 1j
        self.assertEqual(
            _local_csr_change_flags(before, _local_stored_csr_snapshot(mat)), 2
        )


if __name__ == "__main__":
    unittest.main()
