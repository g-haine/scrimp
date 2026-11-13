# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""MPI integration tests for the GetFEM/PETSc algebra bridge.

Run with real GetFEM/PETSc (not the lightweight tests/__init__.py stubs):
    mpiexec -n 1 python -m unittest -v tests.test_mpi_algebra
    mpiexec -n 2 python -m unittest -v tests.test_mpi_algebra
    mpiexec -n 4 python -m unittest -v tests.test_mpi_algebra

The synthetic GetFEM sparse matrix carries PARTIAL contributions on each
rank with GLOBAL row/column numbering, matching GetFEM MPI level 2.
This tests the bridge contract, not the physical finite-element discretization.
"""

import unittest
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
from petsc4py import PETSc

from scrimp.dphs import DPHS
from scrimp.utils import linalg


COMM = PETSc.COMM_WORLD
RANK = COMM.getRank()
SIZE = COMM.getSize()
N = max(8, SIZE * 2)


class FakeGetFEMDistributedMatrix:
    """GetFEM CSC-of-transpose view of a partial GLOBAL-NxN matrix."""

    def __init__(self, count=N, ranks=SIZE):
        # Every process contributes 1/P of the same stencil, including
        # off-rank rows. Summing partial data must recover a single operator.
        self.row_ptr = np.arange(0, 2 * count + 1, 2, dtype=PETSc.IntType)
        self.col_idx = np.empty(2 * count, dtype=PETSc.IntType)
        self.values = np.empty(2 * count, dtype=PETSc.ScalarType)
        for k in range(count):
            self.col_idx[2*k:2*k+2] = (k, (k + 1) % count)
            self.values[2*k:2*k+2] = (2.0 / ranks, -1.0 / ranks)

    def transpose(self):
        # linalg.convert_gmm_to_petsc always transposes a GetFEM copy
        # to obtain this CSR view.
        return None

    def csc_ind(self):
        return self.row_ptr, self.col_idx

    def csc_val(self):
        return self.values


class TestMPIAlgebra(unittest.TestCase):
    def test_distributed_matrix_contributions_are_added_and_not_accumulated(self):
        mat = PETSc.Mat().createAIJ(size=(N, N), comm=COMM)
        mat.setUp()
        fake = FakeGetFEMDistributedMatrix()
        try:
            with patch.object(linalg.gf, "Spmat", return_value=fake):
                for _ in range(2):
                    # Reassembly must REPLACE old entries, not double the
                    # operator across Newton iterations.
                    linalg.convert_gmm_to_petsc(object(), mat)

                    x = mat.createVecRight()
                    y = mat.createVecLeft()
                    x.set(1.0)
                    mat.mult(x, y)
                    lo, hi = y.getOwnershipRange()
                    np.testing.assert_allclose(
                        y.getArray(readonly=True),
                        np.ones(hi - lo),
                        atol=1e-11,
                        rtol=1e-11,
                    )
                    self.assertEqual(mat.getSize(), (N, N))
                    x.destroy()
                    y.destroy()
        finally:
            mat.destroy()

    def test_rhs_inserted_on_owned_rows_only(self):
        rhs = PETSc.Vec().createMPI(N, comm=COMM)
        model = SimpleNamespace(
            nbdof=lambda: N,
            assembly=lambda **kwargs: None,
            rhs=lambda: np.arange(N, dtype=float) + 1.0,
        )
        system = object.__new__(DPHS)
        system.bricks = {}
        system.gf_model = model
        system.rhs = rhs
        try:
            system.assemble_rhs()
            lo, hi = rhs.getOwnershipRange()
            np.testing.assert_allclose(
                rhs.getArray(readonly=True),
                np.arange(lo, hi, dtype=float) + 1.0,
            )
            self.assertEqual(rhs.getSize(), N)
        finally:
            rhs.destroy()

    def test_algebraic_tolerances_preserve_distributed_layout(self):
        residual = PETSc.Vec().createMPI(N, comm=COMM)
        system = object.__new__(DPHS)
        system.F = residual
        system.gf_model = SimpleNamespace(
            nbdof=lambda: N,
            interval_of_variable=lambda name: (2, 1),
        )
        system.ports = {
            "a": SimpleNamespace(
                get_algebraic=lambda: True,
                get_substituted=lambda: True,
                get_flow=lambda: "a",
            ),
        }
        captured = {}
        ts = SimpleNamespace(
            getTolerances=lambda: (1e-5, 1e-7),
            setTolerances=lambda **kwargs: captured.update(kwargs),
        )
        try:
            system.exclude_algebraic_var_from_lte(ts)
            atol = captured["atol"]
            lo, hi = atol.getOwnershipRange()
            expected = np.full(hi - lo, 1e-7, dtype=float)
            if lo <= 2 < hi:
                expected[2 - lo] = np.inf
            np.testing.assert_allclose(
                atol.getArray(readonly=True), expected
            )
            self.assertEqual(atol.getSize(), residual.getSize())
            atol.destroy()
        finally:
            residual.destroy()


if __name__ == "__main__":
    unittest.main()
