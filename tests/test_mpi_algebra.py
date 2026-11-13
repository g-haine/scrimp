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

The synthetic matrices exercise BOTH genuinely distributed contributions
and fully replicated GetFEM matrices with GLOBAL row/column numbering.
The real-GetFEM test checks a physical mass integral with MPI level 2.
"""

import os
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
        # Each process contributes some different rows, except row 0
        # where all processes contribute 1/P. The global matrix is
        # therefore the original stencil and includes overlapping entries.
        self.row_ptr = np.arange(0, 2 * count + 1, 2, dtype=PETSc.IntType)
        self.col_idx = np.empty(2 * count, dtype=PETSc.IntType)
        self.values = np.empty(2 * count, dtype=PETSc.ScalarType)
        for k in range(count):
            self.col_idx[2*k:2*k+2] = (k, (k + 1) % count)
            factor = (1.0 / ranks if k == 0 else
                      (1.0 if k % ranks == RANK else 0.0))
            self.values[2*k:2*k+2] = (2.0 * factor, -1.0 * factor)

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
            with patch.dict(os.environ, {"SCRIMP_GETFEM_MATRIX_MODE": "distributed"}), \
                    patch.object(linalg.gf, "Spmat", return_value=fake):
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

    def test_replicated_matrix_does_not_scale_with_rank_count(self):
        """Replicated GetFEM Python matrices must be inserted just once."""
        mat = PETSc.Mat().createAIJ(size=(N, N), comm=COMM)
        mat.setUp()
        fake = FakeGetFEMDistributedMatrix()
        # Replace the distributed contributions with an exact full copy
        # of the same matrix on every rank.
        for k in range(N):
            fake.values[2*k:2*k+2] = (2.0, -1.0)
        try:
            with patch.dict(os.environ, {"SCRIMP_GETFEM_MATRIX_MODE": "replicated"}), \
                    patch.object(linalg.gf, "Spmat", return_value=fake):
                for _ in range(2):
                    linalg.convert_gmm_to_petsc(object(), mat)
                    x, y = mat.createVecRight(), mat.createVecLeft()
                    try:
                        x.set(1.0)
                        mat.mult(x, y)
                        lo, hi = y.getOwnershipRange()
                        np.testing.assert_allclose(
                            y.getArray(readonly=True),
                            np.ones(hi-lo),
                            atol=1e-11, rtol=1e-11,
                        )
                    finally:
                        x.destroy()
                        y.destroy()
        finally:
            mat.destroy()

    def test_real_getfem_mpi2_mass_matrix_is_rank_independent(self):
        """A real GetFEM model must integrate the constant 1 on [0,1].

        This is much stronger than the synthetic partial-matrix test:
        it detects double-counted GetFEM contributions, including those
        represented with zero-valued sparse slots on some MPI ranks.
        The identity integral is independent of the number of processes.
        """
        gf = linalg.gf
        if int(gf.util_mpi_parallelism_level()) != 2:
            self.skipTest("requires GetFEM MPI level 2")
        mesh = gf.Mesh("cartesian", np.linspace(0.0, 1.0, 17))
        mf = gf.MeshFem(mesh, 1)
        mf.set_fem(gf.Fem("FEM_PK(1,1)"))
        mim = gf.MeshIm(mesh, gf.Integ("IM_GAUSS1D(5)"))
        model = gf.Model("real")
        model.add_fem_variable("u", mf)
        model.add_linear_term(mim, "u*Test_u")
        model.assembly(option="build_matrix")

        n = model.nbdof()
        self.assertEqual(n, mf.nbdof())
        mat = PETSc.Mat().createAIJ(size=(n, n), comm=COMM)
        mat.setUp()
        mat.setOption(PETSc.Mat.Option.NEW_NONZERO_ALLOCATION_ERR, False)
        try:
            with patch.dict(os.environ, {"SCRIMP_GETFEM_MATRIX_MODE": "replicated"}):
                linalg.convert_gmm_to_petsc(model.tangent_matrix(), mat)
            x, y = mat.createVecRight(), mat.createVecLeft()
            try:
                x.set(1.0)
                mat.mult(x, y)
                # Exact integral of 1*1 over the unit interval = 1.
                measured = float(np.real(x.dot(y)))
                self.assertAlmostEqual(measured, 1.0, places=9)
                self.assertEqual(mat.getSize(), (n, n))
            finally:
                x.destroy()
                y.destroy()
        finally:
            mat.destroy()

    def test_real_getfem_nonlinear_tangent_across_newton_iterates(self):
        """A changing GetFEM nonlinear Jacobian must not change MPI mode.

        For R(u)=integral(u**2 * Test_u), the tangent at a constant u=c
        satisfies 1^T J 1 = 2*c on [0,1]. Verify this for two separate
        GetFEM assemblies into the SAME PETSc matrix.
        """
        gf = linalg.gf
        mesh = gf.Mesh("cartesian", np.linspace(0.0, 1.0, 17))
        mf = gf.MeshFem(mesh, 1)
        mf.set_fem(gf.Fem("FEM_PK(1,1)"))
        mim = gf.MeshIm(mesh, gf.Integ("IM_GAUSS1D(5)"))
        model = gf.Model("real")
        model.add_fem_variable("u", mf)
        model.add_nonlinear_term(mim, "u*u*Test_u")
        n = model.nbdof()
        mat = PETSc.Mat().createAIJ(size=(n, n), comm=COMM)
        mat.setUp()
        mat.setOption(PETSc.Mat.Option.NEW_NONZERO_ALLOCATION_ERR, False)
        try:
            with patch.dict(os.environ, {"SCRIMP_GETFEM_MATRIX_MODE": "replicated"}):
                for c in (1.0, 2.0):
                    model.set_variable("u", np.full(n, c, dtype=float))
                    model.assembly(option="build_matrix")
                    linalg.convert_gmm_to_petsc(model.tangent_matrix(), mat)
                    x, y = mat.createVecRight(), mat.createVecLeft()
                    try:
                        x.set(1.0)
                        mat.mult(x, y)
                        self.assertAlmostEqual(
                            float(np.real(x.dot(y))), 2.0*c, places=9
                        )
                    finally:
                        x.destroy()
                        y.destroy()
        finally:
            mat.destroy()

    def test_nonlinear_tangents_reuse_mpi_matrices_without_accumulation(self):
        """Newton reassembly must not replace PETSc MPIAIJ matrix handles.

        This regression also checks that variable sparsity may be inserted
        into the same preallocated matrix and old values do not accumulate.
        The original SWE crash was observed when temporary matrices were
        reclaimed by PETSc.garbage_cleanup in TSMonitor.
        """
        def linear_diagonal(value):
            mat = PETSc.Mat().createAIJ(size=(N, N), nnz=1, comm=COMM)
            mat.setUp()
            diag = mat.createVecLeft()
            try:
                diag.set(value)
                mat.setDiagonal(diag)
                mat.assemble()
            finally:
                diag.destroy()
            return mat

        model = SimpleNamespace(
            nbdof=lambda: N,
            assembly=lambda **kwargs: None,
            tangent_matrix=lambda: None,
        )
        system = object.__new__(DPHS)
        system.gf_model = model
        system.bricks = {}
        system.mass = linear_diagonal(3.0)
        system.stiffness = linear_diagonal(7.0)
        system.nl_mass = PETSc.Mat().create(comm=COMM)
        system.nl_stiffness = PETSc.Mat().create(comm=COMM)
        system.tangent_mass = PETSc.Mat().create(comm=COMM)
        system.tangent_stiffness = PETSc.Mat().create(comm=COMM)
        matrices = (
            system.mass, system.stiffness, system.nl_mass,
            system.nl_stiffness, system.tangent_mass,
            system.tangent_stiffness,
        )
        try:
            for mat in (system.nl_mass, system.nl_stiffness,
                        system.tangent_mass, system.tangent_stiffness):
                system.init_matrix(mat)
            original_mass = system.tangent_mass
            original_stiffness = system.tangent_stiffness
            fake = FakeGetFEMDistributedMatrix()
            with patch.dict(os.environ,
                            {"SCRIMP_GETFEM_MATRIX_MODE": "replicated"}), \
                    patch.object(linalg.gf, "Spmat", return_value=fake):
                for coefficient in (1.0, 2.0, 0.5, 3.0, 1.0):
                    for k in range(N):
                        fake.values[2*k:2*k+2] = (
                            2.0 * coefficient, -coefficient
                        )
                    system.assemble_nl_mass()
                    system.assemble_nl_stiffness()
                    self.assertIs(system.tangent_mass, original_mass)
                    self.assertIs(system.tangent_stiffness, original_stiffness)
                    for mat, linear_value in (
                        (system.tangent_mass, 3.0),
                        (system.tangent_stiffness, 7.0),
                    ):
                        x, y = mat.createVecRight(), mat.createVecLeft()
                        try:
                            x.set(1.0)
                            mat.mult(x, y)
                            lo, hi = y.getOwnershipRange()
                            np.testing.assert_allclose(
                                y.getArray(readonly=True),
                                np.full(hi-lo, linear_value + coefficient),
                                rtol=1e-12, atol=1e-12,
                            )
                        finally:
                            x.destroy()
                            y.destroy()
        finally:
            for mat in matrices:
                mat.destroy()

    def test_init_matrix_preserves_distributed_vector_storage(self):
        """Preallocating the diagonal must never replace a Vec's buffer.

        Regression for DPHS.init_matrix using Vec.setArray(0), which is
        an invalid way to fill the distributed PETSc vector with zero.
        Repeatedly initialize, assemble, then use a PETSc-owned vector.
        """
        system = object.__new__(DPHS)
        system.gf_model = SimpleNamespace(nbdof=lambda: N)
        for _ in range(4):
            matrix = PETSc.Mat().create(comm=COMM)
            try:
                system.init_matrix(matrix)
                self.assertEqual(matrix.getSize(), (N, N))
                diag = matrix.getDiagonal()
                try:
                    lo, hi = diag.getOwnershipRange()
                    self.assertEqual(len(diag.getArray(readonly=True)), hi-lo)
                    np.testing.assert_allclose(diag.getArray(readonly=True), 0.0)
                finally:
                    diag.destroy()
                x, y = matrix.createVecRight(), matrix.createVecLeft()
                try:
                    x.set(1.0)
                    matrix.mult(x, y)
                    np.testing.assert_allclose(y.getArray(readonly=True), 0.0)
                finally:
                    x.destroy()
                    y.destroy()
            finally:
                matrix.destroy()

    def test_distinct_operator_and_preconditioner_get_same_jacobian(self):
        """TS's distinct A/P matrices must both receive J=K+shift*M."""
        def diagonal_matrix(value):
            mat = PETSc.Mat().createAIJ(size=(N, N), nnz=1, comm=COMM)
            mat.setUp()
            diag = mat.createVecLeft()
            diag.set(value)
            mat.setDiagonal(diag)
            mat.assemble()
            diag.destroy()
            return mat

        stiffness = diagonal_matrix(2.0)
        mass = diagonal_matrix(3.0)
        A = diagonal_matrix(0.0)
        P = diagonal_matrix(0.0)
        system = object.__new__(DPHS)
        system.gf_model = SimpleNamespace(set_time=lambda t: None)
        system.lumped_variables = set()
        system.linear_mass = True
        system.linear_stiffness = True
        system.tangent_stiffness = stiffness
        system.tangent_mass = mass
        try:
            system.IJacobian(None, 0.0, None, None, 5.0, A, P)
            for mat in (A, P):
                x, y = mat.createVecRight(), mat.createVecLeft()
                try:
                    x.set(1.0)
                    mat.mult(x, y)
                    lo, hi = y.getOwnershipRange()
                    np.testing.assert_allclose(
                        y.getArray(readonly=True),
                        np.full(hi - lo, 17.0),
                        atol=1e-12,
                    )
                finally:
                    x.destroy()
                    y.destroy()
        finally:
            for mat in (stiffness, mass, A, P):
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
