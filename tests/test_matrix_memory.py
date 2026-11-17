# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""No-PETSc unit checks of the compact matrix-memory algebra.

PETSc MPI/AIJ/MatCopy behavior is separately validated in the native
short acoustic and Saint-Venant comparisons; these tests establish
the exact unchanged finite-difference formula and object reuse.
"""
import importlib.util
from pathlib import Path
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "scrimp" / "utils" / "matrix_memory.py"
SPEC = importlib.util.spec_from_file_location("_scrimp_matrix_memory", SOURCE)
MEMORY = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MEMORY)


class _Mat:
    creations = 0

    def __init__(self, values):
        self.values = np.asarray(values, dtype=np.float64).copy()
        self.destroyed = False
        self.assemble_count = 0
        _Mat.creations += 1

    def destroy(self):
        self.destroyed = True

    def copy(self, result=None, structure=None):
        if result is None:
            return _Mat(self.values)
        result.values = self.values.copy()
        return result

    def axpy(self, alpha, other, structure=None):
        self.values += alpha * other.values

    def scale(self, alpha):
        self.values *= alpha

    def assemble(self):
        self.assemble_count += 1


class _Vec:
    creations = 0

    def __init__(self, values):
        self.values = np.asarray(values, dtype=np.float64).copy()
        _Vec.creations += 1

    def duplicate(self):
        return _Vec(np.zeros_like(self.values))

    def copy(self, result=None):
        if result is None:
            return _Vec(self.values)
        result.values[:] = self.values
        return result

    def axpy(self, alpha, other):
        self.values += alpha * other.values


class _Owner:
    """M(z)=B0+z0*B1+z1*B2, a symmetric affine mass operator."""

    B0 = np.array([[4., 0.], [0., 3.]])
    B1 = np.array([[0.75, 0.1], [0.1, -0.25]])
    B2 = np.array([[0.2, 0.5], [0.5, 0.3]])

    def __init__(self, initial):
        self.current = initial.values.copy()
        self.tangent_mass = _Mat(self.B0)
        self.assembly_count = 0
        self._mass_hessian_fd_probe = None
        self._mass_hessian_fd_matrix = None

    def _sync_lumped_trial_state(self, t, state):
        self.current = state.values.copy()

    def assemble_nl_mass(self):
        self.assembly_count += 1
        self.tangent_mass.values[:] = (
            self.B0 + self.current[0] * self.B1
            + self.current[1] * self.B2
        )


class CompactMatrixMemoryTests(unittest.TestCase):
    def test_linear_tangent_alias_disposes_only_unused_placeholder(self):
        base = _Mat([[1., 0.], [0., 2.]])
        placeholder = _Mat([[0., 0.], [0., 0.]])
        result = MEMORY.alias_readonly_linear_tangent(base, placeholder)
        self.assertIs(result, base)
        self.assertFalse(base.destroyed)
        self.assertTrue(placeholder.destroyed)
        MEMORY.alias_readonly_linear_tangent(base, base)
        self.assertFalse(base.destroyed)

    def test_hessian_reuses_one_matrix_and_one_vec_at_every_callback(self):
        eta = 1e-4
        z = _Vec([2., -1.])
        zd = _Vec([0.4, -0.7])
        owner = _Owner(z)
        reference_perturbation = zd.values[0] * owner.B1 + zd.values[1] * owner.B2
        matrix_creations_before = _Mat.creations
        vec_creations_before = _Vec.creations
        last_workspace = None
        last_probe = None

        for repeat in range(5):
            new_state = _Vec([2. + 0.2 * repeat, -1. + 0.1 * repeat])
            base_jacobian = np.array([[4., 1.], [1., 5.]])
            jacobian = _Mat(base_jacobian)
            state_before = new_state.values.copy()
            MEMORY.add_mass_hessian_fd_reused(
                owner, 0.01 * repeat, new_state, zd, jacobian,
                eta, structure="different",
            )
            np.testing.assert_allclose(
                jacobian.values, base_jacobian + reference_perturbation,
                # Central differences of double-precision values
                # around O(1) matrices with eta=1e-4 incur O(1e-12)
                # roundoff. Bound it tightly without demanding exactness.
                rtol=0., atol=2e-11,
            )
            np.testing.assert_array_equal(new_state.values, state_before)
            np.testing.assert_array_equal(owner.current, state_before)
            np.testing.assert_allclose(
                owner.tangent_mass.values,
                owner.B0 + state_before[0] * owner.B1
                + state_before[1] * owner.B2,
                rtol=0., atol=1e-14,
            )
            self.assertEqual(jacobian.assemble_count, 1)
            self.assertIsNotNone(owner._mass_hessian_fd_matrix)
            self.assertIsNotNone(owner._mass_hessian_fd_probe)
            self.assertEqual(
                owner._mass_hessian_fd_matrix_allocations, 1
            )
            self.assertEqual(
                owner._mass_hessian_fd_probe_allocations, 1
            )
            self.assertEqual(
                owner._mass_hessian_fd_matrix_reuses, repeat
            )
            if last_workspace is not None:
                self.assertIs(owner._mass_hessian_fd_matrix, last_workspace)
                self.assertIs(owner._mass_hessian_fd_probe, last_probe)
            last_workspace = owner._mass_hessian_fd_matrix
            last_probe = owner._mass_hessian_fd_probe

        # Apart from test-local Jacobians, a single persistent matrix
        # was allocated for the Hessian across all five calls.
        self.assertEqual(
            _Mat.creations - matrix_creations_before, 6
        )
        self.assertEqual(
            _Vec.creations - vec_creations_before, 6
        )

    def test_hessian_restores_base_state_on_exception(self):
        z = _Vec([0.5, 0.25])
        zd = _Vec([0.1, 0.2])
        owner = _Owner(z)
        real_assemble = owner.assemble_nl_mass
        failures = {"count": 0}

        def faulty_assemble():
            failures["count"] += 1
            if failures["count"] == 2:
                raise ValueError("synthetic GetFEM assembly failure")
            return real_assemble()

        owner.assemble_nl_mass = faulty_assemble
        with self.assertRaisesRegex(ValueError, "synthetic"):
            MEMORY.add_mass_hessian_fd_reused(
                owner, 0., z, zd, _Mat(np.eye(2)),
                1e-4, structure="different",
            )
        np.testing.assert_array_equal(owner.current, z.values)
        np.testing.assert_allclose(
            owner.tangent_mass.values,
            owner.B0 + z.values[0] * owner.B1
            + z.values[1] * owner.B2,
            rtol=0., atol=1e-14,
        )


if __name__ == "__main__":
    unittest.main()
