# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Accepted pseudo-TS state must always reach GetFEM, regardless of linear flags.

A prior source-brick bug accidentally provided this synchronization through
DPHS.monitor() only when linear_stiffness=False. The fix must never revert
the correct RHS-only source classification or alter Newton tangent assembly.
"""

from contextlib import contextmanager
import unittest

import numpy as np

from scrimp.dphs import DPHS


class MockProfiler:
    def __init__(self):
        self.calls = []

    @contextmanager
    def measure(self, label):
        self.calls.append(label)
        yield


class MockVector:
    def __init__(self, values):
        self.values = np.asarray(values, dtype=np.float64)

    def getArray(self, readonly=False):
        return self.values


class MockScatter:
    def __init__(self):
        self.calls = 0

    def scatter(self, source, destination):
        self.calls += 1
        destination.values[:] = source.values


class MockGetfem:
    def __init__(self, values):
        self.data = np.asarray(values, dtype=np.float64).copy()

    def to_variables(self, values):
        self.data = np.asarray(values, dtype=np.float64).copy()

    def from_variables(self):
        return self.data.copy()


class TestPseudoStateSynchronization(unittest.TestCase):
    def test_accepted_full_state_written_independently_of_flags(self):
        target = np.asarray([0.0047176, 14.3678, 0.0202897, 5.76279])
        for linear_stiffness in (True, False):
            with self.subTest(linear_stiffness=linear_stiffness):
                obj = DPHS.__new__(DPHS)
                obj.linear_mass = True
                obj.linear_stiffness = linear_stiffness
                obj.profiler = MockProfiler()
                scatter = MockScatter()
                obj.Scatter = (scatter, None)
                obj.GlobalVec = MockVector(np.zeros(4))
                obj.gf_model = MockGetfem([0.0, 14.3679, 0.0, 0.0])
                pseudo_solution = MockVector(target)
                obj._sync_accepted_pseudo_initial_state(pseudo_solution)

                np.testing.assert_array_equal(
                    obj.gf_model.from_variables(), target
                )
                np.testing.assert_array_equal(obj.GlobalVec.values, target)
                self.assertEqual(scatter.calls, 1)
                self.assertEqual(
                    obj.profiler.calls, ["initial.petsc.accepted_state_sync"]
                )
                self.assertEqual(obj.linear_stiffness, linear_stiffness)
                # State must be independently owned by GetFEM:
                pseudo_solution.values[:] = -7.
                np.testing.assert_array_equal(
                    obj.gf_model.from_variables(), target
                )


if __name__ == "__main__":
    unittest.main()
