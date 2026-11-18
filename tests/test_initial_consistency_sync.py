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

import ast
import importlib.util
from pathlib import Path
import unittest

import numpy as np

# CI needs neither PETSc nor gmsh. Run the exact pure transfer helper;
# native MPI tests validate the entire DPHS integration.
SRC_ROOT = Path(__file__).resolve().parents[1]
HELPER_SOURCE = SRC_ROOT / "scrimp" / "utils" / "initial_state.py"
SPEC = importlib.util.spec_from_file_location("scrimp_initial_state_isolated", HELPER_SOURCE)
HELPER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(HELPER)


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
                scatter = MockScatter()
                global_vec = MockVector(np.zeros(4))
                gf_model = MockGetfem([0.0, 14.3679, 0.0, 0.0])
                pseudo_solution = MockVector(target)
                # This data transfer has no input or dependency on
                # linear_stiffness, regardless of the model flag value.
                _ = linear_stiffness
                HELPER.publish_accepted_pseudo_state(
                    scatter, pseudo_solution, global_vec, gf_model
                )
                np.testing.assert_array_equal(gf_model.from_variables(), target)
                np.testing.assert_array_equal(global_vec.values, target)
                self.assertEqual(scatter.calls, 1)
                # GetFEM owns an independent state copy:
                pseudo_solution.values[:] = -7.
                np.testing.assert_array_equal(gf_model.from_variables(), target)


    def test_dphs_initialization_actually_calls_transfer(self):
        """Check production wiring without importing native libraries."""
        module = ast.parse(
            (SRC_ROOT / "scrimp" / "dphs.py").read_text(encoding="utf-8")
        )
        dphs = next(
            node for node in module.body
            if isinstance(node, ast.ClassDef) and node.name == "DPHS"
        )
        functions = {
            node.name: node for node in dphs.body
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
        }
        self.assertIn("_sync_accepted_pseudo_initial_state", functions)
        self.assertIn("init_step", functions)
        helper_calls = {
            call.func.id for call in ast.walk(
                functions["_sync_accepted_pseudo_initial_state"]
            )
            if isinstance(call, ast.Call) and isinstance(call.func, ast.Name)
        }
        self.assertIn("publish_accepted_pseudo_state", helper_calls)
        init_calls = {
            call.func.attr for call in ast.walk(functions["init_step"])
            if isinstance(call, ast.Call) and isinstance(call.func, ast.Attribute)
        }
        self.assertIn("_sync_accepted_pseudo_initial_state", init_calls)

if __name__ == "__main__":
    unittest.main()
