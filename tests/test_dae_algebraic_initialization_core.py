# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Source/isolated helper regression; no native GetFEM/PETSc in GitHub CI.

Native MPI 1/2 scientific and lifecycle qualification is a separate
mandatory gate before changing the default initialization strategy.
"""

import ast
import importlib.util
from pathlib import Path
import sys
import types
import unittest
from unittest.mock import patch

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
HELPER = ROOT / "scrimp/utils/dae_initialization.py"
DPHS = ROOT / "scrimp/dphs.py"


class FakeOptions:
    def __init__(self, values=None):
        self.values = values or {}

    def getString(self, key, default=None):
        return str(self.values.get(key, default))


class FakeState:
    def __init__(self, name):
        self.name = name

    def get_name(self):
        return self.name


class FakePort:
    def __init__(self, flow, effort, *, algebraic=False,
                 substituted=False):
        self.flow = flow
        self.effort = effort
        self.algebraic = algebraic
        self.substituted = substituted

    def get_algebraic(self):
        return self.algebraic

    def get_flow(self):
        return self.flow

    def get_effort(self):
        return self.effort

    def get_substituted(self):
        return self.substituted


class FakeGFModel:
    def __init__(self, slots):
        self.slots = slots

    def nbdof(self):
        return max(start+width for start, width in self.slots.values())

    def interval_of_variable(self, name):
        return self.slots[name]


class FakeSystem:
    def __init__(self):
        self.time_scheme = FakeOptions()
        self.linear_mass = True
        self.linear_stiffness = True
        self.lumped_variables = set()
        self.states = {"q": FakeState("q"), "p": FakeState("p")}
        self.ports = {
            "q": FakePort("q", "e_q"),
            "p": FakePort("p", "e_p"),
            "L": FakePort("U_L", "Y_L", algebraic=True),
        }
        self.gf_model = FakeGFModel({
            "q": (0,2), "p": (2,1), "e_q": (3,2),
            "e_p": (5,1), "U_L": (6,1), "Y_L": (7,1),
        })


def helper_source():
    return HELPER.read_text(encoding="utf-8")


def load_lightweight_helper():
    """Exercise unchanged helper functions without installing PETSc."""
    fake_petsc = types.ModuleType("petsc4py")
    fake_petsc.PETSc = types.SimpleNamespace()
    fake_mpi = types.ModuleType("mpi4py")
    fake_mpi.MPI = types.SimpleNamespace()
    with patch.dict(
        sys.modules, {"petsc4py": fake_petsc, "mpi4py": fake_mpi},
    ):
        spec = importlib.util.spec_from_file_location(
            "scrimp_dae_helper_isolated", HELPER
        )
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
    return module


class TestCoreDAEInitialization(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.helper = load_lightweight_helper()

    def test_valid_linear_layout_is_exact_partition(self):
        n, diff, idx = self.helper._intervals(FakeSystem())
        self.assertEqual(n, 8)
        self.assertEqual(diff["q"], (0,2))
        self.assertEqual(diff["p"], (2,1))
        np.testing.assert_array_equal(
            idx, np.array([3,4,5,6,7], dtype=np.int64)
        )

    def test_unsupported_models_refused(self):
        for mode in ("linear_mass", "linear_stiffness"):
            with self.subTest(mode=mode):
                sys = FakeSystem()
                setattr(sys, mode, False)
                with self.assertRaisesRegex(ValueError, "only linear"):
                    self.helper._intervals(sys)
        s = FakeSystem()
        s.lumped_variables.add("a")
        with self.assertRaisesRegex(ValueError, "lumped"):
            self.helper._intervals(s)
        s = FakeSystem()
        s.ports["L"].substituted = True
        with self.assertRaisesRegex(ValueError, "substituted"):
            self.helper._intervals(s)

    def test_unclassified_or_overlapping_variables_are_refused(self):
        s = FakeSystem()
        s.gf_model.slots["Y_L"] = (8, 1)
        with self.assertRaisesRegex(ValueError, "Unclassified"):
            self.helper._intervals(s)
        s = FakeSystem()
        s.gf_model.slots["U_L"] = (5,1)
        with self.assertRaisesRegex(ValueError, "Unclassified|overlap"):
            self.helper._intervals(s)

    def test_separate_feasibility_accuracy_and_pass_limits(self):
        s = FakeSystem()
        self.assertEqual(
            self.helper._options(s),
            (1e-10,1e-8,1e-11,3,1e-13,1e-13),
        )
        for vals in (
            {"init_algebraic_accuracy": "nan"},
            {"init_algebraic_atol": "-1"},
            {"init_algebraic_max_passes": "0"},
            {"init_algebraic_max_passes": "21"},
        ):
            with self.subTest(values=vals):
                s.time_scheme = FakeOptions(vals)
                with self.assertRaises(ValueError):
                    self.helper._options(s)

    def test_candidate_default_preserves_explicit_pseudo_bdf_flow(self):
        tree = ast.parse(DPHS.read_text(encoding="utf-8"))
        cls = next(
            item for item in tree.body
            if isinstance(item, ast.ClassDef) and item.name == "DPHS"
        )
        init = next(
            item for item in cls.body
            if isinstance(item, ast.FunctionDef) and item.name == "init_step"
        )
        source = ast.get_source_segment(
            DPHS.read_text(encoding="utf-8"), init
        )
        self.assertLess(
            source.index('if init_type == "algebraic":'),
            source.index('init_dt = float(self.time_scheme["init_step_dt"])'),
        )
        self.assertIn("self.initialization_report = initialize_linear_algebraic(self)", source)
        self.assertIn('if init_type not in ("pseudo", "bdf"):', source)
        self.assertIn('"init_step_ts_type", "algebraic"', source)
        self.assertIn('if init_type not in ("pseudo", "bdf"):', source)
        self.assertIn('init_step_from_negative_time', source)
        self.assertIn('TS.solve(InitVec)', source)

    def test_uses_native_ifunction_only_and_distributed_submatrix(self):
        src = helper_source()
        for required in (
            "system.IFunction(None, t0, z, zdot, residual)",
            "system.stiffness.createSubMatrix(is_alg, is_alg)",
            "PETSc.KSP.Type.LSQR",
            "dtype=PETSc.IntType",
            "system.Scatter[0].scatter(z, system.GlobalVec)",
            "np.array_equal(",
            "InconsistentDAEInitialData",
            "system.initialization_report = report",
            "if not accepted:",
            "max_passes",
        ):
            self.assertIn(required, src)
        self.assertNotIn("from sandbox", src)
        self.assertNotIn("import scipy", src)
        self.assertNotIn("PETSc.TS()", src)


if __name__ == "__main__":
    unittest.main()
