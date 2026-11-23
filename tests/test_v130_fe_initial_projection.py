# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Source-only contracts for the universal default FE L2 initializer.

The expression translator is extracted from actual production AST.
This suite does not claim a native GetFEM/PETSc, MPI, or physics PASS.
"""

import ast
import math
import numpy as np
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]
PROJECTION = ROOT / "scrimp/utils/fe_initialization.py"
DPHS = ROOT / "scrimp/dphs.py"


def actual_real_coefficient_bridge():
    tree = ast.parse(PROJECTION.read_text(encoding="utf-8"))
    method = next(n for n in tree.body if isinstance(n, ast.FunctionDef)
                  and n.name == "real_getfem_projection_coefficients")
    env = {"np": np}
    exec(compile(ast.Module(body=[method], type_ignores=[]),
                 str(PROJECTION), "exec"), env)
    return env["real_getfem_projection_coefficients"]


def actual_transfer_policy():
    tree = ast.parse(PROJECTION.read_text(encoding="utf-8"))
    method = next(n for n in tree.body if isinstance(n, ast.FunctionDef)
                  and n.name == "fe_projection_transfer_policy")
    env = {}
    exec(compile(ast.Module(body=[method], type_ignores=[]),
                 str(PROJECTION), "exec"), env)
    return env["fe_projection_transfer_policy"]


def actual_rhs_scaler():
    tree = ast.parse(PROJECTION.read_text(encoding="utf-8"))
    method = next(n for n in tree.body if isinstance(n, ast.FunctionDef)
                  and n.name == "scale_getfem_projection_rhs")
    env = {"np": np}
    exec(compile(ast.Module(body=[method], type_ignores=[]),
                 str(PROJECTION), "exec"), env)
    return env["scale_getfem_projection_rhs"]


def actual_expression_translator():
    tree = ast.parse(PROJECTION.read_text(encoding="utf-8"))
    method = next(n for n in tree.body if isinstance(n, ast.FunctionDef)
                  and n.name == "python_expression_to_gwfl")
    env = {"ast": ast, "math": math, "_SCALAR_FUNCTIONS": frozenset((
        "sin", "cos", "tan", "asin", "acos", "atan", "atan2",
        "exp", "log", "log10", "sqrt", "abs", "sqr",
    ))}
    exec(compile(ast.Module(body=[method], type_ignores=[]),
                 str(PROJECTION), "exec"), env)
    return env["python_expression_to_gwfl"]


class GetFEMMPITransferPolicyTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.policy = staticmethod(actual_transfer_policy())

    def test_mpi2_replicated_rejected(self):
        with self.assertRaisesRegex(ValueError, "requires .*distributed"):
            self.policy(2, "replicated")

    def test_distributed_linear_mass_is_coo_eligible(self):
        self.assertTrue(self.policy(1, "distributed"))
        self.assertTrue(self.policy(2, "distributed"))

    def test_serial_compatibility(self):
        self.assertFalse(self.policy(1, "replicated"))

    def test_invalid_modes_and_mpi_size_rejected(self):
        with self.assertRaises(ValueError):
            self.policy(0, "distributed")
        with self.assertRaises(ValueError):
            self.policy(2, "invalid")

    def test_uses_existing_bridge_without_separate_coo_flag(self):
        code = PROJECTION.read_text(encoding="utf-8")
        self.assertIn("coo_eligible = fe_projection_transfer_policy(", code)
        self.assertIn("distributed_coo_eligible=coo_eligible", code)
        self.assertNotIn("SCRIMP_FE_INITIAL_COO", code)


class GetFEMMPIProjectionRHSScaleTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.scale_rhs = staticmethod(actual_rhs_scaler())

    def test_normal_scale_unchanged(self):
        np.testing.assert_array_equal(
            self.scale_rhs(np.array([1., -2., 3.]), 1.0),
            [1., -2., 3.],
        )

    def test_duplicate_rhs_is_multiplied_not_divided(self):
        # Known SCRIMP GetFEM MPI-2 calibration returns 1/P
        # for P-fold duplicated RHS; never invert this factor.
        np.testing.assert_array_equal(
            self.scale_rhs(np.array([2., 4., -6.]), 0.5),
            [1., 2., -3.],
        )

    def test_bad_scale_and_rhs_rejected(self):
        for factor in (0., -1., float("inf"), float("nan")):
            with self.subTest(factor=factor):
                with self.assertRaises(ValueError):
                    self.scale_rhs(np.array([1.]), factor)
        with self.assertRaises(ValueError):
            self.scale_rhs(np.array([float("nan")]), 1.)

    def test_production_uses_the_calibrated_rhs(self):
        source = PROJECTION.read_text(encoding="utf-8")
        self.assertIn("gf_rhs = scale_getfem_projection_rhs(", source)
        self.assertIn("gf_rhs[lo:hi],", source)
        self.assertNotIn("gf_rhs[lo:hi] / float(rhs_scale)", source)


class RealGetFEMCoefficientBridgeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.bridge = staticmethod(actual_real_coefficient_bridge())

    def test_complex_petsc_real_coefficients_become_float64(self):
        value = np.array([1+0j, -2+1e-15j, 0.5-1e-15j],
                         dtype=np.complex128)
        output = self.bridge(value)
        self.assertEqual(output.dtype, np.dtype("float64"))
        self.assertTrue(output.flags.c_contiguous)
        np.testing.assert_array_equal(output, [1., -2., 0.5])

    def test_real_coefficients_stay_exact(self):
        data = np.array([0., 1.5, -7.], dtype=np.float64)
        output = self.bridge(data)
        np.testing.assert_array_equal(output, data)
        self.assertTrue(output.flags.c_contiguous)

    def test_reject_significant_imaginary_component(self):
        with self.assertRaisesRegex(ValueError, "imaginary"):
            self.bridge(np.array([1+1e-3j, 2+0j]))

    def test_reject_invalid_real_projection_coefficients(self):
        for case in (np.array([float("nan")]),
                     np.array([float("inf")]),
                     np.array([], dtype=float),
                     np.array([[1.,2.]]),
                     np.array(["1"], dtype=str)):
            with self.subTest(case=repr(case)):
                with self.assertRaises(ValueError):
                    self.bridge(case)


class UniversalFEInitialExpressionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.translate = staticmethod(actual_expression_translator())

    def test_scalar_affine_and_nonpolynomial_numpy(self):
        u = self.translate("5.*np.exp(-25*((x-0.6)**2+y*y))", 1)
        self.assertIn("exp(", u)
        self.assertIn("pow(", u)
        self.assertIn("X(1)", u)
        self.assertIn("X(2)", u)
        self.assertNotIn("np.", u)
        self.assertNotIn("**", u)

    def test_vector_xy_and_numpy_array(self):
        for s in ("[x+y,2.*y-x]", "np.array([x+y,2.*y-x])"):
            u = self.translate(s, 2)
            self.assertTrue(u.startswith("["))
            self.assertIn(";", u)
            self.assertIn("X(1)", u)
            self.assertIn("X(2)", u)

    def test_any_number_of_vector_components(self):
        u = self.translate("[x,y,z]", 3)
        self.assertEqual(u, "[X(1);X(2);X(3)]")

    def test_math_constants_and_power(self):
        u = self.translate("np.pi*np.sin(x)+math.sqrt(y**2)", 1)
        self.assertIn("sin(", u)
        self.assertIn("sqrt(", u)
        self.assertIn("pow(", u)
        self.assertNotIn("math.", u)

    def test_bad_vector_shape_fails_closed(self):
        with self.assertRaisesRegex(ValueError, "components"):
            self.translate("[x,y]", 3)
        with self.assertRaisesRegex(ValueError, "Expected"):
            self.translate("x", 2)

    def test_unsafe_and_unavailable_expressions_fail_closed(self):
        for expr in (
            "__import__('os').system('true')", "np.random.rand(5)",
            "np.where(x>0,x,0)", "foo(x)", "x[0]", "x^2",
            "np.array(x)", "float('nan')", "1/0 if x>0 else 0",
        ):
            with self.subTest(expression=expr):
                with self.assertRaises((ValueError, SyntaxError)):
                    self.translate(expr, 1)

    def test_production_uses_sparse_petsc_not_dense_numpy(self):
        src = PROJECTION.read_text(encoding="utf-8")
        self.assertIn("convert_gmm_to_petsc(", src)
        self.assertIn("PETSc.Mat().createAIJ(", src)
        self.assertIn('pc.setFactorSolverType("mumps")', src)
        self.assertIn("getOwnershipRange()", src)
        self.assertIn("PETSc.Scatter().toAll(x)", src)
        self.assertIn("real_getfem_projection_coefficients(", src)
        self.assertNotIn("np.linalg.solve(", src)

    def test_native_mpi_pilot_uses_distributed_bridge(self):
        runner = (
            ROOT / "scripts/run-v130-fe-initial-mpi.sh"
        ).read_text(encoding="utf-8")
        self.assertIn(
            "--env SCRIMP_GETFEM_MATRIX_MODE=distributed", runner
        )
        self.assertIn(
            "--env SCRIMP_GETFEM_BRIDGE_DISTRIBUTED=coo", runner
        )
        self.assertNotIn(
            "--env SCRIMP_GETFEM_MATRIX_MODE=replicated", runner
        )
        self.assertNotIn("SCRIMP_FE_INITIAL_COO", runner)
        pilot = (
            ROOT / "sandbox/probe_v130_fe_initial_mpi.py"
        ).read_text(encoding="utf-8")
        self.assertIn('if mode != "distributed":', pilot)
        self.assertIn('"getfem_matrix_mode": mode', pilot)

    def test_production_default_l2_and_explicit_legacy(self):
        tree = ast.parse(DPHS.read_text(encoding="utf-8"))
        cls = next(x for x in tree.body if isinstance(x, ast.ClassDef)
                   and x.name == "DPHS")
        method = next(x for x in cls.body if isinstance(x, ast.FunctionDef)
                      and x.name == "set_initial_value")
        self.assertEqual(len(method.args.kwonlyargs), 2)
        self.assertEqual(
            [arg.arg for arg in method.args.kwonlyargs],
            ["method", "projection_mim"],
        )
        self.assertEqual(
            [v.value for v in method.args.kw_defaults],
            ["l2", None],
        )
        actual = ast.get_source_segment(DPHS.read_text(encoding="utf-8"),
                                        method)
        self.assertIn('if method == "l2":', actual)
        self.assertIn("project_fe_initial_value(", actual)
        self.assertIn('if method not in ("interpolation", "l2"):', actual)
        self.assertIn("if name_variable in self.lumped_variables:", actual)
        self.assertNotIn(
            'raise ValueError("L2 FE projection does not apply to lumped states")',
            actual
        )


if __name__ == "__main__":
    unittest.main()
