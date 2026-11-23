# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Source-only regression for all production GetFEM -> PETSc COO paths.

This tests the actual production AST without importing GetFEM/PETSc.
Native MPI1/2 nonlinear Jacobian and dynamic-graph tests remain
separate mandatory physical gates.
"""

import ast
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]
DPHS = ROOT / "scrimp/dphs.py"
LINALG = ROOT / "scrimp/utils/linalg.py"


class NonlinearCOOAssemblySourceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.dphs = DPHS.read_text(encoding="utf-8")
        cls.linalg = LINALG.read_text(encoding="utf-8")
        cls.dp_ast = ast.parse(cls.dphs)
        cls.la_ast = ast.parse(cls.linalg)

    def test_all_four_production_getfem_matrix_bridges_are_coo_eligible(self):
        cls = next(n for n in self.dp_ast.body
                   if isinstance(n, ast.ClassDef) and n.name == "DPHS")
        names = (
            "assemble_mass", "assemble_stiffness",
            "assemble_nl_mass", "assemble_nl_stiffness"
        )
        for name in names:
            method = next(n for n in cls.body
                          if isinstance(n, ast.FunctionDef) and n.name == name)
            calls = [n for n in ast.walk(method)
                     if isinstance(n, ast.Call)
                     and isinstance(n.func, ast.Name)
                     and n.func.id == "convert_gmm_to_petsc"]
            with self.subTest(method=name):
                self.assertEqual(len(calls), 1)
                kw = {arg.arg: arg.value for arg in calls[0].keywords}
                self.assertIn("distributed_coo_eligible", kw)
                self.assertIsInstance(kw["distributed_coo_eligible"],
                                      ast.Constant)
                self.assertIs(kw["distributed_coo_eligible"].value, True)

    def test_distributed_is_actual_implicit_default(self):
        self.assertIn(
            '"SCRIMP_GETFEM_MATRIX_MODE", "distributed"', self.linalg
        )
        projection = (
            ROOT / "scrimp/utils/fe_initialization.py"
        ).read_text(encoding="utf-8")
        self.assertIn(
            'os.environ.get("SCRIMP_GETFEM_MATRIX_MODE", "distributed")',
            projection
        )

    def test_bridge_default_is_coo_not_ijv(self):
        fn = next(n for n in self.la_ast.body
                  if isinstance(n, ast.FunctionDef)
                  and n.name == "convert_gmm_to_petsc")
        arg = next(i for i,v in enumerate(fn.args.kwonlyargs)
                   if v.arg == "distributed_coo_eligible")
        self.assertIs(fn.args.kw_defaults[arg].value, True)
        self.assertIn(
            '"SCRIMP_GETFEM_BRIDGE_DISTRIBUTED", "coo"', self.linalg
        )
        self.assertIn('mode == "distributed" and distributed_method == "coo"',
                      self.linalg)

    def test_coo_refreshes_pattern_and_keeps_owned_diagonal(self):
        fn = next(n for n in self.la_ast.body
                  if isinstance(n, ast.FunctionDef)
                  and n.name == "convert_gmm_to_petsc")
        source = ast.get_source_segment(self.linalg, fn)
        self.assertIn("B.setPreallocationCOO(coo_i, coo_j)", source)
        self.assertIn("B.setValuesCOO(", source)
        self.assertLess(source.index("B.setPreallocationCOO("),
                        source.index("B.setValuesCOO("))
        self.assertIn("np.arange(lo, hi, dtype=PETSc.IntType)", source)
        self.assertIn("coo_i = np.ascontiguousarray(", source)
        self.assertIn("coo_j = np.ascontiguousarray(", source)

    def test_nonlinear_tangents_keep_same_mat_handle_and_dynamic_axpy(self):
        cls = next(n for n in self.dp_ast.body
                   if isinstance(n, ast.ClassDef) and n.name == "DPHS")
        for name in ("assemble_nl_mass", "assemble_nl_stiffness"):
            method = next(n for n in cls.body
                          if isinstance(n, ast.FunctionDef) and n.name == name)
            segment = ast.get_source_segment(self.dphs, method)
            with self.subTest(name=name):
                self.assertIn(
                    "PETSc.Mat.Structure.DIFFERENT_NONZERO_PATTERN", segment
                )
                self.assertIn(".zeroEntries()", segment)
                self.assertIn(".axpy(", segment)
                self.assertNotIn("self.nl_mass = self.nl_mass.copy()", segment)
                self.assertNotIn(
                    "self.nl_stiffness = self.nl_stiffness.copy()", segment
                )


if __name__ == "__main__":
    unittest.main()
