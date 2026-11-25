# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
"""Source-only unit tests of SCRIMP's Nédélec FEM aliases.

These use the *actual FEM class AST* plus a small GetFEM-constructor stub.
They test descriptor construction and fail-closed behavior without loading
PETSc or native GetFEM. They are not native Maxwell acceptance tests.
"""
import ast
import logging
from pathlib import Path
import unittest

SOURCE = Path(__file__).resolve().parents[1] / "scrimp" / "fem.py"


class StubMesh:
    def __init__(self, dim):
        self.dimension = dim

    def dim(self):
        return self.dimension


class StubMeshFem:
    def __init__(self, mesh, qdim):
        self.mesh, self.dimension = mesh, qdim
        self.descriptor = None

    def set_fem(self, descriptor):
        self.descriptor = descriptor

    def display(self):
        return None

    def qdim(self):
        return self.dimension

    def nbdof(self):
        return 10


class FakeGetFEM:
    MeshFem = StubMeshFem

    def __init__(self):
        self.unavailable = set()
        self.constructed = []

    def Fem(self, descriptor):
        self.constructed.append(descriptor)
        if descriptor in self.unavailable:
            raise AssertionError("GetFEM family unavailable")
        return descriptor


class ScrimpNedelecAliasesTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.source = SOURCE.read_text(encoding="utf-8")
        tree = ast.parse(cls.source)
        definitions = [n for n in tree.body
                       if isinstance(n, ast.ClassDef) and n.name == "FEM"]
        assert len(definitions) == 1
        cls._class_ast = definitions[0]

    def setUp(self):
        self.native = FakeGetFEM()
        module = ast.Module(body=[self._class_ast], type_ignores=[])
        ast.fix_missing_locations(module)
        namespace = {
            "gf": self.native, "logging": logging, "rank": 1,
        }
        exec(compile(module, str(SOURCE), "exec"), namespace)
        self.FEM = namespace["FEM"]

    def make(self, family, order, dim=3):
        fem = self.FEM("E", order, family)
        fem.set_mesh(StubMesh(dim))
        fem.set_dim(dim)
        return fem

    def test_both_patched_families_and_long_names(self):
        descriptors = {
            ("NED1", 0, 3): "FEM_NEDELEC1(3,0)",
            ("NED1", 1, 3): "FEM_NEDELEC1(3,1)",
            ("NEDELEC1", 2, 2): "FEM_NEDELEC1(2,2)",
            ("NED2", 1, 3): "FEM_NEDELEC2(3,1)",
            ("NED2", 2, 3): "FEM_NEDELEC2(3,2)",
            ("NEDELEC2", 1, 2): "FEM_NEDELEC2(2,1)",
        }
        for (family, order, dim), expected in descriptors.items():
            with self.subTest(family=family, order=order, dim=dim):
                fem = self.make(family, order, dim)
                fem.set_fem()
                self.assertTrue(fem.get_isSet())
                self.assertEqual(fem.get_type(), family)
                self.assertEqual(fem.get_order(), order)
                self.assertEqual(fem.get_fem().descriptor, expected)

    def test_existing_RT_and_BDM_shortcuts_unchanged(self):
        for family, k, expected in (
            ("RT", 0, "FEM_RTK(3,0)"),
            ("BDM", 1, "FEM_BDMK(3,1)"),
            ("CG", 1, "FEM_PK(3,1)"),
            ("DG", 1, "FEM_PK_DISCONTINUOUS(3,1)"),
            ("FEM_NEDELEC(3)", -1, "FEM_NEDELEC(3)"),
        ):
            with self.subTest(family=family):
                fem = self.make(family, k)
                fem.set_fem()
                self.assertEqual(fem.get_fem().descriptor, expected)

    def test_NED2_explicit_k1_k2_only(self):
        for order in (-1, 0, 3, 4, 1.5):
            with self.subTest(k=order):
                fem = self.make("NED2", order)
                with self.assertRaisesRegex(ValueError, "k=1 or k=2"):
                    fem.set_fem()
                self.assertFalse(fem.get_isSet())

    def test_NED1_nonnegative_integer_only(self):
        for order in (-1, -2, 1.5):
            with self.subTest(k=order):
                fem = self.make("NED1", order)
                with self.assertRaisesRegex(ValueError, "k>=0"):
                    fem.set_fem()
                self.assertFalse(fem.get_isSet())

    def test_simplices_dimension_two_or_three_only(self):
        for family, k in (("NED1", 0), ("NED2", 1)):
            for dim in (1, 4):
                with self.subTest(family=family, dim=dim):
                    fem = self.make(family, k, dim)
                    with self.assertRaisesRegex(ValueError, "dimension 2/3"):
                        fem.set_fem()

    def test_missing_native_patch_is_never_converted_to_CG(self):
        for family, k, desc in (
            ("NED1", 1, "FEM_NEDELEC1(3,1)"),
            ("NED2", 1, "FEM_NEDELEC2(3,1)"),
            ("NEDELEC2", 2, "FEM_NEDELEC2(3,2)"),
        ):
            with self.subTest(family=family):
                fem = self.make(family, k)
                self.native.unavailable.add(desc)
                with self.assertRaisesRegex(RuntimeError, "unavailable in GetFEM"):
                    fem.set_fem()
                self.assertEqual(fem.get_type(), family)
                self.assertEqual(fem.get_order(), k)
                self.assertFalse(fem.get_isSet())
                self.assertNotIn("FEM_PK(3,1)", self.native.constructed)
                self.native.unavailable.clear()


if __name__ == "__main__":
    unittest.main()
