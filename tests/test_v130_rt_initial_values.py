# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Actual DPHS.set_initial_value extracted from source, without native PETSc.

Exercise the production method with fake GetFEM array interfaces. These
tests validate coefficient ordering and index filtering; native RT0
MeshFem interpolation and PETSc assembly still require Docker.
"""

import ast
import logging
from pathlib import Path
import unittest

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
DPHS_SOURCE = ROOT / "scrimp" / "dphs.py"


def production_initial_value_function():
    tree = ast.parse(DPHS_SOURCE.read_text(encoding="utf-8"))
    cls = next(n for n in tree.body
               if isinstance(n, ast.ClassDef) and n.name == "DPHS")
    method = next(n for n in cls.body
                  if isinstance(n, ast.FunctionDef)
                  and n.name == "set_initial_value")
    code = compile(ast.Module(body=[method], type_ignores=[]),
                   str(DPHS_SOURCE), "exec")
    environment = {"np": np, "logging": logging, "rank": 1}
    exec(code, environment)
    return environment["set_initial_value"]


class MockMeshFem:
    def __init__(self, values, dofs):
        self.values = np.asarray(values)
        self.indices = np.asarray(dofs)

    def eval(self, expression, gbl, loc):
        return self.values.copy()

    def nbdof(self):
        return self.values.size

    def basic_dof_on_region(self, region):
        assert region == 2
        return self.indices

    def qdim(self):
        return 2


class MockPort:
    def __init__(self, mesh_fem, region=2):
        self.fem=mesh_fem
        self.region=region

    def get_flow(self):
        return "q"

    def get_effort(self):
        return "q"

    def get_region(self):
        return self.region

    def get_fem(self):
        return self.fem


class MockModel:
    def __init__(self, mesh_fem):
        self.fem=mesh_fem
        self.initial=None

    def mesh_fem_of_variable(self, name):
        assert name == "q"
        return self.fem

    def set_variable(self, name, vector):
        assert name == "q"
        self.initial=np.asarray(vector).copy()


class MockDPHS:
    def __init__(self, fem, region=2):
        self.initial_value_set={"q":False}
        self.lumped_variables=set()
        self.ports={"q":MockPort(fem,region)}
        self.gf_model=MockModel(fem)

    def set_from_vector(self, name, x):
        self.gf_model.set_variable(name,x)


class RestrictedInitialValueTests(unittest.TestCase):
    def setUp(self):
        self.set_initial_value=production_initial_value_function()

    def test_rt0_odd_filtered_flux_dofs_preserved_exactly(self):
        # Native observed: global RT0 nbdof=4794,
        # wave-region filtered DOFs=3111, qdim=2.
        all_coefficients=np.arange(4794,dtype=float).reshape(2,2397)
        filtered=np.arange(3111,dtype=int)
        self.assertEqual(filtered.size % 2, 1)
        sim=MockDPHS(MockMeshFem(all_coefficients,filtered))
        self.set_initial_value(sim, "q", "[0.,0.]")
        actual=sim.gf_model.initial
        np.testing.assert_array_equal(actual,all_coefficients.ravel()[filtered])
        self.assertEqual(actual.shape,(3111,))
        self.assertTrue(sim.initial_value_set["q"])

    def test_zero_rt0_field_is_not_shortened_or_rejected(self):
        values=np.zeros(4794,dtype=float)
        sim=MockDPHS(MockMeshFem(values,np.arange(3111)))
        self.set_initial_value(sim,"q","[0.,0.]")
        self.assertEqual(sim.gf_model.initial.shape,(3111,))
        self.assertTrue(np.all(sim.gf_model.initial==0))

    def test_vector_cg_interleaved_order_not_changed(self):
        values=np.arange(24,dtype=float).reshape(2,12)
        dofs=np.array([1,5,7,10,11,14,18,23])
        sim=MockDPHS(MockMeshFem(values,dofs))
        self.set_initial_value(sim,"q","[1.,2.]")
        np.testing.assert_array_equal(
            sim.gf_model.initial,values.reshape(-1)[dofs]
        )
        self.assertEqual(sim.gf_model.initial.size,len(dofs))

    def test_scalar_field_region_order_is_unchanged(self):
        values=np.array([1.,2.,3.,4.,5.])
        dofs=np.array([4,0,2],dtype=int)
        sim=MockDPHS(MockMeshFem(values,dofs))
        self.set_initial_value(sim,"q","1.")
        np.testing.assert_array_equal(sim.gf_model.initial,[5.,1.,3.])

    def test_full_region_does_not_modify_existing_evaluation_shape(self):
        values=np.arange(12).reshape(2,6)
        sim=MockDPHS(MockMeshFem(values,[]),region=None)
        self.set_initial_value(sim,"q","[1.,2.]")
        np.testing.assert_array_equal(sim.gf_model.initial,values)
        self.assertEqual(sim.gf_model.initial.shape,(2,6))

    def test_malformed_getfem_index_is_rejected_before_mutation(self):
        for indices in ([-1,2], [4794], [1.25,2.0]):
            with self.subTest(indices=indices):
                sim=MockDPHS(
                    MockMeshFem(np.zeros(4794),indices)
                )
                with self.assertRaisesRegex(ValueError,"DOF indices"):
                    self.set_initial_value(sim,"q","[0.,0.]")
                self.assertIsNone(sim.gf_model.initial)
                self.assertFalse(sim.initial_value_set["q"])

    def test_missing_variable_or_mismatched_fem_evaluation_fails(self):
        sim=MockDPHS(MockMeshFem(np.zeros(8),[0,1]))
        with self.assertRaises(ValueError):
            self.set_initial_value(sim,"missing","1.")
        # Invariant: evaluator gives exactly fem.nbdof coefficients.
        class BrokenFEM(MockMeshFem):
            def nbdof(self):
                return 9
        sim=MockDPHS(BrokenFEM(np.zeros(8),[0,1]))
        with self.assertRaisesRegex(ValueError,"expected 9"):
            self.set_initial_value(sim,"q","[0.,0.]")


if __name__ == "__main__":
    unittest.main()
