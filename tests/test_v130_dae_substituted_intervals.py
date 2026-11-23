# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Test actual production _intervals without PETSc or GetFEM.

These tests prove coefficient classification only; native DAE
feasibility and physics remain separate gates.
"""
import ast
from pathlib import Path
import unittest
import numpy as np

SOURCE = Path(__file__).resolve().parents[1] / "scrimp/utils/dae_initialization.py"

def actual_intervals():
    tree = ast.parse(SOURCE.read_text(encoding="utf-8"))
    fn = next(n for n in tree.body if isinstance(n, ast.FunctionDef)
              and n.name == "_intervals")
    env = {"np": np}
    exec(compile(ast.Module(body=[fn], type_ignores=[]), str(SOURCE), "exec"), env)
    return env["_intervals"]

class State:
    def __init__(self, name): self.name = name
    def get_name(self): return self.name

class Port:
    def __init__(self, flow, effort, algebraic, substituted):
        self.flow, self.effort = flow, effort
        self.algebraic, self.substituted = algebraic, substituted
    def get_flow(self): return self.flow
    def get_effort(self): return self.effort
    def get_algebraic(self): return self.algebraic
    def get_substituted(self): return self.substituted

class Model:
    def __init__(self, intervals): self.intervals = intervals
    def nbdof(self): return max(start + size for start, size in self.intervals.values())
    def interval_of_variable(self, name): return self.intervals[name]

class System:
    def __init__(self, intervals, states, ports, *, linear=True, lumped=False):
        self.gf_model = Model(intervals)
        self.states = {x: State(x) for x in states}
        self.ports = {str(i): x for i,x in enumerate(ports)}
        self.linear_mass = self.linear_stiffness = linear
        self.lumped_variables = {"scalar"} if lumped else set()

class SubstitutedDAEIntervals(unittest.TestCase):
    @classmethod
    def setUpClass(cls): cls.classify = staticmethod(actual_intervals())

    def test_substituted_dynamic_does_not_add_effort(self):
        s = System({"q":(0,5), "Y":(5,2), "U":(7,2)}, ["q"],
                   [Port("q","e_q",False,True), Port("Y","U",True,False)])
        n,d,a = self.classify(s)
        self.assertEqual((n,d), (9,{"q":(0,5)}))
        np.testing.assert_array_equal(a,np.arange(5,9))

    def test_substituted_algebraic_identical_names_and_odd_intrinsic_dofs(self):
        s = System({"T":(0,3), "e_Q":(3,7), "Y":(10,2), "U":(12,2)}, ["T"],
                   [Port("T","T",False,True), Port("e_Q","e_Q",True,True),
                    Port("Y","U",True,False)])
        n,d,a = self.classify(s)
        self.assertEqual((n,d), (14,{"T":(0,3)}))
        np.testing.assert_array_equal(a,np.arange(3,14))

    def test_substituted_algebraic_distinct_non_registered_effort(self):
        s = System({"p":(0,3), "f_Q":(3,5)}, ["p"],
                   [Port("p","e_p",False,True), Port("f_Q","e_Q",True,True)])
        _,_,a = self.classify(s)
        np.testing.assert_array_equal(a,np.arange(3,8))

    def test_unsubstituted_dynamic_effort_still_algebraic(self):
        s = System({"q":(0,2),"e_q":(2,2)}, ["q"],
                   [Port("q","e_q",False,False)])
        _,_,a = self.classify(s)
        np.testing.assert_array_equal(a,[2,3])

    def test_repeated_algebraic_variable_refused(self):
        s = System({"q":(0,2),"Y":(2,2),"U":(4,2)}, ["q"],
                   [Port("Y","U",True,False),Port("Y","unused",True,True)])
        with self.assertRaisesRegex(ValueError,"Repeated"): self.classify(s)

    def test_unclassified_getfem_dofs_refused(self):
        s = System({"q":(0,2),"Y":(2,2),"U":(4,2),"orphan":(6,1)}, ["q"],
                   [Port("Y","U",True,False)])
        with self.assertRaisesRegex(ValueError,"Unclassified"): self.classify(s)

    def test_overlapping_getfem_dofs_refused(self):
        s = System({"q":(0,2),"Y":(1,2),"U":(4,2)}, ["q"],
                   [Port("Y","U",True,False)])
        with self.assertRaisesRegex(ValueError,"overlap"): self.classify(s)

    def test_lumped_guard_unchanged(self):
        s = System({"q":(0,2),"Y":(2,2),"U":(4,2)}, ["q"],
                   [Port("Y","U",True,False)],lumped=True)
        with self.assertRaisesRegex(ValueError,"lumped"): self.classify(s)

    def test_nonlinear_affine_substituted_still_refused(self):
        s = System({"q":(0,2),"e_Q":(2,3)}, ["q"],
                   [Port("q","q",False,True),Port("e_Q","e_Q",True,True)],
                   linear=False)
        with self.assertRaisesRegex(ValueError,"Nonlinear-affine.*substituted"):
            self.classify(s,allow_nonlinear=True)

    def test_nonlinear_guard_and_legacy_unsubstituted_mode(self):
        s = System({"q":(0,2),"Y":(2,2),"U":(4,2)}, ["q"],
                   [Port("Y","U",True,False)],linear=False)
        with self.assertRaisesRegex(ValueError,"linear mass"): self.classify(s)
        n,_,a = self.classify(s,allow_nonlinear=True)
        self.assertEqual(n,6)
        np.testing.assert_array_equal(a,[2,3,4,5])

if __name__ == "__main__":
    unittest.main()
