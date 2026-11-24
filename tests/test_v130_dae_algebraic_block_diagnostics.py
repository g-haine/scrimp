# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Source-only execution of actual DAE diagnostic AST (no PETSc/GetFEM)."""

import ast
from pathlib import Path
from types import SimpleNamespace
import unittest

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT/"scrimp/utils/dae_initialization.py"


def actual_diagnostic():
    tree=ast.parse(SOURCE.read_text(encoding="utf-8"))
    fn=next(x for x in tree.body
            if isinstance(x,ast.FunctionDef)
            and x.name=="_algebraic_residual_by_variable")
    env={"np":np,"MPI":SimpleNamespace(SUM="sum")}
    exec(compile(ast.Module(body=[fn],type_ignores=[]),
                 str(SOURCE),"exec"),env)
    return env["_algebraic_residual_by_variable"]


class DAEBlocksTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.fn=staticmethod(actual_diagnostic())

    def test_native_block_norms_partition_without_differential_rows(self):
        class Port:
            def __init__(self, name, flow, effort, algebraic, substituted,region):
                self.data=name,flow,effort,algebraic,substituted,region
            def get_name(self): return self.data[0]
            def get_flow(self): return self.data[1]
            def get_effort(self): return self.data[2]
            def get_algebraic(self): return self.data[3]
            def get_substituted(self): return self.data[4]
            def get_region(self): return self.data[5]
        intervals={"p":(0,2),"Y_bnd":(2,2),"U_bnd":(4,1),"e_Q":(5,2)}
        ports={
            "boundary":Port("Boundary","Y_bnd","U_bnd",True,False,20),
            "heat":Port("Heat flux","e_Q","e_Q",False,True,1),
        }
        system=SimpleNamespace(
            ports=ports,
            gf_model=SimpleNamespace(interval_of_variable=lambda n:intervals[n]),
        )
        class Residual:
            def getOwnershipRange(self): return (0,7)
            def getArray(self,readonly=True):
                return np.array([100.,100.,3.,4.,12.,0.,5.])
        class Comm:
            def allreduce(self,val,op=None): return val
        r=self.fn(system,Residual(),Comm())
        self.assertEqual(set(r),{"Y_bnd","U_bnd"})
        self.assertEqual(r["Y_bnd"]["residual_l2"],5.)
        self.assertEqual(r["U_bnd"]["residual_l2"],12.)
        self.assertEqual(r["Y_bnd"]["region"],20)
        self.assertAlmostEqual(
            sum(x["residual_l2"]**2 for x in r.values()),
            13.**2,
        )

    def test_algebraic_residual_is_reported_only_after_true_failure(self):
        source=SOURCE.read_text(encoding="utf-8")
        fn=ast.parse(source)
        initializer=next(x for x in fn.body
                         if isinstance(x,ast.FunctionDef)
                         and x.name=="initialize_linear_algebraic")
        body=ast.get_source_segment(source,initializer)
        self.assertIn('if final_norm > tolerance:',body)
        self.assertIn("_algebraic_residual_by_variable(",body)
        self.assertIn('"algebraic_residual_by_variable"',body)
        self.assertIn('"dominant_algebraic_variables"',body)
        self.assertIn("system.initialization_report = report",body)
        self.assertIn("system.gf_model.to_variables(initial.copy())",body)
        self.assertIn("raise InconsistentDAEInitialData(",body)

    def test_original_heat_wave_has_explicit_algebraic_gate_and_bad_outer_trace(self):
        source=(ROOT/"examples/heat_wave.py").read_text(encoding="utf-8")
        self.assertIn('init_step_ts_type="algebraic"',source)
        self.assertIn('hw.set_control("Boundary", "0.")',source)
        self.assertIn(
            'hw.set_initial_value("p", "5.*np.exp(-25*((x-0.6)*(x-0.6)+y*y))")',
            source,
        )
        self.assertGreater(5*np.exp(-4),0.09)


if __name__=="__main__":
    unittest.main()
