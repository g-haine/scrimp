# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Host-only GetFEM-compatible fake test of public Port.compute sampling.

This executes the actual AST of Port.compute and its result getters:
no fake sampling formulas replacing source code, no PETSc/GetFEM import,
no MPI claim. True native integrated physics needs a Docker host.
"""

import ast
import logging
from pathlib import Path
from types import SimpleNamespace
import time
import unittest

import numpy as np


ROOT=Path(__file__).resolve().parents[1]
PORT=ROOT/"scrimp"/"port.py"


def source_port_class():
    tree=ast.parse(PORT.read_text(encoding="utf-8"))
    original=next(
        node for node in tree.body
        if isinstance(node,ast.ClassDef) and node.name=="Port"
    )
    requested={"compute","get_power_sampling","get_power"}
    methods=[
        node for node in original.body if isinstance(node,ast.FunctionDef)
        and node.name in requested
    ]
    assert len(methods)==len(requested)
    module=ast.Module(
        body=[ast.ClassDef(
            name="Port",bases=[],keywords=[],body=methods,decorator_list=[]
        )],type_ignores=[]
    )

    def fake_asm(_generic,_method,_order,_form,_region,model):
        if model.fail:
            raise FloatingPointError("Injected native GetFEM error")
        return float(model.state[0]**2+7.*model.time)

    namespace={
        "np":np,"logging":logging,"time":time,"rank":1,
        "gf":SimpleNamespace(Model=object,asm=fake_asm),
        "Domain":object,
    }
    exec(compile(ast.fix_missing_locations(module),str(PORT),"exec"),
         namespace)
    return namespace["Port"]


class FakeGetFEM:
    def __init__(self):
        self.time=None
        self.state=None
        self.snapshots=[]
        self.fail=False

    def set_time(self,t):
        self.time=float(t)

    def to_variables(self,z):
        self.state=np.array(z,copy=True)
        self.snapshots.append((self.time,self.state.copy()))


class GenericPortSamplingTests(unittest.TestCase):
    def setUp(self):
        cls=source_port_class()
        for name,value in (
            ("get_algebraic",True),("get_substituted",False),
            ("get_region",None),("get_kind","scalar-field"),
            ("get_flow","flow"),("get_effort","effort"),
            ("get_mesh_id",0),("get_name","synthetic-port"),
        ):
            setattr(cls,name,lambda self,value=value:value)
        cls.set_is_computed=lambda self:setattr(self,"_Port__is_computed",True)
        self.port=cls()
        self.port._Port__is_computed=False
        self.port._Port__power_sampling="unknown"
        self.port._Port__power=[]
        self.model=FakeGetFEM()
        self.domain=SimpleNamespace(_int_method={0:"actual-method"})
        self.solution={
            "t":[0.,.1,.3],
            "z":[np.array([1.]),np.array([3.]),np.array([9.])],
        }

    def test_default_generic_instantaneous_states_and_times(self):
        self.port.compute(self.solution,self.model,self.domain)
        self.assertEqual(self.port.get_power_sampling(),"instantaneous")
        np.testing.assert_allclose(self.port.get_power(),[1.,9.7,83.1])
        self.assertEqual(
            [x[0] for x in self.model.snapshots],[0.,.1,.3]
        )
        np.testing.assert_allclose(
            [x[1][0] for x in self.model.snapshots],[1.,3.,9.]
        )

    def test_explicit_cn_step_preserves_previous_native_state_convention(self):
        self.port.compute(
            self.solution,self.model,self.domain,sampling="cn_step"
        )
        self.assertEqual(self.port.get_power_sampling(),"cn_step")
        np.testing.assert_allclose(self.port.get_power(),[1.,4.7,38.1])
        np.testing.assert_allclose(
            [x[1][0] for x in self.model.snapshots],[1.,2.,6.]
        )
        self.assertEqual(
            [x[0] for x in self.model.snapshots],[0.,.1,.3]
        )
        # The time is t[k], not midpoint physical time.
        self.assertAlmostEqual(self.port.get_power()[1],4.+7.*.1)

    def test_cannot_apply_crank_nicolson_by_accident(self):
        self.port.compute(
            self.solution,self.model,self.domain,sampling="instantaneous"
        )
        nodal=np.array(self.port.get_power())
        self.port.compute(
            self.solution,self.model,self.domain,sampling="cn_step"
        )
        stage=np.array(self.port.get_power())
        self.assertGreater(np.max(abs(stage-nodal)),40.)

    def test_reject_unrecognized_sampling_without_changing_old_result(self):
        self.port.compute(self.solution,self.model,self.domain)
        saved=self.port.get_power()
        with self.assertRaisesRegex(ValueError,"Unsupported port-power sampling"):
            self.port.compute(
                self.solution,self.model,self.domain,sampling="bdf_step"
            )
        self.assertEqual(self.port.get_power_sampling(),"instantaneous")
        self.assertEqual(saved,self.port.get_power())

    def test_invalid_time_history_fails(self):
        for t in ([0.,.1,.1],[0.,.2,.1],[0.,float("nan"),.3],
                  [0.,.1],[]):
            with self.subTest(times=t),self.assertRaises(ValueError):
                self.port.compute(
                    {"t":t,"z":self.solution["z"]},
                    self.model,self.domain
                )

    def test_assembly_failure_does_not_leave_old_power_marked_computed(self):
        self.port.compute(self.solution,self.model,self.domain)
        self.model.fail=True
        with self.assertRaisesRegex(FloatingPointError,"Injected native"):
            self.port.compute(self.solution,self.model,self.domain)
        self.assertEqual(self.port.get_power_sampling(),"unknown")
        with self.assertRaises(AssertionError):
            self.port.get_power()

    def test_original_dphs_dispatch_uses_sampling_not_embedded_cn(self):
        source=(ROOT/"scrimp"/"dphs.py").read_text(encoding="utf-8")
        self.assertIn(
            'sampling = "cn_step" if ts_type == "cn" else "instantaneous"',
            source
        )
        self.assertIn(
            "port.compute(\n                self.solution, self.gf_model, self.domain, sampling=sampling",
            source
        )
        self.assertIn('sampling = port.get_power_sampling()',source)
        self.assertIn('"instantaneous_trapezoid"',source)
        self.assertNotIn(
            'convention = "cn_step" if scheme == "cn"',
            source
        )


if __name__=="__main__":
    unittest.main()
