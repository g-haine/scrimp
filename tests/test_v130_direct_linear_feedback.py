# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Fail-closed linear-feedback classification: actual production methods.

These tests execute the real DPHS.set_control(), add_brick() and
set_linear_flags() methods using a fake GetFEM assembly recorder.
Native distributed GetFEM/PETSc testing remains a separate gate.
"""

import ast
import logging
from pathlib import Path
import re
import unittest


ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/"scrimp"/"dphs.py"


class StubBrick:
    def __init__(self, name, form, regions, linear=True, dt=False,
                 position="constitutive", mesh_id=0, explicit=False):
        self.name=name
        self.form=form
        self.regions=regions
        self.linear=linear
        self.dt=dt
        self.position=position
        self.mesh_id=mesh_id
        self.explicit=explicit
        self.ids=[]

    def get_name(self): return self.name
    def get_form(self): return self.form
    def get_regions(self): return self.regions
    def get_linear(self): return self.linear
    def get_dt(self): return self.dt
    def get_position(self): return self.position
    def get_explicit(self): return self.explicit
    def get_mesh_id(self): return self.mesh_id
    def add_id_brick_to_list(self, identifier): self.ids.append(identifier)


def actual_dpHs_methods():
    astree=ast.parse(SOURCE.read_text(encoding="utf-8"))
    helpers=[
        n for n in astree.body
        if isinstance(n,ast.FunctionDef) and n.name in {
            "_control_expression_is_direct_linear_feedback",
            "_control_expression_is_prescribed",
            "_all_bricks_are_linear",
        }
    ]
    cls=next(n for n in astree.body
             if isinstance(n,ast.ClassDef) and n.name=="DPHS")
    methods=[
        n for n in cls.body if isinstance(n,ast.FunctionDef)
        and n.name in {"set_control","add_brick","set_linear_flags"}
    ]
    assert len(helpers)==3 and len(methods)==3
    env={
        "re":re, "logging":logging, "rank":1,
        "Brick":StubBrick,
        "_PRESCRIBED_CONTROL_SYMBOLS":frozenset({
            "t","x","y","z","sin","cos","exp",
        }),
    }
    exec(compile(ast.Module(body=helpers+methods,type_ignores=[]),
                 str(SOURCE),"exec"),env)
    return env


class StubControl:
    def __init__(self, name="Interface Heat", variable="U_T"):
        self.name=name
        self.variable=variable
    def get_kind(self): return "scalar-field"
    def get_name_control(self): return self.variable
    def get_region(self): return 10
    def get_mesh_id(self): return 0
    def get_name(self): return self.name
    def get_isSet(self): return True


class StubPort:
    def __init__(self, flow,effort):
        self.flow=flow
        self.effort=effort
    def get_flow(self): return self.flow
    def get_effort(self): return self.effort


class RecordingGetFEM:
    def __init__(self): self.calls=[]
    def add_linear_term(self, mim, form, region):
        self.calls.append(("linear",form,region))
        return len(self.calls)
    def add_nonlinear_term(self, mim, form, region):
        self.calls.append(("nonlinear",form,region))
        return len(self.calls)
    def add_source_term_brick(self, mim, name, form, region):
        self.calls.append(("source",name,form,region))
        return len(self.calls)


class Recorder:
    def __init__(self):
        self.controls={
            "Interface Heat":StubControl(),
            "Interface Wave":StubControl("Interface Wave","U_w"),
        }
        self.ports={
            "heat_feedback":StubPort("Y_T","U_T"),
            "wave_feedback":StubPort("Y_w","U_w"),
        }
        self.bricks={}
        self.gf_model=RecordingGetFEM()
        self.domain=type("Domain",(),{"_int_method":["mock-mim"]})()
        self.linear_mass=True
        self.linear_stiffness=True


class DirectLinearControlTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.methods=actual_dpHs_methods()

    def model(self):
        return Recorder()

    def add_control(self, control="Interface Heat", expression="Y_w"):
        sim=self.model()
        # Execute the real DPHS.add_brick method, rather than only
        # recording its arguments: this checks GetFEM routing.
        sim.add_brick=lambda brick: self.methods["add_brick"](sim,brick)
        self.methods["set_control"](sim,control,expression)
        return sim

    def test_only_bare_known_model_symbols_are_exactly_linear(self):
        helper=self.methods["_control_expression_is_direct_linear_feedback"]
        vars={"Y_w","Y_T","p"}
        for text in ("Y_w","+Y_T","-Y_T"," - Y_w ","p"):
            with self.subTest(text=text):
                self.assertTrue(helper(text,vars))
        for text in (
            "Y_w*Y_w","Y_w*Y_T","Y_w+Y_T",
            "2*Y_w","sin(Y_w)","Y_w(t)","-macro",
            "0.","sin(t)","Y_w+1.","(Y_w)","Y_w/2",
            "Y_T**2","Y_w, Y_T","",
        ):
            with self.subTest(text=text):
                self.assertFalse(helper(text,vars))
        self.assertFalse(helper(None,vars))

    def test_actual_add_brick_routes_bare_feedback_to_linear_tangent(self):
        for control,expr in (
            ("Interface Heat","Y_w"),
            ("Interface Heat"," + Y_w "),
            ("Interface Wave","Y_T"),
            ("Interface Wave","-Y_T"),
        ):
            with self.subTest(control=control,expr=expr):
                model=self.add_control(control,expr)
                name=("U_T" if control=="Interface Heat" else "U_w")
                self.assertEqual(
                    [item[0] for item in model.gf_model.calls],
                    ["linear","linear"],
                )
                self.assertEqual(
                    model.gf_model.calls[-1][1],
                    f"({expr})*Test_{name}",
                )
                self.assertEqual(model.gf_model.calls[-1][2],10)
                self.assertTrue(
                    model.bricks[name+"_feedback"].get_linear()
                )
                self.methods["set_linear_flags"](model)
                self.assertTrue(model.linear_mass)
                self.assertTrue(model.linear_stiffness)
                self.assertTrue(self.methods["_all_bricks_are_linear"](model.bricks))

    def test_nonlinear_feedback_keeps_real_getfem_newton_tangent(self):
        for expr in ("Y_w*Y_w","sin(Y_w)","Y_w*Y_T","opaque_macro"):
            with self.subTest(expr=expr):
                model=self.add_control(expression=expr)
                self.assertEqual(
                    [item[0] for item in model.gf_model.calls],
                    ["linear","nonlinear"],
                )
                self.assertFalse(model.bricks["U_T_feedback"].get_linear())
                self.methods["set_linear_flags"](model)
                self.assertFalse(model.linear_stiffness)
                self.assertFalse(self.methods["_all_bricks_are_linear"](model.bricks))

    def test_open_loop_source_remains_source_not_linear_feedback(self):
        model=self.add_control(expression="0.")
        self.assertEqual(
            [item[0] for item in model.gf_model.calls],
            ["linear","source"],
        )
        self.methods["set_linear_flags"](model)
        self.assertTrue(model.linear_stiffness)
        self.assertTrue(self.methods["_all_bricks_are_linear"](model.bricks))
        self.assertTrue(
            getattr(model.bricks["U_T_source"],
                    "_scrimp_prescribed_control_source",False)
        )

    def test_no_rewrite_of_original_fe_form_or_published_code(self):
        source=SOURCE.read_text(encoding="utf-8")
        self.assertIn(
            'f"({expression}){times}Test_{u}"',source
        )
        self.assertIn(
            'linear=direct_linear_feedback',source
        )
        self.assertNotIn('position="source",\n'
                         '                    linear=direct_linear_feedback',source)


if __name__=="__main__":
    unittest.main()
