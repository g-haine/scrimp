"""Standalone source and numerical unit regressions for public v1.3 APIs.

No GetFEM/PETSc import is required for these deliberately isolated
negative tests. Native MPI physical tests remain separately necessary.
"""
import ast
import importlib.util
import logging
from pathlib import Path
import tempfile
import types
import unittest

import numpy as np


ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"scrimp/utils/example_output.py"
SAMPLING=ROOT/"scrimp/utils/power_sampling.py"
HAMILTONIAN=ROOT/"scrimp/hamiltonian.py"
DPHS=ROOT/"scrimp/dphs.py"


def load_standalone(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def hamiltonian_with_mock_getfem():
    """Compile only the real classes using a tiny, explicit GF mock."""
    src=ast.parse(HAMILTONIAN.read_text(encoding="utf-8"))
    classes=[node for node in src.body
             if isinstance(node,ast.ClassDef)
             and node.name in ("Term","Hamiltonian")]
    model_type=type("Model",(),{})
    namespace={
        "np":np,"logging":logging,"time":__import__("time"),
        "rank":1,
        "gf":types.SimpleNamespace(
            Model=model_type,
            asm=lambda mode,mim,order,expr,region,model:
                float(model.current_time+model.current_state[0])
        ),
        "Domain":type("Domain",(),{}),
    }
    program=ast.Module(body=classes,type_ignores=[])
    exec(compile(program,str(HAMILTONIAN),"exec"),namespace)
    return namespace["Hamiltonian"],namespace["Term"],model_type


class UserOutputPowerHamiltonianTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.paths=load_standalone("scrimp_example_output_test",OUT)
        cls.power=load_standalone("scrimp_power_sampling_test",SAMPLING)
        cls.core=DPHS.read_text(encoding="utf-8")
        ast.parse(cls.core)
        ast.parse(HAMILTONIAN.read_text(encoding="utf-8"))

    def test_normal_output_root_and_typed_files_are_deterministic(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=self.paths.output_root("heat_wave",working_dir=tmp)
            self.assertEqual(root,Path(tmp)/"outputs"/"heat_wave")
            # Resolving alone must not write anything.
            self.assertFalse(root.exists())
            figure=self.paths.output_file(root,"Hamiltonian.png","png")
            self.assertEqual(figure,root/"png"/"Hamiltonian.png")
            self.assertTrue(figure.parent.is_dir())
            self.assertEqual(
                self.paths.output_file(root,"Hamiltonian.svg","svg"),
                root/"svg"/"Hamiltonian.svg"
            )
            self.assertEqual(
                self.paths.output_file(root,"T.pvd","mesh",create=False),
                root/"mesh"/"T.pvd"
            )
            self.assertEqual(
                self.paths.output_root("ignored",output_dir=root),root
            )

    def test_output_traversal_and_bad_kind_rejected(self):
        for name in ("../bad","/bad","x/y","..",""):
            with self.subTest(name=name):
                with self.assertRaises(ValueError):
                    self.paths.output_root(name)
        root=Path("/tmp/scrimp-output-negative-control")
        with self.assertRaises(ValueError):
            self.paths.output_file(root,"../leak.png","png",create=False)
        with self.assertRaises(ValueError):
            self.paths.output_file(root,"data.png","random",create=False)
        with self.assertRaises(ValueError):
            self.paths.output_file(root,"wrong.svg","png",create=False)

    def test_cn_consecutive_true_TS_steps_only(self):
        p=self.power.select_port_power_sampling
        t=np.linspace(0.,.08,9)
        self.assertEqual(
            p("cn",t,accepted_steps=8,dt=.01,adapt_type="none"),
            "cn_step"
        )
        self.assertEqual(
            p("cn",t[::2],accepted_steps=8,dt=.01,adapt_type="none"),
            "instantaneous"
        )
        self.assertEqual(
            p("cn",t,accepted_steps=None,dt=.01,adapt_type="none"),
            "instantaneous"
        )
        self.assertEqual(
            p("bdf",t,accepted_steps=8,dt=.01,adapt_type="none"),
            "instantaneous"
        )
        with self.assertRaisesRegex(ValueError,"consecutive"):
            p("cn",t[::2],accepted_steps=8,dt=.01,
              adapt_type="none",requested="cn_step")

    def test_cn_time_grid_guard_and_nonmonotone(self):
        p=self.power.select_port_power_sampling
        t=np.array([0.,.01,.021,.03])
        self.assertEqual(
            p("cn",t,accepted_steps=3,dt=.01,adapt_type="none"),
            "instantaneous"
        )
        self.assertEqual(
            p("cn",t,accepted_steps=3,adapt_type="basic"),
            "cn_step"
        )
        with self.assertRaises(ValueError):
            p("cn",[0.,.01,.01],accepted_steps=2)

    def test_hamiltonian_native_class_recompute_is_idempotent_and_time_aware(self):
        Hamiltonian,Term,Model=hamiltonian_with_mock_getfem()
        h=Hamiltonian("Energy")
        t=Term("physical time-dependent term","t+u",[1])
        h.add_term(t)
        model=Model()
        model.current_time=-1.
        model.current_state=None
        model.set_time=lambda time:setattr(model,"current_time",time)
        model.to_variables=lambda state:setattr(model,"current_state",state)
        domain=types.SimpleNamespace(_int_method=[object()])
        states={"t":[0.,.1,.2],"z":[[1.],[2.],[3.]]}
        h.compute(states,model,domain)
        self.assertEqual(t.get_values(),[1.,2.1,3.2])
        h.compute(states,model,domain)
        self.assertEqual(t.get_values(),[1.,2.1,3.2])
        # Invalidate the cached result on addition of new term.
        h.add_term(Term("extra","u",[1]))
        self.assertFalse(h.get_is_computed())
        h.compute(states,model,domain)
        self.assertEqual(len(h[1].get_values()),3)
        h.invalidate()
        self.assertFalse(h.get_is_computed())
        h.compute(states,model,domain)
        self.assertEqual(t.get_values(),[1.,2.1,3.2])

    def test_invalid_hamiltonian_history_does_not_corrupt_accepted_values(self):
        Hamiltonian,Term,Model=hamiltonian_with_mock_getfem()
        h=Hamiltonian("Energy")
        term=Term("E","u",[1])
        h.add_term(term)
        model=Model()
        model.set_time=lambda t:setattr(model,"current_time",t)
        model.to_variables=lambda state:setattr(model,"current_state",state)
        domain=types.SimpleNamespace(_int_method=[object()])
        h.compute({"t":[0.,.01],"z":[[1.],[2.]]},model,domain)
        with self.assertRaises(ValueError):
            h.compute({"t":[0.,.01,.01],"z":[[1.],[2.],[3.]]},
                      model,domain)
        self.assertEqual(term.get_values(),[1.,2.])
        self.assertFalse(h.get_is_computed())

    def test_production_source_uses_new_core_api_with_legacy_fallback(self):
        for literal in (
            "def set_output_name(",
            "def get_output_path(",
            "def set_output_directory(",
            "source = Path(outputs_path) / \"mesh\" / (domain.get_name()+\".msh\")",
            "self.hamiltonian.invalidate()",
            "self.powers_computed = False",
            "select_port_power_sampling(",
            "accepted_steps=getattr(self,\"last_ts_steps\",None)",
            "CN powers: stored states are not all consecutive accepted",
            "def plot_Hamiltonian(",
            'self.get_output_path("png",stem+".png",create=True)',
            'self.get_output_path("svg",stem+".svg",create=True)',
            'self.get_output_path("npz",stem+".npz",create=True)',
            "np.savez_compressed(npz,**traces)",
            "self.plot_powers(ax,HamTot=HamTot)",
            'path = os.path.join(outputs_path, "pv", name_variable)',
            'path = os.path.join(self._output_root, "mesh", name_variable)',
        ):
            with self.subTest(literal=literal):
                self.assertIn(literal,self.core)


if __name__=="__main__":
    unittest.main()
