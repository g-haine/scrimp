# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Non-native contract: canonical examples directly create coherent outputs.

No claim of GetFEM physical acceptance is made by these source tests.
The nine-case native acceptance is the separate user-output smoke.
"""
import ast
from pathlib import Path
import unittest

ROOT=Path(__file__).resolve().parents[1]
MODELS={
    "mass_spring_damper.py":"mass_spring_damper",
    "nonlinear_pendulum.py":"nonlinear_pendulum_",
    "heat.py":"heat",
    "heat_wave.py":"heat_wave",
    "shallow_water.py":"shallow_water",
    "wave.py":"wave_2d",
    "wave_dissipative.py":"wave_dissipative",
    "wave_coenergy.py":"wave_coenergy",
}


class UserExampleOutputLayoutTests(unittest.TestCase):
    def test_all_eight_canonical_public_examples_set_one_system_root(self):
        for name,system in MODELS.items():
            with self.subTest(script=name):
                source=(ROOT/"examples"/name).read_text(encoding="utf-8")
                ast.parse(source,filename=name)
                self.assertIn("set_output_name(",source)
                self.assertIn(system,source)
                self.assertIn("plot_Hamiltonian(",source)
                self.assertNotIn('Path.cwd()/"outputs"/"png"',source)

    def test_every_canonical_domain_is_scoped_to_its_output_root(self):
        # The package now imports without eager global outputs/dirs.
        # A single forgotten S.Domain(...) creates a legacy
        # outputs/mesh sibling and breaks the user's file layout.
        for filename in MODELS:
            src=(ROOT/"examples"/filename).read_text(encoding="utf-8")
            tree=ast.parse(src,filename=filename)
            for node in ast.walk(tree):
                if not isinstance(node,ast.Call):
                    continue
                fn=node.func
                if not (isinstance(fn,ast.Attribute)
                        and isinstance(fn.value,ast.Name)
                        and fn.value.id=="S" and fn.attr=="Domain"):
                    continue
                names={keyword.arg for keyword in node.keywords}
                with self.subTest(example=filename,line=node.lineno):
                    self.assertIn(
                        "mesh_output_dir",names,
                        f"{filename}:Domain at line {node.lineno} "
                        "would generate a global outputs/mesh sibling",
                    )

    def test_native_smoke_rejects_any_extra_global_output_sibling(self):
        src=(ROOT/"scripts/run-v130-example-user-output-smoke.sh").read_text()
        self.assertIn('output_parent=work/"outputs"',src)
        self.assertIn('siblings!=[system]',src)
        self.assertIn('unexpected global output siblings',src)
        self.assertIn('len(msh)!=1',src)
        self.assertIn('"no_generic_legacy_output_siblings":True',src)

    def test_heat_wave_user_can_see_real_subsystem_transfer(self):
        src=(ROOT/"examples/heat_wave.py").read_text(encoding="utf-8")
        for needle in (
            "def render_heat_wave_energy_ports(",
            'get_assembly_region()',
            'system.get_quantity(expr,region=region,CN=consecutive)',
            "select_port_power_sampling(",
            'consecutive=(sampling=="cn_step")',
            'np.diff(total)+dt*(d_heat[1:]-p_physical[1:])',
            'H_heat=heat,H_wave=wave,H_total=total',
            'heat_wave_subsystem_energies.png',
            'heat_wave_interface_powers.png',
            'heat_wave_energy_work.png',
            'heat_wave_observables',
            'hw.export_to_pv("T")',
            'hw.export_to_pv("p")',
            'hw.heat_wave_energy_port_visuals = render_heat_wave_energy_ports(',
        ):
            with self.subTest(needle=needle):
                self.assertIn(needle,src)
        self.assertNotIn("from sandbox",src)
        self.assertNotIn("if not physics_only: pass",src)

    def test_normal_examples_do_not_require_sandbox_figures(self):
        src=(ROOT/"scripts/run-v130-example-user-output-smoke.sh").read_text()
        for script in (
            "examples.mass_spring_damper",
            "examples.nonlinear_pendulum",
            "examples.heat",
            "examples.heat_wave",
            "examples.shallow_water",
            "examples.wave",
            "examples.wave_dissipative",
            "examples.wave_coenergy",
        ):
            self.assertIn(script,src)
        self.assertIn("run_case heat_wave_mpi2 heat_wave 2",src)
        self.assertIn('root/"mesh"',src)
        self.assertIn('root/"png"',src)
        self.assertIn('root/"svg"',src)
        self.assertIn('root/"npz"',src)
        self.assertIn('ROOT',src.upper())  # User workdir persists, not image
        self.assertIn("scripts/check-v130-published-immutable.py",src)
        self.assertIn("PASS_V130_ALL_NINE_REAL_USER_EXAMPLE_OUTPUT_CASES",src)

    def test_renderer_manifest_paths_are_directories_not_dicts(self):
        # Regression for the former heat.py metadata-shadowing bug:
        # render_heat_visuals assigned metadata=Path(.../json), then
        # metadata={...}, so metadata/"manifest.json" raised TypeError
        # AFTER a successful genuine GetFEM run and figure export.
        heat=(ROOT/"examples/heat.py").read_text(encoding="utf-8")
        self.assertIn('metadata_dir=folder/"json"',heat)
        self.assertIn('(metadata_dir/"manifest.json").write_text(',heat)
        self.assertIn('metadata={',heat)
        for name in MODELS:
            source=(ROOT/"examples"/name).read_text(encoding="utf-8")
            tree=ast.parse(source,filename=name)
            for func in (node for node in ast.walk(tree)
                         if isinstance(node,(ast.FunctionDef,ast.AsyncFunctionDef))):
                # Check each function's local dictionaries: a dict must
                # never be re-used as the left side of pathlib '/'.
                dict_names=set()
                for node in ast.walk(func):
                    if isinstance(node,ast.Assign) and isinstance(node.value,ast.Dict):
                        dict_names.update(
                            target.id for target in node.targets
                            if isinstance(target,ast.Name)
                        )
                for node in ast.walk(func):
                    if (isinstance(node,ast.BinOp)
                            and isinstance(node.op,ast.Div)
                            and isinstance(node.left,ast.Name)):
                        self.assertNotIn(
                            node.left.id,dict_names,
                            f"{name}:{func.name}: {node.left.id} is a dict, not an output folder",
                        )

    def test_published_examples_not_touched_by_output_migration(self):
        src=(ROOT/"scripts/run-v130-example-user-output-smoke.sh").read_text()
        self.assertNotIn("examples/published:",src)
        self.assertIn("scripts/check-v130-published-immutable.py",src)


if __name__=="__main__":
    unittest.main()
