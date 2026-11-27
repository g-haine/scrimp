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

    def test_published_examples_not_touched_by_output_migration(self):
        src=(ROOT/"scripts/run-v130-example-user-output-smoke.sh").read_text()
        self.assertNotIn("examples/published:",src)
        self.assertIn("scripts/check-v130-published-immutable.py",src)


if __name__=="__main__":
    unittest.main()
