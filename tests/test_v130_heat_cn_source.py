# SCRIMP v1.3.0 — source-only contract for upgraded canonical heat example.
# A native PETSc/GetFEM physical PASS requires the separate Docker pilot.
import ast
from pathlib import Path
import unittest

ROOT=Path(__file__).resolve().parents[1]
HEAT=ROOT/"examples/heat.py"
BASELINE=ROOT/"sandbox/probe_v130_heat_physical.py"
BASELINE_SCRIPT=ROOT/"scripts/run-v130-heat-example-baseline.sh"


class HeatCNSourceContract(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.source=HEAT.read_text(encoding="utf-8")
        cls.tree=ast.parse(cls.source,filename=str(HEAT))

    def test_public_heat_function_defaults_to_cn(self):
        fn=next(node for node in self.tree.body
                if isinstance(node,ast.FunctionDef) and node.name=="heat_eq")
        defaults={
            a.arg: d.value if isinstance(d,ast.Constant) else None
            for a,d in zip(fn.args.kwonlyargs,fn.args.kw_defaults)
        }
        self.assertEqual(defaults["time_stepper"],"cn")
        self.assertIsNone(defaults["initializer"])
        self.assertIsNone(defaults["init_profile"])
        self.assertIn('init_profile = "strict" if time_stepper == "cn" else "default"',
                      self.source)
        self.assertFalse(defaults["showcase"])
        self.assertIn('initializer = "algebraic" if time_stepper == "cn" else "pseudo"',
                      self.source)
        self.assertIn('scheme["ts_bdf_order"] = 2',self.source)
        self.assertIn('init_step_ts_type=initializer',self.source)
        self.assertIn('ts_adapt_type="none"',self.source)

    def test_original_weak_forms_controls_and_initial_data_retained(self):
        for fragment in (
            'S.FEM(states[0].get_name(), 1, FEM="DG")',
            'S.FEM(ports[0].get_name(), 2, FEM="CG")',
            'S.Brick("M_T", "T*rho*Test_T", [1], dt=True',
            'S.Brick("D", "-Div(J_Q)*Test_T", [1]',
            'S.Brick("CR_Q", "f_Q.Lambda.Test_J_Q", [1])',
            'expressions = ["1.", "1.", "1.", "0.2*T"]',
            'heat.set_initial_value("T", "1. + 2.*np.exp(-50*',
            'S.Term("L^2-norm", "0.5*T*rho*T", [1])',
        ):
            with self.subTest(fragment=fragment):
                self.assertIn(fragment,self.source)

    def test_actual_native_cn_power_controls_and_visuals(self):
        for fragment in (
            'system.get_quantity(',
            '"f_Q.Lambda.f_Q", region=1, CN=True',
            '"Boundary control (bottom)"',
            '"Boundary control (right)"',
            '"Boundary control (top)"',
            '"Boundary control (left)"',
            'port.get_power_sampling() != "cn_step"',
            'heat.export_to_pv("T", t="All")',
            'render_heat_visuals(',
            'heat_temperature_snapshots.png',
            'heat_temperature_evolution.gif',
            'heat_energy_and_ports.png',
            'no_spatial_or_time_interpolation',
            'parser.set_defaults(showcase=True)',
        ):
            with self.subTest(fragment=fragment):
                self.assertIn(fragment,self.source)

    def test_strict_init_options_are_applied_inside_heat_scheme(self):
        # Setting PETSc.Options() before DPHS() is ineffective because
        # DPHS.__init__ resets that global PETSc option database.
        for fragment in (
            'if init_profile == "strict":',
            'init_algebraic_max_passes=8',
            'init_algebraic_accuracy=1e-13',
            'init_algebraic_ksp_rtol=1e-14',
            'init_algebraic_ksp_atol=1e-14',
            'heat.set_time_scheme(**scheme)',
            '"--init-profile",choices=("default","strict")',
            '"strict" if args.time_stepper == "cn" else "default"',
            'init_profile=effective_init_profile',
            '"init_profile_applied"',
        ):
            with self.subTest(fragment=fragment):
                self.assertIn(fragment,self.source)
        self.assertLess(
            self.source.index('init_algebraic_max_passes=8'),
            self.source.index('heat.set_time_scheme(**scheme)'),
        )
        pilot=(ROOT/"scripts/run-v130-heat-cn-pilot.sh").read_text(
            encoding="utf-8")
        self.assertIn('--init-profile "$profile"',pilot)
        self.assertIn("Strict PETSc initializer options NOT APPLIED",pilot)
        self.assertNotIn("probe_v130_heat_cn_init_precision.py",pilot)

    def test_historical_bdf2_pseudo_reproducibility_is_explicit(self):
        self.assertIn(
            'heat_eq(time_stepper="bdf", initializer="pseudo")',
            BASELINE.read_text(encoding="utf-8")
        )
        runner=BASELINE_SCRIPT.read_text(encoding="utf-8")
        self.assertIn("--time-stepper bdf --initializer pseudo --no-showcase",
                      runner)


if __name__=="__main__":
    unittest.main()
