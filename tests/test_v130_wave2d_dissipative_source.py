# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# See LICENSE for license information.

"""Source-only physics contract for the 2A -> 2B dissipative wave transition.

This test intentionally does not claim real GetFEM/PETSc qualification.
It protects the additive five-brick resistive port from silently
changing the anisotropic conservative membrane, signs or FE spaces.
"""

import ast
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]
A = ROOT/"examples"/"wave.py"
B = ROOT/"examples"/"wave_dissipative.py"
PILOT = ROOT/"scripts"/"run-v130-wave2d-dissipative-pilot.sh"


def wave_function(path):
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    functions = [n for n in tree.body
                 if isinstance(n, ast.FunctionDef) and n.name=="wave_eq"]
    if len(functions) != 1:
        raise AssertionError("One public wave_eq() function required")
    return functions[0]


def weak_bricks(node):
    found = {}
    for call in ast.walk(node):
        if not (isinstance(call, ast.Call) and
                isinstance(call.func, ast.Attribute) and
                call.func.attr=="Brick"):
            continue
        if len(call.args) < 3:
            raise AssertionError("Expected explicit name, GWFL and regions")
        name = ast.literal_eval(call.args[0])
        found[name] = (
            ast.literal_eval(call.args[1]),
            tuple(ast.literal_eval(call.args[2])),
            tuple(sorted((kw.arg, ast.literal_eval(kw.value))
                         for kw in call.keywords)),
        )
    return found


class Wave2BDissipativeSourceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.source2a = A.read_text(encoding="utf-8")
        cls.source2b = B.read_text(encoding="utf-8")
        cls.func2a = wave_function(A)
        cls.func2b = wave_function(B)
        cls.base_bricks = weak_bricks(cls.func2a)
        cls.damped_bricks = weak_bricks(cls.func2b)

    def test_same_original_2a_dirac_and_constitutive_bricks(self):
        self.assertEqual(len(self.base_bricks), 20)
        self.assertEqual(len(self.damped_bricks), 25)
        for name, term in self.base_bricks.items():
            with self.subTest(brick=name):
                self.assertEqual(self.damped_bricks[name], term)

    def test_exactly_five_extra_resistive_bricks(self):
        added = {name:self.damped_bricks[name]
                 for name in self.damped_bricks.keys()-self.base_bricks.keys()}
        self.assertEqual(set(added),
                         {"M_r", "I_r", "-I_r^T", "-M_e_r", "CR_r"})
        self.assertEqual(added["M_r"][0], "f_r*Test_f_r")
        self.assertEqual(added["I_r"][0], "e_r*Test_p")
        self.assertEqual(added["-I_r^T"][0], "-e_p*Test_f_r")
        self.assertEqual(added["-M_e_r"][0], "-e_r*Test_e_r")
        self.assertEqual(added["CR_r"][0], "nu*f_r*Test_e_r")

    def test_same_original_anisotropic_material_and_four_boundary_regions(self):
        for key in (
            '("bottom", 10, "B")',
            '("right", 11, "R")',
            '("top", 12, "T")',
            '("left", 13, "L")',
            '"L": 2., "l": 1., "h": mesh_size',
            '"[[5+x,x*y],[x*y,2+y]]"',
            '"3-x", "p"',
            'position="flow" if is_left else "effort"',
            'S.FEM("q", 1, "DG")',
            'S.FEM("p", 2, "CG")',
            '"0.2*sin(4*pi*t)*sin(pi*y)"',
        ):
            with self.subTest(physics=key):
                self.assertIn(key, self.source2a)
                self.assertIn(key, self.source2b)

    def test_resistive_port_nonnegative_viscosity_and_matching_cg2(self):
        self.assertIn('"Damping", "f_r", "e_r", "scalar-field"', self.source2b)
        self.assertIn('dissipative=True, region=1', self.source2b)
        self.assertIn('S.FEM("Damping", 2, "CG")', self.source2b)
        self.assertIn('f"{damping_strength}*0.5*(2.0-x)"', self.source2b)
        self.assertIn('damping_strength >= 0.', self.source2b)
        self.assertIn('np.isfinite(damping_strength)', self.source2b)
        self.assertNotIn('init_step_ts_type="pseudo"', self.source2b)

    def test_initial_compatibility_and_same_fixed_CN(self):
        for term in (
            'system.set_initial_value("q", "[0., 0.]")',
            'system.set_initial_value("p", "0.")',
            'init_step_ts_type="algebraic", init_step=True',
            'ts_type="cn", ts_adapt_type="none"',
            'ts_exact_final_time="matchstep"',
            'ksp_type="gmres", pc_type="lu"',
            'dt=time_step, dt_save=time_step',
            '"0.5*q.T.q"',
            '"0.5*p*p/rho"',
        ):
            with self.subTest(model=term):
                self.assertIn(term, self.source2a)
                self.assertIn(term, self.source2b)

    def test_independent_original_GetFEM_volume_dissipation(self):
        for term in (
            'system.ports["Damping"]',
            'resistive.get_power_sampling() != "cn_step"',
            '"nu*f_r*f_r"',
            '"(f_r+e_p)*(f_r+e_p)"',
            '"(e_r-nu*f_r)*(e_r-nu*f_r)"',
            'system.gf_model.to_variables(mean)',
            '0.5*(z[k-1]+z[k]) if k else z[k]',
            'np.diff(t) * dissipated_power[1:]',
            'step_deltaH - step_work + step_loss',
            'np.cumsum(step_loss)',
            'H-H[0]-W+L',
            'independent_volume_quadrature_matches_resistive_port',
        ):
            with self.subTest(audit=term):
                self.assertIn(term,self.source2b)
        self.assertNotIn(
            '"resistive_effort_constitutive_law": bool(', self.source2b
        )

    def test_raw_native_quantities_never_invented(self):
        # Match the actual returned observations schema, not the
        # incidental source spelling (dict(H=H) vs {"H": H}).
        tree = ast.parse(self.source2b)
        functions = [
            node for node in tree.body
            if isinstance(node, ast.FunctionDef)
            and node.name == "measure_wave_2d"
        ]
        self.assertEqual(len(functions), 1)
        returns = [
            node.value for node in ast.walk(functions[0])
            if isinstance(node, ast.Return)
            and isinstance(node.value, ast.Call)
            and isinstance(node.value.func, ast.Name)
            and node.value.func.id == "dict"
        ]
        self.assertEqual(len(returns), 1)
        self.assertEqual(returns[0].args, [])
        self.assertTrue(all(k.arg is not None for k in returns[0].keywords))
        returned_keys = [k.arg for k in returns[0].keywords]
        self.assertEqual(len(returned_keys), len(set(returned_keys)))
        for key in (
            "t", "z", "H", "signed_work",
            "cumulative_loss", "dissipated_power",
            "independently_assembled_volume_loss",
            "resistive_flow_constraint",
            "resistive_effort_constraint",
            "energy_step_defects",
        ):
            with self.subTest(quantity=key):
                self.assertIn(key, returned_keys)
        self.assertIn("state_SHA256", self.source2b)
        self.assertIn('for name in ("q", "p", "e_q", "e_p", "f_r", "e_r")',
                      self.source2b)

    def test_2b_sphinx_is_literal_native_source_and_accurate_physics(self):
        document=(ROOT/"docs"/"source"/"examples"/"wave.rst").read_text(
            encoding="utf-8"
        ).split("2B — Adding a resistive port: heterogeneous damping")[1]
        self.assertEqual(document.count(
            ".. literalinclude:: ../../../examples/wave_dissipative.py"
        ),8)
        for k in range(1,9):
            start=f":start-after: # TUTORIAL 2B-{k:02d}:"
            end=f":end-before: # TUTORIAL 2B-{k+1:02d}:"
            with self.subTest(marker=k):
                self.assertIn(start,document)
                self.assertIn(end,document)
        for term in ("f_r=-e_p", "e_r=\\nu f_r",
                     "nu*f_r*f_r", "H", "four",
                     "Port.compute()", "weak projection"):
            with self.subTest(term=term):
                self.assertIn(term,document)
        self.assertNotIn('init_step_ts_type="pseudo"',document)

    def test_pilot_is_only_original_mumps_model_and_immutable_published(self):
        script=PILOT.read_text(encoding="utf-8")
        self.assertIn('for ranks in 1 2; do',script)
        self.assertIn('--solver mumps',script)
        self.assertIn('--damping-strength 1',script)
        self.assertIn('python tests/test_v130_wave2d_dissipative_source.py',
                      script)
        self.assertIn('scripts/check-v130-published-immutable.py',script)
        self.assertNotIn('--solver superlu_dist',script)
        self.assertNotIn('--solver superlu ',script)
        self.assertIn('dissipated_power',script)
        self.assertIn('independently_assembled_volume_loss',script)


if __name__ == "__main__":
    unittest.main()
