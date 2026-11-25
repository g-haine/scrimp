# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
"""Static contract for genuine 2C coenergy and the visualization pipeline.

Native GetFEM/PETSc scientific PASS is ONLY reported by the pilot,
never inferred from these tests.
"""
import ast
from pathlib import Path
import re
import unittest


ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/"examples"/"wave_coenergy.py"
DOC=ROOT/"docs"/"source"/"examples"/"wave_coenergy.rst"
RUNNER=ROOT/"scripts"/"run-v130-wave2d-coenergy-pilot.sh"
RENDERER=ROOT/"scripts"/"render-v130-wave2d-native.py"


class WaveCoenergySourceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.source=SOURCE.read_text(encoding="utf-8")
        cls.doc=DOC.read_text(encoding="utf-8")
        cls.runner=RUNNER.read_text(encoding="utf-8")
        cls.renderer=RENDERER.read_text(encoding="utf-8")
        cls.ast=ast.parse(cls.source)
        cls.ast_renderer=ast.parse(cls.renderer)

    def test_true_compliance_is_algebraic_inverse_of_2a_tensor(self):
        variables={}
        for node in self.ast.body:
            if isinstance(node,ast.Assign) and len(node.targets)==1:
                target=node.targets[0]
                if isinstance(target,ast.Name) and target.id in ("DET_T","TINV"):
                    expr=ast.Expression(body=node.value)
                    variables[target.id]=eval(
                        compile(expr,str(SOURCE),"eval"),{"__builtins__":{}},variables
                    )
        self.assertEqual(set(variables),{"DET_T","TINV"})
        expr=variables["TINV"]
        for x,y in ((0.,0.),(.4,.2),(-.8,.3),(.7,-.6)):
            inv=eval(expr,{"__builtins__":{}},{"x":x,"y":y})
            t=((5+x,x*y),(x*y,2+y))
            product=[[sum(t[i][k]*inv[k][j] for k in range(2))
                      for j in range(2)] for i in range(2)]
            for i in range(2):
                for j in range(2):
                    self.assertAlmostEqual(product[i][j],float(i==j),places=13)
        self.assertIn('S.Brick("M_q", "q.Tinv.Test_q"',self.source)
        self.assertIn('S.Term("Potential energy", "0.5*q.Tinv.q"',self.source)

    def test_two_region_coenergy_and_strictly_positive_resistive_inverse(self):
        for snippet in (
            'S.Domain("Concentric"',
            'S.CoState("e_q", "Stress", states[0], substituted=True)',
            'S.CoState("e_p", "Velocity", states[1], substituted=True)',
            'S.Port("Damping", "e_r", "e_r", "scalar-field", substituted=True, region=1)',
            'S.Brick("M_r", "e_r/nu*Test_e_r", [1], position="flow")',
            'S.Brick("I_r", "e_r*Test_p", [1], position="effort")',
            'S.Brick("-I_r^T", "-p*Test_e_r", [1], position="effort")',
            '"0.2+10*(0.36-(x*x+y*y))"',
        ):
            with self.subTest(snippet=snippet):
                self.assertIn(snippet,self.source)
        for x,y in ((0,0),(.59,0),(.0,.59)):
            self.assertGreater(0.2+10*(.36-x*x-y*y),0)
        self.assertNotIn('init_step_ts_type="pseudo"',self.source)

    def test_closed_loop_feedback_is_passive_and_both_losses_are_measured(self):
        self.assertIn('expressions = [f"{FEEDBACK_GAIN}*Y"]',self.source)
        self.assertIn('step_boundary = dt*P[1:]',self.source)
        self.assertIn('np.max(np.abs(feedback-P))',self.source)
        self.assertIn('"positive_nontrivial_feedback_loss"',self.source)
        self.assertIn('f"{FEEDBACK_GAIN}*Y*Y"',self.source)
        self.assertIn('"e_r*e_r/nu"',self.source)
        self.assertIn('"-p*e_r"',self.source)
        self.assertIn('"nu*p*p"',self.source)
        self.assertIn('"cumulative_CN_dissipative_energy_balance"',self.source)
        self.assertIn('"independent_weak_internal_loss_matches"',self.source)
        self.assertIn('"independent_feedback_quadrature_matches"',self.source)
        self.assertIn("resistive_weak_residual_L2",self.source)
        self.assertNotIn("resistive.get_power()",self.source)

    def test_port_orientation_is_original_grad_grad_not_mixed_2a(self):
        # In the Grad-Grad coenergy assembly, the flow brick is negated
        # by DPHS.add_brick; together with -B^T this implies Y=-p.
        # Therefore H'=-integral(U*Y)-D_r. U=+kY is passive.
        self.assertIn(
            'S.Brick("M_Y", "Y*Test_Y", [20], position="flow")',
            self.source,
        )
        self.assertIn(
            'S.Brick("-B^T", "-p*Test_Y", [20], position="effort")',
            self.source,
        )
        self.assertIn('expressions = [f"{FEEDBACK_GAIN}*Y"]',self.source)
        self.assertIn('step_boundary = dt*P[1:]',self.source)
        self.assertIn('"independent_feedback_quadrature_matches"',self.source)
        self.assertNotIn('step_boundary = -dt*P[1:]',self.source)

    def test_exact_dae_and_ten_fixed_cn_steps(self):
        for s in (
            'init_step_ts_type="algebraic", init_step=True',
            'ts_type="cn", ts_adapt_type="none"',
            'ts_exact_final_time="matchstep"',
            'time_step=0.008',
            'final_time=0.08',
            'pc_type="lu"',
            "system.last_ts_rejections==0",
            'wave.set_initial_value("q", "[0., 0.]")',
        ):
            with self.subTest(term=s):
                self.assertIn(s,self.source)
        self.assertNotIn("wave.plot_Hamiltonian(save_figure=True)",self.source)

    def test_genuine_native_fe_vtk_and_archive_preserved(self):
        self.assertIn("np.savez_compressed(",self.source)
        self.assertIn('for name in ("q","p","e_r"):',self.source)
        for term in ("--check --export --output /artifacts",
                     "mumps_mpi1","mumps_mpi2",
                     "rank_traces_agree","energy_step_defects",
                     "summary.json","zipfile.ZipFile",
                     "check-v130-published-immutable.py",
                     "render-v130-wave2d-native.py"):
            with self.subTest(term=term):
                self.assertIn(term,self.runner)

    def test_visuals_are_native_not_synthetic(self):
        for term in ("vtk_binary","vtk_snapshot","field_history",
                     "VTK FE field / mesh size mismatch",
                     "np.load(io.BytesIO(data.read(",
                     "field_snapshots",".svg",
                     "p_evolution.gif","energy_balance",
                     "energy_and_ports.csv","provenance.json",
                     "Refusing unqualified selected native physics case"):
            with self.subTest(term=term):
                self.assertIn(term,self.renderer)

    def test_sphinx_includes_same_canonical_source(self):
        self.assertEqual(self.doc.count(".. literalinclude:: ../../../examples/wave_coenergy.py"),8)
        for i in range(1,9):
            self.assertIn(f"# TUTORIAL 2C-{i:02d}:",self.source)
            self.assertIn(f":start-after: # TUTORIAL 2C-{i:02d}:",self.doc)
        self.assertIn("T^{-1}",self.doc)
        self.assertIn("p_evolution.gif",self.doc)
        self.assertNotIn("Wave_CoEnergy_Hamiltonian.png",self.doc)


if __name__=="__main__":
    unittest.main()
