# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Source-only preflight of the strict native Heat–Wave physical pilot.

This does NOT qualify GetFEM/PETSc or the physical model: it prevents
a failed Python heredoc or a changed scientific gate wasting MPI time.
"""
import ast
from pathlib import Path
import unittest

ROOT=Path(__file__).resolve().parents[1]
RUNNER=ROOT/"scripts/run-v130-heat-wave-physical-pilot.sh"
EXAMPLE=ROOT/"examples/heat_wave.py"
DIAG=ROOT/"sandbox/measure_v130_heat_wave_physics.py"
PROBE=ROOT/"sandbox/probe_v130_main_2d_examples.py"


class HeatWavePilotSourceTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.runner=RUNNER.read_text(encoding="utf-8")
        cls.model=EXAMPLE.read_text(encoding="utf-8")
        cls.diagnostic=DIAG.read_text(encoding="utf-8")
        cls.probe=PROBE.read_text(encoding="utf-8")

    def test_embedded_python_compiles(self):
        marker = "python3 - \"$out\" <<'PY'\n"
        self.assertEqual(self.runner.count(marker),1)
        verifier=self.runner.split(marker,1)[1].split("\nPY\n",1)[0]
        ast.parse(verifier, filename=str(RUNNER))
        self.assertIn('raise SystemExit(0 if report["status"].startswith("PASS_") else 1)',
                      verifier)

    def test_true_physics_gates_and_explicit_nonqualification(self):
        for token in (
            "FAIL_OR_UNQUALIFIED_HEAT_WAVE_PHYSICAL",
            "PASS_HEAT_WAVE_PHYSICAL_MUMPS_MPI1",
            "CN_heat_balance",
            "CN_wave_balance",
            "CN_total_independent_balance",
            "oriented_interface_power_cancels",
            "outer_homogeneous_boundary_power_zero",
            "nonnegative_nontrivial_heat_dissipation",
            "with zipfile.ZipFile",
            '"native_exit_code":rc',
            "MPI2_qualified",
            "full_horizon_qualified",
        ):
            self.assertIn(token,self.runner)
        self.assertIn("CN_per_step_heat_balance_error",self.diagnostic)
        self.assertIn("CN_per_step_wave_balance_error",self.diagnostic)
        self.assertIn("CN_cumulative_total_balance_error",self.diagnostic)
        self.assertIn("collect_heat_wave_physics(",self.probe)

    def test_cn_physical_balance_is_the_actual_executable_expression(self):
        """Check the three-step CN identity by AST, not text formatting."""
        marker = "python3 - \"$out\" <<'PY'\n"
        verifier = self.runner.split(marker, 1)[1].split("\nPY\n", 1)[0]
        tree = ast.parse(verifier)
        values = {}
        for node in ast.walk(tree):
            if isinstance(node, ast.Assign):
                for target in node.targets:
                    if isinstance(target, ast.Name) and target.id in ("dH", "dt", "res"):
                        values[target.id] = node.value
        self.assertEqual(set(values), {"dH", "dt", "res"})
        for name, expression in (
            ("dH", "np.diff(H)"),
            ("dt", "np.diff(t)"),
            ("res", "dH + dt*(d-interface-exterior)"),
        ):
            with self.subTest(name=name):
                expected = ast.parse(expression, mode="eval").body
                self.assertEqual(
                    ast.dump(values[name], include_attributes=False),
                    ast.dump(expected, include_attributes=False),
                )

    def test_tutorial_links_qualified_interface_assembly_policy(self):
        tutorial = ROOT/"docs/source/examples/heat_wave.rst"
        decision = ROOT/"docs/V130_INTERFACE_PV_PRODUCTION_DECISION_20261124.md"
        self.assertTrue(decision.is_file())
        text = tutorial.read_text(encoding="utf-8")
        self.assertIn("Automatic assembly of conforming interfaces", text)
        self.assertIn("V130_INTERFACE_PV_PRODUCTION_DECISION_20261124.md", text)
        self.assertIn("correct subdomain side", text)
        self.assertIn("U_T = Y_w", text)
        self.assertIn("U_w = Y_T", text)

    def test_real_shipped_model_and_readonly_native_mount(self):
        for token in (
            'S.FEM("Heat flux", k, "RT")',
            'S.FEM("q", k, "RT")',
            'hw.set_control("Interface Heat", "Y_w")',
            'hw.set_control("Interface Wave", "Y_T")',
            'init_step_ts_type="algebraic"',
            'ts_type="cn"',
            'ts_adapt_type="none"',
        ):
            self.assertIn(token,self.model)
        for token in (
            'scripts/check-v130-published-immutable.py',
            'scripts/run-v130-heat-wave-physical-pilot.sh',
        ):
            if token == 'scripts/run-v130-heat-wave-physical-pilot.sh':
                self.assertTrue(RUNNER.is_file())
            else:
                self.assertIn(token,self.runner)
        self.assertIn('mpiexec --oversubscribe --bind-to none -n 1',self.runner)
        self.assertIn('python -u /probe-2d.py --example heat_wave --backend mumps',
                      self.runner)
        self.assertIn('examples/heat_wave.py:/workspace/scrimp/examples/heat_wave.py:ro',
                      self.runner)


if __name__=="__main__":
    unittest.main()
