# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Source-only scientific regressions for the first user tutorial.

The actual observe_wave_1d function is extracted from its AST: no
GetFEM/PETSc imports or synthetic claims of native qualification.
The fake time series is ONLY for host logic and negative controls.
"""

import ast
from pathlib import Path
import types
import unittest

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
EXAMPLE = ROOT / "examples" / "wave_1d.py"
DOC = ROOT / "docs" / "source" / "started.rst"
EXAMPLE_SOURCE = EXAMPLE.read_text(encoding="utf-8")
DOC_SOURCE = DOC.read_text(encoding="utf-8")


def original_observer():
    module = ast.parse(EXAMPLE_SOURCE, filename=str(EXAMPLE))
    method = next(
        n for n in module.body
        if isinstance(n, ast.FunctionDef) and n.name == "observe_wave_1d"
    )
    extract = ast.fix_missing_locations(ast.Module(
        body=[method], type_ignores=[],
    ))
    namespace = {"np": np}
    exec(compile(extract, str(EXAMPLE), "exec"), namespace)
    return namespace["observe_wave_1d"]


class Port:
    def __init__(self, powers, sampling="cn_step"):
        self.powers = powers
        self.sampling = sampling

    def get_power(self):
        return self.powers

    def get_power_sampling(self):
        return self.sampling


class System:
    def __init__(self):
        self.solution = {
            "t": np.array([0., 0.004, 0.008]),
            "z": np.zeros((3, 4)),
        }
        self.H = np.array([1., 1.001, 1.003])
        self.ports = {
            "Boundary control (left)": Port([0., .25, .5]),
            "Boundary control (right)": Port([0., 0., 0.]),
        }
        self.powers_computed = True
        self.last_ts_converged_reason = 2
        self.last_ts_steps = 2
        self.last_ts_rejections = 0
        self.last_ts_snes_failures = 0
        self.initialization_report = {
            "status": "PASS", "differential_initial_exact": True,
        }

    def get_Hamiltonian(self):
        return self.H


class FirstTutorialTests(unittest.TestCase):
    def setUp(self):
        self.observe = original_observer()

    def test_python_and_document_syntax_sources_parse(self):
        compile(EXAMPLE_SOURCE, str(EXAMPLE), "exec")
        self.assertIn(".. code-block:: bash", DOC_SOURCE)
        self.assertIn("examples/wave_1d.py", DOC_SOURCE)

    def test_both_outward_neumann_controls_are_positive_weak_bricks(self):
        self.assertIn('"U_L*Test_p"', EXAMPLE_SOURCE)
        self.assertIn('"U_R*Test_p"', EXAMPLE_SOURCE)
        self.assertNotIn('"-U_L*Test_p"', EXAMPLE_SOURCE)

    def test_observations_are_physical_velocity_and_negative_flow_side(self):
        for name in ("L", "R"):
            self.assertIn(f'"-Y_{name}*Test_Y_{name}"', EXAMPLE_SOURCE)
            self.assertIn(f'"-e_p*Test_Y_{name}"', EXAMPLE_SOURCE)
        self.assertIn('"-Y_L * Test_Y_L"', DOC_SOURCE)
        self.assertIn('"-Y_R * Test_Y_R"', DOC_SOURCE)

    def test_strict_initial_and_fixed_CN_discretization(self):
        for value in ('init_step_ts_type="algebraic"', 'ts_type="cn"',
                      'ts_adapt_type="none"', '"4*x*(1-x)"'):
            self.assertIn(value, EXAMPLE_SOURCE)

    def test_actual_observer_accepts_correct_step_work(self):
        result = self.observe(System())
        self.assertTrue(all(result["gates"].values()))
        self.assertLess(float(np.max(np.abs(result["energy_defect"]))),1e-14)

    def test_actual_observer_rejects_endpoint_trapezoid_instead_of_CN(self):
        sys = System()
        sys.H[-1] = 1.002
        self.assertFalse(self.observe(sys)["gates"]["stepwise_signed_energy_balance"])

    def test_actual_observer_rejects_wrong_boundary_power_sign(self):
        sys = System()
        sys.ports["Boundary control (left)"].powers = [0., -.25, -.5]
        self.assertFalse(self.observe(sys)["gates"]["cumulative_signed_energy_balance"])

    def test_actual_observer_rejects_missing_differential_preservation(self):
        sys = System()
        sys.initialization_report["differential_initial_exact"] = False
        self.assertFalse(self.observe(sys)["gates"][
            "strict_DAE_init_without_changing_dynamic_data"
        ])

    def test_actual_observer_rejects_non_CN_sampling(self):
        sys = System()
        sys.ports["Boundary control (left)"].sampling = "instantaneous"
        with self.assertRaisesRegex(RuntimeError, "requires actual CN"):
            self.observe(sys)

    def test_actual_observer_rejects_nonfinite_trajectory(self):
        sys = System()
        sys.solution["z"][1,0] = np.nan
        with self.assertRaisesRegex(RuntimeError, "Nonfinite"):
            self.observe(sys)

    def test_public_solvers_honest_about_parallel_limitations(self):
        self.assertIn('choices=("mumps", "superlu", "superlu_dist")',
                      EXAMPLE_SOURCE)
        self.assertIn('args.solver == "superlu" and comm.size != 1',
                      EXAMPLE_SOURCE)
        self.assertIn('SuperLU_DIST', DOC_SOURCE)


if __name__ == "__main__":
    unittest.main()
