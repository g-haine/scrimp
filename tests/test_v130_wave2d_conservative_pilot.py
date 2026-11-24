# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Source-only scientific guards for the canonical conservative 2D wave.

No GetFEM or PETSc imports, no manufactured 5-backend PASS. Validate
that q.T.q really CONTRACTS the parameter tensor T, that all
four mixed boundary algebraic equations use physical outward signs,
that initial fields and the left velocity input are compatible
at t=0, and that this first pilot remains independent of its
future damped and co-energy teaching examples.
"""

import ast
import math
from pathlib import Path
import unittest

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
WAVE = ROOT/"examples"/"wave.py"
DOC = ROOT/"docs"/"source"/"examples"/"wave.rst"
DAMPED = ROOT/"examples"/"wave_dissipative.py"
COENERGY = ROOT/"examples"/"wave_coenergy.py"


def source(path):
    return path.read_text(encoding="utf-8")


class Mixed2DWavePilot(unittest.TestCase):
    def test_2a_is_canonical_and_2b_2c_remain_separate(self):
        for path in (WAVE, DAMPED, COENERGY):
            with self.subTest(path=path):
                ast.parse(source(path), filename=str(path))
        self.assertIn('def wave_eq(', source(WAVE))
        self.assertIn('init_step_ts_type="algebraic"', source(WAVE))
        self.assertNotIn('init_step_ts_type="pseudo"', source(WAVE))
        # The other two examples are NOT silently promoted before
        # their distinct damping and co-energy scientific gates pass.
        self.assertIn('init_step_ts_type="pseudo"', source(DAMPED))
        self.assertIn('init_step_ts_type="pseudo"', source(COENERGY))

    def test_real_hamiltonian_contains_tensor_T(self):
        text = source(WAVE)
        self.assertIn(
            'S.Term(\n        "Potential energy (q dot T q)", "0.5*q.T.q", [1],',
            text,
        )
        self.assertIn(
            '"CR_q", "q.T.Test_e_q", [1]',
            text,
        )
        self.assertIn(
            '"T", "Anisotropic membrane tension", "tensor-field",',
            text,
        )
        self.assertNotIn('"0.5*q.q"', text)

    def test_original_tensor_positive_definite_over_rectangle(self):
        # The unchanged original symmetric T(x,y) has det(T)>0.
        for x in np.linspace(0., 2., 13):
            for y in np.linspace(0., 1., 11):
                T = np.array([[5+x, x*y], [x*y, 2+y]])
                self.assertGreater(np.linalg.eigvalsh(T).min(), 0.)
                self.assertGreater(3-x, 0.)

    def test_neumann_observation_and_dirichlet_velocity_signs(self):
        text = source(WAVE)
        for suffix in ("B", "R", "T"):
            self.assertIn(
                f'"-Y_{suffix}*Test_Y_{suffix}"',
                text,
            )
            self.assertIn(
                f'"-e_p*Test_Y_{suffix}"',
                text,
            )
            self.assertIn(f'"U_{suffix}*Test_p"', text)
        # Left boundary: U_L is physically velocity, Y_L force,
        # both contributing +U_L*Y_L to physical boundary work.
        for expr in (
            '"-U_L*Test_Y_L"', '"-e_p*Test_Y_L"',
            '"Y_L*Test_p"',
        ):
            self.assertIn(expr, text)

    def test_zero_initial_data_compatible_with_zero_control_at_t0(self):
        text = source(WAVE)
        self.assertIn('system.set_initial_value("q", "[0., 0.]")', text)
        self.assertIn('system.set_initial_value("p", "0.")', text)
        self.assertIn('"0.2*sin(4*pi*t)*sin(pi*y)"', text)
        for y in np.linspace(0., 1., 11):
            self.assertAlmostEqual(
                0.2*math.sin(4*math.pi*0.)*math.sin(math.pi*y),
                0.,
            )
        self.assertGreater(
            0.2*math.sin(4*math.pi*.04)*math.sin(math.pi*.5),
            1e-2,
        )

    def test_native_requirements_are_fail_closed(self):
        text=source(WAVE)
        for needle in (
            'ts_type="cn"', 'ts_adapt_type="none"',
            'ts_exact_final_time="matchstep"',
            'init_step=True', '"pc_type"', '"cn_step"',
            'system.last_ts_steps == len(t) - 1',
            'original_differential_initial_state_preserved',
            'nontrivial_signed_left_boundary_work',
            'zero_power_from_three_unforced_Neumann_boundaries',
            'each_CN_step_mechanical_energy_balance',
            'cumulative_CN_mechanical_energy_balance',
            '"0.5*q.T.q"', '"0.5*p*p/rho"',
        ):
            if needle == '"pc_type"':
                needle = 'pc_type="lu"'
            with self.subTest(needle=needle):
                self.assertIn(needle,text)

    def test_sphinx_conservative_uses_actual_source_in_seven_steps(self):
        doc = source(DOC)
        code = source(WAVE)
        self.assertIn("normal", doc.lower())
        self.assertIn("Neumann", doc)
        self.assertIn("wave_coenergy", doc)
        self.assertIn("wave_dissipative", doc)
        self.assertIn(r"q\cdot Tq", doc)
        self.assertEqual(doc.count(".. literalinclude:: ../../../examples/wave.py"), 8)
        self.assertEqual(doc.count("   :dedent: 4"), 7)
        self.assertIn("# TUTORIAL 2A-00:", code)
        self.assertIn("# TUTORIAL 2A-08:", code)
        self.assertIn("   :start-after: # TUTORIAL 2A-00:", doc)
        self.assertIn("   :end-before: # TUTORIAL 2A-08:", doc)
        self.assertIn("mesh_size = 0.2", doc)
        self.assertIn("final_time = 0.08", doc)
        self.assertIn("time_step = 0.008", doc)
        for j in range(1, 8):
            token=f"# TUTORIAL 2A-0{j}:"
            with self.subTest(step=j):
                self.assertIn(token, code)
                self.assertIn(f"   :start-after: {token}", doc)
        # The Sphinx code is included by reference, not copied
        # from another file whose arguments can drift.
        self.assertNotIn('init_step_ts_type="pseudo"', doc[:doc.index(
            "2B — Adding damping to the DPHS"
        )])
        self.assertIn("Adding damping", doc)

    def test_five_lu_acceptance_runner_uses_canonical_source(self):
        runner=source(ROOT/"scripts"/"run-v130-wave2d-conservative-acceptance.sh")
        self.assertIn("mumps:1 mumps:2 superlu:1 superlu_dist:1 superlu_dist:2", runner)
        self.assertIn('python -u -m examples.wave', runner)
        self.assertIn('"z"', runner)
        self.assertIn("power_left", runner)
        self.assertIn("e_q", runner)
        self.assertIn("e_p", runner)
        self.assertIn("PASS_WAVE2D_CONSERVATIVE_FIVE_LU_MPI_CN_PORTS_VTK", runner)
        self.assertNotIn("probe_v130_wave2d_conservative.py", runner)

    def test_qualified_models_are_not_modified_as_hidden_prerequisite(self):
        text = source(WAVE)
        self.assertNotIn("examples.published", text)
        self.assertNotIn("wave_dissipative", text)
        self.assertNotIn("wave_coenergy", text)
        self.assertIn("Brugnoli, Haine, Matignon (2022)", text)
        self.assertIn("10.4208/ijnam2023-1005", text)


if __name__ == "__main__":
    unittest.main()
