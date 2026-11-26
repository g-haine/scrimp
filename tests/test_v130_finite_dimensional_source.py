"""Source-only SCRIMP v1.3 finite-dimensional pedagogy/physical-model checks.

Never pretend this qualifies a GetFEM/PETSc execution: the companion
original numerical pilot is scripts/run-v130-mass-spring-physical-pilot.sh.
"""
import ast
from pathlib import Path
import unittest

ROOT=Path(__file__).resolve().parents[1]


class FiniteDimensionalV130Source(unittest.TestCase):
    def test_pedagogic_documentation_order(self):
        text=(ROOT/"docs/source/examples.rst").read_text(encoding="utf-8")
        order=[
            "wave","wave_coenergy","maxwell_3d",
            "mass_spring_damper","nonlinear_pendulum",
            "heat","heat_wave","shallow_water",
        ]
        entries=[
            line.strip().removeprefix("examples/")
            for line in text.splitlines()
            if line.strip().startswith("examples/")
        ]
        self.assertEqual(entries,order)
        self.assertIn(":doc:`first 1D vibrating-string tutorial <started>`",text)
        self.assertEqual(len(entries),len(set(entries)))

    def test_original_mass_spring_model_not_rewritten(self):
        path=ROOT/"examples/mass_spring_damper.py"
        source=path.read_text(encoding="utf-8")
        ast.parse(source,filename=str(path))
        for fragment in (
            "mass = 1.5",
            "stiffness = 5.0",
            "damping = 0.25",
            'sys.gf_model.add_macro("force", "0.2*sin(2*pi*t)")',
            'S.State("q", "Displacement", "lumped")',
            'S.State("p", "Linear momentum", "lumped")',
            'sys.add_lumped_port("q")',
            'sys.add_lumped_port("p")',
            'S.Brick("J_qp", "e_p*Test_q/ell"',
            'S.Brick("J_pq", "-e_q*Test_p/ell"',
            '"-damp_c*e_p*Test_p/ell"',
            '"(spring_k*q-e_q)*Test_e_q/ell"',
            '"(p/osc_mass-e_p)*Test_e_p/ell"',
            'sys.set_initial_value("q", 0.4)',
            'sys.set_initial_value("p", 0.0)',
            'init_step_ts_type="pseudo"',
            'ts_type="bdf", ts_bdf_order=2',
        ):
            with self.subTest(fragment=fragment):
                self.assertIn(fragment,source)

    def test_original_model_native_audit_uses_ODE_not_a_surrogate(self):
        probe=(ROOT/"sandbox/probe_v130_mass_spring_physics.py"
               ).read_text(encoding="utf-8")
        ast.parse(probe)
        for fragment in (
            "from examples.mass_spring_damper import mass_spring_damper",
            "mass_spring_damper(",
            'system.get_solution("q")',
            'system.get_solution("p")',
            'system.get_solution("e_q")',
            'system.get_solution("e_p")',
            "system.get_Hamiltonian()",
            "ode_q_defect",
            "ode_p_defect",
            "bdf2_power_defect",
            "rk4_reference",
            "BDF2_k_ge_2_no_unverified_startup",
        ):
            with self.subTest(fragment=fragment):
                self.assertIn(fragment,probe)
        self.assertNotIn("S.DPHS(",probe)
        self.assertNotIn("S.Brick(",probe)
        runner=(ROOT/"scripts/run-v130-mass-spring-physical-pilot.sh"
                ).read_text(encoding="utf-8")
        for fragment in (
            'python -u /probe-msd.py',
            'mpiexec --oversubscribe --bind-to none -n "$mpi"',
            'for mpi in 1 2',
            "check-v130-published-immutable.py",
            "MPI1_MPI2_all_DAE_and_physical_observable_parity",
        ):
            with self.subTest(fragment=fragment):
                self.assertIn(fragment,runner)


if __name__=="__main__":
    unittest.main()
