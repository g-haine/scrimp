"""Source-only contract for full-sine, lumped nonlinear pendulum pilot.

This is NOT native GetFEM/PETSc acceptance. The user-host driver
scripts/run-v130-nonlinear-pendulum-physical-pilot.sh is required.
"""
import ast
from pathlib import Path
import unittest

ROOT=Path(__file__).resolve().parents[1]


class NonlinearPendulumSourceV130(unittest.TestCase):
    def test_original_nonlinear_physics_and_explicit_IC_preserved(self):
        path=ROOT/"examples/nonlinear_pendulum.py"
        text=path.read_text(encoding="utf-8")
        ast.parse(text,filename=str(path))
        for fragment in (
            "gravity = 9.81",
            "rotational_damping = 0.08",
            'S.State("theta", "Angular position (rad)", "lumped")',
            'S.State("p", "Angular momentum", "lumped")',
            'pend.add_lumped_port("theta")',
            'pend.add_lumped_port("p")',
            'S.Brick("J_theta_p", "e_p*Test_theta/ell"',
            'S.Brick("J_p_theta", "-e_theta*Test_p/ell"',
            '"-damp*e_p*Test_p/ell"',
            '"gravity_torque*sin(theta)*Test_e_theta/ell"',
            'S.Brick("CR_momentum", "(p/inertia-e_p)*Test_e_p/ell"',
            'pend.set_initial_value("theta", theta0)',
            'pend.set_initial_value("p", p0)',
            "mass*gravity*length*np.sin(theta0)",
            'pend.gf_model.set_variable("e_p", np.array([p0/inertia]))',
            'initialization="exact"',
            '"init_step": False',
            '"init_step_ts_type": "pseudo", "init_step": True',
            'ts_type="bdf", ts_bdf_order=2',
        ):
            with self.subTest(fragment=fragment):
                self.assertIn(fragment,text)

    def test_native_probe_does_not_replace_nonlinear_model(self):
        path=ROOT/"sandbox/probe_v130_nonlinear_pendulum_physics.py"
        probe=path.read_text(encoding="utf-8")
        ast.parse(probe,filename=str(path))
        for fragment in (
            "from examples.nonlinear_pendulum import nonlinear_pendulum",
            "system=nonlinear_pendulum(",
            'system.get_solution("theta")',
            'system.get_solution("p")',
            'system.get_solution("e_theta")',
            'system.get_solution("e_p")',
            'system.get_Hamiltonian()',
            "nonlinear_step_power=(",
            "etheta[2:]*D2theta+ep[2:]*D2p",
            "G*np.sin(theta)",
            "rk4_reference",
            "PETSc_BDF2_startup_internal_stage_not_claimed",
            "no_pseudo=(args.initialization==\"exact\")",
        ):
            with self.subTest(fragment=fragment):
                self.assertIn(fragment,probe)
        self.assertNotIn("S.DPHS(",probe)
        self.assertNotIn("S.Brick(",probe)
        # A quadratic H(z)-difference identity is INVALID for nonlinear sin.
        self.assertNotIn("energy(2*",probe)
        self.assertNotIn("H[2:]-H[1:-1]",probe)

    def test_native_driver_checks_mpi_and_force(self):
        source=(ROOT/"scripts/run-v130-nonlinear-pendulum-physical-pilot.sh"
                ).read_text(encoding="utf-8")
        for fragment in (
            'for case in "free:0" "forced:0.1"',
            'for mpi in 1 2',
            "--initialization exact",
            "check-v130-published-immutable.py",
            "BDF2_nonquadratic_coenergy_power_defect",
            "parity_all_four_original_lumped_DAE_variables",
        ):
            with self.subTest(fragment=fragment):
                self.assertIn(fragment,source)


if __name__=="__main__":
    unittest.main()
