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
        # The only supported path initializes the exact nonlinear efforts.
        self.assertNotIn("legacy_pseudo",text)
        self.assertNotIn("init_step_ts_type",text)
        self.assertNotIn("initialization=",text)
        self.assertNotIn('ts_type="bdf"',text)
        self.assertNotIn("BDF2",text)
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
            'init_step=False,',
            'ts_type="cn",',
        ):
            with self.subTest(fragment=fragment):
                self.assertIn(fragment,text)

    def test_native_probe_does_not_replace_nonlinear_model(self):
        path=ROOT/"sandbox/probe_v130_nonlinear_pendulum_physics.py"
        probe=path.read_text(encoding="utf-8")
        ast.parse(probe,filename=str(path))
        for fragment in (
            "from examples.nonlinear_pendulum import nonlinear_pendulum",
            "system = nonlinear_pendulum(",
            'system.get_solution("theta")',
            'system.get_solution("p")',
            'system.get_solution("e_theta")',
            'system.get_solution("e_p")',
            'system.get_Hamiltonian()',
            "cn_power_defect = cn_coenergy_rate+cn_dissipation-cn_input_power",
            "mean_etheta = .5*(etheta[1:]+etheta[:-1])",
            "mean_ep = .5*(ep[1:]+ep[:-1])",
            "cn_nonlinear_chain_rule_defect = cn_H_rate-cn_coenergy_rate",
            "CN_coenergy_power_identity",
            "G*np.sin(theta)",
            "rk4_reference",
            '"PETSc_TSCN_endpoint_trapezoidal"',
            'system.time_scheme.getString("init_step","true")!="false"',
        ):
            with self.subTest(fragment=fragment):
                self.assertIn(fragment,probe)
        self.assertNotIn("S.DPHS(",probe)
        self.assertNotIn("S.Brick(",probe)
        self.assertNotIn("legacy_pseudo",probe)
        self.assertNotIn("--initialization",probe)
        # A quadratic H(z)-difference identity is INVALID for nonlinear sin.
        self.assertNotIn("energy(2*",probe)
        self.assertNotIn("H[2:]-H[1:-1]",probe)
        self.assertNotIn("BDF2",probe)

    def test_canonical_example_generates_own_CN_figures(self):
        source=(ROOT/"examples/nonlinear_pendulum.py"
                ).read_text(encoding="utf-8")
        ast.parse(source)
        for fragment in (
            'if plot:',
            'pend.pendulum_visuals = render_pendulum_visuals(',
            'def render_pendulum_visuals(',
            'system.get_solution("theta")',
            'system.get_solution("p")',
            'system.get_solution("e_theta")',
            'system.get_solution("e_p")',
            'system.get_Hamiltonian()',
            'pendulum_trajectory.png',
            'pendulum_phase_portrait.png',
            'pendulum_energy_and_power.png',
            'pendulum_motion.gif',
            'cn_chain_rule_defect=cn_H_rate-cn_coenergy_rate',
            'nonquadratic_H_CN_chain_rule_defect_not_assumed_zero',
            'nonlinear_pendulum(plot=True)',
        ):
            with self.subTest(fragment=fragment):
                self.assertIn(fragment,source)
        self.assertNotIn("BDF2",source)

    def test_native_driver_checks_mpi_and_force(self):
        source=(ROOT/"scripts/run-v130-nonlinear-pendulum-physical-pilot.sh"
                ).read_text(encoding="utf-8")
        self.assertNotIn("--initialization",source)
        self.assertNotIn("BDF2",source)
        for fragment in (
            'for case in "free:0" "forced:0.1"',
            'for mpi in 1 2',
            '--solver mumps',
            "check-v130-published-immutable.py",
            "CN_coenergy_power_defect",
            "parity_all_four_original_lumped_DAE_variables",
        ):
            with self.subTest(fragment=fragment):
                self.assertIn(fragment,source)


if __name__=="__main__":
    unittest.main()
