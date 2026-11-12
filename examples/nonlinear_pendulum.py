# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Nonlinear damped pendulum, modeled as a lumped port-Hamiltonian system.

States are the angle theta (radians) and angular momentum p = I*theta_dot.
For a point mass m at distance L from the pivot:
    I = m*L**2
    H = p**2/(2*I) + m*g*L*(1-cos(theta))
    e_theta = m*g*L*sin(theta),  e_p = p/I
    theta_dot = e_p
    p_dot = -e_theta - damping*e_p + tau(t)

This example keeps the full nonlinear sine restoring torque. It does
not use the small-angle approximation sin(theta) = theta.
"""

import numpy as np
import scrimp as S


def nonlinear_pendulum(final_time=12.0, time_step=0.01,
                       initial_angle=1.2, torque_amplitude=0.0,
                       plot=False):
    """Integrate a nonlinear (possibly forced) pendulum with damping."""
    mass = 1.0
    length = 1.0
    gravity = 9.81
    rotational_damping = 0.08
    inertia = mass * length**2
    interval_length = 1.0

    pend = S.DPHS("real")
    pend.set_domain(S.Domain("Interval", {"L": interval_length, "h": 1.0}))

    angle = S.State("theta", "Angular position (rad)", "lumped")
    momentum = S.State("p", "Angular momentum", "lumped")
    for state in (angle, momentum):
        pend.add_state(state)
    pend.add_costate(S.CoState("e_theta", "Restoring torque", angle))
    pend.add_costate(S.CoState("e_p", "Angular velocity", momentum))
    pend.add_lumped_port("theta")
    pend.add_lumped_port("p")

    for name, value in {
        "ell": interval_length,
        "inertia": inertia,
        "gravity_torque": mass*gravity*length,
        "damp": rotational_damping,
    }.items():
        pend.gf_model.add_macro(name, str(float(value)))
    pend.gf_model.add_macro(
        "torque", f"{float(torque_amplitude)}*sin(2*pi*0.3*t)"
    )

    bricks = [
        S.Brick("M_theta", "theta*Test_theta/ell", [1],
                dt=True, position="flow"),
        S.Brick("M_p", "p*Test_p/ell", [1],
                dt=True, position="flow"),

        # Skew-symmetric canonical coupling.
        S.Brick("J_theta_p", "e_p*Test_theta/ell", [1],
                position="effort"),
        S.Brick("J_p_theta", "-e_theta*Test_p/ell", [1],
                position="effort"),
        S.Brick("Viscous_damping", "-damp*e_p*Test_p/ell", [1],
                position="effort"),
        S.Brick("External_torque", "torque*Test_p/ell", [1],
                position="effort", linear=False),

        # Co-energy constitutive relations.
        S.Brick("CR_theta_mass", "-e_theta*Test_e_theta/ell", [1]),
        S.Brick("CR_theta_nonlinear",
                "gravity_torque*sin(theta)*Test_e_theta/ell",
                [1], linear=False),
        S.Brick("CR_momentum", "(p/inertia-e_p)*Test_e_p/ell", [1]),
    ]
    for brick in bricks:
        pend.add_brick(brick)

    pend.set_initial_value("theta", float(initial_angle))
    pend.set_initial_value("p", 0.0)

    pend.set_time_scheme(
        ts_type="bdf", ts_bdf_order=2,
        ts_exact_final_time="matchstep",
        ts_adapt_type="none",
        ksp_type="gmres", pc_type="lu",
        pc_factor_mat_solver_type="superlu",
        t_0=0.0, t_f=float(final_time),
        dt=float(time_step), dt_save=float(time_step),
        init_step=True,
    )
    pend.solve()

    pend.hamiltonian.set_name("Nonlinear pendulum mechanical energy")
    pend.hamiltonian.add_term(
        S.Term("Gravitational potential",
               "gravity_torque*(1-cos(theta))/ell", [1])
    )
    pend.hamiltonian.add_term(
        S.Term("Kinetic", "p*p/(2*inertia*ell)", [1])
    )

    if plot:
        import matplotlib.pyplot as plt

        t = np.asarray(pend.solution["t"], dtype=float)
        angles = np.asarray(pend.get_solution("theta"), dtype=float).ravel()
        pvals = np.asarray(pend.get_solution("p"), dtype=float).ravel()
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7, 7))
        ax1.plot(t, angles, label="theta (rad)")
        ax1.set(xlabel="Time", ylabel="Angle (rad)")
        ax1.grid(True)
        ax2.plot(angles, pvals/inertia)
        ax2.set(xlabel="Angle (rad)", ylabel="Angular velocity")
        ax2.grid(True)
        fig.tight_layout()
        pend.plot_Hamiltonian(with_powers=False)
        plt.show()

    return pend


if __name__ == "__main__":
    nonlinear_pendulum(plot=True)
