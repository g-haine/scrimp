# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Viscous 1-D shallow water in a translating rigid tank (coupled PDE--ODE).

Experimental sandbox model based on SWE_Tank.pdf and swe_diss_tank.ipynb.
Requires the global scalar port extension in this branch of SCRIMP.

Fluid conservative states (per unit length):
    a_h = b*h                       (fluid cross-section)
    mom = rho*a_h*(u + v)           (absolute fluid momentum)
Lumped states:
    P = M_T*v + integral(mom) dx    (total momentum)
    pos                              (tank position)
Algebraic scalar:
    v = (P - integral(mom) dx) / M_T (tank velocity)

Hamiltonian:
    H = integral(mom**2/(2*rho*a_h) + rho*g*a_h**2/(2*b)) dx
        + M_T*v**2/2.

Using the conservative momentum removes the state-dependent mass matrix
of the original alpha_2 formulation. The modulated fluid Dirac operator is
    J12 e_m = -partial_x(a_h*e_m),
    J21 e_h = -a_h*partial_x(e_h),
    J22 e_m = -(partial_x(mom*e_m) + mom*partial_x(e_m)),
and the dissipation is -integral(rho*nu*a_h*|partial_x(e_m)|**2).
Weak forms drop the boundary power flux (closed / impermeable tank).
See README_lumped_tank.md for the remaining numerical caveats.
"""

import numpy as np
import scrimp as S
import getfem as gf

# GetFEM trace verbosity is unrelated to its unconditional C++ MPI
# "Assembly time" messages (compiled with GETFEM_PARA_LEVEL > 0).
gf.util_trace_level(0)


def sloshing_tank(
    mesh_size=0.05, final_time=0.2, plot=False,
    time_step=0.002, save_step=None, adaptive=True,
):
    """Build and integrate the translating-tank model; return a SCRIMP DPHS.

    Set adaptive=False for temporal convergence studies with fixed time steps.
    The final time is matched exactly; all runs save the final state.
    """
    if mesh_size <= 0 or final_time <= 0 or time_step <= 0:
        raise ValueError("mesh_size, final_time and time_step must be positive")
    if save_step is None:
        save_step = 0.01 if adaptive else time_step
    if save_step <= 0:
        raise ValueError("save_step must be positive")
    length = 2.0
    width = 1.0
    density = 1.0
    gravity = 9.81
    viscosity = 0.04
    tank_mass = 5.0

    tank = S.DPHS("real")
    tank.set_domain(S.Domain("Interval", {"L": length, "h": mesh_size}))
    # Physical geometry used by the optional MP4 renderer; not display scaling.
    tank.sloshing_visualization = {"length": length, "width": width}

    a_h = S.State("a_h", "Fluid cross-sectional area", "scalar-field")
    mom = S.State("mom", "Fluid absolute linear momentum", "scalar-field")
    total_p = S.State("P", "Total fluid plus tank momentum", "lumped")
    position = S.State("pos", "Tank position", "lumped")
    states = [a_h, mom, total_p, position]
    costates = [
        S.CoState("e_h", "Fluid area co-energy", a_h),
        S.CoState("e_m", "Relative fluid velocity", mom),
        S.CoState("e_P", "Tank velocity", total_p),
        S.CoState("e_pos", "Tank position co-energy", position),
    ]
    for state in states:
        tank.add_state(state)
    for costate in costates:
        tank.add_costate(costate)

    # An algebraic substituted port contributes only the unknown v, with
    # exactly ONE global degree of freedom, no elementwise DG(0) duplicates.
    tank.add_port(
        S.Port("tank_velocity", "v", "unused", "lumped",
               algebraic=True, substituted=True)
    )

    # The effort variables of dynamic ports share their state's discretization.
    tank.add_FEM(S.FEM("a_h", 2, "CG"))
    tank.add_FEM(S.FEM("mom", 1, "CG"))
    tank.add_lumped_port("P")
    tank.add_lumped_port("pos")
    tank.add_lumped_port("tank_velocity")

    # GetFEM macros are constants used in GWFL expressions below.
    for symbol, value in {
        "ell": length,
        "width": width,
        "rho": density,
        "grav": gravity,
        "nu": viscosity,
        "Mt": tank_mass,
    }.items():
        tank.gf_model.add_macro(symbol, str(float(value)))
    tank.gf_model.add_macro("force", "0.1*sin(2*pi*t)")

    bricks = [
        # Linear time derivatives (no state-dependent mass matrices).
        S.Brick("M_area", "a_h*Test_a_h", [1], dt=True, position="flow"),
        S.Brick("M_momentum", "mom*Test_mom", [1], dt=True,
                position="flow"),
        # The integral of Test_P/ell (likewise Test_pos/ell) equals 1.
        S.Brick("M_total_P", "P*Test_P/ell", [1], dt=True,
                position="flow"),
        S.Brick("M_position", "pos*Test_pos/ell", [1], dt=True,
                position="flow"),

        # Fluid skew-symmetric Dirac structure, with no boundary power input.
        S.Brick("J_area_mom", "a_h*e_m*Grad(Test_a_h)", [1],
                linear=False, position="effort"),
        S.Brick("J_mom_area", "-a_h*Grad(e_h)*Test_mom", [1],
                linear=False, position="effort"),
        S.Brick("J_mom_mom",
                "mom*e_m*Grad(Test_mom)-mom*Grad(e_m)*Test_mom",
                [1], linear=False, position="effort"),
        S.Brick("Viscous_dissipation",
                "-rho*nu*a_h*Grad(e_m)*Grad(Test_mom)",
                [1], linear=False, position="effort"),

        # Canonical lumped pair (P,pos), with external mechanical force.
        S.Brick("J_P_pos", "-e_pos*Test_P/ell", [1],
                position="effort"),
        S.Brick("J_pos_P", "e_P*Test_pos/ell", [1],
                position="effort"),
        S.Brick("Applied_force", "force*Test_P/ell", [1],
                linear=False, position="effort"),

        # Fluid effort e_h = rho*g*a_h/width - mom^2/(2*rho*a_h^2).
        S.Brick("CR_h_mass", "-e_h*Test_e_h", [1]),
        S.Brick("CR_h_gravity", "rho*grav*a_h/width*Test_e_h", [1]),
        S.Brick("CR_h_kinetic",
                "-0.5*mom*mom/(rho*a_h*a_h)*Test_e_h",
                [1], linear=False),

        # Fluid effort e_m = mom/(rho*a_h) - v = u.
        S.Brick("CR_m_mass", "-e_m*Test_e_m", [1]),
        S.Brick("CR_m_momentum", "mom/(rho*a_h)*Test_e_m", [1],
                linear=False),
        S.Brick("CR_m_velocity", "-v*Test_e_m", [1]),

        # Lumped efforts e_P = v and e_pos = 0.
        S.Brick("CR_P", "(v-e_P)*Test_e_P/ell", [1]),
        S.Brick("CR_position", "-e_pos*Test_e_pos/ell", [1]),

        # Global constraint: Mt*v + int(mom) dx - P = 0.
        # Integration over the *whole* interval assembles one scalar row.
        S.Brick("Tank_velocity_constraint",
                "(mom+(Mt*v-P)/ell)*Test_v", [1]),
    ]
    for brick in bricks:
        tank.add_brick(brick)

    # All four differential states have consistent initial data.
    # The initial fluid velocity and tank velocity are zero.
    tank.set_initial_value("a_h",
                           "2.0+0.06*np.exp(-12.0*(x-0.8)**2)")
    tank.set_initial_value("mom", "0.0")
    tank.set_initial_value("P", 0.0)
    tank.set_initial_value("pos", 0.0)
    tank.gf_model.set_variable("v", np.array([0.0]))

    time_options = dict(
        ts_type="bdf",
        ts_bdf_order=2,
        # Avoid time overshoot: previously t_f=0.2 stopped at 0.206996
        # and SCRIMP subsequently discarded the last step.
        ts_exact_final_time="matchstep",
        ksp_type="gmres",
        pc_type="lu",
        pc_factor_mat_solver_type="superlu",
        t_0=0.0,
        t_f=float(final_time),
        dt=float(time_step),
        dt_save=float(save_step),
        init_step=True,
    )
    if not adaptive:
        time_options["ts_adapt_type"] = "none"
    tank.set_time_scheme(**time_options)

    tank.solve()

    tank.hamiltonian.set_name("Total fluid and tank mechanical energy")
    tank.hamiltonian.add_term(
        S.Term("Fluid kinetic", "mom*mom/(2*rho*a_h)", [1])
    )
    tank.hamiltonian.add_term(
        S.Term("Fluid gravitational", "rho*grav*a_h*a_h/(2*width)", [1])
    )
    tank.hamiltonian.add_term(
        S.Term("Tank kinetic", "Mt*v*v/(2*ell)", [1])
    )

    # Check the algebraic global constraint and total momentum balance.
    times = np.asarray(tank.solution["t"], dtype=float)
    p_vals = np.asarray(tank.get_solution("P"), dtype=float).ravel()
    v_vals = np.asarray(tank.get_solution("v"), dtype=float).ravel()
    fluid_momentum = np.asarray(tank.get_quantity("mom"), dtype=float).ravel()
    defect = tank_mass*v_vals + fluid_momentum - p_vals
    constraint_error = float(np.max(np.abs(defect)))
    print("Max |Mt*v + integral(mom) - P| =", constraint_error)
    print("Last stored time =", times[-1], "(target =", final_time, ")")
    force_values = 0.1*np.sin(2*np.pi*times)
    fluid_mass = np.asarray(tank.get_quantity("a_h"), dtype=float).ravel()
    mass_error = float(np.max(np.abs(fluid_mass-fluid_mass[0])))
    print("Max |fluid mass - initial mass| =", mass_error)

    # BDF2 may add numerical dissipation; inspect the resulting defect.
    viscous_diss = np.asarray(
        tank.get_quantity("rho*nu*a_h*(Grad_e_m.Grad_e_m)"),
        dtype=float,
    ).ravel()
    energy = np.asarray(tank.get_Hamiltonian(), dtype=float).ravel()

    if len(times) > 1:
        impulse = np.r_[0.0, np.cumsum(
            0.5*(force_values[1:]+force_values[:-1])*np.diff(times)
        )]
        impulse_error = float(
            np.max(np.abs(p_vals-p_vals[0]-impulse))
        )
        # Exact impulse for the specific force(t)=0.1*sin(2*pi*t).
        exact_impulse = 0.1/(2*np.pi)*(1-np.cos(2*np.pi*times))
        exact_impulse_error = float(
            np.max(np.abs(p_vals-p_vals[0]-exact_impulse))
        )
        print("Max |P-P(0)-trapezoidal impulse| =", impulse_error)
        print("Max |P-P(0)-exact impulse| =", exact_impulse_error)
        power = force_values*v_vals - viscous_diss
        expected_increment = np.r_[0.0, np.cumsum(
            0.5*(power[1:]+power[:-1])*np.diff(times)
        )]
        energy_error = float(
            np.max(np.abs(energy-energy[0]-expected_increment))
        )
        print("Max |H-H(0)-integral(force*v-dissipation)| =", energy_error)
    else:
        impulse_error = exact_impulse_error = energy_error = float("nan")
        print("At least two saved time steps are needed for work/energy checks.")

    tank.sloshing_metrics = {
        "requested_final_time": float(final_time),
        "last_saved_time": float(times[-1]),
        "mesh_size": float(mesh_size),
        "time_step": float(time_step),
        "adaptive": bool(adaptive),
        "constraint_error": constraint_error,
        "mass_error": mass_error,
        "impulse_trapezoid_error": impulse_error,
        "impulse_exact_error": exact_impulse_error,
        "energy_balance_error": energy_error,
        "energy_final_change": float(energy[-1]-energy[0]),
    }

    if plot:
        # The gravitational baseline (~40 for these parameters) hides both
        # kinetic energies on an absolute-H plot. Plot their changes instead.
        import matplotlib.pyplot as plt

        fluid_ke = np.asarray(
            tank.get_quantity("mom*mom/(2*rho*a_h)"), dtype=float
        ).ravel()
        tank_ke = 0.5*tank_mass*v_vals*v_vals

        fig, ax = plt.subplots()
        ax.plot(times, fluid_ke, label="Fluid kinetic")
        ax.plot(times, tank_ke, label="Tank kinetic")
        ax.set(
            title="Kinetic energies (separate from gravitational baseline)",
            xlabel="Time t", ylabel="Energy",
        )
        ax.grid(True)
        ax.legend()
        fig.tight_layout()

        fig, ax = plt.subplots()
        ax.plot(times, energy-energy[0], label="H(t) - H(0)")
        if len(times) > 1:
            ax.plot(times, expected_increment, "--",
                    label="Integrated work - viscous dissipation")
            ax.plot(times, energy-energy[0]-expected_increment, ":",
                    label="Balance defect")
        ax.set(
            title="Total mechanical energy increment and power balance",
            xlabel="Time t", ylabel="Energy relative to initial state",
        )
        ax.grid(True)
        ax.legend()
        fig.tight_layout()
        plt.show()
    return tank


if __name__ == "__main__":
    sloshing_tank(plot=True)
