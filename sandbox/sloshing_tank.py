# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
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


def sloshing_tank(mesh_size=0.05, final_time=0.2, plot=False):
    """Build and integrate the translating-tank model; return a SCRIMP DPHS."""
    length = 2.0
    width = 1.0
    density = 1.0
    gravity = 9.81
    viscosity = 0.04
    tank_mass = 5.0

    tank = S.DPHS("real")
    tank.set_domain(S.Domain("Interval", {"L": length, "h": mesh_size}))

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

    tank.set_time_scheme(
        ts_type="bdf",
        ts_bdf_order=2,
        ksp_type="gmres",
        pc_type="lu",
        pc_factor_mat_solver_type="superlu",
        t_0=0.0,
        t_f=float(final_time),
        dt=0.002,
        dt_save=0.01,
        init_step=True,
    )

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
    print("Max |Mt*v + integral(mom) - P| =", np.max(np.abs(defect)))
    force_values = 0.1*np.sin(2*np.pi*times)
    if len(times) > 1:
        impulse = np.r_[0.0, np.cumsum(
            0.5*(force_values[1:]+force_values[:-1])*np.diff(times)
        )]
        print("Max |P-P(0)-integral(force)| =",
              np.max(np.abs(p_vals-p_vals[0]-impulse)))

    if plot:
        tank.plot_Hamiltonian()
    return tank


if __name__ == "__main__":
    sloshing_tank(plot=True)
