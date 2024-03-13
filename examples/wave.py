# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2024 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""
- file:             examples/wave.py
- authors:          Giuseppe Ferraro, Ghislain Haine
- date:             22 nov. 2022
- brief:            wave equations
"""

import scrimp as S
from itertools import zip_longest

def wave_eq():
    """A structure-preserving discretization of the wave equation with boundary control

    Formulation DAE (energy/co-energy), Grad-Grad, Mixed boundary condition on the Rectangle
    """

    # Init the distributed port-Hamiltonian system
    wave = S.DPHS("real")

    # Set the domain (using the built-in geometry `Rectangle`)
    # Omega = 1, Gamma_Bottom = 10, Gamma_Right = 11, Gamma_Top = 12, Gamma_Left = 13
    wave.set_domain(S.Domain("Rectangle", {"L": 2.0, "l": 1.0, "h": 0.1}))

    ## Define the variables and their discretizations

    states = [
        S.State("q", "Strain", "vector-field"),
        S.State("p", "Linear momentum", "scalar-field"),
    ]
    costates = [
        S.CoState("e_q", "Stress", states[0]),
        S.CoState("e_p", "Velocity", states[1]),
    ]
    ports = [
        S.Port("Damping", "f_r", "e_r", "scalar-field"),
    ]
    params = [
        S.Parameter("T", "Young's modulus", "tensor-field", "[[5+x,x*y],[x*y,2+y]]", "q"),
        S.Parameter("rho", "Mass density", "scalar-field", "3-x", "p"),
        S.Parameter("nu", "viscosity", "scalar-field", "0.5*(2.0-x)", "Damping"),
    ]

    control_ports = [
        S.Control_Port(
            "Boundary control (bottom)",
            "U_B",
            "Normal force",
            "Y_B",
            "Velocity trace",
            "scalar-field",
            region=10,
        ),
        S.Control_Port(
            "Boundary control (right)",
            "U_R",
            "Normal force",
            "Y_R",
            "Velocity trace",
            "scalar-field",
            region=11,
        ),
        S.Control_Port(
            "Boundary control (top)",
            "U_T",
            "Normal force",
            "Y_T",
            "Velocity trace",
            "scalar-field",
            region=12,
        ),
        S.Control_Port(
            "Boundary control (left)",
            "U_L",
            "Velocity trace",
            "Y_L",
            "Normal force",
            "scalar-field",
            region=13,
            position="flow",
        ),
    ]

    # Define the Finite Elements Method of each port
    FEMs = [
        S.FEM(states[0].get_name(), 1, "DG"),
        S.FEM(states[1].get_name(), 2, "CG"),
        S.FEM(ports[0].get_name(), 1, "DG"),
        S.FEM(control_ports[0].get_name(), 1, "DG"),
        S.FEM(control_ports[1].get_name(), 1, "DG"),
        S.FEM(control_ports[2].get_name(), 1, "DG"),
        S.FEM(control_ports[3].get_name(), 1, "DG"),
    ]

    for state, costate, param, fem, port, control_port in zip_longest(
        states, costates, params, FEMs, ports, control_ports
    ):
        if state is not None:
            # Add a state
            wave.add_state(state)
        if costate is not None:
            # Add its co-state
            wave.add_costate(costate)
        if fem is not None:
            # Add a Finite Element Method to the `port`
            wave.add_FEM(fem)
        if param is not None:
            # Add a (possibly space-varying) parameter to the `port`
            wave.add_parameter(param)
        if port is not None:
            # Add a resistive `port`
            wave.add_port(port)
        if control_port is not None:
            # Add a control `port` on the bottom part of the boundary (Neumann, thus position='effort' - default)
            wave.add_control_port(control_port)

    ## Set Hamiltonian
    wave.hamiltonian.set_name("Mechanical energy")

    terms = [
        S.Term("Potential energy", "0.5*q.T.q", [1]),
        S.Term("Kinetic energy", "0.5*p*p/rho", [1]),
    ]

    for term in terms:
        # Set the Hamiltonian (can be done later, even after solve)
        wave.hamiltonian.add_term(term)

    ## Define the pHs via `Brick` == non-zero block matrices == variational terms
    bricks = [
        # Add the mass matrices from the left-hand side: the `flow` part of the Dirac structure
        S.Brick("M_q", "q.Test_q", [1], dt=True, position="flow"),
        S.Brick("M_p", "p*Test_p", [1], dt=True, position="flow"),
        S.Brick("M_r", "f_r*Test_f_r", [1], position="flow"),
        S.Brick("M_Y_B", "Y_B*Test_Y_B", [10], position="flow"),
        S.Brick("M_Y_R", "Y_R*Test_Y_R", [11], position="flow"),
        S.Brick("M_Y_T", "Y_T*Test_Y_T", [12], position="flow"),
        # The Dirichlet term is applied via Lagrange multiplier == the colocated output
        S.Brick("M_Y_L", "U_L*Test_Y_L", [13], position="flow"),
        # Add the matrices from the right-hand side: the `effort` part of the Dirac structure
        S.Brick("D", "Grad(e_p).Test_q", [1], position="effort"),
        S.Brick("-D^T", "-e_q.Grad(Test_p)", [1], position="effort"),
        S.Brick("I_r", "e_r*Test_p", [1], position="effort"),
        S.Brick("B_B", "U_B*Test_p", [10], position="effort"),
        S.Brick("B_R", "U_R*Test_p", [11], position="effort"),
        S.Brick("B_T", "U_T*Test_p", [12], position="effort"),
        # The Dirichlet term is applied via Lagrange multiplier == the colocated output
        S.Brick("B_L", "Y_L*Test_p", [13], position="effort"),
        S.Brick("-I_r^T", "-e_p*Test_f_r", [1], position="effort"),
        S.Brick("C_B", "-e_p*Test_Y_B", [10], position="effort"),
        S.Brick("C_R", "-e_p*Test_Y_R", [11], position="effort"),
        S.Brick("C_T", "-e_p*Test_Y_T", [12], position="effort"),
        S.Brick("C_L", "-e_p*Test_Y_L", [13], position="effort"),
        ## Define the constitutive relations
        # Hooke's law under implicit form `- M_e_q e_q + CR_q q = 0`
        S.Brick("-M_e_q", "-e_q.Test_e_q", [1]),
        S.Brick("CR_q", "q.T.Test_e_q", [1]),
        # Linear momentum definition under implicit form `- M_e_p e_p + CR_p p = 0`
        S.Brick("-M_e_p", "-e_p*Test_e_p", [1]),
        S.Brick("CR_p", "p/rho*Test_e_p", [1]),
        # Linear viscous fluid damping `- M_e_r e_r + CR_r f_r = 0`
        S.Brick("-M_e_r", "-e_r*Test_e_r", [1]),
        S.Brick("CR_r", "nu*f_r*Test_e_r", [1]),
    ]

    for brick in bricks:
        wave.add_brick(brick)

    ## Initialize the problem
    t_f = 5.0
    expressions = ["0.", "0.", "0.", f"0.1*sin(4.*t)*sin(4*pi*y)*exp(-10.*pow((0.5*{t_f}-t),2))"]

    for control_port, expression in zip(control_ports, expressions):
        # Set the control functions: it automatically constructs the related `Brick`s such that `- M_u u + f(t) = 0`
        wave.set_control(control_port.get_name(), expression)

    # Set the initial data
    wave.set_initial_value("q", "[0., 0.]")
    wave.set_initial_value("p", "3**(-20*((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5)))")

    ## Solve in time

    # Define the time scheme
    wave.set_time_scheme(ts_type="cn",
                         t_f=t_f, 
                         dt_save=0.01,
                         )

    # Solve
    wave.solve()

    ## Post-processing

    # Plot the Hamiltonian with the power supplied at the boundary
    wave.plot_Hamiltonian(save_figure=True)

    # Export solutions for ParaView
    # wave.export_to_pv("q")
    # wave.export_to_pv("p")

    return wave  # For consol use

if __name__ == "__main__":
    wave = wave_eq()
    
