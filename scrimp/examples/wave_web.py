# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2022 Ghislain Haine
#
# See the LICENSE file in the root directory for license information.
#
# github: https://github.com/g-haine/scrimp

"""
- file:             examples/wave.py
- author:           Ghislain Haine
- date:             22 nov. 2022
- last modified:    12 dec. 2022
- brief:            wave system
"""
from scrimp import *
from scrimp.utils.mesh import set_verbose_gf
from itertools import zip_longest


def wave_eq():
    """
    A structure-preserving discretization of the wave equation with boundary control

    Formulation DAE (energy/co-energy), Grad-Grad, Mixed boundary condition on the Rectangle
    """

    set_verbose_gf(0)

    # Init the distributed port-Hamiltonian system
    wave = DPHS("real")

    # Set the domain (using the built-in geometry `Rectangle`)
    # Omega = 1, Gamma_Bottom = 10, Gamma_Right = 11, Gamma_Top = 12, Gamma_Left = 13
    wave.set_domain(Domain("Interval", {"L": 1.0, "h": 0.01}))

    ## Define the variables and their discretizations

    states = [
        State("q", "Strain", "scalar-field"),
        State("p", "Linear momentum", "scalar-field"),
    ]
    costates = [
        CoState("e_q", "Stress", states[0]),
        CoState("e_p", "Velocity", states[1]),
    ]
    ports = []
    params = [
        Parameter("T", "Young's modulus", "scalar-field", "1", "q"),
        Parameter("rho", "Mass density", "scalar-field", "1 + x*(1-x)", "p"),
    ]

    control_ports = [
        Control_Port(
            "Boundary control (left)",
            "U_L",
            "Velocity",
            "Y_L",
            "Normal force",
            "scalar-field",
            region=10,
        ),
        Control_Port(
            "Boundary control (right)",
            "U_R",
            "Velocity",
            "Y_R",
            "Normal force",
            "scalar-field",
            region=11,
        ),
    ]

    FEMs = [
        # name of the variable: (is the same of states, ports and controls ports), order, FEM
        FEM(states[0].get_name(), 2),
        FEM(states[1].get_name(), 1),
        FEM(control_ports[0].get_name(), 1),
        FEM(control_ports[1].get_name(), 1),
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
        Term("Kinetic energy", "0.5*p*p/rho", [1]),
        Term("Potential energy", "0.5*q*T*q", [1]),
    ]

    for term in terms:
        # Set the Hamiltonian (can be done later, even after solve)
        wave.hamiltonian.add_term(term)

    ## Define the Dirac structure via getfem `brick` = non-zero block matrix
    bricks = [
        Brick("M_q", "q.Test_q", [1], dt=True, position="flow"),
        Brick("M_p", "p*Test_p", [1], dt=True, position="flow"),
        Brick("M_Y_L", "Y_L*Test_Y_L", [10], position="flow"),
        Brick("M_Y_R", "Y_R*Test_Y_R", [11], position="flow"),
        Brick("D", "Grad(e_p).Test_q", [1], position="effort"),
        Brick("-D^T", "-e_q.Grad(Test_p)", [1], position="effort"),
        Brick("B_L", "U_L*Test_p", [10], position="effort"),
        Brick("B_R", "U_R*Test_p", [11], position="effort"),
        Brick("-B_L^T", "-e_p*Test_Y_L", [10], position="effort"),
        Brick("-B_R^T", "-e_p*Test_Y_R", [11], position="effort"),
        Brick("-M_e_q", "-e_q.Test_e_q", [1]),
        Brick("CR_q", "q.T*Test_e_q", [1]),
        Brick("-M_e_p", "-e_p*Test_e_p", [1]),
        Brick("CR_p", "p/rho*Test_e_p", [1]),
    ]

    for brick in bricks:
        wave.add_brick(brick)

    ## Initialize the problem
    expressions = ["sin(2*pi*t)", "0."]

    for control_port, expression in zip(control_ports, expressions):
        # Set the control functions (automatic construction of bricks such that -M_u u + f(t) = 0)
        wave.set_control(control_port.get_name(), expression)

    wave.set_initial_value("q", "10.")
    wave.set_initial_value("p", "np.exp(-50.*(x-0.5)*(x-0.5))")

    wave.solve()

    # Plot the Hamiltonian with the power supplied at the boundary
    wave.plot_Hamiltonian(save_figure=True)

    # Export solutions for ParaView
    wave.export_to_pv("q")
    wave.export_to_pv("p")

    # Plot the matrices representing the Dirac structure
    # wave.spy_Dirac()

    return wave  # For consol use


if __name__ == "__main__":
    wave = wave_eq()
