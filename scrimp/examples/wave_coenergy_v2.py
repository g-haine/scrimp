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


def wave_coenergy():
    """
    A structure-preserving discretization of the wave equation with boundary control

    Formulation co-energy, Grad-Grad, output feedback law at the boundary, damping on a subdomain
    """

    set_verbose_gf(0)

    # Init the distributed port-Hamiltonian system
    wave = DPHS("real")

    # Set the domain (using the built-in geometry `Rectangle`)
    # Omega = 1, Gamma_Bottom = 10, Gamma_Right = 11, Gamma_Top = 12, Gamma_Left = 13
    wave.set_domain(Domain("Concentric", {"R": 1.0, "r": 0.6, "h": 0.15}))

    ## Define the variables and their discretizations

    states = [
        State("q", "Stress", "vector-field"),
        State("p", "Velocity", "scalar-field"),
    ]
    costates = [
        CoState("e_q", "Stress", states[0], substituted=True),
        CoState("e_p", "Velocity", states[1], substituted=True),
    ]
    ports = [
        Port("Damping", "e_r", "e_r", "scalar-field", substituted=True, region=1),
    ]
    params = [
        Parameter(
            "Tinv",
            "Young's modulus inverse",
            "tensor-field",
            "[[5+x,x*y],[x*y,2+y]]",
            "q",
        ),
        Parameter("rho", "Mass density", "scalar-field", "3-x", "p"),
        Parameter(
            "nu",
            "Viscosity",
            "scalar-field",
            "10*(0.36-(x*x+y*y))",
            ports[0].get_name(),
        ),
    ]

    control_ports = [
        Control_Port(
            "Boundary control",
            "U",
            "Normal force",
            "Y",
            "Velocity trace",
            "scalar-field",
            region=20,
        ),
    ]

    FEMs = [
        # name of the variable: (is the same of states, ports and controls ports), order, FEM
        FEM(states[0].get_name(), 1, "DG"),
        FEM(states[1].get_name(), 2, "CG"),
        FEM(ports[0].get_name(), 1, "DG"),
        FEM(control_ports[0].get_name(), 1, "DG"),
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
        Term("Potential energy", "0.5*q.Tinv.q", [1, 2]),
        Term("Kinetic energy", "0.5*p*p*rho", [1, 2]),
    ]

    for term in terms:
        # Set the Hamiltonian (can be done later, even after solve)
        wave.hamiltonian.add_term(term)

    ## Define the Dirac structure via getfem `brick` = non-zero block matrix
    bricks = [
        # Add the mass matrices from the left-hand side: the `flow` part of the Dirac structure
        Brick("M_q", "q.Tinv.Test_q", [1, 2], dt=True, position="flow"),
        Brick("M_p", "p*rho*Test_p", [1, 2], dt=True, position="flow"),
        Brick("M_r", "e_r/nu*Test_e_r", [1], position="flow"),
        Brick("M_Y", "Y*Test_Y", [20], position="flow"),
        # Add the matrices from the right-hand side: the `effort` part of the Dirac structure
        Brick("D", "Grad(p).Test_q", [1, 2], position="effort"),
        Brick("-D^T", "-q.Grad(Test_p)", [1, 2], position="effort"),
        Brick("I_r", "e_r*Test_p", [1], position="effort"),
        Brick("B", "U*Test_p", [20], position="effort"),
        Brick("-I_r^T", "-p*Test_e_r", [1], position="effort"),
        Brick("-B^T", "-p*Test_Y", [20], position="effort"),
    ]

    for brick in bricks:
        wave.add_brick(brick)

    ## Initialize the problem
    expressions = ["0.005*Y"]

    for control_port, expression in zip(control_ports, expressions):
        # Set the control functions (automatic construction of bricks such that -M_u u + f(t) = 0)
        wave.set_control(control_port.get_name(), expression)

    # Set the initial data
    wave.set_initial_value("q", "[0., 0.]")
    wave.set_initial_value("p", "2.72**(-20*((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5)))")

    ## Solve in time

    # Define the time scheme (default ts_type='cn', t_f=1, dt=0.01, etc.)
    wave.set_time_scheme(t_f=1.0, dt=0.01, dt_save=0.01)

    # Solve
    wave.solve()

    ## Post-processing

    # Plot the Hamiltonian with the power supplied at the boundary
    wave.plot_Hamiltonian()

    # Export solutions for ParaView
    wave.export_to_pv("q")
    wave.export_to_pv("p")

    # Plot the matrices representing the Dirac structure
    # wave.spy_Dirac()

    return wave  # For consol use


if __name__ == "__main__":
    wave = wave_coenergy()
