# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2025 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""
- file:             examples/wave.py
- authors:          Giuseppe Ferraro, Ghislain Haine
- date:             22 nov. 2022
- brief:            2D wave equations
"""

# Import scrimp
import scrimp as S

def wave_eq():
    """A structure-preserving discretization of the wave equation with mixed boundary control

    Formulation DAE (energy/co-energy), Grad-Grad, Mixed boundary condition on the Rectangle
    Undamped case.
    """
    
    # Init the distributed port-Hamiltonian system
    wave = S.DPHS("real")

    # Set the domain (using the built-in geometry `Rectangle`)
    # Labels: Omega = 1, Gamma_Bottom = 10, Gamma_Right = 11, Gamma_Top = 12, Gamma_Left = 13
    rectangle = S.Domain("Rectangle", {"L": 2.0, "l": 1.0, "h": 0.1})

    # And add it to the dphs
    wave.set_domain(rectangle)

    # Define the variables
    states = [
        S.State("q", "Strain", "vector-field"),
        S.State("p", "Linear momentum", "scalar-field"),
    ]
    costates = [
        S.CoState("e_q", "Stress", states[0]),
        S.CoState("e_p", "Velocity", states[1]),
    ]

    # Add them to the dphs
    for state in states:
        wave.add_state(state)
    for costate in costates:
        wave.add_costate(costate)

    # Define the control ports
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

    # Add them to the dphs
    for ctrl_port in control_ports:
        wave.add_control_port(ctrl_port)

    # Define the Finite Elements Method of each port
    FEMs = [
        S.FEM(states[0].get_name(), 1, "DG"),
        S.FEM(states[1].get_name(), 2, "CG"),
        S.FEM(control_ports[0].get_name(), 1, "DG"),
        S.FEM(control_ports[1].get_name(), 1, "DG"),
        S.FEM(control_ports[2].get_name(), 1, "DG"),
        S.FEM(control_ports[3].get_name(), 1, "DG"),
    ]

    # Add them to the dphs
    for FEM in FEMs:
        wave.add_FEM(FEM)

    # Define physical parameters
    parameters = [
        S.Parameter("T", "Young's modulus", "tensor-field", "[[5+x,x*y],[x*y,2+y]]", "q"),
        S.Parameter("rho", "Mass density", "scalar-field", "3-x", "p"),
    ]

    # Add them to the dphs
    for parameter in parameters:
        wave.add_parameter(parameter)

    # Define the pHs via `Brick` == non-zero block matrices == variational terms
    bricks = [
        ## Define the Dirac structure
        # Define the mass matrices from the left-hand side: the `flow` part of the Dirac structure
        S.Brick("M_q", "q.Test_q", [1], dt=True, position="flow"),
        S.Brick("M_p", "p*Test_p", [1], dt=True, position="flow"),
        S.Brick("M_Y_B", "Y_B*Test_Y_B", [10], position="flow"),
        S.Brick("M_Y_R", "Y_R*Test_Y_R", [11], position="flow"),
        S.Brick("M_Y_T", "Y_T*Test_Y_T", [12], position="flow"),
        # The Dirichlet term is applied via Lagrange multiplier == the colocated output
        S.Brick("M_Y_L", "U_L*Test_Y_L", [13], position="flow"),
        # Define the matrices from the right-hand side: the `effort` part of the Dirac structure
        S.Brick("D", "Grad(e_p).Test_q", [1], position="effort"),
        S.Brick("-D^T", "-e_q.Grad(Test_p)", [1], position="effort"),
        S.Brick("B_B", "U_B*Test_p", [10], position="effort"),
        S.Brick("B_R", "U_R*Test_p", [11], position="effort"),
        S.Brick("B_T", "U_T*Test_p", [12], position="effort"),
        # The Dirichlet term is applied via Lagrange multiplier == the colocated output
        S.Brick("B_L", "Y_L*Test_p", [13], position="effort"),
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
    ]

    # Add all these `Bricks` to the dphs
    for brick in bricks:
        wave.add_brick(brick)

    ## Initialize the problem
    # The controls expression, ordered as the control_ports
    t_f = 5.
    expressions = ["0.", "0.", "0.", f"0.1*sin(4.*t)*sin(4*pi*y)*exp(-10.*pow((0.5*{t_f}-t),2))"]

    # Add each expression to its control_port
    for control_port, expression in zip(control_ports, expressions):
        # Set the control functions: it automatically constructs the related `Brick`s such that `- M_u u + f(t) = 0`
        wave.set_control(control_port.get_name(), expression)

    # Set the initial data
    q_0 = "[0., 0.]"
    wave.set_initial_value("q", q_0)
    p_0 = "3**(-20*((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5)))"
    wave.set_initial_value("p", p_0)

    ## Solve in time
    # Define the time scheme ("cn" is Crank-Nicolson)
    wave.set_time_scheme(ts_type="cn",
                         t_f=t_f,
                         dt_save=0.01,
                         )

    # Solve
    wave.solve()

    ## Post-processing
    # Set Hamiltonian's name
    wave.hamiltonian.set_name("Mechanical energy")
    # Define each Hamiltonian Term
    terms = [
        S.Term("Potential energy", "0.5*q.T.q", [1]),
        S.Term("Kinetic energy", "0.5*p*p/rho", [1]),
    ]
    # Add them to the Hamiltonian
    for term in terms:
        wave.hamiltonian.add_term(term)

    # Plot the Hamiltonian and save the output
    wave.plot_Hamiltonian(save_figure=True, filename="Hamiltonian_Wave_2D_Conservative.png")
    
    return wave

if __name__ == "__main__":
    wave = wave_eq()
    
