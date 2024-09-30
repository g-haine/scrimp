# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2024 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""
- file:             examples/heat.py
- authors:          Giuseppe Ferraro, Ghislain Haine
- date:             22 nov. 2022
- brief:            2D heat equation with Lyapunov Hamiltonian
"""

# Import scrimp
import scrimp as S

def heat_eq():
    """A structure-preserving discretization of the heat equation with mixed boundary control

    Formulation with substitution of the co-state,
    Lyapunov L^2 functional, Div-Div, Mixed boundary condition on the Rectangle (including
    impedance-like absorbing boundary condition).
    """

    # Init the distributed port-Hamiltonian system
    heat = S.DPHS("real")

    # Set the domain (using the built-in geometry `Rectangle`)
    # Labels: Omega = 1, Gamma_Bottom = 10, Gamma_Right = 11, Gamma_Top = 12, Gamma_Left = 13
    heat.set_domain(S.Domain("Rectangle", {"L": 2.0, "l": 1.0, "h": 0.1}))

    # Define the variables and their discretizations and add them to the dphs
    states = [
        S.State("T", "Temperature", "scalar-field"),
    ]
    costates = [
        # Substituted=True indicates that only one variable has to be discretized on this port
        S.CoState("T", "Temperature", states[0], substituted=True)
    ]

    ports = [
        S.Port("Heat flux", "f_Q", "J_Q", "vector-field"),
    ]

    control_ports = [
        S.Control_Port(
            "Boundary control (bottom)",
            "U_B",
            "Temperature",
            "Y_B",
            "- Normal heat flux",
            "scalar-field",
            region=10,
            position="effort",
        ),
        S.Control_Port(
            "Boundary control (right)",
            "U_R",
            "Temperature",
            "Y_R",
            "- Normal heat flux",
            "scalar-field",
            region=11,
            position="effort",
        ),
        S.Control_Port(
            "Boundary control (top)",
            "U_T",
            "Temperature",
            "Y_T",
            "- Normal heat flux",
            "scalar-field",
            region=12,
            position="effort",
        ),
        S.Control_Port(
            "Boundary control (left)",
            "U_L",
            "- Normal heat flux",
            "Y_L",
            "Temperature",
            "scalar-field",
            region=13,
            position="flow",
        ),
    ]

    for state in states:
        heat.add_state(state)
    for costate in costates:
        heat.add_costate(costate)
    for port in ports:
        heat.add_port(port)
    for ctrl_port in control_ports:
        heat.add_control_port(ctrl_port)

    FEMs = [
        S.FEM(states[0].get_name(), 1, FEM="DG"),
        S.FEM(ports[0].get_name(), 2, FEM="CG"),
        S.FEM(control_ports[0].get_name(), 1, FEM="DG"),
        S.FEM(control_ports[1].get_name(), 1, FEM="DG"),
        S.FEM(control_ports[2].get_name(), 1, FEM="DG"),
        S.FEM(control_ports[3].get_name(), 1, FEM="DG"),
    ]
    for FEM in FEMs:
        heat.add_FEM(FEM)

    # Define the physical parameters
    parameters = [
        S.Parameter("rho", "Mass density times heat capacity", "scalar-field", "3.", "T"),
        S.Parameter(
            "Lambda",
            "Heat conductivity",
            "tensor-field",
            "[[1e-2,0.],[0.,1e-2]]",
            "Heat flux",
        ),
    ]
    # Add them to the dphs
    for parameter in parameters:
        heat.add_parameter(parameter)

    # Define the Dirac structure and the constitutive relations block matrices as `Brick`
    bricks = [
        # Add the mass matrices from the left-hand side: the `flow` part of the Dirac structure
        S.Brick("M_T", "T*rho*Test_T", [1], dt=True, position="flow"),
        S.Brick("M_Q", "f_Q.Test_f_Q", [1], position="flow"),
        S.Brick("M_Y_B", "Y_B*Test_Y_B", [10], position="flow"),
        S.Brick("M_Y_R", "Y_R*Test_Y_R", [11], position="flow"),
        S.Brick("M_Y_T", "Y_T*Test_Y_T", [12], position="flow"),
        # Normal trace is imposed by Lagrange multiplier on the left side == the collocated output
        S.Brick("M_Y_L", "U_L*Test_Y_L", [13], position="flow"),
        # Add the matrices from the right-hand side: the `effort` part of the Dirac structure
        S.Brick("D", "-Div(J_Q)*Test_T", [1], position="effort"),
        S.Brick("-D^T", "T*Div(Test_f_Q)", [1], position="effort"),
        S.Brick("B_B", "-U_B*Test_f_Q.Normal", [10], position="effort"),
        S.Brick("B_R", "-U_R*Test_f_Q.Normal", [11], position="effort"),
        S.Brick("B_T", "-U_T*Test_f_Q.Normal", [12], position="effort"),
        # Normal trace is imposed by Lagrange multiplier on the left side == the collocated output
        S.Brick("B_L", "-Y_L*Test_f_Q.Normal", [13], position="effort"),
        S.Brick("C_B", "J_Q.Normal*Test_Y_B", [10], position="effort"),
        S.Brick("C_R", "J_Q.Normal*Test_Y_R", [11], position="effort"),
        S.Brick("C_T", "J_Q.Normal*Test_Y_T", [12], position="effort"),
        S.Brick("C_L", "J_Q.Normal*Test_Y_L", [13], position="effort"),
        ## Define the constitutive relations as getfem `brick`
        # Fourier's law under implicit form - M_e_Q e_Q + CR_Q Q = 0
        S.Brick("-M_J_Q", "-J_Q.Test_J_Q", [1]),
        S.Brick("CR_Q", "f_Q.Lambda.Test_J_Q", [1]),
    ]
    for brick in bricks:
        heat.add_brick(brick)

    # Initialize the problem
    expressions = ["1.", "1.", "1.", "0.2*T"]

    for control_port, expression in zip(control_ports, expressions):
        # Set the control functions (automatic construction of bricks such that -M_u u + f(t) = 0)
        heat.set_control(control_port.get_name(), expression)

    # Set the initial data
    heat.set_initial_value("T", "1. + 2.*np.exp(-50*((x-1)*(x-1)+(y-0.5)*(y-0.5))**2)")

    ## Solve in time
    # Define the time scheme ("bdf" is backward differentiation formula)
    heat.set_time_scheme(t_f=5.,
                         ts_type="bdf", 
                         ts_bdf_order=2, 
                         dt=0.01,
                         )

    # Solve
    heat.solve()

    ## Post-processing
    # Set Hamiltonian name
    heat.hamiltonian.set_name("Lyapunov formulation")
    # Define the term
    terms = [
        S.Term("L^2-norm", "0.5*T*rho*T", [1]),
    ]
    # Add them to the Hamiltonian
    for term in terms:
        heat.hamiltonian.add_term(term)

    # Plot the Hamiltonian
    heat.plot_Hamiltonian(save_figure=True, filename="Hamiltonian_Heat_2D.png")

    return heat

if __name__ == "__main__":
    heat = heat_eq()
