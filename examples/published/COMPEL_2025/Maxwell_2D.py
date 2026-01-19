# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""
- file:             examples/published/COMPEL_2025/Maxwell_2D.py
- authors:          Ghislain Haine
- date:             10 jan. 2024
- brief:            2D Maxwell equations in co-energy formulation
"""

import scrimp as S
from itertools import zip_longest

def Maxwell_2D():
    """A structure-preserving discretization of the 2D Maxwell equations with boundary control

    Waveguide experiment.

    Formulation co-energy, Grad-Grad, output feedback law at the boundary, damping on a subdomain
    
    This script can be used to reproduce the full order model used in the examples shown in the paper:
    @article{Gouzien_2025,
        title={{Port-Hamiltonian reduced order modelling of the 2D Maxwell equations}},
        volume={},
        ISSN={},
        DOI={},
        number={},
        journal={},
        publisher={},
        author={Mattéo Gouzien, Charles Poussot-Vassal, Ghislain Haine, Denis Matignon},
        year={2025},
        pages={}
    }
    
    Example of use: to run the 2D Maxwell problem:
        `python Maxwell_2D.py`
    
    See the file `scrimp.yml` in the folder of this script to know the exact version of each library used for the results obtained in the paper.
    
    Returns:
        the DPHS object
    """

    # Init the distributed port-Hamiltonian system
    M = S.DPHS("real")
    
    # Useful macros for the weak forms
    M.gf_model.add_macro('div(v)', 'Trace(Grad(v))')
    M.gf_model.add_macro('Rot', '[[0,1],[-1,0]]')
    M.gf_model.add_macro('Tangent', '(Rot*Normal)')
    M.gf_model.add_macro('Curl2D(v)', 'div(Rot*v)')
    M.gf_model.add_macro('GradPerp(v)', 'Rot.Grad(v)')
    M.gf_model.add_macro('Gyro(v)', 'Curl2D(v)*Rot')

    # Set the domain (using the built-in geometry `Rectangle`)
    # Omega = 1, Gamma_Bottom = 10, Gamma_Right = 11, Gamma_Top = 12, Gamma_Left = 13
    M.set_domain(S.Domain("Rectangle", {"L": 2.0, "l": 0.1, "h": 0.01}))

    ## Define the variables and their discretizations
    eps = 8.85e-3
    mu = 1.25e3
    sigma = 1.0e-4
    from math import sqrt
    Zinf = sqrt(mu/eps)
    states = [
        S.State("E", "Electric field", "vector-field"),
        S.State("H", "Magnetic field", "scalar-field"),
    ]
    costates = [
        S.CoState("E", "Electric field", states[0], substituted=True),
        S.CoState("H", "Magnetic field", states[1], substituted=True),
    ]
    params = [
        S.Parameter(
            "eps",
            "permitivitty",
            "tensor-field",
            f"[[{eps},0],[0,{eps}]]",
            "E",
        ),
        S.Parameter(
            "sigma", 
            "conductivity", 
            "tensor-field", 
            f"[[{sigma},0],[0,{sigma}]]", 
            "E"),
        S.Parameter(
            "mu",
            "permeability",
            "scalar-field",
            f"{mu}",
            "H",
        ),
        S.Parameter(
            "zw",
            "impedance",
            "scalar-field",
            "0.",
            "H",
        ),
        S.Parameter(
            "zinf",
            "impedance",
            "scalar-field",
            f"{Zinf}",
            "H",
        ),
    ]
    
    ports = []

    control_ports = [
        S.Control_Port(
            "Boundary control (B)",
            "U_B",
            "...",
            "Y_B",
            "...",
            "scalar-field",
            region=10,
            position="effort",
        ),
        S.Control_Port(
            "Boundary control (R)",
            "U_R",
            "...",
            "Y_R",
            "...",
            "scalar-field",
            region=11,
            position="effort",
        ),
        S.Control_Port(
            "Boundary control (T)",
            "U_T",
            "...",
            "Y_T",
            "...",
            "scalar-field",
            region=12,
            position="effort",
        ),
        S.Control_Port(
            "Boundary control (L)",
            "U_L",
            "...",
            "Y_L",
            "...",
            "scalar-field",
            region=13,
            position="effort",
        ),
    ]

    FEMs = [
        # name of the variable: (is the same of states, ports and controls ports), order, FEM
        S.FEM("E", 1, "DG"),
        S.FEM("H", 2, "CG"),
        S.FEM("Boundary control (B)", 1, "DG"),
        S.FEM("Boundary control (R)", 1, "DG"),
        S.FEM("Boundary control (T)", 1, "DG"),
        S.FEM("Boundary control (L)", 1, "DG"),
    ]

    for state, costate, param, fem, port, control_port in zip_longest(
        states, costates, params, FEMs, ports, control_ports
    ):
        if state is not None:
            # Add a state
            M.add_state(state)
        if costate is not None:
            # Add its co-state
            M.add_costate(costate)
        if fem is not None:
            # Add a Finite Element Method to the `port`
            M.add_FEM(fem)
        if param is not None:
            # Add a (possibly space-varying) parameter to the `port`
            M.add_parameter(param)
        if port is not None:
            # Add a resistive `port`
            M.add_port(port)
        if control_port is not None:
            # Add a control `port` on the bottom part of the boundary (Neumann, thus position='effort' - default)
            M.add_control_port(control_port)

    ## Set Hamiltonian
    M.hamiltonian.set_name("Electromagnetic energy")

    terms = [
        S.Term("Electric energy", "0.5*E.eps.E", [1]),
        S.Term("Magnetic energy", "0.5*H*mu*H", [1]),
    ]

    for term in terms:
        # Set the Hamiltonian (can be done later, even after solve)
        M.hamiltonian.add_term(term)

    ## Define the Dirac structure via getfem `brick` = non-zero block matrix
    bricks = [
        # Add the mass matrices from the left-hand side: the `flow` part of the Dirac structure
        S.Brick("M_q", "E.eps.Test_E", [1], dt=True, position="flow"),
        S.Brick("M_p", "H*mu*Test_H", [1], dt=True, position="flow"),
        S.Brick("M_Y_B", "Y_B*Test_Y_B", [10], position="flow"),
        S.Brick("M_Y_R", "Y_R*Test_Y_R", [11], position="flow"),
        S.Brick("M_Y_T", "Y_T*Test_Y_T", [12], position="flow"),
        S.Brick("M_Y_L", "Y_L*Test_Y_L", [13], position="flow"),
        # Add the matrices from the right-hand side: the `effort` part of the Dirac structure
        S.Brick("-D^T", "GradPerp(H).Test_E", [1], position="effort"),
        S.Brick("D", "-E.GradPerp(Test_H)", [1], position="effort"),
        # Joule
        S.Brick("R", "-E.sigma.Test_E", [1], position="effort"),
        # IBC
        S.Brick("IBC", "-H*zinf*Test_H", [11,13], position="effort"),
        # Controls and observations
        S.Brick("B_B", "U_B*Test_H", [10], position="effort"),
        S.Brick("-B_B^T", "-H*Test_Y_B", [10], position="effort"),
        S.Brick("B_R", "U_R*Test_H", [11], position="effort"),
        S.Brick("-B_R^T", "-H*Test_Y_R", [11], position="effort"),
        S.Brick("B_T", "U_T*Test_H", [12], position="effort"),
        S.Brick("-B_T^T", "-H*Test_Y_T", [12], position="effort"),
        S.Brick("B_L", "U_L*Test_H", [13], position="effort"),
        S.Brick("-B_L^T", "-H*Test_Y_L", [13], position="effort"),
    ]

    for brick in bricks:
        M.add_brick(brick)

    ## Initialize the problem
    expressions = ["0.", 
                   "0.",
                   "0.",
                   "exp(-(pow(t-0.5,2)/0.05))"]

    for control_port, expression in zip(control_ports, expressions):
        # Set the control functions (automatic construction of bricks such that -M_u u + f(t) = 0)
        M.set_control(control_port.get_name(), expression)

    # Set the initial data
    M.set_initial_value("E", "[0., 0.]")
    M.set_initial_value("H", "0.")

    ## Solve in time

    # Define the time scheme
    M.set_time_scheme(
                     t_f=10.0, 
                     dt=0.01, 
                     dt_save=0.01,
                     )

    # Solve
    M.solve()

    ## Post-processing

    # Plot the Hamiltonian with the power supplied at the boundary
    M.plot_Hamiltonian(with_powers=False)

    # Export solutions for ParaView
    M.export_to_pv("E")
    M.export_to_pv("H")
    
    import numpy as np
    t = np.array(M.solution["t"])
    err_B = M.get_quantity("pow(Y_B+H,2)",region=10)
    err_R = M.get_quantity("pow(Y_R+H,2)",region=11)
    err_T = M.get_quantity("pow(Y_T+H,2)",region=12)
    err_L = M.get_quantity("pow(Y_L+H,2)",region=13)
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=[8, 5])
    ax = fig.add_subplot(111)
    ax.plot(t,err_B,t,err_R,t,err_T,t,err_L)
    plt.show()
    
    # Export Matrices
    from scrimp.utils.linalg import convert_PETSc_to_scipy
    import scipy.io as sio
    
    to_exp = dict()
    
    E_dofs = M.gf_model.interval_of_variable("E")
    E_range = [E_dofs[0]+k for k in range(E_dofs[1])]
    to_exp['E_range'] = E_range
    H_dofs = M.gf_model.interval_of_variable("H")
    H_range = [H_dofs[0]+k for k in range(H_dofs[1])]
    to_exp['H_range'] = H_range
    U_B_dofs = M.gf_model.interval_of_variable("U_B")
    U_B_range = [U_B_dofs[0]+k for k in range(U_B_dofs[1])]
    to_exp['U_B_range'] = U_B_range
    Y_B_dofs = M.gf_model.interval_of_variable("Y_B")
    Y_B_range = [Y_B_dofs[0]+k for k in range(Y_B_dofs[1])]
    to_exp['Y_B_range'] = Y_B_range
    U_R_dofs = M.gf_model.interval_of_variable("U_R")
    U_R_range = [U_R_dofs[0]+k for k in range(U_R_dofs[1])]
    to_exp['U_R_range'] = U_R_range
    Y_R_dofs = M.gf_model.interval_of_variable("Y_R")
    Y_R_range = [Y_R_dofs[0]+k for k in range(Y_R_dofs[1])]
    to_exp['Y_R_range'] = Y_R_range
    U_T_dofs = M.gf_model.interval_of_variable("U_T")
    U_T_range = [U_T_dofs[0]+k for k in range(U_T_dofs[1])]
    to_exp['U_T_range'] = U_T_range
    Y_T_dofs = M.gf_model.interval_of_variable("Y_T")
    Y_T_range = [Y_T_dofs[0]+k for k in range(Y_T_dofs[1])]
    to_exp['Y_T_range'] = Y_T_range
    U_L_dofs = M.gf_model.interval_of_variable("U_L")
    U_L_range = [U_L_dofs[0]+k for k in range(U_L_dofs[1])]
    to_exp['U_L_range'] = U_L_range
    Y_L_dofs = M.gf_model.interval_of_variable("Y_L")
    Y_L_range = [Y_L_dofs[0]+k for k in range(Y_L_dofs[1])]
    to_exp['Y_L_range'] = Y_L_range
    
    E_mat = convert_PETSc_to_scipy(M.tangent_mass)
    to_exp['E_mat'] = E_mat
    JR_mat = convert_PETSc_to_scipy(M.tangent_stiffness)
    to_exp['JR_mat'] = JR_mat
    
    sio.savemat('Numelec_HIGH.mat', {'NUMELEC': to_exp})
    
    return M

if __name__ == "__main__":
    Maxwell = Maxwell_2D()
