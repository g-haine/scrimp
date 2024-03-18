# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2024 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""
- file:             examples/shallow_water.py
- authors:          Ghislain Haine
- date:             22 nov. 2022
- brief:            inviscid shallow water equations
"""

# Import scrimp
import scrimp as S

def shallow_water_eq():
    """A structure-preserving discretization of the inviscid shallow-water equation

    Formulation Grad-Grad, homogeneous boundary condition, on a tank
    """

    # Init the distributed port-Hamiltonian system
    swe = S.DPHS("real")
    
    # Set the domain (using the built-in geometry `Rectangle`)
    # Labels: Omega = 1, Gamma_Bottom = 10, Gamma_Right = 11, Gamma_Top = 12, Gamma_Left = 13
    swe.set_domain(S.Domain("Rectangle", {"L": 2.0, "l": 0.5, "h": 0.1}))
    
    # Define the states and costates
    states = [
        S.State("h", "Fluid height", "scalar-field"),
        S.State("p", "Linear momentum", "vector-field"),
    ]
    costates = [
        S.CoState("e_h", "Pressure", states[0]),
        S.CoState("e_p", "Velocity", states[1]),
    ]
    
    # Add them to the dphs
    for state in states:
        swe.add_state(state)
    for costate in costates:
        swe.add_costate(costate)

    # Define the control ports
    control_ports = [
        S.Control_Port(
            "Boundary control 0",
            "U_0",
            "Normal velocity",
            "Y_0",
            "Fluid height",
            "scalar-field",
            region=10,
            position="effort",
        ),
        S.Control_Port(
            "Boundary control 1",
            "U_1",
            "Normal velocity",
            "Y_1",
            "Fluid height",
            "scalar-field",
            region=11,
            position="effort",
        ),
        S.Control_Port(
            "Boundary control 2",
            "U_2",
            "Normal velocity",
            "Y_2",
            "Fluid height",
            "scalar-field",
            region=12,
            position="effort",
        ),
        S.Control_Port(
            "Boundary control 3",
            "U_3",
            "Normal velocity",
            "Y_3",
            "Fluid height",
            "scalar-field",
            region=13,
            position="effort",
        ),
    ]
    
    # Add them to the dphs
    for ctrl_port in control_ports:
        swe.add_control_port(ctrl_port)

    # Define the Finite Elements Method of each port
    k = 1
    FEMs = [
        S.FEM(states[0].get_name(), k+1, FEM="CG"),
        S.FEM(states[1].get_name(), k, FEM="DG"),
        S.FEM(control_ports[0].get_name(), k, "DG"),
        S.FEM(control_ports[1].get_name(), k, "DG"),
        S.FEM(control_ports[2].get_name(), k, "DG"),
        S.FEM(control_ports[3].get_name(), k, "DG"),
    ]
    
    # Add them to the dphs
    for FEM in FEMs:
        swe.add_FEM(FEM)
    
    # Define physical parameters
    rho = 1000.
    g = 10.
    parameters = [
        S.Parameter("rho", "Mass density", "scalar-field", f"{rho}", "h"),
        S.Parameter("g", "Gravity", "scalar-field", f"{g}", "h"),
    ]
    
    # Add them to the dphs
    for parameter in parameters:
        swe.add_parameter(parameter)

    # Define the pHs via `Brick` == non-zero block matrices == variational terms
    # Some macros for the sake of readability
    swe.gf_model.add_macro('div(v)', 'Trace(Grad(v))')
    swe.gf_model.add_macro('Rot', '[[0,1],[-1,0]]')
    swe.gf_model.add_macro('Curl2D(v)', 'div(Rot*v)')
    swe.gf_model.add_macro('Gyro(v)', 'Curl2D(v)*Rot')
    bricks = [
        # Define the mass matrices of the left-hand side of the "Dirac structure" (position="flow")
        S.Brick("M_h", "h * Test_h", [1], dt=True, position="flow"),
        S.Brick("M_p", "h * p . Test_p", [1], dt=True, linear=False, position="flow"),
        S.Brick("M_Y_0", "Y_0 * Test_Y_0", [10], position="flow"),
        S.Brick("M_Y_1", "Y_1 * Test_Y_1", [11], position="flow"),
        S.Brick("M_Y_2", "Y_2 * Test_Y_2", [12], position="flow"),
        S.Brick("M_Y_3", "Y_3 * Test_Y_3", [13], position="flow"),
        
        # Define the first line of the right-hand side of the "Dirac structure" (position="effort")
        S.Brick("-D^T", "h * e_p . Grad(Test_h)", [1], linear=False, position="effort"),
        # with the boundary control
        S.Brick("B_0", "- U_0 * Test_h", [10], position="effort"),
        S.Brick("B_1", "- U_1 * Test_h", [11], position="effort"),
        S.Brick("B_2", "- U_2 * Test_h", [12], position="effort"),
        S.Brick("B_3", "- U_3 * Test_h", [13], position="effort"),
        # Define the second line of the right-hand side of the "Dirac structure" (position="effort")
        S.Brick("D", "- Grad(e_h) . Test_p * h", [1], linear=False, position="effort"),
        # with the gyroscopic term (beware that "Curl" is not available in the GWFL of getfem)
        S.Brick("G", "(Gyro(p) * e_p) . Test_p", [1], linear=False, explicit=True, position="effort"),
        # Define the third line of the right-hand side of the "Dirac structure" (position="effort")
        S.Brick("C_0", "- e_h * Test_Y_0", [10], position="effort"),
        S.Brick("C_1", "- e_h * Test_Y_1", [11], position="effort"),
        S.Brick("C_2", "- e_h * Test_Y_2", [12], position="effort"),
        S.Brick("C_3", "- e_h * Test_Y_3", [13], position="effort"),
        
        ## Define the constitutive relations (position="constitutive", the default value)
        # For e_h: first the mass matrix WITH A MINUS because we want an implicit formulation 0 = - M e_h + F(h)
        S.Brick("-M_e_h", "- e_h * Test_e_h", [1]),
        # second the linear part as a linear brick
        S.Brick("Q_h", "rho * g * h * Test_e_h", [1]),
        # third the non-linear part as a non-linear brick (linear=False)
        S.Brick("P_h", "0.5 * (p . p) / rho * Test_e_h", [1], linear=False),
        # For e_p: first the mass matrix WITH A MINUS because we want an implicit formulation 0 = - M e_p + F(p)
        S.Brick("-M_e_p", "- e_p . Test_e_p", [1]),
        # second the LINEAR brick
        S.Brick("Q_p", "p / rho . Test_e_p", [1]),
    ]
    for brick in bricks:
        swe.add_brick(brick)

    ## Initialize the problem
    swe.set_control("Boundary control 0", "0.")
    swe.set_control("Boundary control 1", "0.")
    swe.set_control("Boundary control 2", "0.")
    swe.set_control("Boundary control 3", "0.")

    # Set the initial data
    swe.set_initial_value("h", "3. - (np.sign(x-0.5)+1)/3.")
    swe.set_initial_value("p", f"[ 0., 0.]")

    ## Solve in time
    # Define the time scheme
    swe.set_time_scheme(
        ts_type="bdf",
        ts_bdf_order=4,
        t_f=0.5,
        dt=0.0001,
        dt_save=0.01,
    )

    # Solve the system in time
    swe.solve()

    ## Post-processing
    # Set Hamiltonian's name
    swe.hamiltonian.set_name("Mechanical energy")
    # Define Hamiltonian terms
    terms = [
        S.Term("Kinetic energy", "0.5*h*p.p/rho", [1]),
        S.Term("Potential energy", "0.5*rho*g*h*h", [1]),
    ]
    # Add them to the Hamiltonian
    for term in terms:
        swe.hamiltonian.add_term(term)
    # Plot the Hamiltonian
    swe.plot_Hamiltonian(save_figure=True, filename="Hamiltonian_Inviscid_Shallow_Water_2D.png")
    
    # Export for ParaView
    swe.export_to_pv("h")

    return swe

if __name__ == "__main__":
    swe = shallow_water_eq()
    
