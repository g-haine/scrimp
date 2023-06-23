# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2022 Ghislain Haine
#
# See the LICENSE file in the root directory for license information.
#
# github: https://github.com/g-haine/scrimp

"""
- file:             examples/swe.py
- author:           Ghislain Haine
- date:             22 nov. 2022
- last modified:    12 dec. 2022
- brief:            swe system
"""
from scrimp import *
from scrimp.utils.mesh import set_verbose_gf
from itertools import zip_longest


def shallow_water():
    """
    A structure-preserving discretization of the shallow-water equation with boundary control

    !TO DO: add Navier-Stokes dissipation (and in particular vorticity)

    Formulation DAE, Grad-Grad, uniform boundary condition on the Disk
    """

    set_verbose_gf(0)

    # Init the distributed port-Hamiltonian system
    swe = DPHS("real")
    
    ## Macros
    swe.gf_model.add_macro('div(v)', 'Trace(Grad(v))')
    swe.gf_model.add_macro('D(v)', 'Sym(Grad(v))')
    swe.gf_model.add_macro('Rot', '[[0,1],[-1,0]]')
    swe.gf_model.add_macro('Tangent', '(Rot*Normal)')
    swe.gf_model.add_macro('Curl2D(v)', 'div(Rot*v)')
    swe.gf_model.add_macro('Gyro(v)', 'Curl2D(v)*Rot')

    # Set the domain (using the built-in geometry `Rectangle`)
    # Omega = 1, Gamma_Bottom = 10, Gamma_Right = 11, Gamma_Top = 12, Gamma_Left = 13
    swe.set_domain(Domain("Disk", {"R": 1, "h": 0.15}))
    swe.domain.get_mesh()[0].region_merge(10,11)
    swe.domain.get_mesh()[0].region_merge(10,12)
    swe.domain.get_mesh()[0].region_merge(10,13)
    swe.domain.get_mesh()[0].delete_region([11,12,13])
    
    ## Define the variables and their discretizations

    states = [
        State("h", "Fluid height", "scalar-field"),
        State("p", "Linear momentum", "vector-field"),
    ]
    costates = [
        CoState("e_h", "Pressure", states[0]),
        CoState("e_p", "Velocity", states[1]),
    ]
    ports = []
    rho = 1025.0
    g = 9.81
    nu = 0.001
    params = [
        Parameter("rho", "Mass density", "scalar-field", f"{rho}", "h"),
        Parameter("g", "Gravity", "scalar-field", f"{g}", "h"),
        Parameter("nu", "Viscosity", "scalar-field", f"{nu}", "h"),
    ]

    control_ports = [
        Control_Port(
            "Boundary control 0",
            "U_0",
            "Normal velocity",
            "Y_0",
            "Fluid height",
            "scalar-field",
            region=10,
            position="effort",
        ),
        Control_Port(
            "Boundary control 1",
            "U_1",
            "Velocity",
            "Y_1",
            "~Normal derivative",
            "vector-field",
            region=10,
            position="effort",
        ),
    ]

    FEMs = [
        # name of the variable: (is the same of states, ports and controls ports), order, FEM
        FEM(states[0].get_name(), 2, FEM="CG"),
        FEM(states[1].get_name(), 1, FEM="CG"),
        FEM(control_ports[0].get_name(), 1, "CG"),
        FEM(control_ports[1].get_name(), 1, "CG"),
    ]

    for state, costate, param, fem, port, control_port in zip_longest(
        states, costates, params, FEMs, ports, control_ports
    ):
        if state is not None:
            # Add a state
            swe.add_state(state)
        if costate is not None:
            # Add its co-state
            swe.add_costate(costate)
        if fem is not None:
            # Add a Finite Element Method to the `port`
            swe.add_FEM(fem)
        if param is not None:
            # Add a (possibly space-varying) parameter to the `port`
            swe.add_parameter(param)
        if port is not None:
            # Add a resistive `port`
            swe.add_port(port)
        if control_port is not None:
            # Add a control `port` on the bottom part of the boundary (Neumann, thus position='effort' - default)
            swe.add_control_port(control_port)

    ## Set Hamiltonian
    swe.hamiltonian.set_name("Mechanical energy")

    terms = [
        Term("Kinetic energy", "0.5*h*p.p/rho", [1]),
        Term("Potential energy", "0.5*rho*g*h*h", [1]),
    ]

    for term in terms:
        # Set the Hamiltonian (can be done later, even after solve)
        swe.hamiltonian.add_term(term)

    ## Define the Dirac structure via getfem `brick` = non-zero block matrix
    bricks = [
        # Define the mass matrices of the left-hand side of the "Dirac structure" (position="flow")
        Brick("M_h", "h * Test_h", [1], dt=True, position="flow"),
        Brick("M_p", "h * p . Test_p", [1], dt=True, linear=False, position="flow"),
        Brick("M_Y_0", "Y_0 * Test_Y_0", [10], position="flow"),
        Brick("M_Y_1", "U_1 . Test_Y_1", [10], position="flow"),
        
        # Define the first line of the right-hand side of the "Dirac structure" (position="effort")
        Brick("-D^T", "h * e_p . Grad(Test_h)", [1], linear=False, position="effort"),
        # with the boundary control
        Brick("B_0", "- U_0 * Test_h", [10], position="effort"),
        # Define the second line of the right-hand side of the "Dirac structure" (position="effort")
        Brick("D", "- Grad(e_h) . Test_p * h", [1], linear=False, position="effort"),
        # with the gyroscopic term (beware that "Curl" is not available in the GWFL of getfem)
        Brick("G", "h * (Gyro(p) * e_p) . Test_p", [1], linear=False, 
               explicit=True, position="effort"),
        Brick("D(v)", "- 2 * nu * h * D(e_p) : D(Test_p)", [1], linear=False, 
               explicit=False, position="effort"),
        Brick("div(v)", "- 2 * nu * h * div(e_p) * div(Test_p)", [1], linear=False, 
               explicit=False, position="effort"),
        Brick("B_1", "Y_1 . Test_p", [10], position="effort"),
        # Define the third line of the right-hand side of the "Dirac structure" (position="effort")
        Brick("C_0", "e_h * Test_Y_0", [10], position="effort"),
        Brick("C_1", "- e_p . Test_Y_1", [10], position="effort"),
        # Define the constitutive relations (position="constitutive", the default value)
        # For e_h: first the mass matrix WITH A MINUS because we want an implicit formulation 0 = - M e_h + F(h)
        Brick("-M_e_h", "- e_h * Test_e_h", [1]),
        # second the linear part as a linear brick
        Brick("CR_h_lin", "rho * g * h * Test_e_h", [1]),
        # third the non-linear part as a non-linear brick (linear=False)
        Brick("CR_h_nl", "0.5 * (p . p) / rho * Test_e_h", [1], linear=False, 
              explicit=True),
        # For e_p: first the mass matrix WITH A MINUS because we want an implicit formulation 0 = - M e_p + F(p)
        Brick("-M_e_p", "- e_p . Test_e_p", [1]),
        # second the non-linear brick (linear=False)
        Brick("CR_p", "p / rho . Test_e_p", [1]),
    ]

    for brick in bricks:
        swe.add_brick(brick)

    ## Initialize the problem
    swe.set_control("Boundary control 0", "0.")
    swe.set_control("Boundary control 1", "0.*Normal + 2.*Tangent")

    # Set the initial data
    swe.set_initial_value(
        "h", "5."
    )
    swe.set_initial_value(
        "p",
        f"[ -2*{rho}*np.sin(0.5*np.pi*np.sqrt(x*x+y*y))*np.sin(np.arctan2(y,x)), 2*{rho}*np.sin(0.5*np.pi*np.sqrt(x*x+y*y))*np.cos(np.arctan2(y,x))]",
    )

    ## Solve in time

    # Define the time scheme (default ts_type='cn', t_f=1, dt=0.01, etc.)
    swe.set_time_scheme(
        ts_type="bdf",
        bdf_orde=6,
        ksp_type="gmres",
        pc_type="lu", # pc_factor_mat_solver_type='mumps',
        t_0=0.0,
        t_f=0.5,
        dt=0.0001,
        dt_save=0.01,
        ts_adapt_dt_min=0.00001,
        init_step=True,
    )

    # Solve the system in time
    swe.solve()

    # Plot the Hamiltonian
    swe.plot_Hamiltonian(with_powers=True)

    # Saving solutions for ParaView post-processing
    swe.export_to_pv("h")
    swe.export_to_pv("e_p")

    return swe  # For consol use

if __name__ == "__main__":
    swe = shallow_water()
