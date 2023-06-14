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

    # Set the domain (using the built-in geometry `Rectangle`)
    # Omega = 1, Gamma_Bottom = 10, Gamma_Right = 11, Gamma_Top = 12, Gamma_Left = 13
    swe.set_domain(Domain("Disk", {"R": 1, "h": 0.15}))

    ## Define the variables and their discretizations

    states = [
        State("h", "Fluid height", "scalar-field"),
        State("p", "Linear momentum", "vector-field"),
    ]
    costates = [
        CoState("e_h", "Pressure", states[0]),
        CoState("e_p", "Rate of flow", states[1]),
    ]
    ports = []
    params = [
        Parameter("rho", "Mass density", "scalar-field", "1000.", "h"),
        Parameter("g", "Gravity", "scalar-field", "10.", "h"),
    ]

    control_ports = [
        Control_Port(
            "Boundary control",
            "U",
            "Normal rate of flow",
            "Y",
            "Fluid height",
            "scalar-field",
            region=10,
            position="effort",
        ),
    ]

    FEMs = [
        # name of the variable: (is the same of states, ports and controls ports), order, FEM
        FEM(states[0].get_name(), 2, FEM="CG"),
        FEM(states[1].get_name(), 1, FEM="DG"),
        FEM(control_ports[0].get_name(), 2, "CG"),
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
        # Define the mass matrices of the left-hand side of the 'Dirac structure' (position='flow')
        Brick("M_h", "h*Test_h", [1], dt=True, position="flow"),
        Brick("M_p", "p.Test_p", [1], dt=True, position="flow"),
        Brick("M_Y", "Y*Test_Y", [10], position="flow"),
        # Define the first line of the right-hand side of the 'Dirac structure' (position='effort')
        Brick("-D^T", "e_p.Grad(Test_h)", [1], position="effort"),
        # with the boundary control
        Brick("B", "-U*Test_h", [10], position="effort"),
        # Define the second line of the right-hand side of the 'Dirac structure' (position='effort')
        Brick("D", "-Grad(e_h).Test_p", [1], position="effort"),
        # with the gyroscopic term (beware that 'Curl' is not available in the GWFL of getfem)
        Brick(
            "G",
            "Trace([[0,1],[-1,0]]*Grad(p/rho))*rho*[[0, 1],[-1, 0]]*e_p/h.Test_p",
            [1],
            linear=False,
            position="effort",
        ),
        # Define the third line of the right-hand side of the 'Dirac structure' (position='effort')
        Brick("C", "e_h*Test_Y", [10], position="effort"),
        # Define the constitutive relations (position='constitutive', the default value)
        # For e_h: first the mass matrix WITH A MINUS because we want an implicit formulation 0 = - M e_h + F(h)
        Brick("-M_e_h", "-e_h*Test_e_h", [1]),
        # second the linear part as a linear brick
        Brick("CR_h_lin", "rho*g*h*Test_e_h", [1]),
        # third the non-linear part as a non-linear brick (linear=False)
        Brick("CR_h_nl", "0.5*p.p*Test_e_h/rho", [1], linear=False),
        # For e_p: first the mass matrix WITH A MINUS because we want an implicit formulation 0 = - M e_p + F(p)
        Brick("-M_e_p", "-e_p.Test_e_p", [1]),
        # second the non-linear brick (linear=False)
        Brick("CR_p", "h*p.Test_e_p/rho", [1], linear=False),
    ]

    for brick in bricks:
        swe.add_brick(brick)

    ## Initialize the problem
    expressions = ["0."]

    for control_port, expression in zip(control_ports, expressions):
        # Set the control functions (automatic construction of bricks such that -M_u u + f(t) = 0)
        swe.set_control(control_port.get_name(), expression)

    # Set the initial data
    swe.set_initial_value(
        "h", "10. + 0.1*np.exp(-50*((x-0.5)*(x-0.5) + (y-0.3)*(y-0.3)))"
    )
    swe.set_initial_value(
        "p",
        "[ - 10.*np.sin(np.pi*np.sqrt(x*x+y*y))*np.sin(np.arctan2(y,x)), 10.*np.sin(np.pi*np.sqrt(x*x+y*y))*np.cos(np.arctan2(y,x))]",
    )

    ## Solve in time

    # Define the time scheme (default ts_type='cn', t_f=1, dt=0.01, etc.)
    swe.set_time_scheme(
        ts_type="cn",
        ksp_type="preonly",
        pc_type="lu",  # pc_factor_mat_solver_type='mumps',
        t_0=0.0,
        t_f=0.5,
        dt=0.01,
        dt_save=0.01,
        init_step=True,
    )

    # Solve the system in time
    swe.solve()

    # Plot the Hamiltonian (with_powers=True will add the time-integration of the product f*e for each algebraic port (f,e))
    swe.plot_Hamiltonian(with_powers=True)

    # Saving solutions for ParaView post-processing
    swe.export_to_pv("h")
    swe.export_to_pv("p")

    return swe  # For consol use


if __name__ == "__main__":
    swe = shallow_water()
