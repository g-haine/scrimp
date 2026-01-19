# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file in the root directory for license information.
#
# github: https://github.com/g-haine/scrimp

from scrimp import *
from itertools import zip_longest


def heat_eq():
    set_verbose_gf(0)

    # Init the distributed port-Hamiltonian system
    heat = DPHS("real")

    # Set the domain (using the built-in geometry `Rectangle`)
    # Define the variables and their discretizations
    heat.set_domain(Domain("Rectangle", {"L": 2.0, "l": 1.0, "h": 0.15}))

    # Define State/s`)
    states = [
        State("T", "Temp", "scalar-field", None, 0)
    ]

    # Define Co-State/s`)
    costates = [
        CoState("T", "Temp", states[0], True)
    ]

    # Define Ports/s`)
    ports = [
        Port("heat flux", "f_Q", "e_Q", "vector-field", 0, True, False, None)
    ]

    # Define parameters/s`)
    parameters = [
        Parameter("rho", "mass density", "scalar-field", "1.", "T"),
        Parameter("Lambda", "heat conductivity", "tensor-field",
                  "[[1.,0.],[0.,1.]]", "heat flux")
    ]

    # Define control_ports/s)
    control_ports = [
        Control_Port("Boundary control (bottom)", "U_B", "TEMP", "Y_B",
                     "normal heat flux", "scalar-field", 10, "effort", 0),

        Control_Port("Boundary control (right)", "U_R", "TEMP", "Y_R",
                     "normal heat flux", "scalar-field", 11, "effort", 0),

        Control_Port("Boundary control (top)", "U_T", "TEMP", "Y_T",
                     "normal heat flux", "scalar-field", 12, "effort", 0),

        Control_Port("Boundary control (left)", "U_L", "normal heat flux",
                     "Y_L", "TEMP", "scalar-field", 13, "flow", 0)
    ]

    # Define FEM/s`)
    FEMs = [
        FEM("T", 1, "CG"),
        FEM("heat flux", 2, "CG"),
        FEM("Boundary control (bottom)", 1, "CG"),
        FEM("Boundary control (right)", 1, "CG"),
        FEM("Boundary control (top)", 1, "CG"),
        FEM("Boundary control (left)", 1, "CG")
    ]

    for state, costate, param, fem, port, control_port in zip_longest(
        states, costates, parameters, FEMs, ports, control_ports
    ):
        # Add a state
        if state is not None:
            heat.add_state(state)

    # Add its co-state
        if costate is not None:
            heat.add_costate(costate)

    # Add a Finite Element Method to the `port`
        if fem is not None:
            heat.add_FEM(fem)

    # Add a (possibly space-varying) parameter to the `port`
        if param is not None:
            heat.add_parameter(param)

    # Add a resistive `port`
        if port is not None:
            heat.add_port(port)

    # Add a control `port` on the boundary (Neumann, thus position='effort' - default)
        if control_port is not None:
            heat.add_control_port(control_port)

    # Set Hamiltonian
    heat.hamiltonian.set_name("Lyapunov formulation")

    # Define Term/s`)
    terms = [
        Term("L^2-norm", "0.5*T*rho*T", [1], 0)
    ]

    for term in terms:
        heat.hamiltonian.add_term(term)

    # Define Bricks/s`)
    bricks = [
        Brick("M_T", "T*rho*Test_T", [1], True, True, "flow", 0),
        Brick("M_Q", "f_Q.Test_f_Q", [1], True, False, "flow", 0),
        Brick("M_Y_B", "Y_B*Test_Y_B", [10], True, False, "flow", 0),
        Brick("M_Y_R", "Y_R*Test_Y_R", [11], True, False, "flow", 0),
        Brick("M_Y_T", "Y_T*Test_Y_T", [12], True, False, "flow", 0),
        Brick("M_Y_L", "U_L*Test_Y_L", [13], True, False, "flow", 0),
        Brick("D", "-Div(e_Q)*Test_T", [1], True, False, "effort", 0),
        Brick("-D^T", "T*Div(Test_f_Q)", [1], True, False, "effort", 0),
        Brick("B_B", "-U_B*Test_f_Q.Normal", [10], True, False, "effort", 0),
        Brick("B_R", "-U_R*Test_f_Q.Normal", [11], True, False, "effort", 0),
        Brick("B_T", "-U_T*Test_f_Q.Normal", [12], True, False, "effort", 0),
        Brick("B_L", "-Y_L*Test_f_Q.Normal", [13], True, False, "effort", 0),
        Brick("C_B", "e_Q.Normal*Test_Y_B", [10], True, False, "effort", 0),
        Brick("C_R", "e_Q.Normal*Test_Y_R", [11], True, False, "effort", 0),
        Brick("C_T", "e_Q.Normal*Test_Y_T", [12], True, False, "effort", 0),
        Brick("C_L", "e_Q.Normal*Test_Y_L", [13], True, False, "effort", 0),
        Brick("-M_e_Q", "-e_Q.Test_e_Q", [1],
              True, False, "constitutive", 0),
        Brick("CR_Q", "f_Q.Lambda.Test_e_Q", [
              1], True, False, "constitutive", 0)
    ]

    for brick in bricks:
        heat.add_brick(brick)

    # Define expression/s`)
    expressions = ["t", "t", "t", "0.2"]

    for control_port, expression in zip(control_ports, expressions):
        heat.set_control(control_port.get_name(), expression)

    # Define initial_value/s`)

    heat.set_initial_value("T", "np.exp(-50*((x-1)*(x-1)+(y-0.5)*(y-0.5))**2)")

    # Solve in time
    heat.set_time_scheme(t_f=1, dt=0.01, pc_type="jacobi")

    # Solve
    heat.solve()

    # Plot the Hamiltonian with the power supplied at the boundary
    heat.plot_Hamiltonian(save_figure=True)

    return heat


if __name__ == "__main__":
    heat = heat_eq()
