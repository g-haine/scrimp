# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2022 Ghislain Haine
#
# See the LICENSE file in the root directory for license information.
#
# github: https://github.com/g-haine/scrimp

from scrimp import *
from scrimp.utils.mesh import set_verbose_gf
from itertools import zip_longest


def fsdj_eq():
    set_verbose_gf(0)

    # Init the distributed port-Hamiltonian system
    fsdj = DPHS("real")

    # Set the domain (using the built-in geometry `Rectangle`)
    fsdj.set_domain(Domain("Rectangle",{"L":1,"l":2,"h":3}))    ## Define the variables and their discretizations

    # Define State/s`)
    states = [
    State("fgds","df","scalar-field",df,0)
    ]

    # Define Co-State/s`)
    costates = [
    CoState("2","rfd",states[0],False)
    ]

    # Define Ports/s`)
    ports = [
    Port("p1","dfh","h","scalar-field",0,True,False,234)
    ]

    # Define parameters/s`)
    parameters = [
    Parameter("sda","sda","scalar-field","dfs","p1")
    ]

    # Define control_ports/s)
    control_ports = [
    Control_Port("Boundary control (bottom)","U_B","","Y_B","","scalar-field",10,"effort",0),
    Control_Port("Boundary control (right)","U_R","","Y_R","","scalar-field",11,"effort",0),
    Control_Port("Boundary control (top)","U_T","","Y_T","","scalar-field",12,"effort",0),
    Control_Port("Boundary control (left)","U_L","","Y_L","","scalar-field",13,"effort",0)
    ]

    # Define FEM/s`)
    FEMs = [
    FEM("fgds",,"CG"),
    FEM("p1",,"CG"),
    FEM("Boundary control (bottom)",,"CG"),
    FEM("Boundary control (right)",,"CG"),
    FEM("Boundary control (top)",,"CG"),
    FEM("Boundary control (left)",,"CG")
    ]

    for state, costate, param, fem, port, control_port in zip_longest(
            states, costates, parameters, FEMs, ports, control_ports
        ):
            # Add a state
            if state is not None:
                fsdj.add_state(state)

        # Add its co-state
            if costate is not None:
                fsdj.add_costate(costate)

        # Add a Finite Element Method to the `port`
            if fem is not None:
                fsdj.add_FEM(fem)

        # Add a (possibly space-varying) parameter to the `port`
            if param is not None:
                fsdj.add_parameter(param)

        # Add a resistive `port`
            if port is not None:
                fsdj.add_port(port)

        # Add a control `port` on the boundary (Neumann, thus position='effort' - default)
            if control_port is not None:
                fsdj.add_control_port(control_port)

    ## Set Hamiltonian
    fsdj.hamiltonian.set_name("")

    # Define Term/s`)
    terms = [
    Term("","",[],0)
    ]

    for term in terms:
        fsdj.hamiltonian.add_term(term)


    # Define Bricks/s`)
    bricks = [
    Brick("","",[],True,False,"constitutive",0)
    ]

    for brick in bricks:
        fsdj.add_brick(brick)


    # Define expression/s`)
    expressions = ["","",""]

    for control_port, expression in zip(control_ports, expressions):
        fsdj.set_control(control_port.get_name(), expression)


    # Define initial_value/s
    fsdj.set_initial_value("","")

    # Solve in time
    fsdj.set_time_scheme(ts_equation_type=,ts_type="",ksp_type="",pc_type="",pc_factor_mat_solver_type=,t_0=,t_f=,dt=,dt_save=,ts_adapt_dt_max=,init_step=)

    # Solve
    fsdj.solve()

    # Plot the Hamiltonian with the power supplied at the boundary
    fsdj.plot_Hamiltonian(save_figure=True)

    return fsdj

if __name__ == "__main__":
    fsdj = fsdj_eq()