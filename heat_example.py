# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2022 Ghislain Haine
#
# See the LICENSE file in the root directory for license information.
#
# github: https://github.com/g-haine/scrimp

from scrimp import *
from itertools import zip_longest


def heat_eq():
    # Init the distributed port-Hamiltonian system
    heat = DPHS("real")

    # Set the domain (using the built-in geometry `Rectangle`)
    heat.set_domain(Domain("Rectangle",{"L":2,"l":1,"h":0.1}))    ## Define the variables and their discretizations

    # Define State/s`)
    states = [
    State("T","Temperature","scalar-field",None,0)
    ]

    # Define Co-State/s`)
    costates = [
    CoState("T","Temperature",states[0],True)
    ]

    # Define Ports/s`)
    ports = [
    Port("Heat flux","f_Q","J_Q","vector-field",0,True,False,True,None)
    ]

    # Define parameters/s`)
    parameters = [
    Parameter("rho","Mass density","scalar-field","1","T"),
    Parameter("Lambda","Heat conductivity","tensor-field","[[1,0],[0,1]]","Heat flux"),
    Parameter("cv","Heat Capacity","scalar-field","1","T")
    ]

    # Define control_ports/s)
    control_ports = [
    Control_Port("B.C. (bottom)","U_B","Normal heat flux","Y_B","Temperature","scalar-field",10,"effort",0),
    Control_Port("B.C. (right)","U_R","Normal heat flux","Y_R","Temperature","scalar-field",11,"effort",0),
    Control_Port("B.C. (top)","U_T","Normal heat flux","Y_T","Temperature","scalar-field",12,"effort",0),
    Control_Port("B.C. (left)","U_L","Normal heat flux","Y_L","Temperature","scalar-field",13,"effort",0)
    ]

    # Define FEM/s`)
    FEMs = [
    FEM("T",2,"CG"),
    FEM("Heat flux",1,"CG"),
    FEM("B.C. (bottom)",1,"CG"),
    FEM("B.C. (right)",1,"CG"),
    FEM("B.C. (top)",1,"CG"),
    FEM("B.C. (left)",1,"CG")
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

    ## Set Hamiltonian
    heat.hamiltonian.set_name("Lyapunov formulation")

    # Define Term/s`)
    terms = [
    Term("L^2-norm","0.5*T*rho*T",[1],0)
    ]

    for term in terms:
        heat.hamiltonian.add_term(term)


    # Define Bricks/s`)
    bricks = [
    Brick("P_T","rho*cv*T*Test_T",[1],True,True,"flow",0),
    Brick("M_Q","f_Q.Test_f_Q",[1],True,False,"flow",0),
    Brick("M_Y_B","Y_B*Test_Y_B",[10],True,False,"flow",0),
    Brick("M_Y_R","Y_R*Test_Y_R",[11],True,False,"flow",0),
    Brick("M_Y_T","Y_T*Test_Y_T",[12],True,False,"flow",0),
    Brick("M_Y_L","Y_L*Test_Y_L",[13],True,False,"flow",0),
    Brick("G","J_Q.Grad(Test_T)",[1],True,False,"effort",0),
    Brick("-G^T","-Grad(T).Test_f_Q",[1],True,False,"effort",0),
    Brick("B_B","U_B*Test_T",[10],True,False,"effort",0),
    Brick("B_R","U_R*Test_T",[11],True,False,"effort",0),
    Brick("B_T","U_T*Test_T",[12],True,False,"effort",0),
    Brick("B_L","U_L*Test_T",[13],True,False,"effort",0),
    Brick("-B_B^T","-T*Test_Y_B",[10],True,False,"effort",0),
    Brick("-B_R^T","-T*Test_Y_R",[11],True,False,"effort",0),
    Brick("-B_T^T","-T*Test_Y_T",[12],True,False,"effort",0),
    Brick("-B_L^T","-T*Test_Y_L",[13],True,False,"effort",0),
    Brick("-M_Q","-J_Q.Test_J_Q",[1],True,False,"constitutive",0),
    Brick("L_Q","f_Q.Lambda.Test_J_Q",[1],True,False,"constitutive",0)
    ]

    for brick in bricks:
        heat.add_brick(brick)


    # Define expression/s`)
    expressions = ["-5","-7","3","3"]

    for control_port, expression in zip(control_ports, expressions):
        heat.set_control(control_port.get_name(), expression)


    # Define initial_value/s
    heat.set_initial_value("T","x*x+y*y+3*x-5*y")

    # Solve in time
    heat.set_time_scheme(t_0=0,t_f=1,dt=0.01,ts_type="bdf",ts_bdf_order=2,ts_adapt_dt_min=1e-6,dt_save=0.01)

    # Solve
    heat.solve()

    # Plot the Hamiltonian with the power supplied at the boundary
    heat.plot_Hamiltonian(save_figure=True)

    return heat

if __name__ == "__main__":
    heat = heat_eq()