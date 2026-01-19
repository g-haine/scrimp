# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""
- file:             examples/published/LHMNC_2024/heat.py
- authors:          Giuseppe Ferraro, Michel Fournié, Ghislain Haine
- date:             19 jan. 2024
- brief:            2D example of the heat equation
"""

from scrimp import *

def heat():
    """
    A structure-preserving discretization of the heat equation with heat flux 
    boundary controls, as a port-Hamiltonian system.
    
    Discretisation use the grad-grad formulation, thus prescribed explicitly the
    opposite of the normal heat flux.
    
    The L^\infty([0,t_f],L^2(\Omega)) error between the exact solution and the 
    computed solution is given.
    
    This script can be used to reproduce the examples shown in the paper:
    @article{Ferraro_2024,
        title={{Simulation and control of interactions in multi-physics, a Python package for port-Hamiltonian systems}},
        volume={58},
        ISSN={2405-8963},
        DOI={10.1016/j.ifacol.2024.08.267},
        number={6},
        journal={IFAC-PapersOnLine},
        publisher={Elsevier BV},
        author={Ferraro, Giuseppe and Fournié, Michel and Haine, Ghislain},
        year={2024},
        pages={119--124}
    }
    
    Example of use: to run the heat problem on 2 cores:
        `mpirun -n 2 python heat.py`
    
    See the file `scrimp.yml` in the folder of this script to know the exact version of each library used for the results obtained in the paper.
    
    Returns:
        the DPHS object
    """

    heat = DPHS("real")
    
    heat.set_domain( 
        Domain("Rectangle", {"L": 2.0,"l": 1.0, "h": 0.1}) )
    
    states = [ State("T", "Temperature", "scalar-field") ]
    costates = [ CoState("T", "Temperature", states[0], substituted=True) ]
    
    ports = [ Port("Heat flux", "f_Q", "J_Q", "vector-field") ]
    
    control_ports = [
        Control_Port("B.C.(bottom)", "U_B", "Normal heat flux", "Y_B", "Temperature", "scalar-field", 10),
        Control_Port("B.C.(right)", "U_R", "Normal heat flux", "Y_R", "Temperature", "scalar-field", 11),
        Control_Port("B.C.(top)", "U_T", "Normal heat flux", "Y_T", "Temperature", "scalar-field", 12),
        Control_Port("B.C.(left)", "U_L", "Normal heat flux", "Y_L", "Temperature", "scalar-field", 13),
        ]
    
    FEMs = [
        FEM("T", 2, "CG"),
        FEM("Heat flux", 1, "DG"),
        FEM("B.C.(bottom)", 1, "CG"),
        FEM("B.C.(right)", 1, "CG"),
        FEM("B.C.(top)", 1, "CG"),
        FEM("B.C.(left)", 1, "CG"),
        ]
        
    parameters = [
        Parameter("rho", "(Mass density)x(Heat capacity)", "scalar-field", "1.", "T"),
        Parameter("cv", "heat capacity", "scalar-field", "1.", "T"),
        Parameter("lambda", "heat conductivity", "tensor-field", "[[1.,0.],[0.,1.]]", "Heat flux"),
        ]
    
    for (state,costate) in zip(states,costates):
        heat.add_state(state)
        heat.add_costate(costate)
    for port in ports:
        heat.add_port(port)
    for control_port in control_ports:
        heat.add_control_port(control_port)
    for fem in FEMs:
        heat.add_FEM(fem)
    for param in parameters:
        heat.add_parameter(param)
    
    heat.hamiltonian.set_name("Lyapunov formulation")
    terms = [ Term("L^2-norm", "0.5*rho*cv*T*T", [1], 0) ]
    for term in terms:
         heat.hamiltonian.add_term(term)

    bricks = [
        Brick("P_T", "rho*cv*T*Test_T", [1], dt=True, position="flow"),
        Brick("M_Q", "f_Q.Test_f_Q", [1], position="flow"),
        Brick("M_Y_B", "Y_B*Test_Y_B", [10], position="flow"),
        Brick("M_Y_R", "Y_R*Test_Y_R", [11], position="flow"),
        Brick("M_Y_T", "Y_T*Test_Y_T", [12], position="flow"),
        Brick("M_Y_L", "Y_L*Test_Y_L", [13], position="flow"),
        
        Brick("G", "J_Q.Grad(Test_T)", [1], position="effort"),
        Brick("-G^T", "-Grad(T).Test_f_Q", [1], position="effort"),
        
        Brick("B_B", "U_B*Test_T", [10], position="effort"),
        Brick("B_R", "U_R*Test_T", [11], position="effort"),
        Brick("B_T", "U_T*Test_T", [12], position="effort"),
        Brick("B_L", "U_L*Test_T", [13], position="effort"),
        Brick("-B_B^T", "-T*Test_Y_B", [10], position="effort"),
        Brick("-B_R^T", "-T*Test_Y_R", [11], position="effort"),
        Brick("-B_T^T", "-T*Test_Y_T", [12], position="effort"),
        Brick("-B_L^T", "-T*Test_Y_L", [13], position="effort"),
        
        Brick("-M_Q","-J_Q.Test_J_Q",[1],position="constitutive"),
        Brick("L_Q","f_Q.lambda.Test_J_Q",[1],position="constitutive")
        ]
    for brick in bricks:
         heat.add_brick(brick)
         
    expressions = ["-5", "-7", "3", "3"]
    for control_port, expression in zip(control_ports, expressions):
         heat.set_control(control_port.get_name(), expression)
    
    heat.set_initial_value("T", "x*x+y*y+3*x-5*y")
    
    heat.set_time_scheme(t_f=1., dt=0.01, init_step=False)
    heat.solve()
    
    heat.plot_Hamiltonian()
    err_L2 = heat.get_quantity("pow(T-4*t-x*x-y*y-3*x+5*y,2)")
    from math import sqrt
    print("Maximal L^2-error over time:", sqrt(max(err_L2)))
    
    return heat

if __name__ == "__main__":
    heat = heat()
