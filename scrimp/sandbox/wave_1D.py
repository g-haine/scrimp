# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2023 Ghislain Haine
#
# See the LICENSE file in the root directory for license information.
#
# github: https://github.com/g-haine/scrimp

"""
- file:             sandbox/wave_1D.py
- author:           Ghislain Haine
- date:             14 dec. 2022
- last modified:    14 dec. 2022
- brief:            1D wave equation
"""

# import tracemalloc
# tracemalloc.start()

from scrimp import *
from scrimp.utils.mesh import set_verbose_gf
set_verbose_gf(0)
    
wave = DPHS("real")

domain = Domain("Interval", {"L": 1., "h": 0.01})
wave.set_domain(domain)

alpha_q = State("q", "Strain", "scalar-field")
alpha_p = State("p", "Linear momentum", "scalar-field")
wave.add_state(alpha_q)
wave.add_state(alpha_p)

e_q = CoState("e_q", "Stress", alpha_q)
e_p = CoState("e_p", "Velocity", alpha_p)
wave.add_costate(e_q)
wave.add_costate(e_p)

left_end = Control_Port("Boundary control (left)", "U_L", "Normal force", "Y_L", "Velocity", "scalar-field", region=10)
right_end = Control_Port("Boundary control (right)", "U_R", "Normal force", "Y_R", "Velocity", "scalar-field", region=11)
wave.add_control_port(left_end)
wave.add_control_port(right_end)

wave.hamiltonian.set_name("Mechanical energy")
terms = [
    Term("Potential energy", "0.5*q*T*q", [1]),
    Term("Kinetic energy", "0.5*p*p/rho", [1]),
]

for term in terms:
    wave.hamiltonian.add_term(term)

V_q = FEM("q", 1)
V_p = FEM("p", 2)
V_L = FEM("Boundary control (left)", 1)
V_R = FEM("Boundary control (right)", 1)
wave.add_FEM(V_q)
wave.add_FEM(V_p)
wave.add_FEM(V_L)
wave.add_FEM(V_R)

T = Parameter("T", "Young\'s modulus", "scalar-field", "3", "q")
rho = Parameter("rho", "Mass density", "scalar-field", "2", "p")
wave.add_parameter(T)
wave.add_parameter(rho)

bricks = [
    # M matrix, on the flow side
    Brick("M_q", "q * Test_q", [1], dt=True, position="flow"),
    Brick("M_p", "p * Test_p", [1], dt=True, position="flow"),
    Brick("M_Y_L", "Y_L * Test_Y_L", [10], position="flow"),
    Brick("M_Y_R", "Y_R * Test_Y_R", [11], position="flow"),
    
    # J matrix, on the effort side
    Brick("D", "Grad(e_p) * Test_q", [1], position="effort"),

    Brick("-D^T", "-e_q * Grad(Test_p)", [1], position="effort"),
    Brick("B_L", "-U_L * Test_p", [10], position="effort"),
    Brick("B_R", "U_R * Test_p", [11], position="effort"),

    Brick("-B_L^T", "e_p * Test_Y_L", [10], position="effort"),
    Brick("-B_R^T", "-e_p * Test_Y_R", [11], position="effort"),
    
    # Constitutive relations
    Brick("-M_e_q", "-e_q * Test_e_q", [1]),
    Brick("CR_q", "q*T * Test_e_q", [1]),

    Brick("-M_e_p", "-e_p * Test_e_p", [1]),
    Brick("CR_p", "p/rho * Test_e_p", [1]),
    ]

for brick in bricks:
    wave.add_brick(brick)

## Initialize the problem
expression_left = "sin(2*pi*t)"
expression_right = "0."
wave.set_control("Boundary control (left)", expression_left)
wave.set_control("Boundary control (right)", expression_right)

q_init = "np.exp(-50*(x-0.5)*(x-0.5))"
p_init = "0."
wave.set_initial_value("q", q_init)
wave.set_initial_value("p", p_init)



# snapshot1 = tracemalloc.take_snapshot()


wave.set_time_scheme()#ts_type="cn")
wave.solve()

wave.plot_Hamiltonian()



# snapshot2 = tracemalloc.take_snapshot()



# top_stats = snapshot2.compare_to(snapshot1, 'lineno')

# print("[ Top 50 differences ]")
# for stat in top_stats[:50]:
#     print(stat)
