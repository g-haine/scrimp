# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2023 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""
- file:             sandbox/wave_1D.py
- authors:          Giuseppe Ferraro, Ghislain Haine
- date:             14 dec. 2022
- brief:            1D wave equation
"""

import scrimp as S
    
wave = S.DPHS("real")

domain = S.Domain("Interval", {"L": 1., "h": 0.01})
wave.set_domain(domain)

alpha_q = S.State("q", "Strain", "scalar-field")
alpha_p = S.State("p", "Linear momentum", "scalar-field")
wave.add_state(alpha_q)
wave.add_state(alpha_p)

e_q = S.CoState("e_q", "Stress", alpha_q)
e_p = S.CoState("e_p", "Velocity", alpha_p)
wave.add_costate(e_q)
wave.add_costate(e_p)

left_end = S.Control_Port("Boundary control (left)", "U_L", "Normal force", "Y_L", "Velocity", "scalar-field", region=10)
right_end = S.Control_Port("Boundary control (right)", "U_R", "Normal force", "Y_R", "Velocity", "scalar-field", region=11)
wave.add_control_port(left_end)
wave.add_control_port(right_end)

V_q = S.FEM("q", 2)
V_p = S.FEM("p", 1)
V_L = S.FEM("Boundary control (left)", 1)
V_R = S.FEM("Boundary control (right)", 1)
wave.add_FEM(V_q)
wave.add_FEM(V_p)
wave.add_FEM(V_L)
wave.add_FEM(V_R)

T = S.Parameter("T", "Young's modulus", "scalar-field", "1", "q")
rho = S.Parameter("rho", "Mass density", "scalar-field", "1 + x*(1-x)", "p")
wave.add_parameter(T)
wave.add_parameter(rho)

bricks = [
    # M matrix, on the flow side
    S.Brick("M_q", "q * Test_q", [1], dt=True, position="flow"),
    S.Brick("M_p", "p * Test_p", [1], dt=True, position="flow"),
    S.Brick("M_Y_L", "Y_L * Test_Y_L", [10], position="flow"),
    S.Brick("M_Y_R", "Y_R * Test_Y_R", [11], position="flow"),
    
    # J matrix, on the effort side
    S.Brick("D", "Grad(e_p) * Test_q", [1], position="effort"),

    S.Brick("-D^T", "-e_q * Grad(Test_p)", [1], position="effort"),
    S.Brick("B_L", "-U_L * Test_p", [10], position="effort"),
    S.Brick("B_R", "U_R * Test_p", [11], position="effort"),

    S.Brick("-B_L^T", "e_p * Test_Y_L", [10], position="effort"),
    S.Brick("-B_R^T", "-e_p * Test_Y_R", [11], position="effort"),
    
    # Constitutive relations
    S.Brick("-M_e_q", "-e_q * Test_e_q", [1]),
    S.Brick("CR_q", "q*T * Test_e_q", [1]),

    S.Brick("-M_e_p", "-e_p * Test_e_p", [1]),
    S.Brick("CR_p", "p/rho * Test_e_p", [1]),
    ]

for brick in bricks:
    wave.add_brick(brick)

expression_left = "-sin(2*pi*t)"
expression_right = "0."
wave.set_control("Boundary control (left)", expression_left)
wave.set_control("Boundary control (right)", expression_right)

q_init = "2.*np.exp(-50.*(x-0.5)*(x-0.5))"
p_init = "0."
wave.set_initial_value("q", q_init)
wave.set_initial_value("p", p_init)

wave.solve()

wave.hamiltonian.set_name("Mechanical energy")
terms = [
    S.Term("Kinetic energy", "0.5*p*p/rho", [1]),
    S.Term("Potential energy", "0.5*q*T*q", [1]),
]
for term in terms:
    wave.hamiltonian.add_term(term)

wave.plot_Hamiltonian()

