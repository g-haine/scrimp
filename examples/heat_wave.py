# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2024 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""
- file:             sandbox/heat_wave.py
- authors:          Giuseppe Ferraro, Ghislain Haine
- date:             15 dec. 2022
- brief:            a 2D coupled heat-wave system
"""

import scrimp as S

hw = S.DPHS("real")

hw.set_domain(S.Domain("Concentric", {"R": 1., "r": 0.6, "h": 0.1}))

T = S.State("T", "Temperature", "scalar-field", region=1)
hw.add_state(T)
hw.add_costate(S.CoState("T", "Temperature", T, substituted=True))
V_T = S.FEM("T", 2)
hw.add_FEM(V_T)

Flux_Q = S.Port("Heat flux", "f_Q", "e_Q", "vector-field", dissipative=True, substituted=True, region=1)
hw.add_port(Flux_Q)
V_Q = S.FEM("Heat flux", 1)
hw.add_FEM(V_Q)

hw.add_control_port(S.Control_Port("Interface Heat", "U_T", "Heat flux", 
                                   "Y_T", "Temperature", "scalar-field", 
                                   region=10, position="effort"))
V_int_T = S.FEM("Interface Heat", 1)
hw.add_FEM(V_int_T)

p = S.State("p", "Velocity", "scalar-field", region=2)
hw.add_state(p)
hw.add_costate(S.CoState("p", "Velocity", p, substituted=True))
V_p = S.FEM("p", 1)
hw.add_FEM(V_p)

q = S.State("q", "Stress", "vector-field", region=2)
hw.add_state(q)
hw.add_costate(S.CoState("q", "Stress", q, substituted=True))
V_q = S.FEM("q", 2)
hw.add_FEM(V_q)

hw.add_control_port(S.Control_Port("Interface Wave", "U_w", "Velocity", 
                                   "Y_w", "Velocity", "scalar-field", 
                                   region=10, position="effort"))
V_int_w = S.FEM("Interface Wave", 1)
hw.add_FEM(V_int_w)

hw.add_control_port(S.Control_Port("Boundary Wave", "U_w_bnd", "Velocity", 
                                   "Y_w_bnd", "Velocity", "scalar-field", 
                                   region=20, position="flow"))
V_bnd_w = S.FEM("Boundary Wave", 1)
hw.add_FEM(V_bnd_w)

hw.add_parameter(S.Parameter("rho_CV", "Mass density times heat capacity", "scalar-field", "1.", "T"))
hw.add_parameter(S.Parameter("Lambda_inv", "Heat conductivity", "tensor-field", "[[1., 0.],[0., 1.]]", "Heat flux"))

hw.add_parameter(S.Parameter("rho", "Mass density", "scalar-field", "1.", "p"))
hw.add_parameter(S.Parameter("Kappa_inv", "Young\"s modulus", "tensor-field", "[[1., 0.],[0., 1.]]", "q"))

# === Heat grad-grad

hw.add_brick(S.Brick("M_T", "T*rho_CV*Test_T", [1], dt=True, position="flow"))
hw.add_brick(S.Brick("M_Q", "f_Q.Lambda_inv.Test_f_Q", [1], position="flow"))

hw.add_brick(S.Brick("M_Y_T", "Y_T*Test_Y_T", [10], position="flow"))

hw.add_brick(S.Brick("D_T", "f_Q.Grad(Test_T)", [1], position="effort"))
hw.add_brick(S.Brick("B_T", "-U_T*Test_T", [10], position="effort"))

hw.add_brick(S.Brick("-D_T^T", "-Grad(T).Test_f_Q", [1], position="effort"))

hw.add_brick(S.Brick("-B_T^T", "T*Test_Y_T", [10], position="effort"))

# === Wave div-div

hw.add_brick(S.Brick("M_p", "p*rho*Test_p", [2], dt=True, position="flow"))
hw.add_brick(S.Brick("M_q", "q.Kappa_inv.Test_q", [2], dt=True, position="flow"))

hw.add_brick(S.Brick("M_Y_w", "Y_w*Test_Y_w", [10], position="flow"))
hw.add_brick(S.Brick("M_Y_w_bnd", "U_w_bnd*Test_Y_w_bnd", [20], position="flow"))

hw.add_brick(S.Brick("D_w", "Div(q)*Test_p", [2], position="effort"))

hw.add_brick(S.Brick("-D_w^T", "-p*Div(Test_q)", [2], position="effort"))
hw.add_brick(S.Brick("B_w", "-U_w*Test_q.Normal", [10], position="effort"))
hw.add_brick(S.Brick("B_w_bnd", "-Y_w_bnd*Test_q.Normal", [20], position="effort"))

hw.add_brick(S.Brick("-B_w^T", "q.Normal*Test_Y_w", [10], position="effort"))
hw.add_brick(S.Brick("-B_w_bnd^T", "q.Normal*Test_Y_w_bnd", [20], position="effort"))

# === Gyrator

hw.set_control("Interface Heat", "Y_w") # Care must be taken, the Normal vector has an opposite sign
hw.set_control("Interface Wave", "Y_T")
hw.set_control("Boundary Wave", "0.")

hw.set_initial_value("T", "2.*np.exp(-25*((x-0.6)*(x-0.6)+y*y))")
hw.set_initial_value("p", "2.*np.exp(-25*((x-0.6)*(x-0.6)+y*y))")
hw.set_initial_value("q", "[0.,0.]")

hw.set_time_scheme(t_0=0., t_f=2., dt=0.01)

hw.solve()

hw.hamiltonian.add_term(S.Term("Lyapunov heat", "0.5*rho_CV*T*T", [1]))
hw.hamiltonian.add_term(S.Term("Kinetic energy", "0.5*rho*p*p", [2]))
hw.hamiltonian.add_term(S.Term("Potential energy", "0.5*q.Kappa_inv.q", [2]))

hw.plot_Hamiltonian()

hw.export_to_pv("T")
hw.export_to_pv("p")

# hw.spy_Dirac()
