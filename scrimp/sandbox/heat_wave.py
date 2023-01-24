# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2023 Ghislain Haine
#
# See the LICENSE file in the root directory for license information.
#
# github: https://github.com/g-haine/scrimp

"""
- file:             sandbox/heat_wave.py
- author:           Ghislain Haine
- date:             15 dec. 2022
- last modified:    15 dec. 2022
- brief:            a 2D coupled heat-wave system

!TO DO: Correct Lagrange multiplier assignement
"""

from scrimp import set_verbose_gf
set_verbose_gf(0)

from scrimp.dpHs import dpHs
    
hw = dpHs('real')

hw.set_domain('Concentric', {'R': 1., 'r': 0.6, 'h': 0.1})

hw.add_state('T', 'Temperature', 'scalar-field', region=1)
hw.add_costate('T', 'Temperature', 'T', substituted=True)
hw.add_FEM('T', 1, FEM='CG')

hw.add_port('Q', 'J_Q', 'J_Q', 'vector-field', 
            algebraic=True, substituted=True, region=1)
hw.add_FEM('Q', 1, FEM='CG')

hw.add_control_port('Interface Heat', 'U_T', 'Heat flux', 'Y_T', 'Temperature', 'scalar-field', 
                    region=10, position='flow')
hw.add_FEM('Interface Heat', 1, FEM='CG')

hw.add_state('p', 'Velocity', 'scalar-field', region=2)
hw.add_costate('p', 'Velocity', 'p', substituted=True)
hw.add_FEM('p', 1, FEM='CG')

hw.add_state('q', 'Stress', 'vector-field', region=2)
hw.add_costate('q', 'Stress', 'q', substituted=True)
hw.add_FEM('q', 1, FEM='CG')

hw.add_control_port('Interface Wave', 'U_w', 'Velocity', 'Y_w', 'Velocity', 'scalar-field', 
                    region=10, position='effort')
hw.add_FEM('Interface Wave', 1, FEM='CG')

hw.add_control_port('Boundary Wave', 'U_w_bnd', 'Velocity', 'Y_w_bnd', 'Velocity', 'scalar-field', 
                    region=20, position='effort')
hw.add_FEM('Boundary Wave', 1, FEM='CG')

hw.add_parameter('rho_CV', 'Mass density times heat capacity', 'scalar-field', '1.', 'T')
hw.add_parameter('Lambda_inv', 'Heat conductivity', 'tensor-field', '[[1., 0.],[0., 1.]]', 'Q')

hw.add_parameter('rho', 'Mass density', 'scalar-field', '1.', 'T')
hw.add_parameter('Kappa_inv', 'Young\'s modulus', 'tensor-field', '[[1., 0.],[0., 1.]]', 'q')

hw.add_brick('M_T', 'T*rho_CV*Test_T', [1], dt=True, position='flow')
hw.add_brick('M_Q', 'J_Q.Lambda_inv.Test_J_Q', [1], position='flow')
hw.add_brick('M_Y_T', 'U_T*Test_Y_T', [10], position='flow') # Y_T = Lagrange multiplier
hw.add_brick('M_p', 'p*rho*Test_p', [2], dt=True, position='flow')
hw.add_brick('M_q', 'q.Kappa_inv.Test_q', [2], dt=True, position='flow')
hw.add_brick('M_Y_w', 'Y_w*Test_Y_w', [10], position='flow')
hw.add_brick('M_Y_w_bnd', 'Y_w_bnd*Test_Y_w_bnd', [20], position='flow')

hw.add_brick('D_T', '-Div(J_Q)*Test_T', [1], position='effort')

hw.add_brick('-D_T^T', 'T*Div(Test_J_Q)', [1], position='effort')
hw.add_brick('B_T', '-Y_T*Test_J_Q.Normal', [10], position='effort') # Y_T = Lagrange multiplier

hw.add_brick('-B_T^T', 'J_Q.Normal*Test_Y_T', [10], position='effort')

hw.add_brick('D_w', 'Div(q)*Test_p', [2], position='effort')

hw.add_brick('-D_w^T', '-p*Div(Test_q)', [2], position='effort')
hw.add_brick('B_w', 'U_w*Test_q.Normal', [10], position='effort')
hw.add_brick('B_w_bnd', 'U_w_bnd*Test_q.Normal', [20], position='effort')

hw.add_brick('-B_w^T', '-q.Normal*Test_Y_w', [10], position='effort')
hw.add_brick('-B_w_bnd^T', '-q.Normal*Test_Y_w_bnd', [20], position='effort')

hw.set_control('Interface Heat', 'Y_w')
hw.set_control('Interface Wave', '-Y_T')
hw.set_control('Boundary Wave', '0.')

hw.set_initial_value('T', '5.*np.exp(-20*((x-0.6)*(x-0.6)+y*y))')
hw.set_initial_value('p', '5.*np.exp(-20*((x-0.6)*(x-0.6)+y*y))')
hw.set_initial_value('q', '[0.,0.]')

hw.set_time_scheme(ts_type='cn',
                   t_0=0., t_f=0.5, dt=0.01)

hw.solve()

hw.set_Hamiltonian_term('Lyapunov heat', '0.5*rho_CV*T*T', [1])
hw.set_Hamiltonian_term('Kinetic energy', '0.5*rho*p*p', [2])
hw.set_Hamiltonian_term('Potential energy', '0.5*q.Kappa_inv.q', [2])

hw.plot_Hamiltonian()

hw.spy_Dirac()

hw.export_to_pv('T')
hw.export_to_pv('p')
