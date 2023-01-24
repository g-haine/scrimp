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

from scrimp.dpHs import dpHs
    
wave = dpHs('real')

wave.set_domain('Interval', {'L': 1., 'h': 0.01})

wave.add_state('q', 'Strain', 'scalar-field')
wave.add_costate('e_q', 'Stress', 'q')

wave.add_state('p', 'Linear momentum', 'scalar-field')
wave.add_costate('e_p', 'velocity', 'p')

wave.add_control_port('Boundary control (left)', 'U_L', 'Normal force', 'Y_L', 'Velocity', 'scalar-field', region=10)
wave.add_control_port('Boundary control (right)', 'U_R', 'Normal force', 'Y_R', 'Velocity', 'scalar-field', region=11)

wave.add_FEM('q', 2)
wave.add_FEM('p', 1)
wave.add_FEM('Boundary control (left)', 1)
wave.add_FEM('Boundary control (right)', 1)

wave.add_brick('M_q', 'q*Test_q', [1], dt=True, position='flow')
wave.add_brick('M_p', 'p*Test_p', [1], dt=True, position='flow')
wave.add_brick('M_Y_L', 'Y_L*Test_Y_L', [10], position='flow')
wave.add_brick('M_Y_R', 'Y_R*Test_Y_R', [11], position='flow')

wave.add_brick('D', 'Grad(e_p)*Test_q', [1], position='effort')

wave.add_brick('-D^T', '-e_q*Grad(Test_p)', [1], position='effort')
wave.add_brick('B_L', 'U_L*Test_p', [10], position='effort')
wave.add_brick('B_R', 'U_R*Test_p', [11], position='effort')

wave.add_brick('-B_L^T', '-e_p*Test_Y_L', [10], position='effort')
wave.add_brick('-B_R^T', '-e_p*Test_Y_R', [11], position='effort')

wave.add_parameter('T', 'Young\'s modulus', 'tensor-field', '1', 'q')
wave.add_parameter('rho', 'Mass density', 'scalar-field', '1 + x*(1-x)', 'p')

wave.add_brick('-M_e_q', '-e_q*Test_e_q', [1])
wave.add_brick('CR_q', 'q*T*Test_e_q', [1])

wave.add_brick('-M_e_p', '-e_p*Test_e_p', [1])
wave.add_brick('CR_p', 'p/rho*Test_e_p', [1])

wave.set_control('Boundary control (left)', 'sin(2*pi*t)')
wave.set_control('Boundary control (right)', '0.')

wave.set_initial_value('q', '0.')
wave.set_initial_value('p', 'np.exp(-50.*(x-0.5)*(x-0.5))')

wave.solve()

wave.set_Hamiltonian_term('Kinetic energy', '0.5*p*p/rho', [1])
wave.set_Hamiltonian_term('Potential energy', '0.5*q*T*q', [1])

wave.plot_Hamiltonian()
