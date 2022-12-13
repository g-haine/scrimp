# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2022 Ghislain Haine
#
# See the LICENSE file in the root directory for license information.
#
# github: https://github.com/g-haine/scrimp

"""
- file:             sandbox/Burgers.py
- author:           Ghislain Haine
- date:             06 dec. 2022
- last modified:    13 dec. 2022
- brief:            Burgers / 1D Euler equation
"""

def Burgers():

    from scrimp import set_verbose_gf
    set_verbose_gf(0)
    
    from scrimp.dpHs import dpHs
    
    burgers = dpHs('real')
    
    burgers.set_domain('Interval', {'L': 1., 'h': 0.01})
    
    burgers.add_state('rho', 'density', 'scalar-field')
    burgers.add_costate('h', 'enthalpy', 'rho', substituted=False)
    burgers.add_FEM('rho', 1, FEM='CG')
    
    burgers.add_state('v', 'Velocity', 'scalar-field')
    burgers.add_costate('p', 'Linear momentum', 'v', substituted=False)
    burgers.add_FEM('v', 1, FEM='CG')
    
    burgers.add_control_port('Boundary control (left)', 'U_L', 'description_control', 'Y_L', 'description_observation', 'scalar-field', region=10)
    burgers.add_FEM('Boundary control (left)', 1, FEM='CG')
    burgers.add_control_port('Boundary control (right)', 'U_R', 'description_control', 'Y_R', 'description_observation', 'scalar-field', region=11)
    burgers.add_FEM('Boundary control (right)', 1, FEM='CG')
    
    burgers.add_port('Dissipation', 'f_d', 'e_d', 'scalar-field')
    burgers.add_FEM('Dissipation', 1, FEM='CG')
    burgers.add_parameter('nu', 'Viscosity', 'scalar-field', '0.01', 'Dissipation')
    
    burgers.add_brick('M_rho', 'rho*Test_rho', [1], dt=True, position='flow')
    burgers.add_brick('M_v', 'v*Test_v', [1], dt=True, position='flow')
    burgers.add_brick('M_f_d', 'f_d*Test_f_d', [1], position='flow')
    burgers.add_brick('M_Y_L', 'Y_L*Test_Y_L', [10], position='flow')
    burgers.add_brick('M_Y_R', 'Y_R*Test_Y_R', [11], position='flow')
    
    burgers.add_brick('D', 'p*Grad(Test_rho)', [1], position='effort')
    burgers.add_brick('B_L', '-U_L*Test_rho', [10], position='effort')
    burgers.add_brick('B_R', '-U_R*Test_rho', [11], position='effort')
    
    burgers.add_brick('-D^T', '-Grad(h)*Test_v', [1], position='effort')
    burgers.add_brick('R', 'e_d*Test_v', [1], position='effort')
    
    burgers.add_brick('-R^T', '-p*Test_f_d', [1], position='effort')
    
    burgers.add_brick('C_L', 'h*Test_Y_L', [10], position='effort')
    burgers.add_brick('C_R', 'h*Test_Y_R', [11], position='effort')
    
    burgers.add_brick('M_h', '-h*Test_h', [1])
    burgers.add_brick('CR_h', '0.5*v*v*Test_h', [1], linear=False)
    
    burgers.add_brick('M_p', '-p*Test_p', [1])
    burgers.add_brick('CR_p', 'rho*v*Test_p', [1], linear=False)
    
    burgers.add_brick('M_e_d', '-e_d*Test_e_d', [1])
    burgers.add_brick('CR_d', '-nu*Grad(f_d/rho).Grad(Test_e_d)', [1], linear=False)
    burgers.add_brick('CR_d_b', 'nu*Grad(f_d/rho).Normal*Test_e_d', [10,11], linear=False)
    
    burgers.set_control('Boundary control (left)', '0.')
    burgers.set_control('Boundary control (right)', '0.')
    
    burgers.set_initial_value('rho', '10.')
    burgers.set_initial_value('v', 'np.exp(-50.*(x-0.5)*(x-0.5))')
    
    burgers.set_time_scheme(ts_type='bdf', bdf_order=2, ksp_type='gmres', 
                            pc_type='jacobi', #pc_factor_mat_solver_type='superlu',
                            t_0=0., t_f=0.17, dt=0.01, dt_save=0.01, init_step=True)
    burgers.solve()
    
    burgers.set_Hamiltonian_term('Kinetic energy', '0.5*rho*v*v', [1])
    
    burgers.plot_Hamiltonian(with_powers=True)
    
    burgers.export_to_pv('rho')
    burgers.export_to_pv('v')
    
    return burgers
    
if __name__ == "__main__": 
    burgers = Burgers()
    