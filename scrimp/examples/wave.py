# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2022 Ghislain Haine
#
# See the LICENSE file in the root directory for license information.
#
# github: https://github.com/g-haine/scrimp

"""
- file:             examples/wave.py
- author:           Ghislain Haine
- date:             22 nov. 2022
- last modified:    12 dec. 2022
- brief:            wave system
"""

def wave():
    """
    A structure-preserving discretization of the wave equation with boundary control
    
    Formulation DAE (energy/co-energy), Grad-Grad, Mixed boundary condition on the Rectangle
    """
    
    from scrimp import set_verbose_gf
    set_verbose_gf(0)
    
    from scrimp.dpHs import dpHs
    
    # Init the distributed port-Hamiltonian system
    wave = dpHs('real')
    
    # Set the domain (using the built-in geometry `Rectangle`)
    # Omega = 1, Gamma_Bottom = 10, Gamma_Right = 11, Gamma_Top = 12, Gamma_Left = 13
    wave.set_domain('Rectangle', {'L': 2., 'l': 1., 'h': 0.15})
    
    ## Define the variables and their discretizations
    
    # Add a state
    wave.add_state('q', 'Strain', 'vector-field')
    # Add its co-state
    wave.add_costate('e_q', 'Stress', 'q')
    # Add a Finite Element Method to the `port`
    wave.add_FEM('q', 1, FEM='DG')
    # Add a (possibly space-varying) parameter to the `port`
    wave.add_parameter('T', 'Young\'s modulus', 'tensor-field', '[[5+x,x*y],[x*y,2+y]]', 'q')
    
    # Add a state
    wave.add_state('p', 'Linear momentum', 'scalar-field')
    # Add its co-state
    wave.add_costate('e_p', 'Velocity', 'p')
    # Add a Finite Element Method to the `port`
    wave.add_FEM('p', 2, FEM='CG')
    # Add a (possibly space-varying) parameter to the `port`
    wave.add_parameter('rho', 'Mass density', 'scalar-field', '3-x', 'p')
    
    # Add a resistive `port`
    wave.add_port('Damping', 'f_r', 'e_r', 'scalar-field')
    # Add a FEM on it
    wave.add_FEM('Damping', 1, FEM='DG')
    # Attach a damping parameter to it
    wave.add_parameter('nu', 'viscosity', 'scalar-field', '0.05*x', 'Damping')
    
    # Add a control `port` on the bottom part of the boundary (Neumann, thus position='effort' - default)
    wave.add_control_port('Boundary control (bottom)', 'U_B', 'Normal force', 'Y_B', 'Velocity trace', 'scalar-field', region=10)
    wave.add_FEM('Boundary control (bottom)', 1, FEM='DG')
    # Add a control `port` on the right part of the boundary (Neumann, thus position='effort' - default)
    wave.add_control_port('Boundary control (right)', 'U_R', 'Normal force', 'Y_R', 'Velocity trace', 'scalar-field', region=11)
    wave.add_FEM('Boundary control (right)', 1, FEM='DG')
    # Add a control `port` on the top part of the boundary (Neumann, thus position='effort' - default)
    wave.add_control_port('Boundary control (top)', 'U_T', 'Normal force', 'Y_T', 'Velocity trace', 'scalar-field', region=12)
    wave.add_FEM('Boundary control (top)', 1, FEM='DG')
    # Add a control `port` on the left part of the boundary (Dirichlet, thus position='flow')
    wave.add_control_port('Boundary control (left)', 'U_L', 'Velocity trace', 'Y_L', 'Normal force', 'scalar-field', region=13, position='flow')
    wave.add_FEM('Boundary control (left)', 1, FEM='DG')
    
    # Set the Hamiltonian (can be done later, even after solve)
    wave.set_Hamiltonian_term('Potential energy', '0.5*q.T.q', [1])
    wave.set_Hamiltonian_term('Kinetic energy', '0.5*p*p/rho', [1])
    wave.set_Hamiltonian_name('Mechanical energy')
    
    ## Define the Dirac structure via getfem `brick` = non-zero block matrix
    
    # Add the mass matrices from the left-hand side: the `flow` part of the Dirac structure
    wave.add_brick('M_q', 'q.Test_q', [1], dt=True, position='flow')
    wave.add_brick('M_p', 'p*Test_p', [1], dt=True, position='flow')
    wave.add_brick('M_r', 'f_r*Test_f_r', [1], position='flow')
    wave.add_brick('M_Y_B', 'Y_B*Test_Y_B', [10], position='flow')
    wave.add_brick('M_Y_R', 'Y_R*Test_Y_R', [11], position='flow')
    wave.add_brick('M_Y_T', 'Y_T*Test_Y_T', [12], position='flow')
    # The Dirichlet term is applied via Lagrange multiplier == the colocated output
    wave.add_brick('M_Y_L', 'U_L*Test_Y_L', [13], position='flow')
    
    # Add the matrices from the right-hand side: the `effort` part of the Dirac structure
    wave.add_brick('D', 'Grad(e_p).Test_q', [1], position='effort')
    
    wave.add_brick('-D^T', '-e_q.Grad(Test_p)', [1], position='effort')
    wave.add_brick('I_r', 'e_r*Test_p', [1], position='effort')
    wave.add_brick('B_B', 'U_B*Test_p', [10], position='effort')
    wave.add_brick('B_R', 'U_R*Test_p', [11], position='effort')
    wave.add_brick('B_T', 'U_T*Test_p', [12], position='effort')
    # The Dirichlet term is applied via Lagrange multiplier == the colocated output
    wave.add_brick('B_L', 'Y_L*Test_p', [13], position='effort')
    
    wave.add_brick('-I_r^T', '-e_p*Test_f_r', [1], position='effort')
    
    wave.add_brick('C_B', '-e_p*Test_Y_B', [10], position='effort')
    wave.add_brick('C_R', '-e_p*Test_Y_R', [11], position='effort')
    wave.add_brick('C_T', '-e_p*Test_Y_T', [12], position='effort')
    wave.add_brick('C_L', '-e_p*Test_Y_L', [13], position='effort')
    
    ## Define the constitutive relations as getfem `brick`
    
    # Hooke's law under implicit form - M_e_q e_q + CR_q q = 0
    wave.add_brick('-M_e_q', '-e_q.Test_e_q', [1])
    wave.add_brick('CR_q', 'q.T.Test_e_q', [1])
    
    # Linear momentum definition under implicit form - M_e_p e_p + CR_p p = 0
    wave.add_brick('-M_e_p', '-e_p*Test_e_p', [1])
    wave.add_brick('CR_p', 'p/rho*Test_e_p', [1])
    
    # Linear viscous fluid damping - M_e_r e_r + CR_r f_r = 0
    wave.add_brick('-M_e_r', '-e_r*Test_e_r', [1])
    wave.add_brick('CR_r', 'nu*f_r*Test_e_r', [1])
    
    ## Initialize the problem
    
    # Set the control functions (automatic construction of bricks such that -M_u u + f(t) = 0)
    wave.set_control('Boundary control (bottom)', '0.')
    wave.set_control('Boundary control (right)', '0.')
    wave.set_control('Boundary control (top)', '0.')
    wave.set_control('Boundary control (left)', '0.1*sin(2.*t)*sin(4*pi*y)')
    
    # Set the initial data
    wave.set_initial_value('q', '[0., 0.]')
    wave.set_initial_value('p', '2.72**(-20*((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5)))')
    
    ## Solve in time
    
    # Define the time scheme (default ts_type='cn', t_f=1, dt=0.01, etc.)
    wave.set_time_scheme(t_f=2., dt=0.01, dt_save=0.01)
    
    # Solve
    wave.solve()
    
    ## Post-processing
    
    # Plot the Hamiltonian with the power supplied at the boundary
    wave.plot_Hamiltonian(with_powers=True)
    
    # Export solutions for ParaView
    wave.export_to_pv('q')
    wave.export_to_pv('p')
    
    # Plot the matrices representing the Dirac structure
    wave.spy_Dirac()
    
    return wave # For consol use

if __name__ == "__main__": 
    wave = wave()