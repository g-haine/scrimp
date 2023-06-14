# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2023 Ghislain Haine
#
# See the LICENSE file in the root directory for license information.
#
# github: https://github.com/g-haine/scrimp

"""
- file:             examples/heat.py
- author:           Ghislain Haine
- date:             05 dec. 2022
- last modified:    24 jan. 2023
- brief:            heat system
"""

def heat():
    """
    A structure-preserving discretization of the heat equation with boundary control
    
    Formulation co-energy (constitutive relation from variational derivative of H is substituted),
    Lyapunov L^2 functional, Div-Div, Mixed boundary condition on the Rectangle
    """
    
    from scrimp import set_verbose_gf
    set_verbose_gf(0)
    
    from scrimp import dpHs
    
    # Init the distributed port-Hamiltonian system
    heat = dpHs('real')
    
    # Set the domain (using the built-in geometry `Rectangle`)
    # Omega = 1, Gamma_Bottom = 10, Gamma_Right = 11, Gamma_Top = 12, Gamma_Left = 13
    heat.set_domain('Rectangle', {'L': 2., 'l': 1., 'h': 0.15})
    
    ## Define the variables and their discretizations
    
    # Add a state
    heat.add_state('T', 'Temperature', 'scalar-field')
    # Add its co-state
    heat.add_costate('T', 'Temperature', 'T', substituted=True)
    # Add a Finite Element Method to the `port`
    heat.add_FEM('T', 1, FEM='CG')
    # Add a (possibly space-varying) parameter to the `port`
    heat.add_parameter('rho', 'Mass density times heat capacity', 'scalar-field', '1.', 'T')
    
    # Add an algebraic port
    heat.add_port('Heat flux', 'f_Q', 'e_Q', 'vector-field')
    # Attach a FEM to it
    heat.add_FEM('Heat flux', 2, FEM='CG')
    # Add a (possibly space-varying) parameter to the `port`
    heat.add_parameter('Lambda', 'Heat conductivity', 'tensor-field', '[[1.,0.],[0.,1.]]', 'Heat flux')
    
    # Add a control port
    heat.add_control_port('Boundary control (bottom)', 'U_B', 'Temperature', 'Y_B', 'Normal heat flux', 'scalar-field', region=10, position='effort')
    heat.add_FEM('Boundary control (bottom)', 1, FEM='CG')
    heat.add_control_port('Boundary control (right)', 'U_R', 'Temperature', 'Y_R', 'Normal heat flux', 'scalar-field', region=11, position='effort')
    heat.add_FEM('Boundary control (right)', 1, FEM='CG')
    heat.add_control_port('Boundary control (top)', 'U_T', 'Temperature', 'Y_T', 'Normal heat flux', 'scalar-field', region=12, position='effort')
    heat.add_FEM('Boundary control (top)', 1, FEM='CG')
    heat.add_control_port('Boundary control (left)', 'U_L', 'Normal heat flux', 'Y_L', 'Temperature', 'scalar-field', region=13, position='flow')
    heat.add_FEM('Boundary control (left)', 1, FEM='CG')
    
    ## Define the Dirac structure via getfem `brick` = non-zero block matrix
    
    # Add the mass matrices from the left-hand side: the `flow` part of the Dirac structure
    heat.add_brick('M_T', 'T*rho*Test_T', [1], dt=True, position='flow')
    heat.add_brick('M_Q', 'f_Q.Test_f_Q', [1], position='flow')
    heat.add_brick('M_Y_B', 'Y_B*Test_Y_B', [10], position='flow')
    heat.add_brick('M_Y_R', 'Y_R*Test_Y_R', [11], position='flow')
    heat.add_brick('M_Y_T', 'Y_T*Test_Y_T', [12], position='flow')
    # Normal trace is imposed by Lagrange multiplier on the left side == the collocated output
    heat.add_brick('M_Y_L', 'U_L*Test_Y_L', [13], position='flow')
    
    # Add the matrices from the right-hand side: the `effort` part of the Dirac structure
    heat.add_brick('D', '-Div(e_Q)*Test_T', [1], position='effort')
    heat.add_brick('-D^T', 'T*Div(Test_f_Q)', [1], position='effort')
    heat.add_brick('B_B', '-U_B*Test_f_Q.Normal', [10], position='effort')
    heat.add_brick('B_R', '-U_R*Test_f_Q.Normal', [11], position='effort')
    heat.add_brick('B_T', '-U_T*Test_f_Q.Normal', [12], position='effort')
    # Normal trace is imposed by Lagrange multiplier on the left side == the collocated output
    heat.add_brick('B_L', '-Y_L*Test_f_Q.Normal', [13], position='effort')
    
    heat.add_brick('C_B', 'e_Q.Normal*Test_Y_B', [10], position='effort')
    heat.add_brick('C_R', 'e_Q.Normal*Test_Y_R', [11], position='effort')
    heat.add_brick('C_T', 'e_Q.Normal*Test_Y_T', [12], position='effort')
    # Normal trace is imposed by Lagrange multiplier on the left side == the collocated output
    heat.add_brick('C_L', 'e_Q.Normal*Test_Y_L', [13], position='effort')
    
    ## Define the constitutive relations as getfem `brick`
    
    # Fourier's law under implicit form - M_e_Q e_Q + CR_Q Q = 0
    heat.add_brick('-M_e_Q', '-e_Q.Test_e_Q', [1])
    heat.add_brick('CR_Q', 'f_Q.Lambda.Test_e_Q', [1])
    
    ## Initialize the problem
    
    # Set the control functions (automatic construction of bricks such that -M u + f(t) = 0)
    heat.set_control('Boundary control (bottom)', '0')
    heat.set_control('Boundary control (right)', '0')
    heat.set_control('Boundary control (top)', '0')
    heat.set_control('Boundary control (left)', '-0.5')
    
    # Set the initial data
    heat.set_initial_value('T', 'np.exp(-50*((x-1)*(x-1)+(y-0.5)*(y-0.5))**2)')
    
    ## Solve in time
    
    # Define the time scheme (default ts_type='cn', t_f=1, dt=0.01, etc.)
    heat.set_time_scheme(t_f=1, dt=0.01)
    
    # Solve
    heat.solve()
    
    ## Post-processing
    
    # Set the Hamiltonian
    heat.set_Hamiltonian_term('L^2-norm', '0.5*T*rho*T', [1])
    heat.set_Hamiltonian_name('Lyapunov formulation')
    
    # Plot the Hamiltonian with the power supplied at the boundary
    heat.plot_Hamiltonian()
    
    # Export solutions for ParaView
    heat.export_to_pv('T')
    heat.export_to_pv('e_Q')
    
    return heat # For consol use

if __name__ == "__main__": 
    heat = heat()
