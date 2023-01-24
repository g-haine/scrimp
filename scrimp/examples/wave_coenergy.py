# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2023 Ghislain Haine
#
# See the LICENSE file in the root directory for license information.
#
# github: https://github.com/g-haine/scrimp

"""
- file:             examples/wave_coenergy.py
- author:           Ghislain Haine
- date:             15 dec. 2022
- last modified:    16 dec. 2022
- brief:            wave system
"""

def wave_coenergy():
    """
    A structure-preserving discretization of the wave equation with boundary control
    
    Formulation co-energy, Grad-Grad, output feedback law at the boundary, damping on a subdomain
    """
    
    from scrimp import set_verbose_gf
    set_verbose_gf(0)
    
    from scrimp.dpHs import dpHs
    
    # Init the distributed port-Hamiltonian system
    wave = dpHs('real')
    
    # Set the domain (using the built-in geometry `Concentric`)
    # Omega_in = 1, Omega_out = 2, Interface = 10, Gamma = 20
    wave.set_domain('Concentric', {'R': 1., 'r': 0.6, 'h': 0.15})
    
    ## Define the variables and their discretizations
    
    # Add a state
    wave.add_state('q', 'Stress', 'vector-field')
    # Add its co-state
    wave.add_costate('q', 'Stress', 'q', substituted=True)
    # Add a Finite Element Method to the `port`
    wave.add_FEM('q', 1, FEM='DG')
    # Add a (possibly space-varying) parameter to the `port`
    wave.add_parameter('Tinv', 'Young\'s modulus inverse', 'tensor-field', '[[5+x,x*y],[x*y,2+y]]', 'q')
    
    # Add a state
    wave.add_state('p', 'Velocity', 'scalar-field')
    # Add its co-state
    wave.add_costate('p', 'Velocity', 'p', substituted=True)
    # Add a Finite Element Method to the `port`
    wave.add_FEM('p', 2, FEM='CG')
    # Add a (possibly space-varying) parameter to the `port`
    wave.add_parameter('rho', 'Mass density', 'scalar-field', '3-x', 'p')
    
    # Add a resistive `port`
    wave.add_port('Damping', 'e_r', 'e_r', 'scalar-field', substituted=True, region=1)
    # Add a FEM on it
    wave.add_FEM('Damping', 1, FEM='DG')
    # Attach a damping parameter to it
    wave.add_parameter('nu', 'Viscosity', 'scalar-field', '10*(0.36-(x*x+y*y))', 'Damping')
    
    # Add a control `port` on the bottom part of the boundary (Neumann, thus position='effort' - default)
    wave.add_control_port('Boundary control', 'U', 'Normal force', 'Y', 'Velocity trace', 'scalar-field', region=20)
    wave.add_FEM('Boundary control', 1, FEM='DG')
    
    # Set the Hamiltonian (can be done later, even after solve)
    wave.set_Hamiltonian_term('Potential energy', '0.5*q.Tinv.q', [1,2])
    wave.set_Hamiltonian_term('Kinetic energy', '0.5*p*p*rho', [1,2])
    wave.set_Hamiltonian_name('Mechanical energy')
    
    ## Define the Dirac structure via getfem `brick` = non-zero block matrix
    
    # Add the mass matrices from the left-hand side: the `flow` part of the Dirac structure
    wave.add_brick('M_q', 'q.Tinv.Test_q', [1,2], dt=True, position='flow')
    wave.add_brick('M_p', 'p*rho*Test_p', [1,2], dt=True, position='flow')
    wave.add_brick('M_r', 'e_r/nu*Test_e_r', [1], position='flow')
    wave.add_brick('M_Y', 'Y*Test_Y', [20], position='flow')
    
    # Add the matrices from the right-hand side: the `effort` part of the Dirac structure
    wave.add_brick('D', 'Grad(p).Test_q', [1,2], position='effort')
    
    wave.add_brick('-D^T', '-q.Grad(Test_p)', [1,2], position='effort')
    wave.add_brick('I_r', 'e_r*Test_p', [1], position='effort')
    wave.add_brick('B', 'U*Test_p', [20], position='effort')
    
    wave.add_brick('-I_r^T', '-p*Test_e_r', [1], position='effort')
    
    wave.add_brick('-B^T', '-p*Test_Y', [20], position='effort')
    
    ## Initialize the problem
    
    # Set the control functions (automatic construction of bricks such that -M_u u + f(t) = 0)
    # Output feedback: be aware that Y is the -y "written on paper", care must be taken for the sign!
    wave.set_control('Boundary control', '0.005*Y')
    
    # Set the initial data
    wave.set_initial_value('q', '[0., 0.]')
    wave.set_initial_value('p', '2.72**(-20*((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5)))')
    
    ## Solve in time
    
    # Define the time scheme (default ts_type='cn', t_f=1, dt=0.01, etc.)
    wave.set_time_scheme(t_f=1., dt=0.01, dt_save=0.01)
    
    # Solve
    wave.solve()
    
    ## Post-processing
    
    # Plot the Hamiltonian with the power supplied at the boundary
    wave.plot_Hamiltonian()
    
    # Export solutions for ParaView
    wave.export_to_pv('q')
    wave.export_to_pv('p')
    
    # Plot the matrices representing the Dirac structure
    # wave.spy_Dirac()
    
    return wave # For consol use

if __name__ == "__main__": 
    wave = wave_coenergy()
