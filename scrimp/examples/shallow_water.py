# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2023 Ghislain Haine
#
# See the LICENSE file in the root directory for license information.
#
# github: https://github.com/g-haine/scrimp

"""
- file:             examples/shallow_water.py
- author:           Ghislain Haine
- date:             07 dec. 2022
- last modified:    12 dec. 2022
- brief:            shallow water system
"""

def shallow_water():
    """
    A structure-preserving discretization of the shallow-water equation with boundary control
    
    !TO DO: add Navier-Stokes dissipation (and in particular vorticity)
    
    Formulation DAE, Grad-Grad, uniform boundary condition on the Disk
    """

    from scrimp import set_verbose_gf
    set_verbose_gf(0) # Select the verbosity level of getfem from 0 (quiet) to 3
    
    from scrimp import dpHs
    
    swe = dpHs('real')
    
    # Set the domain
    # Omega = 1, Gamma = 10
    swe.set_domain('Disk', {'R': 1, 'h': 0.15})
    
    # Add the first state and co-state and associate a FEM to the resulting port
    swe.add_state('h', 'Fluid height', 'scalar-field')
    swe.add_costate('e_h', 'Pressure', 'h')
    swe.add_FEM('h', 2, FEM='CG')

    # Add the second state and co-state and associate a FEM to the resulting port
    swe.add_state('p', 'Linear momentum', 'vector-field')
    swe.add_costate('e_p', 'Rate of flow', 'p')
    swe.add_FEM('p', 1, FEM='DG')
    
    # Define two (possibly space-varying) parameters
    swe.add_parameter('rho', 'Mass density', 'scalar-field', '1000.', 'h')
    swe.add_parameter('g', 'Gravity', 'scalar-field', '10.', 'h')
    
    # Define the two terms of the Hamiltonian and set its name
    swe.set_Hamiltonian_term('Kinetic energy', '0.5*h*p.p/rho', [1])
    swe.set_Hamiltonian_term('Potential energy', '0.5*rho*g*h*h', [1])
    swe.set_Hamiltonian_name('Mechanical energy')
    
    # Add a control port on region 10 of the domain (boundary)
    swe.add_control_port('Boundary control', 'U', 'Normal rate of flow', 'Y', 'Fluid height', 'scalar-field', region=10, position='effort')
    swe.add_FEM('Boundary control', 2, FEM='CG')
    
    # Define the mass matrices of the left-hand side of the 'Dirac structure' (position='flow')
    swe.add_brick('M_h', 'h*Test_h', [1], dt=True, position='flow')
    swe.add_brick('M_p', 'p.Test_p', [1], dt=True, position='flow')
    swe.add_brick('M_Y', 'Y*Test_Y', [10], position='flow')
    
    # Define the first line of the right-hand side of the 'Dirac structure' (position='effort')
    swe.add_brick('-D^T', 'e_p.Grad(Test_h)', [1], position='effort')
    # with the boundary control
    swe.add_brick('B', '-U*Test_h', [10], position='effort')
    
    # Define the second line of the right-hand side of the 'Dirac structure' (position='effort')
    swe.add_brick('D', '-Grad(e_h).Test_p', [1], position='effort')
    # with the gyroscopic term (beware that 'Curl' is not available in the GWFL of getfem)
    swe.add_brick('G', 'Trace([[0,1],[-1,0]]*Grad(p/rho))*rho*[[0, 1],[-1, 0]]*e_p/h.Test_p', 
                  [1], linear=False, position='effort')

    # Define the third line of the right-hand side of the 'Dirac structure' (position='effort')
    swe.add_brick('C', 'e_h*Test_Y', [10], position='effort')
    
    # Define the constitutive relations (position='constitutive', the default value)
    # For e_h: first the mass matrix WITH A MINUS because we want an implicit formulation 0 = - M e_h + F(h)
    swe.add_brick('-M_e_h', '-e_h*Test_e_h', [1])
    # second the linear part as a linear brick
    swe.add_brick('CR_h_lin', 'rho*g*h*Test_e_h', [1])
    # third the non-linear part as a non-linear brick (linear=False)
    swe.add_brick('CR_h_nl', '0.5*p.p*Test_e_h/rho', [1], linear=False)
    
    # For e_p: first the mass matrix WITH A MINUS because we want an implicit formulation 0 = - M e_p + F(p)
    swe.add_brick('-M_e_p', '-e_p.Test_e_p', [1])
    # second the non-linear brick (linear=False)
    swe.add_brick('CR_p', 'h*p.Test_e_p/rho', [1], linear=False)
    
    # Finally we set the control function (time-dependent)
    swe.set_control('Boundary control', '0.')
    # and the initial values (space-varying), we can also call numpy with 'np', and use macros for polar coordinates
    swe.set_initial_value('h', '10. + 0.1*np.exp(-50*((x-0.5)*(x-0.5) + (y-0.3)*(y-0.3)))')
    swe.set_initial_value('p', '[ - 10.*np.sin(np.pi*np.sqrt(x*x+y*y))*np.sin(np.arctan2(y,x)), 10.*np.sin(np.pi*np.sqrt(x*x+y*y))*np.cos(np.arctan2(y,x))]')
    
    # We define the time scheme to use (this step may be skipped, a default time-experience is then used)
    swe.set_time_scheme(ts_type='cn', ksp_type='preonly', 
                        pc_type='lu', #pc_factor_mat_solver_type='mumps',
                        t_0=0., t_f=0.5, dt=0.01, dt_save=0.01, init_step=True)
    
    # Solve the system in time
    swe.solve()
    
    # Plot the Hamiltonian (with_powers=True will add the time-integration of the product f*e for each algebraic port (f,e))
    swe.plot_Hamiltonian(with_powers=True)
    
    # Saving solutions for ParaView post-processing
    swe.export_to_pv('h')
    swe.export_to_pv('p')

if __name__ == "__main__": 
    shallow_water()
