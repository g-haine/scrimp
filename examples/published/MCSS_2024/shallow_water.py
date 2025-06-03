# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2025 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""
- file:             examples/published/MCSS_2023/shallow_water_grad.py
- authors:          Ghislain Haine, Flavio Luiz Cardoso-Ribeiro
- date:             20 sep. 2023
- brief:            2D shallow water equation as port-Hamiltonian system
"""

import scrimp as S
from itertools import zip_longest

def shallow_water(experiment=0, formulation="grad"):
    """
    A structure-preserving discretization of the rotational and dissipative shallow 
    water equation with boundary controls as port-Hamiltonian system.
    
    Discretisation use the grad-grad or the div-div formulation.
    
    This script can be used to reproduce the examples shown in the paper:
    @article{Cardoso_Ribeiro_2024,
        title={{Rotational shallow water equations with viscous damping and boundary control: structure-preserving spatial discretization}},
        ISSN={1435-568X},
        DOI={10.1007/s00498-024-00404-6},
        journal={Mathematics of Control, Signals, and Systems},
        publisher={Springer Science and Business Media LLC},
        author={Cardoso-Ribeiro, Flávio Luiz and Haine, Ghislain and Lefèvre, Laurent and Matignon, Denis},
        year={2024}
    }
    
    Example of use: to run experiment 1 with `grad-grad` formulation, write:
        `python shallow_water.py 1 grad`
    
    Outputs may be traced in ParaView with the `PV_trace_experiment.py` file in this folder.
    
    See the file `scrimp.yml` in the folder of this script to know the exact version of each library used for the results obtained in the paper.
    
    Args:
        - experiment (int): the experiment to reproduce:
            * 0: inviscid case with homogeneous boundary control
            * 1: viscous case with homogeneous boundary control
            * 2: emptying tank by the right end, i.e. normal control prescribed
            * 3: rotating tank, i.e. tangent control prescribed
        - formulation (string): the formulation, i.e., where the first integration by parts is done
            * grad: the first integration by parts is done on the mass preservation equation
            * div: the first integration by parts is done on the linear momentum equation
    Returns:
        the DPHS object
    """

    if experiment not in [0,1,2,3]:
        raise ValueError(f'Unknown experiment: {experiment}')
    if formulation not in ["grad","div"]:
        raise ValueError(f'Unknown formulation: {formulation}')
    
    # Experiment parameters
    dx = 0.05               # Mesh size parameter
    rho = 1.0               # Mass density
    g = 0.01                # Gravity constant
    mu = 0.001              # Viscosity
    order = 2
    FEM_h = ["CG", order]       # FEM for the h-type variables
    FEM_p = ["CG", order]       # FEM for the p-type variables
    FEM_b = ["CG", order]       # FEM for the boundary velocity controls and its colocated observations (both tangent and normal components)
    # PETSc time-stepper options
    ts_type = "bdf"
    ts_bdf_order = 2 # Not used if ts_type != `bdf`
    ts_arkimex_type = "a2" # Not used if ts_type != `arkimex`
    ksp_type = "preonly"
    pc_type = "lu"
    pc_factor_mat_solver_type = "superlu"
    t_0 = 0.
    t_f = 5.
    dt = 0.0001
    dt_save = 0.01
    ts_adapt_dt_min = 0.001*dt
    ts_adapt_dt_max = 50.*dt
    init_step = True
    init_step_nb_iter = 2 # Not used if init_step = False
    init_step_ts_type = "pseudo" # Not used if init_step = False
    init_step_dt = 0.1*ts_adapt_dt_min # Not used if init_step = False
    int_scheme = "Explicit" # `Explicit`, `Implicit`, or `CN` for the time integration of powers
    CN = True # Power computation using Crank-Nicolson
    if experiment==0:
        geometry = "Rectangle"  # Geometry of the domain
        L = 2.                  # Length of the rectangular tank
        # FEM of the h-type variables increased without dissipation and grad
        # Otherwise p-type increased without dissipation and div
        if formulation=="grad":
            FEM_h = ["CG", order+1]
        elif formulation=="div":
            FEM_p = ["CG", order+1]
        # Initial data
        k = 0.25
        c = 0.5
        s = 5.
        h_init = f"50. - 0.5*{s} * (np.sign(x-{k}*{L} + {c}*{dx}) - 1) \
                   - 0.25*(np.sign(x-{k}*{L} + {c}*{dx}) + 1) \
                   * (np.sign(x-{k}*{L} - {c}*{dx}) - 1) \
                   * ((-{s}/(2*{c}*{dx}))*x + {s}/2 + {s}*{k}*{L}/(2*{c}*{dx}))"
        p_init = "[ 0., 0. ]"
        # Controls
        U_n = "0."              # Normal value
    elif experiment==1:
        geometry = "Rectangle"  # Geometry of the domain
        L = 2.                  # Length of the rectangular tank
        # Initial data
        k = 0.25
        c = 0.5
        s = 5.
        h_init = f"50. - 0.5*{s} * (np.sign(x-{k}*{L} + {c}*{dx}) - 1) \
                   - 0.25*(np.sign(x-{k}*{L} + {c}*{dx}) + 1) \
                   * (np.sign(x-{k}*{L} - {c}*{dx}) - 1) \
                   * ((-{s}/(2*{c}*{dx}))*x + {s}/2 + {s}*{k}*{L}/(2*{c}*{dx}))"
        p_init = "[ 0., 0. ]"
        # Controls
        U_n = "0."              # Normal value
        U_t = "0."              # Tangent value
    elif experiment==2:
        geometry = "Rectangle"  # Geometry of the domain
        L = 2.                  # Length of the rectangular tank
        # Initial data
        h_init = "50."
        p_init = "[ 0., 0. ]"
        # Controls
        U_n = f"0.1 * min(t,1) * y * ({L}/4-y) * (sign(x-{L})+1)"   # Normal value
        U_t = "0."                                                  # Tangent value
        # Long emptying, double final time
        t_f *= 2.
    elif experiment==3:
        geometry = "Disk"       # Geometry of the domain
        L = 1.                  # Radius of the circular tank
        dx *= 2.                # Need less discretisation
        # Initial data
        h_init = "50."
        c = 0.2
        p_init = f"[ -{c}*np.sin(0.5*np.pi*np.sqrt(x*x+y*y))*np.sin(np.arctan2(y,x)), {c}*np.sin(0.5*np.pi*np.sqrt(x*x+y*y))*np.cos(np.arctan2(y,x)) ]"
        # Controls
        U_n = "0."      # Normal value
        U_t = f"{c} * min( 1, exp(- pow((max(90,t)-100) / (500-max(90,t)),2) ))"    # Tangent value
        # Long time behavior, 100 * final time
        t_f *= 100.
    
    # Init the distributed port-Hamiltonian system
    swe = S.DPHS("real")
    
    # Useful macros for the weak forms
    swe.gf_model.add_macro('div(v)', 'Trace(Grad(v))')
    swe.gf_model.add_macro('D(v)', 'Sym(Grad(v))')
    swe.gf_model.add_macro('Rot', '[[0,1],[-1,0]]')
    swe.gf_model.add_macro('Tangent', '(Rot*Normal)')
    swe.gf_model.add_macro('Curl2D(v)', 'div(Rot*v)')
    swe.gf_model.add_macro('Gyro(v)', 'Curl2D(v)*Rot')

    # Set the domain
    if experiment in [0,1,2]:
        domain = S.Domain(geometry, {"L": L, "l": L/4., "h": dx})
        swe.set_domain(domain)
        # Merging boundaries to avoid multiples boundary ports
        swe.domain.get_mesh()[0].region_merge(10,11)
        swe.domain.get_mesh()[0].region_merge(10,12)
        swe.domain.get_mesh()[0].region_merge(10,13)
        # Deleting merged boundaries
        swe.domain.get_mesh()[0].delete_region([11,12,13])
    else:
        domain = S.Domain(geometry, {"R": L, "h": dx})
        swe.set_domain(domain)
    
    # Define the variables
    states = [
        S.State("h", "Fluid height", "scalar-field"),
        S.State("p", "Linear momentum", "vector-field"),
    ]
    costates = [
        S.CoState("e_h", "Pressure", states[0]),
        S.CoState("e_p", "Velocity", states[1]),
    ]
    # no dissipative ports is to be found, the 'J-R' formulation is used
    ports = []
    # Parameters
    params = [
        S.Parameter("rho", "Mass density", "scalar-field", f"{rho}", "h"),
        S.Parameter("g", "Gravity", "scalar-field", f"{g}", "h"),
    ]
    if experiment in [1,2,3]:
        params.append(S.Parameter("mu", "Viscosity", "scalar-field", f"{mu}", "h"))

    # In experiment 0, the only control is the normal velocity,
    # in other experiment in `grad` formulation, it is imposed by the other control port (which is a vector-field port)
    # Normal control port does not exist in the viscous `div` formulation
    Control_normal = S.Control_Port(
                        r"Control $u^0$",
                        "U_n",
                        "Normal velocity",
                        "Y_n",
                        "Fluid height",
                        "scalar-field",
                        region=10,
                        position="effort",
                    )
    Control_vector = S.Control_Port(
                            r"Control $\mathbf{u}$",
                            "U",
                            "Velocity",
                            "Y",
                            "~Normal derivative",
                            "vector-field",
                            region=10,
                            position="effort",
                        )

    control_ports = []
    if formulation=="grad":
        control_ports.append(Control_normal)
        if experiment in [1,2,3]:
            control_ports.append(Control_vector)
    if formulation=="div":
        if experiment==0:
            control_ports.append(Control_normal)
        if experiment in [1,2,3]:
            control_ports.append(Control_vector)

    # Discretization of each port
    FEMs = [
        S.FEM(states[0].get_name(), FEM_h[1], FEM=FEM_h[0]),
        S.FEM(states[1].get_name(), FEM_p[1], FEM=FEM_p[0]),
    ]
    for control_port in control_ports:
        FEMs.append(S.FEM(control_port.get_name(), FEM_b[1], FEM=FEM_b[0]))

    # Add all objects inside the DPHS
    for state, costate, param, fem, port, control_port in zip_longest(
        states, costates, params, FEMs, ports, control_ports
    ):
        if state is not None:
            # Add a state
            swe.add_state(state)
        if costate is not None:
            # Add its co-state
            swe.add_costate(costate)
        if fem is not None:
            # Add a Finite Element Method to the `port`
            swe.add_FEM(fem)
        if param is not None:
            # Add a (possibly space-varying) parameter to the `port`
            swe.add_parameter(param)
        if port is not None:
            # Add a resistive `port`
            swe.add_port(port)
        if control_port is not None:
            # Add a control `port`
            swe.add_control_port(control_port)

    # Define the Dirac structure via getfem `brick` = non-zero block matrix
    bricks = [
        # Define the mass matrices of the left-hand side of the "Dirac structure"
        S.Brick("M_h", "h * Test_h", [1], dt=True, position="flow"),
        S.Brick("M_p[h]", "h * p . Test_p", [1], dt=True, linear=False, position="flow"),
        # with the gyroscopic term
        S.Brick("G[h,p]", "h * (Gyro(p) * e_p) . Test_p", [1], linear=False, explicit=False, position="effort"),
        
        # Define the constitutive relations
        # For e_h: first the mass matrix WITH A MINUS because we want an implicit formulation: 0 = - M e_h + F(h)
        S.Brick("-M_h", "- e_h * Test_e_h", [1]),
        # second the linear part as a linear brick
        S.Brick("P_h", "rho * g * h * Test_e_h", [1]),
        # third the non-linear part as a non-linear brick
        S.Brick("K_h[p]", "0.5 * (p . p) / rho * Test_e_h", [1], linear=False, explicit=False),
        
        # For e_p: first the mass matrix WITH A MINUS because we want an implicit formulation: 0 = - M e_p + F(p)
        S.Brick("-M_p[h]", "- h * e_p . Test_e_p", [1], linear=False, explicit=False),
        # second the non-linear brick
        S.Brick("K_p[h]", "h * p / rho . Test_e_p", [1], linear=False, explicit=False),
    ]
    
    Brick_Du = S.Brick("-2 mu R_Grad[h]", "- 2 * mu * h * D(e_p) : D(Test_p)", [1], linear=False, explicit=False, position="effort")
    Brick_divu = S.Brick("-2 mu R_div[h]", "- 2 * mu * h * div(e_p) * div(Test_p)", [1], linear=False, explicit=False, position="effort")
    Brick_Lagrange_multiplier = S.Brick("B_b[h]", "h * Y . Test_p", [10], linear=False, explicit=False, position="effort")
    Brick_constraint_mass = S.Brick("M_b[h]", "h * U . Test_Y", [10], linear=False, explicit=False, position="flow")
    Brick_constraint = S.Brick("B_b[h]^T", "h * e_p . Test_Y", [10], linear=False, explicit=False, position="effort")
    
    # The `grad` formulation always needs the following bricks
    if formulation=="grad":
        bricks.append(S.Brick("M_n[h]", "h * Y_n * Test_Y_n", [10], linear=False, explicit=False, position="flow"))
        bricks.append(S.Brick("D[h]", "h * e_p . Grad(Test_h)", [1], linear=False, explicit=False, position="effort"))
        bricks.append(S.Brick("-B_n[h]", "- h * U_n * Test_h", [10], linear=False, explicit=False, position="effort"))
        bricks.append(S.Brick("-D[h]^T", "- h * Grad(e_h) . Test_p", [1], linear=False, explicit=False, position="effort"))
        bricks.append(S.Brick("-B_n[h]^T", "- h * e_h * Test_Y_n", [10], linear=False, explicit=False, position="effort"))
        if experiment in [1,2,3]:
            bricks.append(Brick_Du)
            bricks.append(Brick_divu)
            bricks.append(Brick_Lagrange_multiplier)
            bricks.append(Brick_constraint_mass)
            bricks.append(Brick_constraint)
    # The `div` formulation
    if formulation=="div":
        bricks.append(S.Brick("D[h]", "- div(h * e_p) * Test_h", [1], linear=False, explicit=False, position="effort"))
        bricks.append(S.Brick("-D[h]^T", "e_h * div(h * Test_p)", [1], linear=False, explicit=False, position="effort"))
        if experiment==0:
            bricks.append(S.Brick("B_n[h]", "h * Y_n * Test_p . Normal", [10], linear=False, explicit=False, position="effort"))
            bricks.append(S.Brick("M_n[h]", "h * U_n * Test_Y_n", [10], linear=False, explicit=False, position="flow"))
            bricks.append(S.Brick("B_n[h]^T", "h * (e_p . Normal) * Test_Y_n", [10], linear=False, explicit=False, position="effort"))
        if experiment in [1,2,3]:
            bricks.append(Brick_Du)
            bricks.append(Brick_divu)
            bricks.append(Brick_Lagrange_multiplier)
            bricks.append(Brick_constraint_mass)
            bricks.append(Brick_constraint)
            
    # Add all bricks inside the DPHS
    for brick in bricks:
        swe.add_brick(brick)

    # Set the control
    if formulation=="grad":
        swe.set_control(control_ports[0].get_name(), "("+U_n+")")
        if experiment in [1,2,3]:
            swe.set_control(control_ports[1].get_name(), "(("+U_n+")*Normal + ("+U_t+")*Tangent)")
    if formulation=="div":
        if experiment==0:
            swe.set_control(control_ports[0].get_name(), "("+U_n+")")
        if experiment in [1,2,3]:
            swe.set_control(control_ports[0].get_name(), "(("+U_n+")*Normal + ("+U_t+")*Tangent)")

    # Set the initial data
    swe.set_initial_value("h", h_init)
    swe.set_initial_value("p", p_init)

    # Define the time scheme
    swe.set_time_scheme(
        ts_type=ts_type,
        ts_bdf_order=ts_bdf_order,
        ts_arkimex_type=ts_arkimex_type,
        ksp_type=ksp_type,
        pc_type=pc_type,
        pc_factor_mat_solver_type=pc_factor_mat_solver_type,
        t_0=t_0,
        t_f=t_f,
        dt=dt,
        dt_save=dt_save,
        ts_adapt_dt_min=ts_adapt_dt_min,
        ts_adapt_dt_max=ts_adapt_dt_max,
        init_step=init_step,
        init_step_nb_iter=init_step_nb_iter,
        init_step_ts_type=init_step_ts_type,
        init_step_dt=init_step_dt,
    )
    
    # Solve in time
    swe.solve()

    # Saving solutions for ParaView post-processing
    from petsc4py import PETSc
    comm = PETSc.COMM_WORLD
    rank = comm.getRank()
    import os
    path = os.path.join(os.path.dirname(os.path.realpath(__file__)), os.path.join("outputs", formulation))
    if rank==0:
        if not os.path.isdir(path):
            os.makedirs(path)
    comm.barrier()
    path = os.path.join(path, f"experiment_{experiment}")
    if rank==0:
        if not os.path.isdir(path):
            os.makedirs(path)
    comm.barrier()
    
    swe.export_to_pv("h", path=path)
    swe.export_to_pv("e_p", path=path)

    # Plots
    import matplotlib.pyplot as plt
    plt.ioff()
    import numpy as np
    
    # Set Hamiltonian
    swe.hamiltonian.set_name("Mechanical energy")
    terms = [
        S.Term("Kinetic energy", "0.5 * h * p.p / rho", [1]),
        S.Term("Potential energy", "0.5 * rho * g * h*h", [1]),
    ]
    for term in terms:
        swe.hamiltonian.add_term(term)
    swe.compute_Hamiltonian() # Will compute each term

    # Plot the Hamiltonian using the built-in function
    #swe.plot_Hamiltonian(with_powers=False)
    
    # The computation of powers needs to be consistent with the time scheme.
    # Hence, f.g will be computed as follows: (f_{n+1}+f_n)^T/2 M (g_{n+1}+g_n)/2 using the `CN` keyword
    t = np.array(swe.solution["t"])
    HamTot = np.zeros(t.size)
    TotalEnergy = np.zeros(t.size) # Will contain every energy parts (including integration over time of dissipated and boundary powers)
    Terms = swe.hamiltonian.get_terms()
    for term in Terms:
        HamTot += np.array(term.get_values())
        TotalEnergy += np.array(term.get_values())
    Powers = []
    if experiment in [1,2,3]:
        # The dissipation induced by D(e_p) in the model
        D_diss = np.array(swe.get_quantity("2 * mu * h * D(e_p) : D(e_p)", CN=CN))
        Powers.append(D_diss)
        if int_scheme=="CN":
            int_dt_D_diss = 0.5*(D_diss[0:-1] + D_diss[1:])*(t[1:]-t[0:-1])
        elif int_scheme=="Explicit":
            int_dt_D_diss = D_diss[1:]*(t[1:]-t[0:-1])
        elif int_scheme=="Implicit":
            int_dt_D_diss = D_diss[0:-1]*(t[1:]-t[0:-1])
        else:
            raise ValueError(f'Unknown int_scheme: {int_scheme}')
        int_D_diss = np.array([int_dt_D_diss[0:k].sum() for k in range(len(t))]) # integration over time
        
        # The dissipation induced by div(e_p) in the model
        div_diss = np.array(swe.get_quantity("2 * mu * h * div(e_p) * div(e_p)", CN=CN))
        Powers.append(div_diss)
        if int_scheme=="CN":
            int_dt_div_diss = 0.5*(div_diss[0:-1] + div_diss[1:])*(t[1:]-t[0:-1])
        elif int_scheme=="Explicit":
            int_dt_div_diss = div_diss[1:]*(t[1:]-t[0:-1])
        elif int_scheme=="Implicit":
            int_dt_div_diss = div_diss[0:-1]*(t[1:]-t[0:-1])
        int_div_diss = np.array([int_dt_div_diss[0:k].sum() for k in range(len(t))]) # integration over time
        
        TotalEnergy += int_D_diss + int_div_diss
    
    # The power flowing through the control ports
    Energies = []
    if formulation=="grad":
        power = np.array(swe.get_quantity("-h * U_n * Y_n", CN=CN, region=10))
        Powers.append(power)
        if int_scheme=="CN":
            energy_dt = 0.5*(power[0:-1] + power[1:])*(t[1:]-t[0:-1])
        elif int_scheme=="Explicit":
            energy_dt = power[1:]*(t[1:]-t[0:-1])
        elif int_scheme=="Implicit":
            energy_dt = power[0:-1]*(t[1:]-t[0:-1])
        Energies.append(np.array([energy_dt[0:k].sum() for k in range(len(t))])) # integration over time
        if experiment in [1,2,3]:
            power = np.array(swe.get_quantity("-h * U . Y", CN=CN, region=10))
            Powers.append(power)
            if int_scheme=="CN":
                energy_dt = 0.5*(power[0:-1] + power[1:])*(t[1:]-t[0:-1])
            elif int_scheme=="Explicit":
                energy_dt = power[1:]*(t[1:]-t[0:-1])
            elif int_scheme=="Implicit":
                energy_dt = power[0:-1]*(t[1:]-t[0:-1])
            Energies.append(np.array([energy_dt[0:k].sum() for k in range(len(t))])) # integration over time
    if formulation=="div":
        if experiment==0:
            power = np.array(swe.get_quantity("-h * U_n * Y_n", CN=CN, region=10))
            Powers.append(power)
            if int_scheme=="CN":
                energy_dt = 0.5*(power[0:-1] + power[1:])*(t[1:]-t[0:-1])
            elif int_scheme=="Explicit":
                energy_dt = power[1:]*(t[1:]-t[0:-1])
            elif int_scheme=="Implicit":
                energy_dt = power[0:-1]*(t[1:]-t[0:-1])
            Energies.append(np.array([energy_dt[0:k].sum() for k in range(len(t))])) # integration over time
        if experiment in [1,2,3]:
            power = np.array(swe.get_quantity("-h * U . Y", CN=CN, region=10))
            Powers.append(power)
            if int_scheme=="CN":
                energy_dt = 0.5*(power[0:-1] + power[1:])*(t[1:]-t[0:-1])
            elif int_scheme=="Explicit":
                energy_dt = power[1:]*(t[1:]-t[0:-1])
            elif int_scheme=="Implicit":
                energy_dt = power[0:-1]*(t[1:]-t[0:-1])
            Energies.append(np.array([energy_dt[0:k].sum() for k in range(len(t))])) # integration over time
    
    for energy in Energies:
        TotalEnergy += energy
    
    plt.rcParams['text.usetex'] = True
    plt.rcParams['font.size'] = 32
    figsize = [16,10]
    if rank==0:        
        fig = plt.figure(figsize=figsize, layout="constrained")
        ax = fig.add_subplot(111)
        ax.plot(t, HamTot, "r-")
        fig.legend("Hamiltonian", loc='outside center right')
        plt.xlabel("Time (s)")
        plt.ylabel("Energy (J)")
        fig.suptitle("Evolution of the Hamiltonian")
        ax.grid(axis="both")
        fig.savefig(os.path.join(path,"Hamiltonian.svg"))
        plt.close(fig)
        
        fig = plt.figure(figsize=figsize, layout="constrained")
        ax = fig.add_subplot(111)
        ax.plot(t, (TotalEnergy-TotalEnergy[0])/np.max(TotalEnergy), "r--")
        fig.legend(["Total energy"], loc='outside center right')
        plt.xlabel("Time (s)")
        plt.ylabel("Variation (relative)")
        fig.suptitle("Relative variation of the total energy")
        ax.grid(axis="both")
        fig.savefig(os.path.join(path,"Total_energy_variation.svg"))
        plt.close(fig)
        
        fig = plt.figure(figsize=figsize, layout="constrained")
        fig.suptitle("Evolution of the energies")
        if experiment in [1,2,3]:
            plots = [511,512,513,514,515]
        else:
            plots = [411,412,413,414]
        ax1 = fig.add_subplot(plots[0])
        ax1.plot(t, HamTot, "r-", t, TotalEnergy, "r--")
        ax1.grid(axis="both")
        plt.tick_params('x', labelbottom=False)
        ax2 = fig.add_subplot(plots[1], sharex=ax1)
        ax2.plot(t, Terms[1].get_values(), "b-")
        ax2.grid(axis="both")
        plt.tick_params('x', labelbottom=False)
        ax3 = fig.add_subplot(plots[2], sharex=ax1)
        ax3.plot(t, Terms[0].get_values(), "b--")
        ax3.grid(axis="both")
        plt.tick_params('x', labelbottom=False)
        ax4 = fig.add_subplot(plots[3], sharex=ax1)
        ax4.plot(t, Energies[0], "m-")
        legend = ["Hamiltonian", "Total energy", Terms[1].get_description(), Terms[0].get_description(), r"$\mathcal{S} u^0$"]
        if experiment in [1,2,3] and formulation=="grad":
            ax4.plot(t, Energies[1], "m--")
            legend.append(r"$\mathcal{S} \mathbf{u}$")
        ax4.grid(axis="both")
        if experiment in [1,2,3]:
            plt.tick_params('x', labelbottom=False)
            ax5 = fig.add_subplot(plots[4], sharex=ax1)
            ax5.plot(t, int_D_diss, "g-", t, int_div_diss, "g--")
            legend.append(r"$\mathcal{D}_{Grad}$")
            legend.append(r"$\mathcal{D}_{div}$")
            ax5.grid(axis="both")
        fig.legend(legend, loc='outside center right')
        plt.xlabel("Time (s)")
        fig.supylabel("Energies (J)")
        fig.savefig(os.path.join(path,"Balance_complete.svg"))
        plt.close(fig)
        
        fig = plt.figure(figsize=figsize, layout="constrained")
        fig.suptitle("Evolution of the energies")
        if experiment in [1,2,3]:
            plots = [311,312,313]
        else:
            plots = [211,212]
        ax1 = fig.add_subplot(plots[0])
        ax1.plot(t, HamTot, "r-", t, TotalEnergy, "r--")
        ax1.grid(axis="both")
        plt.tick_params('x', labelbottom=False)
        ax2 = fig.add_subplot(plots[1], sharex=ax1)
        ax2.plot(t, Energies[0], "m-")
        legend = ["Hamiltonian", "Total energy", r"$\mathcal{S} u^0$"]
        if experiment in [1,2,3] and formulation=="grad":
            ax2.plot(t, Energies[1], "m--")
            legend.append(r"$\mathcal{S} \mathbf{u}$")
        ax2.grid(axis="both")
        if experiment in [1,2,3]:
            plt.tick_params('x', labelbottom=False)
            ax3 = fig.add_subplot(plots[2], sharex=ax1)
            ax3.plot(t, int_D_diss, "g-", t, int_div_diss, "g--")
            legend.append(r"$\mathcal{D}_{Grad}$")
            legend.append(r"$\mathcal{D}_{div}$")
            ax3.grid(axis="both")
        fig.legend(legend, loc='outside center right')
        plt.xlabel("Time (s)")
        fig.supylabel("Energies (J)")
        fig.savefig(os.path.join(path,"Balance.svg"))
        plt.close(fig)
        
    # Test if outputs Y are as expected: Y = - e_h * Normal + 2 mu ( Grad(e_p) . Normal + div(e_p) * Normal )
    L2_Error_Y_0 = None
    L1_Error_Y_0 = None
    L2_Error_Y_n = None
    L2_Error_Y_t = None
    L1_Error_Y_n = None
    L1_Error_Y_t = None
    if formulation=="grad" or (formulation=="div" and experiment==0):
        L2_Error_Y_0 = np.array(swe.get_quantity("pow((Y_n + e_h), 2)", region=10))
        L2_Y_0 = np.array(swe.get_quantity("pow((e_h), 2)", region=10))
        L1_Error_Y_0 = np.array(swe.get_quantity("abs(Y_n + e_h)", region=10))
        L1_Y_0 = np.array(swe.get_quantity("abs(e_h)", region=10))
    if (formulation=="grad" or formulation=="div") and (experiment in [1,2,3]):
        L2_Error_Y_n = np.array(swe.get_quantity("pow((Y - 2 * mu * D(e_p) . Normal - 2 * mu * div(e_p) * Normal) . Normal, 2)", region=10))
        L2_Y_n = np.array(swe.get_quantity("pow((2 * mu * D(e_p) . Normal + 2 * mu * div(e_p) * Normal) . Normal, 2)", region=10))
        L2_Error_Y_t = np.array(swe.get_quantity("pow((Y - 2 * mu * D(e_p) . Normal) . Tangent, 2)", region=10))
        L2_Y_t = np.array(swe.get_quantity("pow((2 * mu * D(e_p) . Normal) . Tangent, 2)", region=10))
        L1_Error_Y_n = np.array(swe.get_quantity("abs((Y - 2 * mu * D(e_p) . Normal - 2 * mu * div(e_p) * Normal) . Normal)", region=10))
        L1_Y_n = np.array(swe.get_quantity("abs((2 * mu * D(e_p) . Normal + 2 * mu * div(e_p) * Normal) . Normal)", region=10))
        L1_Error_Y_t = np.array(swe.get_quantity("abs((Y - 2 * mu * D(e_p) . Normal) . Tangent)", region=10))
        L1_Y_t = np.array(swe.get_quantity("abs((2 * mu * D(e_p) . Normal) . Tangent)", region=10))
        
    if rank==0:        
        legend = []
        fig = plt.figure(figsize=figsize, layout="constrained")
        ax = fig.add_subplot(111)
        if L2_Error_Y_0 is not None:
            error_0 = np.sqrt(L2_Error_Y_0)/(1.+np.sqrt(L2_Y_0))
            ax.plot(t, error_0, "k-")
            legend.append(r"$L^2$-error $n$ ($y^0$)")
        if L2_Error_Y_n is not None:
            error_n = np.sqrt(L2_Error_Y_n)/(1.+np.sqrt(L2_Y_n))
            error_t = np.sqrt(L2_Error_Y_t)/(1.+np.sqrt(L2_Y_t))
            ax.plot(t, error_n, "k--", t, error_t, "k:")
            legend.append(r"$L^2$-error $n$ ($\mathbf{y}^\mu \cdot \mathbf{n}$)")
            legend.append(r"$L^2$-error $t$ ($\mathbf{y}^\mu \cdot \mathbf{t}$)")
        fig.legend(legend, loc='outside center right')
        plt.xlabel("Time (s)")
        plt.ylabel(r"$L^2$-error")
        fig.suptitle(r"Relative $L^2$-error: $\left\| Y - Y_{theoretic} \right\|_{L^2(\partial\Omega)} / (1+\left\| Y_{theoretic} \right\|_{L^2(\partial\Omega)})$")
        ax.grid(axis="both")
        fig.savefig(os.path.join(path,"L2_error_on_Y.svg"))
        plt.close(fig)
        legend = []
        
        fig = plt.figure(figsize=figsize, layout="constrained")
        ax = fig.add_subplot(111)
        if L1_Error_Y_0 is not None:
            error_0 = L1_Error_Y_0/(1.+L1_Y_0)
            ax.plot(t, error_0, "k-")
            legend.append(r"$L^1$-error $n$ ($y^0$)")
        if L1_Error_Y_n is not None:
            error_n = L1_Error_Y_n/(1.+L1_Y_n)
            error_t = L1_Error_Y_t/(1.+L1_Y_t)
            ax.plot(t, error_n, "k--", t, error_t, "k:")
            legend.append(r"$L^1$-error $n$ ($\mathbf{y}^\mu \cdot \mathbf{n}$)")
            legend.append(r"$L^1$-error $t$ ($\mathbf{y}^\mu \cdot \mathbf{t}$)")
        fig.legend(legend, loc='outside center right')
        plt.xlabel("Time (s)")
        plt.ylabel(r"$L^1$-error")
        fig.suptitle(r"Relative $L^1$-error: $\left\| Y - Y_{theoretic} \right\|_{L^1(\partial\Omega)} / (1+\left\| Y_{theoretic} \right\|_{L^1(\partial\Omega)})$")
        ax.grid(axis="both")
        fig.savefig(os.path.join(path,"L1_error_on_Y.svg"))
        plt.close(fig)
    
    return swe

if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        experiment = int(sys.argv[1])
        if len(sys.argv)>2:
            formulation = sys.argv[2]
    else:
        experiment = 0
        formulation = "grad"
    swe = shallow_water(experiment,formulation)

