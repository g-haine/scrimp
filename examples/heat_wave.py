# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2025 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""
- file:             sandbox/heat_hw.py
- authors:          Giuseppe Ferraro, Ghislain Haine
- date:             15 dec. 2022
- brief:            a 2D coupled heat-wave system
"""

# Import scrimp
import scrimp as S

def heat_wave_eq(heat_region=1, wave_region=2):
    """A structure-preserving discretization of a coupled heat-wave equation

    Co-energy formulations, heat: div-div, wave: grad-grad, gyrator interconnection
    On the `Concentric` built-in geometry: 1: internal disk, 2: exterior annulus
    
    Args:
        heat_region (int):  the label of the region where the heat equation lies
        wave_region (int):  the label of the region where the wave equation lies
    """
    
    # Init the distributed port-Hamiltonian system
    hw = S.DPHS("real")

    # Set the domain (using the built-in geometry `Concentric`)
    # Labels: Disk = 1, Annulus = 2, Interface = 10, Boundary = 20
    omega = S.Domain("Concentric", {"R": 1.0, "r": 0.6, "h": 0.05})
    
    # And add it to the dphs
    hw.set_domain(omega)
    
    # Define the states and costates, needs the heat and wave region's labels
    states = [
        S.State("T", "Temperature", "scalar-field", region=heat_region),
        S.State("p", "Velocity", "scalar-field", region=wave_region),
        S.State("q", "Stress", "vector-field", region=wave_region),
    ]
    # Use of the `substituted=True` keyword to get the co-energy formulation
    costates = [
        S.CoState("T", "Temperature", states[0], substituted=True),
        S.CoState("p", "Velocity", states[1], substituted=True),
        S.CoState("q", "Stress", states[2], substituted=True),
    ]
    
    # Add them to the dphs
    for state in states:
        hw.add_state(state)
    for costate in costates:
        hw.add_costate(costate)

    # Define the algebraic port
    ports = [
        S.Port("Heat flux", "e_Q", "e_Q", "vector-field", substituted=True, region=heat_region),
    ]
    
    # Add it to the dphs
    for port in ports:
        hw.add_port(port)

    # Define the control ports
    control_ports = [
        S.Control_Port(
            "Interface Heat", 
            "U_T", 
            "Heat flux", 
            "Y_T", 
            "Temperature",
            "scalar-field",
            region=10, 
            position="effort"
        ),
        S.Control_Port(
            "Interface Wave", 
            "U_w", 
            "Velocity", 
            "Y_w", 
            "Velocity", 
            "scalar-field", 
            region=10, 
            position="effort"
        ),
        # This port will be either for the wave or the heat equation
        # It corresponds to the exterior circle of radius R
        S.Control_Port( 
            "Boundary", 
            "U_bnd", 
            "0", 
            "Y_bnd", 
            ".", 
            "scalar-field", 
            region=20, 
            position="flow"
        ),
    ]

    # Add them to the dphs
    for ctrl_port in control_ports:
        hw.add_control_port(ctrl_port)

    # Define the Finite Elements Method of each port
    k = 1
    FEMs = [
        S.FEM("T", k, "DG"),
        S.FEM("Heat flux", k+1, "CG"),
        S.FEM("Interface Heat", k, "DG"),
        S.FEM("p", k+1, "CG"),
        S.FEM("q", k, "DG"),
        S.FEM("Interface Wave", k, "DG"),
        S.FEM("Boundary", k, "DG"),
    ]
    
    # Add them to the dphs
    for FEM in FEMs:
        hw.add_FEM(FEM)
    
    # Define the pHs via `Brick` == non-zero block matrices == variational terms
    # Since we use co-energy formulation, constitutive relations are already taken into
    # account in the mass matrices M_q and M_p
    bricks = [
        # === Heat: div-div
        S.Brick("M_T", "T*Test_T", [heat_region], dt=True, position="flow"),
        S.Brick("M_Q", "e_Q.Test_e_Q", [heat_region], position="flow"),
        S.Brick("M_Y_T", "Y_T*Test_Y_T", [10], position="flow"),
        
        S.Brick("D_T", "-Div(e_Q)*Test_T", [heat_region], position="effort"),
        S.Brick("D_T^T", "T*Div(Test_e_Q)", [heat_region], position="effort"),
        S.Brick("B_T", "U_T*Test_e_Q.Normal", [10], position="effort"),
        S.Brick("B_T^T", "e_Q.Normal*Test_Y_T", [10], position="effort"),

        # === Wave: grad-grad
        S.Brick("M_p", "p*Test_p", [wave_region], dt=True, position="flow"),
        S.Brick("M_q", "q.Test_q", [wave_region], dt=True, position="flow"),
        S.Brick("M_Y_w", "Y_w*Test_Y_w", [10], position="flow"),
        
        S.Brick("D_w", "-q.Grad(Test_p)", [wave_region], position="effort"),
        S.Brick("-D_w^T", "Grad(p).Test_q", [wave_region], position="effort"),
        S.Brick("B_w", "U_w*Test_p", [10], position="effort"),
        S.Brick("B_w^T", "p*Test_Y_w", [10], position="effort"),
    ]
    # === Boundary depends on where is the heat equation / wave equation
    if wave_region==1:
        bricks.append(S.Brick("M_Y_bnd", "Y_bnd*Test_Y_bnd", [20], position="flow"))
        bricks.append(S.Brick("B_bnd", "U_bnd*Test_e_Q.Normal", [20], position="effort"))
        bricks.append(S.Brick("B_bnd^T", "e_Q.Normal*Test_Y_bnd", [20], position="effort"))
    else:
        bricks.append(S.Brick("M_Y_bnd", "U_bnd*Test_Y_bnd", [20], position="flow"))
        bricks.append(S.Brick("B_bnd", "Y_bnd*Test_p", [20], position="effort"))
        bricks.append(S.Brick("B_bnd^T", "p*Test_Y_bnd", [20], position="effort"))
    for brick in bricks:
        hw.add_brick(brick)
    
    # Set the controls
    # === Gyrator interconnection
    hw.set_control("Interface Heat", "Y_w") 
    # CAREFUL: the numerical normal is the same for both sub-domains! Hence the minus sign. 
    hw.set_control("Interface Wave", "-Y_T") 
    # === Dirichlet boundary condition
    hw.set_control("Boundary", "0.")

    # Set the initial data
    hw.set_initial_value("T", "5.*np.exp(-25*((x-0.6)*(x-0.6)+y*y))")
    hw.set_initial_value("p", "5.*np.exp(-25*((x-0.6)*(x-0.6)+y*y))")
    hw.set_initial_value("q", "[0.,0.]")

    ## Solve in time
    # Define the time scheme ("bdf" is backward differentiation formula)
    hw.set_time_scheme(ts_type="bdf",
                       ts_bdf_order=2,
                       t_f=2.,
                       dt=0.001,
                       dt_save=0.05,
                       ksp_type="preonly",
                       pc_type="lu",
                       pc_factor_mat_solver_type="superlu",
                       )

    # Solve
    hw.solve()

    ## Post-processing
    ## Set Hamiltonian's name
    hw.hamiltonian.set_name("Hamiltonian")
    # Define each Hamiltonian Term
    terms = [
        S.Term("Lyapunov heat", "0.5*T*T", [heat_region]),
        S.Term("Kinetic energy", "0.5*p*p", [wave_region]),
        S.Term("Potential energy", "0.5*q.q", [wave_region]),
    ]
    # Add them to the Hamiltonian
    for term in terms:
        hw.hamiltonian.add_term(term)

    # Plot the Hamiltonian and save the output
    hw.plot_Hamiltonian(save_figure=True, filename="Hamiltonian_Heat"+str(heat_region)+"_Wave"+str(wave_region)+"_2D.png")
    
    # Plot the Hamiltonian in log-log scale
    t = hw.solution["t"]
    Hamiltonian = hw.get_Hamiltonian()
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=[8, 5])
    ax = fig.add_subplot(111)
    ax.loglog(t, Hamiltonian)
    ax.grid(axis="both")
    ax.set_xlabel("time t")
    ax.set_ylabel("Hamiltonian")
    ax.set_title("Evolution of the Hamiltonian (log-log)")
    plt.show()
    
    hw.export_to_pv("p")
    hw.export_to_pv("T")
    
    return hw

if __name__ == "__main__":
    
    import sys
    if len(sys.argv)==3:
        heat_region=int(sys.argv[1])
        wave_region=int(sys.argv[2])
    else:
        heat_region=1
        wave_region=2
    hw = heat_wave_eq(heat_region=heat_region, wave_region=wave_region)
    
