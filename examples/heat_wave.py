# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
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
import numpy as np

from scrimp.utils.fe_initialization import enforce_zero_dirichlet_initial_values

def heat_wave_eq(
    heat_region=1, wave_region=2, *, factor_solver=None,
    final_time=2.0, time_step=0.001, save_step=0.05,
    initial_wave_boundary="clamp",
):
    """A physically oriented, energy-dissipative coupled heat-wave example.

    Reference: Haine, Matignon & Monteghetti, Mathematical Reports (2022),
    eqs. (2.1)-(2.2), (3.15), (4.3)-(4.5). This example is NOT the paper's
    numerical experiment: it keeps the original SCRIMP heat div-div and
    wave grad-grad weak forms, choosing compatible mixed FEM spaces.

    The original built-in Concentric geometry has disk=1, annulus=2,
    interface=10, outer boundary=20. CN and the 0-input external
    Dirichlet condition are used to expose the heat dissipation and
    closed interconnection in an interpretable power balance.

    Interface normals always mean each subdomain's outward normal.
    The conservative two-sided conforming interface compiler is
    enabled by default in the v1.3.0 release-preparation source.
    Physical labels, FEM declarations and weak forms are unchanged.
    Ambiguous or unsupported interfaces fail explicitly; use
    SCRIMP_AUTO_INTERFACE=0 only to recover legacy manual assembly
    if your model was already oriented by hand.

    Args:
        heat_region: domain label of the heat equation (1 or 2).
        wave_region: domain label of the wave equation (the other one).
        factor_solver: optional PETSc LU solver (default: SCRIMP MUMPS).
        final_time: physical end time.
        time_step: fixed Crank–Nicolson timestep.
        save_step: accepted trajectory sampling interval.
        initial_wave_boundary: compatible wave velocity initialization
            when the exterior boundary belongs to the wave region.
            "clamp" (default) sets exactly the exterior CG1 degrees
            of freedom to zero, leaving interior L2-projected
            coefficients unchanged. "constrained_l2" computes the
            closest zero-trace field in the physical GetFEM L2 norm
            (MPI1 only). "unconstrained" reproduces the historical,
            generally DAE-incompatible Gaussian for diagnostics only.
            Neither valid method relaxes algebraic DAE constraints.
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

    if {heat_region, wave_region} != {1, 2}:
        raise ValueError("Concentric requires heat/wave regions 1 and 2")

    # The outer boundary is either wave velocity (Dirichlet) or heat
    # temperature (Dirichlet). Its port causality depends on the
    # adjacent subsystem, NOT on the user-selected LU factorization.
    outer_is_wave = (wave_region == 2)
    # On the interface, Normal is the outward normal of the bulk
    # subdomain owning the tested variable. The automatic conforming
    # interface compiler selects that side, independently of whether
    # heat occupies the disk or the annulus. The weak transmission law
    # below therefore needs NO correction based on region numbering.

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
        # Homogeneous Dirichlet: wave velocity=0 (wave annulus)
        # or temperature=0 (heat annulus). Normal traction/flux is
        # the corresponding Lagrange-multiplier output.
        S.Control_Port(
            "Boundary",
            "U_bnd",
            "Velocity" if outer_is_wave else "Temperature",
            "Y_bnd",
            "Normal force" if outer_is_wave else "Normal heat flux",
            "scalar-field",
            region=20,
            position="flow" if outer_is_wave else "effort",
        ),
    ]

    # Add them to the dphs
    for ctrl_port in control_ports:
        hw.add_control_port(ctrl_port)

    # Paper (4.3)-(4.5): DG_k heat scalar, RT_k heat flux,
    # CG_(k+1) interface traces. Our existing wave grad-grad PFEM
    # requires a continuous velocity p: CG_(k+1) rather than the
    # article's div-div DG_k. Wave stress q is now H(div)-conforming
    # RT_k, giving a genuine normal trace.
    # k=0 is the order used in the reference numerical experiments.
    k = 0
    FEMs = [
        S.FEM("T", k, "DG"),
        S.FEM("Heat flux", k, "RT"),
        S.FEM("Interface Heat", k + 1, "CG"),
        S.FEM("p", k + 1, "CG"),
        S.FEM("q", k, "RT"),
        S.FEM("Interface Wave", k + 1, "CG"),
        S.FEM("Boundary", k + 1, "CG"),
    ]
    
    # Add them to the dphs
    for FEM in FEMs:
        hw.add_FEM(FEM)
        # Earlier GetFEM versions could silently fall back to CG
        # on an unavailable intrinsic RT element. Never claim an
        # RT-based model if that fallback occurred.
        if FEM.get_name() in ("Heat flux", "q") and FEM.get_type() != "RT":
            raise RuntimeError(
                "RT finite elements are required for the heat-wave example"
            )

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
        # With e_Q=-grad(T) and the outward HEAT normal, T=U_T
        # enters with a minus in the mixed boundary weak form,
        # for either assignment of heat_region.
        S.Brick(
            "B_T", "-U_T*Test_e_Q.Normal",
            [10], position="effort",
        ),
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
    # Outer homogeneous Dirichlet, consistently with (2.1):
    # - heat exterior: prescribe temperature, observe normal flux;
    # - wave exterior: prescribe velocity, normal traction is the
    #   Lagrange multiplier, following examples/wave.py.
    if not outer_is_wave:
        bricks.append(S.Brick(
            "M_Y_bnd", "Y_bnd*Test_Y_bnd", [20], position="flow"
        ))
        bricks.append(S.Brick(
            "B_bnd", "-U_bnd*Test_e_Q.Normal", [20], position="effort"
        ))
        bricks.append(S.Brick(
            "B_bnd^T", "e_Q.Normal*Test_Y_bnd", [20], position="effort"
        ))
    else:
        bricks.append(S.Brick(
            "M_Y_bnd", "U_bnd*Test_Y_bnd", [20], position="flow"
        ))
        bricks.append(S.Brick(
            "B_bnd", "Y_bnd*Test_p", [20], position="effort"
        ))
        bricks.append(S.Brick(
            "B_bnd^T", "-p*Test_Y_bnd", [20], position="effort"
        ))
    for brick in bricks:
        hw.add_brick(brick)
    
    # Set the controls
    # Physical continuity (2.2): T=p, J_Q.n_heat=q.n_wave, with
    # each Normal outward from its OWN subdomain. The two interface
    # controls impose U_T=Y_w and U_w=Y_T in both region assignments.
    # For the mixed/grad-grad formulation above, the signed interface
    # power is -U_T*Y_T + U_w*Y_w = 0. No disk/annulus sign factor.
    hw.set_control("Interface Heat", "Y_w")
    hw.set_control("Interface Wave", "Y_T")
    # === Dirichlet boundary condition
    hw.set_control("Boundary", "0.")

    # Set the original physical fields with the default true FE L2
    # projection (not a nodal interpolation).
    hw.set_initial_value("T", "5.*np.exp(-25*((x-0.6)*(x-0.6)+y*y))")
    hw.set_initial_value("p", "5.*np.exp(-25*((x-0.6)*(x-0.6)+y*y))")
    hw.set_initial_value("q", "[0.,0.]")

    # For a wave annulus, the Gaussian projected in CG1 has nonzero
    # exterior-boundary coefficients, while homogeneous Dirichlet
    # requires p=0 there. The strict algebraic DAE initializer must
    # REJECT that initial field, not silently alter its dynamic state.
    #
    # Both remedies were previously native-qualified on this exact
    # Concentric Heat-Wave system: (1) clamp only boundary CG1 DOFs,
    # or (2) mass-orthogonal constrained L2 correction. The default
    # is the minimal transparent clamp; both keep T and q unchanged.
    # When heat is outside, there is no exterior wave Dirichlet trace.
    if initial_wave_boundary not in ("clamp", "constrained_l2", "unconstrained"):
        raise ValueError("initial_wave_boundary must be clamp, constrained_l2 "
                         "or unconstrained")
    hw.heat_wave_initialization_report = {
        "method": "not_required_heat_outer" if not outer_is_wave
                  else "unconstrained_historical_diagnostic",
        "boundary_compatibility_enforced": False,
    }
    if outer_is_wave and initial_wave_boundary != "unconstrained":
        original_T = np.asarray(hw.gf_model.variable("T")).reshape(-1).copy()
        original_q = np.asarray(hw.gf_model.variable("q")).reshape(-1).copy()
        report = enforce_zero_dirichlet_initial_values(
            hw, "p", volume_region=wave_region, boundary_region=20,
            method=initial_wave_boundary, factor_solver="mumps",
        )
        if (not np.array_equal(
                np.asarray(hw.gf_model.variable("T")).reshape(-1), original_T
            ) or not np.array_equal(
                np.asarray(hw.gf_model.variable("q")).reshape(-1), original_q
            )):
            raise RuntimeError(
                "Wave boundary compatibility changed the heat/stress initial data"
            )
        report["boundary_compatibility_enforced"] = True
        hw.heat_wave_initialization_report = report

    ## Solve in time
    # CN is the second-order time integrator used for the paper's
    # energy-decay results; fixed accepted steps make the model's
    # dissipation work independently checkable step by step.
    # No explicit solver => SCRIMP's qualified MUMPS default.
    # Other PETSc LU backends remain an opt-in user choice.
    scheme = {}
    if factor_solver is not None:
        scheme["pc_factor_mat_solver_type"] = factor_solver
    hw.set_time_scheme(
        # Use the production constraint-consistent initializer for the
        # index-2 coupled pH-DAE; do not hide an inconsistent Gaussian.
        init_step_ts_type="algebraic",
        ts_type="cn",
        ts_adapt_type="none",
        t_f=float(final_time),
        dt=float(time_step),
        dt_save=float(save_step),
                       ksp_type="preonly",
                       pc_type="lu",
                       **scheme,
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
    
    # A true log-log diagnostic must omit t=0, rather than plotting
    # zero on a logarithmic axis. Keep this separate from the energy
    # identity plot used for the physical qualification.
    import matplotlib.pyplot as plt
    t = np.asarray(hw.solution["t"], dtype=float)
    Hamiltonian = np.asarray(hw.get_Hamiltonian(), dtype=float)
    fig, ax = plt.subplots(figsize=(8, 5))
    positive = (t > 0.) & (Hamiltonian > 0.)
    ax.loglog(t[positive], Hamiltonian[positive])
    ax.grid(axis="both")
    ax.set_xlabel("time t")
    ax.set_ylabel("Hamiltonian")
    ax.set_title("Evolution of the Hamiltonian (log-log)")
    plt.show()
    plt.close(fig)
    
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
    
