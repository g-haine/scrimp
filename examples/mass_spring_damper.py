# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Mass--spring--damper: a finite-dimensional port-Hamiltonian example.

State x=(q,p), with displacement q and momentum p=m*q_dot.
    H(q,p) = k*q**2/2 + p**2/(2*m)
    e_q = dH/dq = k*q,  e_p = dH/dp = p/m
    q_dot = e_p
    p_dot = -e_q - c*e_p + f(t)

Power balance: H_dot = -c*e_p**2 + f*e_p.

Each lumped port has a **single** DOF. The interval used here is only
a normalization domain for GetFEM's weak-form assembly; it is NOT a
spatial discretization of the oscillator.
"""

import numpy as np
from pathlib import Path
import scrimp as S


def mass_spring_damper(
    final_time=8.0, time_step=0.01, plot=False, *,
    factor_solver=None, initialization="exact",
    output_dir=None,
):
    """Simulate the forced, damped oscillator using lumped dynamic ports.

    Returns
    -------
    scrimp.DPHS
        Solved system with Hamiltonian terms and time histories.

    By default, impose the *exact* constitutive co-energies at t=0
    and preserve the prescribed differential state (q0,p0). The
    historical one-step pseudo-time relaxation is reproduced with
    initialization="legacy_pseudo"; it moves the physical initial
    momentum, and is therefore not the qualified default.
    """
    if initialization not in ("exact", "legacy_pseudo"):
        raise ValueError("initialization must be 'exact' or 'legacy_pseudo'")
    mass = 1.5
    stiffness = 5.0
    damping = 0.25
    interval_length = 1.0

    sys = S.DPHS("real")
    sys.set_domain(S.Domain("Interval", {"L": interval_length, "h": 1.0}))

    q = S.State("q", "Displacement", "lumped")
    p = S.State("p", "Linear momentum", "lumped")
    for state in (q, p):
        sys.add_state(state)
    sys.add_costate(S.CoState("e_q", "Spring force", q))
    sys.add_costate(S.CoState("e_p", "Velocity", p))

    # Unlike FEM("q", ...), this registers a single scalar DOF.
    sys.add_lumped_port("q")
    sys.add_lumped_port("p")

    for key, val in {
        "ell": interval_length,
        "osc_mass": mass,
        "spring_k": stiffness,
        "damp_c": damping,
    }.items():
        sys.gf_model.add_macro(key, str(float(val)))
    sys.gf_model.add_macro("force", "0.2*sin(2*pi*t)")

    bricks = [
        S.Brick("M_q", "q*Test_q/ell", [1], dt=True, position="flow"),
        S.Brick("M_p", "p*Test_p/ell", [1], dt=True, position="flow"),
        # Canonical Dirac interconnection J = [[0,1],[-1,0]]
        S.Brick("J_qp", "e_p*Test_q/ell", [1], position="effort"),
        S.Brick("J_pq", "-e_q*Test_p/ell", [1], position="effort"),
        S.Brick("Damping", "-damp_c*e_p*Test_p/ell",
                [1], position="effort"),
        # A time-only source is classified non-linear so that it is
        # assembled on every residual evaluation.
        S.Brick("External_force", "force*Test_p/ell", [1],
                linear=False, position="effort"),
        S.Brick("Spring_law", "(spring_k*q-e_q)*Test_e_q/ell", [1]),
        S.Brick("Momentum_law", "(p/osc_mass-e_p)*Test_e_p/ell", [1]),
    ]
    for brick in bricks:
        sys.add_brick(brick)

    q0, p0 = 0.4, 0.0
    sys.set_initial_value("q", q0)
    sys.set_initial_value("p", p0)
    if initialization == "exact":
        # These are true zero-dimensional GetFEM variables (one DOF
        # each), NOT FE fields to interpolate. The algebraic laws are
        # e_q = dH/dq = k*q and e_p = dH/dp = p/m. Initialize both
        # exactly rather than evolving (q,p) during a pseudo-TS.
        sys.gf_model.set_variable("e_q", np.array([stiffness*q0]))
        sys.gf_model.set_variable("e_p", np.array([p0/mass]))
    # Explicitly expose the mode for physical acceptance diagnostics.
    sys.mass_spring_initialization = initialization
    # An omitted factor_solver leaves the current SCRIMP/PETSc default
    # (MUMPS) untouched. Users can explicitly select superlu (serial),
    # superlu_dist (MPI), or another PETSc-supported LU backend.
    scheme = {}
    if factor_solver is not None:
        scheme["pc_factor_mat_solver_type"] = factor_solver
    # With fully consistent initial differential AND algebraic DOFs
    # no initial TS is necessary. Keep historical pseudo only by
    # explicit opt-in; the physical BDF2 time integrator is unchanged.
    init_options = (
        {"init_step": False}
        if initialization == "exact"
        else {"init_step_ts_type": "pseudo", "init_step": True}
    )
    sys.set_time_scheme(
        ts_type="bdf", ts_bdf_order=2,
        ts_exact_final_time="matchstep",
        ts_adapt_type="none",
        ksp_type="gmres", pc_type="lu",
        t_0=0.0, t_f=float(final_time),
        dt=float(time_step), dt_save=float(time_step),
        **init_options,
        **scheme,
    )
    sys.solve()

    sys.hamiltonian.set_name("Mass--spring--damper energy")
    sys.hamiltonian.add_term(
        S.Term("Spring potential", "spring_k*q*q/(2*ell)", [1])
    )
    sys.hamiltonian.add_term(
        S.Term("Kinetic", "p*p/(2*osc_mass*ell)", [1])
    )

    if plot:
        destination=(
            Path(output_dir) if output_dir is not None
            else Path.cwd()/"outputs"/"showcase"/"mass_spring_damper"
        )
        # Original GetFEM/PETSc accepted states, not probe-generated art.
        sys.mass_spring_visuals = render_mass_spring_visuals(
            sys,destination,mass=mass,stiffness=stiffness,damping=damping
        )
    return sys


def render_mass_spring_visuals(system, folder, *, mass, stiffness, damping):
    """Original lumped-DOF time series, physical power and GIF.

    Only the saved PETSc/GetFEM solution drives every figure and the
    mechanical schematic. This rendering is part of the canonical
    user example, not a research or acceptance helper. It never
    invents an FEM space for a zero-dimensional oscillator.

    BDF2 energy-power identity is evaluated at k>=2 only; the PETSc
    internal startup half-step is not an ordinary BDF1 interval.
    """
    import io
    import json
    from petsc4py import PETSc
    import matplotlib.pyplot as plt
    from PIL import Image

    times=np.asarray(system.solution["t"],dtype=float)
    q=np.asarray(system.get_solution("q"),dtype=float).ravel()
    p=np.asarray(system.get_solution("p"),dtype=float).ravel()
    eq=np.asarray(system.get_solution("e_q"),dtype=float).ravel()
    ep=np.asarray(system.get_solution("e_p"),dtype=float).ravel()
    H=np.asarray(system.get_Hamiltonian(),dtype=float).ravel()
    n=len(times)
    if (n<3 or any(x.shape!=(n,) for x in (q,p,eq,ep,H))
            or not all(np.all(np.isfinite(x))
                       for x in (times,q,p,eq,ep,H))
            or not np.all(np.diff(times)>0.)):
        raise RuntimeError("Missing authentic finite-dimensional GetFEM history")
    if not np.allclose(np.diff(times),np.diff(times)[0],rtol=0.,atol=1e-10):
        raise RuntimeError("Canonical mass-spring BDF2 visualization requires fixed dt")
    dt=float(times[1]-times[0])
    if system.last_ts_converged_reason<=0 or system.last_ts_steps!=n-1:
        raise RuntimeError("Cannot visualize an unconverged lumped simulation")

    force=0.2*np.sin(2.*np.pi*times)
    velocity=p/mass
    dissipation=damping*ep*ep
    input_power=force*ep
    analytic_H=0.5*stiffness*q*q+0.5*p*p/mass

    # Canonical lumped co-energy consistency; check each true saved
    # GetFEM variable before drawing potentially misleading plots.
    coenergy_error=max(float(np.max(np.abs(eq-stiffness*q))),
                       float(np.max(np.abs(ep-velocity))))
    Hamiltonian_error=float(np.max(np.abs(H-analytic_H)))
    if coenergy_error>1e-8 or Hamiltonian_error>1e-8:
        raise RuntimeError(
            "Authentic GetFEM oscillator violates its constitutive law "
            "or Hamiltonian: "+repr((coenergy_error,Hamiltonian_error))
        )

    def energy(position,momentum):
        return .5*stiffness*position*position+momentum*momentum/(2.*mass)

    # Discrete quadratic BDF2 chain-rule identity, not trapezoidal
    # integration of physical work or a claim about PETSc startup.
    a_q,a_p=q[2:],p[2:]
    b_q,b_p=q[1:-1],p[1:-1]
    c_q,c_p=q[:-2],p[:-2]
    bdf2_energy_rate=(
        energy(a_q,a_p)-energy(b_q,b_p)
        +energy(2*a_q-b_q,2*a_p-b_p)
        -energy(2*b_q-c_q,2*b_p-c_p)
        +energy(a_q-2*b_q+c_q,a_p-2*b_p+c_p)
    )/(2.*dt)
    balance_defect=bdf2_energy_rate+dissipation[2:]-input_power[2:]
    max_defect=float(np.max(np.abs(balance_defect)))
    if max_defect>1e-8:
        raise RuntimeError(
            "Original BDF2 oscillator step-power identity failed: "
            +str(max_defect)
        )

    # These scalar quantities exist on every MPI rank. Only the root
    # writes the shared figures; no need to render one copy per rank.
    figures=Path(folder)/"figures"
    if PETSc.COMM_WORLD.getRank()!=0:
        return None
    figures.mkdir(parents=True,exist_ok=True)

    def save_both(fig,stem):
        fig.savefig(figures/(stem+".png"),dpi=180)
        fig.savefig(figures/(stem+".svg"))
        plt.close(fig)

    fig,(ax1,ax2)=plt.subplots(
        2,1,figsize=(10.5,6.5),sharex=True,constrained_layout=True
    )
    ax1.plot(times,q,label="q — displacement")
    ax1.plot(times,ep,label="p/m — velocity")
    ax1.set(ylabel="State / velocity",
            title="Mass–spring–damper · original lumped GetFEM trajectory")
    ax1.legend(loc="best")
    ax2.plot(times,force,label="f(t) — external force")
    ax2.plot(times,p,label="p — momentum")
    ax2.set(xlabel="Time",ylabel="Input / momentum")
    ax2.legend(loc="best")
    for ax in (ax1,ax2):ax.grid(True,alpha=.25)
    save_both(fig,"mass_spring_trajectory")

    fig,ax=plt.subplots(figsize=(8.5,5.8),constrained_layout=True)
    ax.plot(q,ep,color="0.75",linewidth=1.2,zorder=1)
    sc=ax.scatter(q,ep,c=times,cmap="viridis",s=13,zorder=2)
    ax.plot(q[0],ep[0],"o",ms=8,label="t=0")
    ax.plot(q[-1],ep[-1],"s",ms=8,label=f"t={times[-1]:g}")
    fig.colorbar(sc,ax=ax,label="Time")
    ax.set(xlabel="Displacement q",ylabel="Velocity p/m",
           title="Phase portrait · accepted original BDF2 states")
    ax.grid(True,alpha=.25)
    ax.legend()
    save_both(fig,"mass_spring_phase_portrait")

    fig,axes=plt.subplots(
        2,2,figsize=(12.5,7.8),constrained_layout=True
    )
    axes[0,0].plot(times,H,label="H(t) · genuine GetFEM")
    axes[0,0].set(title="Hamiltonian",ylabel="Energy")
    axes[0,0].legend()
    axes[0,1].plot(times,dissipation,label="c (p/m)²")
    axes[0,1].plot(times,input_power,label="f(t) p/m")
    axes[0,1].set(title="Physical power terms",ylabel="Power")
    axes[0,1].legend()
    axes[1,0].plot(times[2:],bdf2_energy_rate,
                   label="Quadratic BDF2 energy rate")
    axes[1,0].plot(times[2:],input_power[2:]-dissipation[2:],
                   "--",label="Input − dissipation")
    axes[1,0].set(
        title="Exact BDF2 step identity (k ≥ 2)",
        ylabel="Energy rate"
    )
    axes[1,0].legend()
    axes[1,1].plot(times[2:],balance_defect,
                   label="BDF2 energy-power defect")
    axes[1,1].set(title="Numerical identity residual",
                  ylabel="Residual")
    axes[1,1].legend()
    for ax in axes.flat:
        ax.set_xlabel("Time")
        ax.grid(True,alpha=.25)
    fig.suptitle("Mass–spring–damper · structure-preserving diagnostics")
    save_both(fig,"mass_spring_energy_and_power")

    # Illustrative mass/spring geometry, explicitly NOT a spatial FEM.
    # Each frame uses one accepted original lumped state. Uniform
    # sample selection does not interpolate time or motion.
    selected=np.unique(np.linspace(
        0,n-1,min(n,36),dtype=int
    ))
    gif_frames=[]
    for i in selected:
        fig,ax=plt.subplots(figsize=(8.2,3.6))
        qi=float(q[i])
        center=1.3+qi
        coil_x=np.linspace(0.25,center-0.13,140)
        coil_y=0.17*np.sin(np.linspace(0,12*np.pi,len(coil_x)))
        coil_y[[0,-1]]=0
        ax.axvline(0.24,color="0.30",linewidth=6)
        ax.plot(coil_x,coil_y,color="0.20",linewidth=2)
        ax.add_patch(plt.Rectangle(
            (center-0.13,-0.25),0.26,0.5,
            facecolor="0.55",edgecolor="0.10"
        ))
        ax.plot([0.1,2.6],[-0.30,-0.30],color="0.55")
        ax.set(xlim=(0,2.6),ylim=(-0.48,0.55),
               title=(f"t={times[i]:.2f}   q={qi:+.4f}   "
                      f"p/m={ep[i]:+.4f}"))
        ax.text(.1,.42,"Schematic only · position = 1.3 + q",
                fontsize=9)
        ax.axis("off")
        fig.tight_layout()
        buffer=io.BytesIO()
        fig.savefig(buffer,format="png",dpi=100)
        plt.close(fig)
        buffer.seek(0)
        with Image.open(buffer) as image:
            gif_frames.append(image.convert("RGB").copy())
    gif_frames[0].save(
        figures/"mass_spring_motion.gif",
        save_all=True,append_images=gif_frames[1:],
        duration=160,loop=0,optimize=False
    )

    history=Path(folder)/"native_lumped_time_history.npz"
    np.savez_compressed(
        history,time=times,q=q,p=p,e_q=eq,e_p=ep,H=H,
        force=force,dissipation=dissipation,input_power=input_power,
        BDF2_step_defect=balance_defect,
    )
    manifest={
        "schema":"scrimp-v130-mass-spring-native-figures-v1",
        "visual_status":"GENERATED_BY_CANONICAL_EXAMPLE_FROM_GETFEM_LUMPED_STATES",
        "original_example":"examples/mass_spring_damper.py",
        "native_GetFEM_saved_states":n,
        "native_original_full_DAE_DOFS":4,
        "time_final":float(times[-1]),
        "no_spatial_FEM_or_time_interpolation":True,
        "motion_is_explicitly_a_schematic":True,
        "original_accepted_states_used_for_GIF":len(selected),
        "max_coenergy_error":coenergy_error,
        "max_Hamiltonian_error":Hamiltonian_error,
        "max_exact_BDF2_power_defect_k_ge_2":max_defect,
        "BDF2_startup_energy_identity_not_claimed":True,
        "figures":[
            "mass_spring_trajectory.png","mass_spring_trajectory.svg",
            "mass_spring_phase_portrait.png","mass_spring_phase_portrait.svg",
            "mass_spring_energy_and_power.png",
            "mass_spring_energy_and_power.svg","mass_spring_motion.gif",
        ],
    }
    (Path(folder)/"manifest.json").write_text(
        json.dumps(manifest,indent=2)+"\n",encoding="utf-8"
    )
    return manifest


if __name__ == "__main__":
    mass_spring_damper(plot=True)
