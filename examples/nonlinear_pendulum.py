# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Nonlinear damped pendulum, modeled as a lumped port-Hamiltonian system.

States are the angle theta (radians) and angular momentum p = I*theta_dot.
For a point mass m at distance L from the pivot, the governing
relationships are::

    I = m*L**2
    H = p**2/(2*I) + m*g*L*(1-cos(theta))
    e_theta = m*g*L*sin(theta),  e_p = p/I
    theta_dot = e_p
    p_dot = -e_theta - damping*e_p + tau(t)

This example keeps the full nonlinear sine restoring torque. It does
not use the small-angle approximation sin(theta) = theta.
"""

import numpy as np
from pathlib import Path
import scrimp as S


def nonlinear_pendulum(final_time=12.0, time_step=0.01,
                       initial_angle=1.2, torque_amplitude=0.0,
                       plot=False, *, factor_solver=None,
                       output_dir=None):
    """Integrate a nonlinear (possibly forced) pendulum with damping.

    Initialize algebraic co-energies from the exact nonlinear sine
    constitutive law, preserving prescribed angle and momentum.
    """
    mass = 1.0
    length = 1.0
    gravity = 9.81
    rotational_damping = 0.08
    inertia = mass * length**2
    interval_length = 1.0

    pend = S.DPHS("real")
    pend.set_domain(S.Domain("Interval", {"L": interval_length, "h": 1.0}))

    angle = S.State("theta", "Angular position (rad)", "lumped")
    momentum = S.State("p", "Angular momentum", "lumped")
    for state in (angle, momentum):
        pend.add_state(state)
    pend.add_costate(S.CoState("e_theta", "Restoring torque", angle))
    pend.add_costate(S.CoState("e_p", "Angular velocity", momentum))
    pend.add_lumped_port("theta")
    pend.add_lumped_port("p")

    for name, value in {
        "ell": interval_length,
        "inertia": inertia,
        "gravity_torque": mass*gravity*length,
        "damp": rotational_damping,
    }.items():
        pend.gf_model.add_macro(name, str(float(value)))
    pend.gf_model.add_macro(
        "torque", f"{float(torque_amplitude)}*sin(2*pi*0.3*t)"
    )

    bricks = [
        S.Brick("M_theta", "theta*Test_theta/ell", [1],
                dt=True, position="flow"),
        S.Brick("M_p", "p*Test_p/ell", [1],
                dt=True, position="flow"),

        # Skew-symmetric canonical coupling.
        S.Brick("J_theta_p", "e_p*Test_theta/ell", [1],
                position="effort"),
        S.Brick("J_p_theta", "-e_theta*Test_p/ell", [1],
                position="effort"),
        S.Brick("Viscous_damping", "-damp*e_p*Test_p/ell", [1],
                position="effort"),
        S.Brick("External_torque", "torque*Test_p/ell", [1],
                position="effort", linear=False),

        # Co-energy constitutive relations.
        S.Brick("CR_theta_mass", "-e_theta*Test_e_theta/ell", [1]),
        S.Brick("CR_theta_nonlinear",
                "gravity_torque*sin(theta)*Test_e_theta/ell",
                [1], linear=False),
        S.Brick("CR_momentum", "(p/inertia-e_p)*Test_e_p/ell", [1]),
    ]
    for brick in bricks:
        pend.add_brick(brick)

    theta0 = float(initial_angle)
    p0 = 0.0
    if not np.isfinite(theta0):
        raise ValueError("initial_angle must be finite")
    pend.set_initial_value("theta", theta0)
    pend.set_initial_value("p", p0)
    # True one-DOF co-energies: full sin(theta0), never linearized.
    # Both constraints are satisfied at t=0 without moving states.
    pend.gf_model.set_variable(
        "e_theta", np.array([mass*gravity*length*np.sin(theta0)])
    )
    pend.gf_model.set_variable("e_p", np.array([p0/inertia]))

    # None leaves the qualified MUMPS default unchanged.
    # Non-default backends are chosen explicitly by example callers.
    scheme = {}
    if factor_solver is not None:
        scheme["pc_factor_mat_solver_type"] = factor_solver
    # No pseudo-time step: the algebraic values are already consistent.
    pend.set_time_scheme(
        ts_type="cn",
        ts_exact_final_time="matchstep",
        ts_adapt_type="none",
        ksp_type="gmres", pc_type="lu",
        t_0=0.0, t_f=float(final_time),
        dt=float(time_step), dt_save=float(time_step),
        # A long nonlinear CN orbit must satisfy the actual
        # algebraic sine law at EVERY saved time, not merely
        # report TS/SNES convergence under loose default
        # relative tolerances. This tightens Newton convergence
        # without changing the original pH equations or CN.
        snes_atol=1e-12, snes_rtol=1e-12, snes_stol=1e-14,
        init_step=False,
        **scheme,
    )
    pend.solve()

    pend.hamiltonian.set_name("Nonlinear pendulum mechanical energy")
    pend.hamiltonian.add_term(
        S.Term("Gravitational potential",
               "gravity_torque*(1-cos(theta))/ell", [1])
    )
    pend.hamiltonian.add_term(
        S.Term("Kinetic", "p*p/(2*inertia*ell)", [1])
    )

    if plot:
        variant = "free" if float(torque_amplitude) == 0.0 else "forced"
        destination = (
            Path(output_dir) if output_dir is not None
            else Path.cwd() / "outputs" / "showcase" / "nonlinear_pendulum"
            / variant
        )
        pend.pendulum_visuals = render_pendulum_visuals(
            pend, destination, inertia=inertia,
            gravitational_torque=mass*gravity*length,
            damping=rotational_damping, length=length,
            torque_amplitude=float(torque_amplitude),
        )
    return pend


def render_pendulum_visuals(
    system, folder, *, inertia, gravitational_torque, damping,
    length, torque_amplitude,
):
    """Render *original* accepted lumped GetFEM states in the example.

    No extra ODE integration, spatial FE simulation, interpolated
    states, or fictitious exact CN Hamiltonian conservation identity.
    Figures are saved so that direct command-line runs work headlessly.
    The motion GIF is explicitly a mechanical schematic.
    """
    import io
    import json

    import matplotlib.pyplot as plt
    from PIL import Image
    from petsc4py import PETSc

    t = np.asarray(system.solution["t"], dtype=float)
    theta = np.asarray(system.get_solution("theta"), dtype=float).ravel()
    p = np.asarray(system.get_solution("p"), dtype=float).ravel()
    e_theta = np.asarray(system.get_solution("e_theta"), dtype=float).ravel()
    e_p = np.asarray(system.get_solution("e_p"), dtype=float).ravel()
    H = np.asarray(system.get_Hamiltonian(), dtype=float).ravel()
    n = len(t)
    if (n < 3 or any(v.shape != (n,) for v in
                     (theta, p, e_theta, e_p, H))
            or not all(np.all(np.isfinite(v)) for v in
                       (t, theta, p, e_theta, e_p, H))
            or not np.all(np.diff(t) > 0)):
        raise RuntimeError("Original nonlinear GetFEM history is incomplete")
    dt = float(t[1] - t[0])
    if not np.allclose(np.diff(t), dt, rtol=0., atol=1e-10):
        raise RuntimeError("Pendulum CN visualization needs fixed dt")
    if (system.time_scheme.getString("ts_type","").lower() != "cn"
            or system.last_ts_converged_reason <= 0
            or system.last_ts_steps != n-1):
        raise RuntimeError("Refusing to render unconverged pendulum")

    # The source force and constitutive laws, not an approximation.
    torque = torque_amplitude * np.sin(2*np.pi*0.3*t)
    dissipation = damping * e_p**2
    input_power = torque * e_p
    analytic_H = (p*p/(2*inertia) +
                  gravitational_torque*(1-np.cos(theta)))
    constitutive_error = max(
        float(np.max(np.abs(e_theta -
                            gravitational_torque*np.sin(theta)))),
        float(np.max(np.abs(e_p - p/inertia))),
    )
    H_error = float(np.max(np.abs(H-analytic_H)))
    if constitutive_error > 1e-8 or H_error > 1e-8:
        # Do not silently weaken the physical accuracy gate.
        # Preserve the authentic four-DOF trajectory for a
        # subsequent scientific diagnosis even if the renderer
        # refuses to draw physically inconsistent curves.
        err_theta=np.abs(e_theta-gravitational_torque*np.sin(theta))
        at=int(np.argmax(err_theta))
        if PETSc.COMM_WORLD.getRank()==0:
            folder=Path(folder)
            folder.mkdir(parents=True,exist_ok=True)
            np.savez_compressed(
                folder/"failed_native_constitutive_trace.npz",
                t=t,theta=theta,p=p,e_theta=e_theta,e_p=e_p,H=H,
                nonlinear_constitutive_defect=(
                    e_theta-gravitational_torque*np.sin(theta)
                ),
            )
        raise RuntimeError(
            "Nonlinear GetFEM constitutive law/Hamiltonian failed: "
            f"max_sine_error={constitutive_error:.12g}, "
            f"max_H_error={H_error:.12g}, "
            f"worst_index={at}, worst_time={t[at]:.12g}, "
            f"worst_theta={theta[at]:.12g}, "
            f"worst_e_theta={e_theta[at]:.12g}; "
            "authentic trajectory preserved in "
            "failed_native_constitutive_trace.npz"
        )

    # PETSc TSCN is the implicit trapezoidal/end-point theta rule,
    # NOT the implicit midpoint rule. Evaluate the nonlinear forces
    # at each accepted end-point, then average for the CN step.
    # Both true nonlinear end-point co-energies are checked above.
    delta_theta=np.diff(theta)/dt
    delta_p=np.diff(p)/dt
    mean_ep=(e_p[1:]+e_p[:-1])/2.
    mean_etheta=(e_theta[1:]+e_theta[:-1])/2.
    mean_torque=(torque[1:]+torque[:-1])/2.
    cn_dissipation=damping*mean_ep**2
    cn_input_power=mean_torque*mean_ep
    cn_ode_theta_defect=delta_theta-mean_ep
    cn_ode_p_defect=delta_p+mean_etheta+damping*mean_ep-mean_torque
    cn_coenergy_rate=mean_etheta*delta_theta+mean_ep*delta_p
    cn_power_defect=cn_coenergy_rate+cn_dissipation-cn_input_power
    max_cn_power_defect=float(np.max(np.abs(cn_power_defect)))
    if (np.max(np.abs(cn_ode_theta_defect))>1e-7
            or np.max(np.abs(cn_ode_p_defect))>1e-7
            or max_cn_power_defect>1e-7):
        raise RuntimeError(
            "Original pendulum trapezoidal CN equation/power identity failed: "
            +str(max_cn_power_defect)
        )

    # Nonquadratic H: average nonlinear e_theta is generally
    # NOT the exact discrete gradient of 1-cos(theta).
    # This exact chain-rule defect is the difference between the
    # observed Hamiltonian rate and the CN co-energy pairing.
    cn_H_rate=np.diff(H)/dt
    cn_chain_rule_defect=cn_H_rate-cn_coenergy_rate
    cn_discrete_work=np.r_[0., np.cumsum(
        dt*(cn_input_power-cn_dissipation)
    )]
    energy_change=H-H[0]
    cn_cumulative_ledger_difference=energy_change-cn_discrete_work
    # The identity below follows from the CN residual plus the
    # exact definition of the nonlinear chain-rule defect. It is NOT
    # a claim that either side vanishes.
    if np.max(np.abs(
        cn_cumulative_ledger_difference[1:]
        -np.cumsum(dt*cn_chain_rule_defect)
    ))>1e-7:
        raise RuntimeError("Nonquadratic CN chain-rule accounting failed")

    # A conventional trapezoid quadrature of CONTINUOUS nodal powers
    # is separate from the genuine CN discrete work above.
    work_rate=input_power-dissipation
    work_trapezoidal=np.zeros(n)
    work_trapezoidal[1:]=np.cumsum(
        .5*(work_rate[:-1]+work_rate[1:])*np.diff(t)
    )
    quadrature_difference=energy_change-work_trapezoidal

    # Scalar histories are present on every rank, but the common
    # user-facing output directory is written only on MPI root.
    if PETSc.COMM_WORLD.getRank() != 0:
        return None
    folder = Path(folder)
    figures = folder/"figures"
    figures.mkdir(parents=True, exist_ok=True)

    def save_figure(fig, name):
        fig.savefig(figures/(name+".png"), dpi=180)
        fig.savefig(figures/(name+".svg"))
        plt.close(fig)

    variant = "free" if torque_amplitude == 0. else "forced"
    fig, axes = plt.subplots(
        2,1, figsize=(10.5,6.7), sharex=True,
        constrained_layout=True,
    )
    axes[0].plot(t,theta,label="Angle θ (rad)")
    axes[0].plot(t,e_p,label="Angular velocity p/I")
    axes[0].set(
        title="Nonlinear pendulum · original GetFEM states · "+variant,
        ylabel="Angle / velocity",
    )
    axes[0].legend(loc="best")
    axes[1].plot(t,p,label="Angular momentum p")
    axes[1].plot(t,torque,label="Applied torque τ(t)")
    axes[1].set(xlabel="Time", ylabel="Momentum / torque")
    axes[1].legend(loc="best")
    for ax in axes: ax.grid(True,alpha=.25)
    save_figure(fig,"pendulum_trajectory")

    fig, ax = plt.subplots(figsize=(8.6,6), constrained_layout=True)
    ax.plot(theta,e_p,color="0.75",zorder=1)
    dots = ax.scatter(theta,e_p,c=t,cmap="viridis",s=9,zorder=2)
    fig.colorbar(dots,ax=ax,label="Time")
    ax.plot(theta[0],e_p[0],"o",ms=8,label="t=0")
    ax.plot(theta[-1],e_p[-1],"s",ms=8,label=f"t={t[-1]:g}")
    ax.set(title="Nonlinear pendulum · original CN phase portrait",
           xlabel="Angle θ (rad)",ylabel="Angular velocity p/I")
    ax.grid(True,alpha=.25)
    ax.legend()
    save_figure(fig,"pendulum_phase_portrait")

    fig, ax = plt.subplots(2,2, figsize=(12.5,8),
                           constrained_layout=True)
    ax[0,0].plot(t,H,label="Native GetFEM Hamiltonian")
    ax[0,0].set(title="Nonquadratic Hamiltonian",ylabel="Energy")
    ax[0,0].legend()
    ax[0,1].plot(t,dissipation,label="Viscous dissipation d·(p/I)²")
    ax[0,1].plot(t,input_power,label="External torque power τ·p/I")
    ax[0,1].set(title="Continuous physical power terms",ylabel="Power")
    ax[0,1].legend()
    step_times=t[1:]
    ax[1,0].plot(step_times,cn_coenergy_rate,
                 label="CN averaged co-energy · Δz/Δt")
    ax[1,0].plot(step_times,cn_input_power-cn_dissipation,
                 "--",label="CN input power − CN dissipation")
    ax[1,0].set(title="CN co-energy power identity (all steps)",
                ylabel="Discrete co-energy rate")
    ax[1,0].legend()
    ax[1,1].plot(t,energy_change,label="H(t)−H(0)")
    ax[1,1].plot(t,cn_discrete_work,"--",
                 label="Cumulative CN port work − dissipation")
    ax[1,1].set(
        title="CN work vs. nonquadratic Hamiltonian",
        ylabel="Energy variation (chain-rule defect ≠ 0)",
    )
    ax[1,1].legend()
    for a in ax.flat:
        a.grid(True,alpha=.25)
        a.set_xlabel("Time")
    fig.suptitle("Nonlinear pendulum · physical vs numerical energy")
    save_figure(fig,"pendulum_energy_and_power")

    # Geometry is illustrative: x=L sin(theta), y=-L cos(theta)
    # at actual accepted time samples only; no spatial FE claim.
    selected = np.unique(np.linspace(
        0,n-1,min(n,36),dtype=int,
    ))
    gif_frames=[]
    for i in selected:
        fig, ax = plt.subplots(figsize=(6.8,5.5))
        x = length*np.sin(theta[i])
        y = -length*np.cos(theta[i])
        ax.plot([0,x],[0,y],color="0.22",linewidth=3)
        ax.plot(x,y,"o",ms=19,color="0.42")
        ax.plot(0,0,"o",ms=6,color="0.10")
        ax.plot([-1.35*length,1.35*length],[-1.2*length,-1.2*length],
                color="0.8",linewidth=1)
        ax.set(xlim=(-1.35*length,1.35*length),
               ylim=(-1.30*length,1.30*length),
               title=(f"t={t[i]:.2f}   θ={theta[i]:+.4f} rad   "
                      f"ω={e_p[i]:+.4f}"))
        ax.text(-1.27*length,1.08*length,
                "Schematic pendulum · accepted GetFEM state",
                fontsize=9)
        ax.set_aspect("equal")
        ax.axis("off")
        fig.tight_layout()
        buffer=io.BytesIO()
        fig.savefig(buffer,format="png",dpi=95)
        plt.close(fig)
        buffer.seek(0)
        with Image.open(buffer) as frame:
            gif_frames.append(frame.convert("RGB").copy())
    gif_frames[0].save(
        figures/"pendulum_motion.gif",
        save_all=True,append_images=gif_frames[1:],
        duration=160,loop=0,optimize=False,
    )

    np.savez_compressed(
        folder/"native_lumped_time_history.npz",
        t=t,theta=theta,p=p,e_theta=e_theta,e_p=e_p,H=H,
        torque=torque,input_power=input_power,
        viscous_dissipation=dissipation,
        CN_coenergy_rate=cn_coenergy_rate,
        CN_power_defect=cn_power_defect,
        CN_ODE_theta_defect=cn_ode_theta_defect,
        CN_ODE_p_defect=cn_ode_p_defect,
        CN_step_input_power=cn_input_power,
        CN_step_dissipation=cn_dissipation,
        CN_chain_rule_defect=cn_chain_rule_defect,
        CN_discrete_work=cn_discrete_work,
        CN_cumulative_ledger_difference=cn_cumulative_ledger_difference,
        trapezoidal_work_minus_dissipation=work_trapezoidal,
        true_Hamiltonian_change=energy_change,
        H_minus_trapezoidal_work=quadrature_difference,
    )
    manifest={
        "schema":"scrimp-v130-pendulum-canonical-figures-v1",
        "status":"EXAMPLE_NATIVE_GETFEM_NONLINEAR_PENDULUM_VISUALS",
        "canonical_example":"examples/nonlinear_pendulum.py",
        "original_saved_GetFEM_DAE_states":n,
        "original_DAE_DOFs":4,
        "torque_amplitude":float(torque_amplitude),
        "variant":variant,
        "GIF_genuine_accepted_frames":len(selected),
        "no_interpolated_or_spatial_FE_states":True,
        "motion_is_schematic":True,
        "time_stepper":"cn",
        "CN_is_trapezoidal_not_implicit_midpoint":True,
        "nonquadratic_H_CN_chain_rule_defect_not_assumed_zero":True,
        "trapezoidal_quadrature_is_not_exact_CN_ledger":True,
        "max_nonlinear_coenergy_constitutive_error":constitutive_error,
        "max_native_Hamiltonian_error":H_error,
        "max_CN_step_coenergy_power_defect":max_cn_power_defect,
        "max_CN_nonlinear_H_chain_rule_defect":float(
            np.max(np.abs(cn_chain_rule_defect))
        ),
        "figures":[
            "pendulum_trajectory.png","pendulum_trajectory.svg",
            "pendulum_phase_portrait.png","pendulum_phase_portrait.svg",
            "pendulum_energy_and_power.png","pendulum_energy_and_power.svg",
            "pendulum_motion.gif",
        ],
    }
    (folder/"manifest.json").write_text(
        json.dumps(manifest,indent=2)+"\n",encoding="utf-8"
    )
    return manifest


if __name__ == "__main__":
    nonlinear_pendulum(plot=True)
