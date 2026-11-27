# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Native physical qualification of the canonical nonlinear pendulum in CN.

Use the actual original GetFEM lumped sine/bricks, no substitute model.
PETSc ts_type=cn is the endpoint trapezoidal Theta=1/2 scheme.
For each saved interval (n-1,n), the exact differential residuals are
  (theta_n-theta_old)/dt = (e_p_n+e_p_old)/2
  (p_n-p_old)/dt = -(e_theta_n+e_theta_old)/2
                           -d*(e_p_n+e_p_old)/2+(tau_n+tau_old)/2.
The co-energy pairing is exact for these CN equations, but not an
exact discrete Hamiltonian difference for H containing 1-cos(theta).
The chain-rule defect is separately measured from actual GetFEM H.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from petsc4py import PETSc

from examples.nonlinear_pendulum import nonlinear_pendulum


G = 9.81
I = 1.0
D = 0.08
FREQUENCY = 0.3


def ode(t, x, amplitude):
    theta, p = x
    return np.array([
        p/I,
        -G*np.sin(theta)-D*(p/I)
        +amplitude*np.sin(2*np.pi*FREQUENCY*t),
    ], dtype=float)


def rk4_reference(t, initial, amplitude):
    """Independent high-accuracy ODE integration (not a SCRIMP model)."""
    out = np.empty((len(t), 2), dtype=float)
    x = np.asarray(initial, dtype=float).copy()
    out[0] = x
    for n in range(1, len(t)):
        dt = float(t[n]-t[n-1])/64.
        time = float(t[n-1])
        for _ in range(64):
            k1 = ode(time, x, amplitude)
            k2 = ode(time+dt/2, x+dt*k1/2, amplitude)
            k3 = ode(time+dt/2, x+dt*k2/2, amplitude)
            k4 = ode(time+dt, x+dt*k3, amplitude)
            x = x + dt*(k1+2*k2+2*k3+k4)/6.
            time += dt
        out[n] = x
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--final-time", type=float, default=.08)
    parser.add_argument("--dt", type=float, default=.01)
    parser.add_argument("--angle", type=float, default=1.2)
    parser.add_argument("--torque", type=float, default=0.)
    parser.add_argument("--solver", default="mumps")
    parser.add_argument("--showcase", action="store_true",
                        help="Figures rendered only by the canonical example")
    parser.add_argument("--reference-atol", type=float, default=3e-3)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    system = nonlinear_pendulum(
        final_time=args.final_time, time_step=args.dt,
        initial_angle=args.angle, torque_amplitude=args.torque,
        plot=args.showcase, factor_solver=args.solver,
    )
    t = np.asarray(system.solution["t"], dtype=float)
    z = np.asarray(system.solution["z"], dtype=float)
    theta = np.asarray(system.get_solution("theta"), dtype=float).ravel()
    p = np.asarray(system.get_solution("p"), dtype=float).ravel()
    etheta = np.asarray(system.get_solution("e_theta"), dtype=float).ravel()
    ep = np.asarray(system.get_solution("e_p"), dtype=float).ravel()
    H = np.asarray(system.get_Hamiltonian(), dtype=float).ravel()
    n = len(t)
    steps = round(args.final_time/args.dt)
    for name, values in (
        ("theta",theta), ("p",p), ("e_theta",etheta),
        ("e_p",ep), ("H",H)
    ):
        if values.shape!=(n,) or not np.all(np.isfinite(values)):
            raise RuntimeError("Invalid native GetFEM "+name+" history")
    if (steps < 1 or n!=steps+1 or z.shape!=(n,4)
            or not np.all(np.isfinite(z))
            or np.max(np.abs(t-np.arange(n)*args.dt))>1e-9):
        raise RuntimeError("Original nonlinear CN full DAE/time history malformed")
    if (system.time_scheme.getString("ts_type","").lower()!="cn"
            or system.time_scheme.getString("init_step","true")!="false"):
        raise RuntimeError("Canonical nonlinear pendulum is not exact-start CN")

    tau = args.torque*np.sin(2*np.pi*FREQUENCY*t)
    analytic_H = .5*p*p/I + G*(1.-np.cos(theta))
    dissipation = D*ep**2
    input_power = tau*ep
    # Endpoint trapezoidal CN, NOT midpoint evaluation of nonlinear sin.
    dtheta = np.diff(theta)/args.dt
    dp = np.diff(p)/args.dt
    mean_ep = .5*(ep[1:]+ep[:-1])
    mean_etheta = .5*(etheta[1:]+etheta[:-1])
    mean_tau = .5*(tau[1:]+tau[:-1])
    cn_theta_defect = dtheta-mean_ep
    cn_p_defect = dp+mean_etheta+D*mean_ep-mean_tau
    cn_coenergy_rate = mean_etheta*dtheta+mean_ep*dp
    cn_dissipation = D*mean_ep**2
    cn_input_power = mean_tau*mean_ep
    cn_power_defect = cn_coenergy_rate+cn_dissipation-cn_input_power
    cn_H_rate = np.diff(H)/args.dt
    cn_nonlinear_chain_rule_defect = cn_H_rate-cn_coenergy_rate
    cn_discrete_work = np.r_[0., np.cumsum(
        args.dt*(cn_input_power-cn_dissipation)
    )]
    cn_cumulative_H_minus_work = H-H[0]-cn_discrete_work
    cn_chain_accounting_defect = (
        cn_cumulative_H_minus_work[1:]
        -np.cumsum(args.dt*cn_nonlinear_chain_rule_defect)
    )

    reference = rk4_reference(t, [args.angle, 0.], args.torque)
    max_ref = float(np.max(np.abs(
        np.column_stack((theta,p))-reference
    )))
    residuals = {
        "constitutive_nonlinear_full_sine":
            float(np.max(np.abs(etheta-G*np.sin(theta)))),
        "constitutive_velocity":
            float(np.max(np.abs(ep-p/I))),
        "native_H_vs_analytic":
            float(np.max(np.abs(H-analytic_H))),
        "CN_ODE_theta_all_steps":float(np.max(np.abs(cn_theta_defect))),
        "CN_ODE_p_all_steps":float(np.max(np.abs(cn_p_defect))),
        "CN_coenergy_power_all_steps":float(np.max(np.abs(cn_power_defect))),
        "CN_nonlinear_H_chain_rule_defect_NOT_assumed_zero":
            float(np.max(np.abs(cn_nonlinear_chain_rule_defect))),
        "CN_cumulative_chain_accounting":
            float(np.max(np.abs(cn_chain_accounting_defect))),
        "RK4_full_sine_ODE_max_theta_or_p_error":max_ref,
        "theta_initial_drift":float(theta[0]-args.angle),
        "p_initial_drift":float(p[0]),
        "e_theta_initial_drift":float(etheta[0]-G*np.sin(args.angle)),
        "e_p_initial_drift":float(ep[0]),
        "max_unforced_H_step_increase":
            float(np.max(np.diff(H))) if args.torque == 0. else None,
    }
    gates = {
        "native_CN_completed":
            system.last_ts_converged_reason > 0
            and system.last_ts_steps==steps
            and system.last_ts_rejections==0
            and system.last_ts_snes_failures==0,
        "exact_prescribed_initial_differential_states":
            max(abs(residuals["theta_initial_drift"]),
                abs(residuals["p_initial_drift"]))<1e-12,
        "exact_initial_nonlinear_constitutive_data":
            max(abs(residuals["e_theta_initial_drift"]),
                abs(residuals["e_p_initial_drift"]))<1e-12,
        "full_sine_constitutive":
            residuals["constitutive_nonlinear_full_sine"]<1e-8,
        "velocity_constitutive":
            residuals["constitutive_velocity"]<1e-8,
        "native_nonquadratic_H":
            residuals["native_H_vs_analytic"]<1e-8,
        "CN_theta_differential_equation":
            residuals["CN_ODE_theta_all_steps"]<1e-7,
        "CN_p_differential_equation":
            residuals["CN_ODE_p_all_steps"]<1e-7,
        "CN_coenergy_power_identity":
            residuals["CN_coenergy_power_all_steps"]<1e-7,
        "CN_nonquadratic_chain_rule_accounting":
            residuals["CN_cumulative_chain_accounting"]<1e-7,
        "independent_full_sine_RK4_reference":
            max_ref<args.reference_atol,
        "nonnegative_H_and_dissipation":
            float(np.min(H))>=-1e-10
            and float(np.min(dissipation))>=-1e-14,
        "unforced_damped_H_bounded_below_initial":
            bool(args.torque!=0 or H[-1]<=H[0]+1e-8),
        "actual_forcing_nonzero_if_requested":
            bool(args.torque==0 or np.max(np.abs(input_power))>1e-9),
    }
    mpi = PETSc.COMM_WORLD.getSize()
    rank = PETSc.COMM_WORLD.getRank()
    args.out.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(
        args.out/f"rank{rank}.npz",
        t=t,z=z,theta=theta,p=p,e_theta=etheta,e_p=ep,H=H,
        torque=tau,input_power=input_power,viscous_dissipation=dissipation,
        RK4_reference_theta=reference[:,0],
        RK4_reference_p=reference[:,1],
        CN_ODE_theta_defect=cn_theta_defect,
        CN_ODE_p_defect=cn_p_defect,
        CN_coenergy_rate=cn_coenergy_rate,
        CN_step_dissipation=cn_dissipation,
        CN_step_input_power=cn_input_power,
        CN_coenergy_power_defect=cn_power_defect,
        CN_nonlinear_chain_rule_defect=cn_nonlinear_chain_rule_defect,
        CN_discrete_work=cn_discrete_work,
        CN_cumulative_H_minus_work=cn_cumulative_H_minus_work,
        CN_chain_accounting_defect=cn_chain_accounting_defect,
    )
    native_figures=getattr(system,"pendulum_visuals",None)
    if args.showcase and rank==0 and native_figures is None:
        raise RuntimeError("Original canonical pendulum did not render figures")
    report = {
        "schema":"scrimp-v130-original-nonlinear-pendulum-CN-physics-v1",
        "status":("PASS_ORIGINAL_NONLINEAR_PENDULUM_CN_NATIVE"
                  if all(gates.values()) else "FAIL_OR_UNQUALIFIED_CN_PENDULUM"),
        "example":"examples/nonlinear_pendulum.py",
        "MPI":mpi,"rank":rank,"solver":args.solver,
        "t_final":args.final_time,"dt":args.dt,
        "original_initial_theta":args.angle,
        "original_initial_p":0.,
        "torque_amplitude":args.torque,
        "original_fully_nonlinear_sine":True,
        "physical_parameters":{
            "m":1.,"length":1.,"gravity":9.81,"inertia":1.,
            "damping":D,"torque_frequency_hz":FREQUENCY,
        },
        "time_scheme":"PETSc_TSCN_endpoint_trapezoidal",
        "initialization":"exact_constitutive",
        "pseudo_TS_executed":False,
        "steps":int(system.last_ts_steps),
        "saved_states":n,
        "DAE_DOFs":int(z.shape[1]),
        "H_initial":float(H[0]),
        "H_final":float(H[-1]),
        "power_identity_scope":
            "CN_endpoint_trapezoidal_coenergy_pairing_all_intervals",
        "nonquadratic_H_chain_rule_defect_NOT_zero":True,
        "RK4_tolerance":args.reference_atol,
        "canonical_showcase":bool(args.showcase),
        "canonical_visuals":native_figures,
        "residuals":residuals,
        "gates":{k:bool(v) for k,v in gates.items()},
        "state_SHA256":hashlib.sha256(z.tobytes()).hexdigest(),
    }
    (args.out/f"rank{rank}.json").write_text(
        json.dumps(report,indent=2)+"\n",encoding="utf-8"
    )
    print("ORIGINAL NONLINEAR PENDULUM CN PHYSICAL REPORT "
          +json.dumps(report,sort_keys=True),flush=True)


if __name__=="__main__":
    main()
