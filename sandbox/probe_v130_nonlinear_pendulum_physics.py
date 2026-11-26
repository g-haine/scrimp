"""Native qualification of ORIGINAL nonlinear 0D pendulum pH-DAE.

No surrogate SCRIMP model: the unmodified original mass, Dirac,
damping and FULL sine constitutive GetFEM bricks are used directly.
The exact algebraic initial data preserve the physical IVP.
For a NONquadratic H, the BDF2 discrete co-energy power identity is
 e_theta[n] D2 theta[n] + e_p[n] D2 p[n]
     + d*e_p[n]**2 - tau(t_n)*e_p[n] == 0
(k>=2). This is NOT the telescoping quadratic H identity.
The actual GetFEM H and an independent RK4 ODE are separately checked.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from petsc4py import PETSc

from examples.nonlinear_pendulum import nonlinear_pendulum


G=9.81
I=1.0
D=.08
FREQUENCY=.3


def ode(t,x,amplitude):
    theta,p=x
    return np.array([
        p/I,
        -G*np.sin(theta)-D*(p/I)
        +amplitude*np.sin(2.*np.pi*FREQUENCY*t)
    ],dtype=float)


def rk4_reference(t,initial,amplitude):
    out=np.empty((len(t),2),dtype=float)
    x=np.asarray(initial,dtype=float).copy()
    out[0]=x
    for n in range(1,len(t)):
        dt=float(t[n]-t[n-1])/64.
        time=float(t[n-1])
        for _ in range(64):
            k1=ode(time,x,amplitude)
            k2=ode(time+dt/2.,x+dt*k1/2.,amplitude)
            k3=ode(time+dt/2.,x+dt*k2/2.,amplitude)
            k4=ode(time+dt,x+dt*k3,amplitude)
            x=x+dt*(k1+2*k2+2*k3+k4)/6.
            time+=dt
        out[n]=x
    return out


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--final-time",type=float,default=.08)
    parser.add_argument("--dt",type=float,default=.01)
    parser.add_argument("--angle",type=float,default=1.2)
    parser.add_argument("--torque",type=float,default=0.)
    parser.add_argument("--solver",default="mumps")
    parser.add_argument("--out",type=Path,required=True)
    args=parser.parse_args()
    system=nonlinear_pendulum(
        final_time=args.final_time,time_step=args.dt,
        initial_angle=args.angle,torque_amplitude=args.torque,
        plot=False,factor_solver=args.solver,
    )

    t=np.asarray(system.solution["t"],dtype=float)
    z=np.asarray(system.solution["z"],dtype=float)
    theta=np.asarray(system.get_solution("theta"),dtype=float).ravel()
    p=np.asarray(system.get_solution("p"),dtype=float).ravel()
    etheta=np.asarray(system.get_solution("e_theta"),dtype=float).ravel()
    ep=np.asarray(system.get_solution("e_p"),dtype=float).ravel()
    H=np.asarray(system.get_Hamiltonian(),dtype=float).ravel()
    n=len(t)
    steps=round(args.final_time/args.dt)
    for key,series in (
        ("theta",theta),("p",p),("e_theta",etheta),("e_p",ep),("H",H),
    ):
        if series.shape!=(n,) or not np.all(np.isfinite(series)):
            raise RuntimeError(f"Invalid original pendulum GetFEM {key} history")
    if (n!=steps+1 or z.shape!=(n,4)
            or not np.all(np.isfinite(z))
            or np.max(np.abs(t-np.arange(n)*args.dt))>1e-10):
        raise RuntimeError("Original PETSc BDF2 time grid/DAE state mismatch")

    tau=args.torque*np.sin(2.*np.pi*FREQUENCY*t)
    analytic_H=.5*p*p/I+G*(1.-np.cos(theta))
    Dphysical=D*(p/I)**2
    Pinput=tau*(p/I)
    D2theta=(3.*theta[2:]-4.*theta[1:-1]+theta[:-2])/(2.*args.dt)
    D2p=(3.*p[2:]-4.*p[1:-1]+p[:-2])/(2.*args.dt)
    ODEtheta=D2theta-ep[2:]
    ODEp=D2p+etheta[2:]+D*ep[2:]-tau[2:]
    # Exact identity for the implicit BDF2 *residual* paired with
    # the true nonlinear co-energies, NOT a false H(z)-difference rule.
    nonlinear_step_power=(
        etheta[2:]*D2theta+ep[2:]*D2p
        +D*ep[2:]**2-tau[2:]*ep[2:]
    )
    # The continuous analytic Hamiltonian is independently tested.
    # BDF2 generally has numerical dissipation; compare H to RK4,
    # not to an invented exact discrete-gradient time rule.
    reference=rk4_reference(t,[args.angle,0.],args.torque)
    max_rk4=float(np.max(np.abs(
        np.column_stack((theta,p))-reference
    )))
    residuals={
        "constitutive_nonlinear_full_sine":
            float(np.max(np.abs(etheta-G*np.sin(theta)))),
        "constitutive_velocity":
            float(np.max(np.abs(ep-p/I))),
        "native_H_vs_analytic":
            float(np.max(np.abs(H-analytic_H))),
        "BDF2_ODE_theta_k_ge_2":float(np.max(np.abs(ODEtheta))),
        "BDF2_ODE_p_k_ge_2":float(np.max(np.abs(ODEp))),
        "nonquadratic_BDF2_coenergy_power_k_ge_2":
            float(np.max(np.abs(nonlinear_step_power))),
        "RK4_original_sine_ODE_max_theta_or_p_error":max_rk4,
        "theta_initial_drift":float(theta[0]-args.angle),
        "p_initial_drift":float(p[0]),
        "e_theta_initial_drift":
            float(etheta[0]-G*np.sin(args.angle)),
        "e_p_initial_drift":float(ep[0]),
    }
    if system.time_scheme.getString("init_step","true")!="false":
        raise RuntimeError("Canonical nonlinear pendulum unexpectedly ran initial TS")

    gates={
        "native_BDF2_completed":
            system.last_ts_converged_reason>0
            and system.last_ts_steps==steps
            and system.last_ts_rejections==0
            and system.last_ts_snes_failures==0,
        "exact_prescribed_initial_differential_states":
            max(abs(residuals["theta_initial_drift"]),
                abs(residuals["p_initial_drift"]))<1e-12,
        "exact_initial_nonlinear_constitutive_data":
            max(abs(residuals["e_theta_initial_drift"]),
                abs(residuals["e_p_initial_drift"]))<1e-12,
        "nonlinear_sine_coenergy_relation":
            residuals["constitutive_nonlinear_full_sine"]<1e-8,
        "lumped_velocity_coenergy_relation":
            residuals["constitutive_velocity"]<1e-8,
        "original_GetFEM_nonquadratic_H":
            residuals["native_H_vs_analytic"]<1e-8,
        "BDF2_theta_differential_equation":
            residuals["BDF2_ODE_theta_k_ge_2"]<1e-7,
        "BDF2_p_differential_equation":
            residuals["BDF2_ODE_p_k_ge_2"]<1e-7,
        "exact_nonquadratic_BDF2_coenergy_power_identity":
            residuals["nonquadratic_BDF2_coenergy_power_k_ge_2"]<1e-7,
        "independent_full_sine_RK4_reference":
            max_rk4<3e-3,
        "nonnegative_H_and_dissipation":
            float(np.min(H))>=-1e-10
            and float(np.min(Dphysical))>=-1e-14,
        "unforced_damped_H_nonincreasing":
            bool(args.torque!=0 or np.max(np.diff(H))<1e-8),
        "actual_force_power_excited_if_requested":
            bool(args.torque==0 or np.max(np.abs(Pinput))>1e-9),
    }
    rank=PETSc.COMM_WORLD.getRank()
    mpi=PETSc.COMM_WORLD.getSize()
    args.out.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(
        args.out/f"rank{rank}.npz",
        t=t,z=z,theta=theta,p=p,e_theta=etheta,e_p=ep,H=H,
        torque=tau,input_power=Pinput,viscous_dissipation=Dphysical,
        RK4_reference_theta=reference[:,0],
        RK4_reference_p=reference[:,1],
        BDF2_ODE_theta_defect=ODEtheta,
        BDF2_ODE_p_defect=ODEp,
        BDF2_nonquadratic_coenergy_power_defect=nonlinear_step_power,
    )
    report={
        "schema":"scrimp-v130-original-nonlinear-pendulum-pH-BDF2-v1",
        "status":("PASS_ORIGINAL_NONLINEAR_PENDULUM_NATIVE_PILOT"
                  if all(gates.values()) else "FAIL_OR_UNQUALIFIED_PENDULUM"),
        "example":"examples/nonlinear_pendulum.py",
        "MPI":mpi,"rank":rank,"solver":args.solver,
        "t_final":args.final_time,"dt":args.dt,
        "original_initial_theta":args.angle,
        "original_initial_p":0.,
        "torque_amplitude":args.torque,
        "original_fully_nonlinear_sine":True,
        "physical_parameters":{"m":1.,"length":1.,"gravity":9.81,
                               "inertia":1.,"damping":.08,
                               "torque_frequency_hz":.3},
        "time_scheme":"original_fixed_BDF2",
        "initialization":"exact_constitutive",
        "pseudo_TS_executed":False,
        "steps":int(system.last_ts_steps),
        "saved_states":n,
        "DAE_DOFs":int(z.shape[1]),
        "H_initial":float(H[0]),
        "H_final":float(H[-1]),
        "power_identity_scope":
            "nonquadratic_coenergy_at_BDF2_D2_z_k_ge_2_not_H_telescoping",
        "PETSc_BDF2_startup_internal_stage_not_claimed":True,
        "residuals":residuals,
        "gates":{k:bool(v) for k,v in gates.items()},
        "state_SHA256":hashlib.sha256(z.tobytes()).hexdigest(),
    }
    (args.out/f"rank{rank}.json").write_text(
        json.dumps(report,indent=2)+"\n",encoding="utf-8"
    )
    print("ORIGINAL NONLINEAR PENDULUM PHYSICAL REPORT "
          +json.dumps(report,sort_keys=True),flush=True)


if __name__=="__main__":
    main()
