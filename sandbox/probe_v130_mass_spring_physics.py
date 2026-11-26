"""Native diagnostic of the ORIGINAL mass_spring_damper.py pH-DAE.

Keep original weak forms, lumped ports, BDF2 time integration, forcing
and parameter values; the exact co-energy start is the accepted default. Check the differential and constitutive relations,
the exact quadratic BDF2 step-energy identity (k>=2), and an
independent high-accuracy RK4 approximation of the true forced ODE.
This is a short qualification pilot, not a modified example.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from petsc4py import PETSc

from examples.mass_spring_damper import mass_spring_damper


MASS=1.5
STIFFNESS=5.0
DAMPING=0.25
INITIAL_Q=0.4
INITIAL_P=0.0


def energy(q, p):
    return 0.5*STIFFNESS*q*q+p*p/(2.*MASS)


def ode(t, x):
    q,p=x
    return np.array([p/MASS, -STIFFNESS*q-(DAMPING/MASS)*p
                     +0.2*np.sin(2.*np.pi*t)], dtype=float)


def rk4_reference(times, initial):
    """Independent RK4 at 32 substeps per actual saved CN/BDF interval."""
    out=np.empty((len(times),2),dtype=float)
    out[0]=initial
    x=np.asarray(initial,dtype=float).copy()
    for n in range(1,len(times)):
        t=float(times[n-1])
        dt=float(times[n]-times[n-1])/32.
        for i in range(32):
            k1=ode(t,x)
            k2=ode(t+dt/2.,x+dt*k1/2.)
            k3=ode(t+dt/2.,x+dt*k2/2.)
            k4=ode(t+dt,x+dt*k3)
            x=x+(dt/6.)*(k1+2.*k2+2.*k3+k4)
            t+=dt
        out[n]=x
    return out


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--final-time",type=float,required=True)
    parser.add_argument("--time-step",type=float,required=True)
    parser.add_argument("--solver",default="mumps")
    parser.add_argument("--showcase",action="store_true",
                        help="Render by canonical example (not by probe)")
    parser.add_argument("--out",type=Path,required=True)
    args=parser.parse_args()

    # Canonical physical example, no parameter/weak-form edits.
    system=mass_spring_damper(
        final_time=args.final_time,
        time_step=args.time_step,
        plot=args.showcase,factor_solver=args.solver,
    )
    t=np.asarray(system.solution["t"],dtype=float)
    z=np.asarray(system.solution["z"],dtype=float)
    q=np.asarray(system.get_solution("q"),dtype=float).ravel()
    p=np.asarray(system.get_solution("p"),dtype=float).ravel()
    eq=np.asarray(system.get_solution("e_q"),dtype=float).ravel()
    ep=np.asarray(system.get_solution("e_p"),dtype=float).ravel()
    H=np.asarray(system.get_Hamiltonian(),dtype=float).ravel()
    n=len(t)
    expected_steps=round(args.final_time/args.time_step)
    if any(arr.shape!=(n,) for arr in (q,p,eq,ep,H)):
        raise RuntimeError(
            f"Original lumped port or Hamiltonian shape failed: "
            f"times={t.shape}, q={q.shape}, p={p.shape}, "
            f"eq={eq.shape}, ep={ep.shape}, H={H.shape}"
        )
    if (z.shape!=(n,4) or n!=expected_steps+1
            or not np.all(np.isfinite(z))
            or not np.all(np.isfinite(H))
            or not np.all(np.diff(t)>0.)):
        raise RuntimeError("Original GetFEM lumped-state/TS trajectory malformed")

    f=0.2*np.sin(2.*np.pi*t)
    vel=p/MASS
    D=DAMPING*vel*vel
    P=f*vel
    physical_H=energy(q,p)

    # Differentiate only original accepted BDF2 states k>=2;
    # PETSc's first BDF2 startup interval is not a BDF1 step.
    h=np.diff(t)
    if not np.allclose(h,args.time_step,rtol=0,atol=1e-10):
        raise RuntimeError("Original BDF2 fixed time step unexpectedly changed")
    derivative_q=(3*q[2:]-4*q[1:-1]+q[:-2])/(2.*args.time_step)
    derivative_p=(3*p[2:]-4*p[1:-1]+p[:-2])/(2.*args.time_step)
    ode_q_defect=derivative_q-ep[2:]
    ode_p_defect=derivative_p+eq[2:]+DAMPING*ep[2:]-f[2:]

    # H(a)-H(b)+H(2a-b)-H(2b-c)+H(a-2b+c) over 2 dt
    # is precisely <grad H(a), (3a-4b+c)/(2 dt)> for
    # quadratic H. No trapezoidal work/quadrature is substituted.
    bdf2_power_rate=(
        physical_H[2:]-physical_H[1:-1]
        +energy(2*q[2:]-q[1:-1],2*p[2:]-p[1:-1])
        -energy(2*q[1:-1]-q[:-2],2*p[1:-1]-p[:-2])
        +energy(q[2:]-2*q[1:-1]+q[:-2],
                p[2:]-2*p[1:-1]+p[:-2])
    )/(2.*args.time_step)
    bdf2_power_defect=bdf2_power_rate + D[2:]-P[2:]
    reference=rk4_reference(t, [INITIAL_Q,INITIAL_P])
    initial=[float(q[0]-INITIAL_Q),float(p[0]-INITIAL_P)]
    max_ref=float(np.max(np.abs(
        np.column_stack((q,p))-reference
    )))
    residuals={
        "coenergy_spring":float(np.max(np.abs(eq-STIFFNESS*q))),
        "coenergy_momentum":float(np.max(np.abs(ep-vel))),
        "Hamiltonian_vs_analytic":float(np.max(np.abs(H-physical_H))),
        "BDF2_ODE_q_k_ge_2":float(np.max(np.abs(ode_q_defect))),
        "BDF2_ODE_p_k_ge_2":float(np.max(np.abs(ode_p_defect))),
        "exact_quadratic_BDF2_step_energy_k_ge_2":
            float(np.max(np.abs(bdf2_power_defect))),
        "RK4_reference_max_q_or_p_error":max_ref,
        "initial_prescribed_q_difference":initial[0],
        "initial_prescribed_p_difference":initial[1],
    }
    # No pseudo-time relaxation in the canonical example. BDF2 startup
    # itself must not be confused with a full BDF1 energy-power step.
    if system.time_scheme.getString("init_step","true") != "false":
        raise RuntimeError("Canonical example unexpectedly ran an initial TS")
    gates={
        "original_TS_converged":
            system.last_ts_converged_reason>0
            and system.last_ts_steps==expected_steps
            and system.last_ts_rejections==0
            and system.last_ts_snes_failures==0,
        "original_numerical_initial_values_preserved":
            max(abs(initial[0]),abs(initial[1]))<1e-8,
        "spring_coenergy_relation":residuals["coenergy_spring"]<1e-8,
        "momentum_coenergy_relation":residuals["coenergy_momentum"]<1e-8,
        "original_GetFEM_Hamiltonian":residuals["Hamiltonian_vs_analytic"]<1e-8,
        "BDF2_differential_q":residuals["BDF2_ODE_q_k_ge_2"]<1e-7,
        "BDF2_differential_p":residuals["BDF2_ODE_p_k_ge_2"]<1e-7,
        "BDF2_exact_step_power":
            residuals["exact_quadratic_BDF2_step_energy_k_ge_2"]<1e-7,
        "independent_forced_ODE_reference":
            max_ref<2e-3,
        "positive_energy_and_dissipation":
            np.all(H>=0.) and np.all(D>=0.),
    }
    rank=PETSc.COMM_WORLD.getRank()
    mpi=PETSc.COMM_WORLD.getSize()
    args.out.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(
        args.out/f"rank{rank}.npz",
        t=t,z=z,q=q,p=p,e_q=eq,e_p=ep,H=H,reference_q=reference[:,0],
        reference_p=reference[:,1],Fourier_or_viscous_dissipation=D,
        external_force=f,external_power=P,
        BDF2_ODE_q_defect=ode_q_defect,
        BDF2_ODE_p_defect=ode_p_defect,
        BDF2_step_power_defect=bdf2_power_defect,
    )
    figures=getattr(system,"mass_spring_visuals",None)
    if args.showcase and PETSc.COMM_WORLD.getRank()==0 and not figures:
        raise RuntimeError("Canonical mass-spring did not publish visuals")
    report={
        "schema":"scrimp-v130-original-mass-spring-bdf2-pH-pilot-v1",
        "status":("PASS_ORIGINAL_LUMPED_ODE_NATIVE_PILOT"
                  if all(gates.values()) else "FAIL_OR_UNQUALIFIED_NATIVE_PILOT"),
        "example":"examples/mass_spring_damper.py",
        "MPI":mpi,"rank":rank,"solver":args.solver,
        "t_final":args.final_time,"dt":args.time_step,
        "steps":system.last_ts_steps,"saved_states":n,
        "full_DAE_DOFs":int(z.shape[1]),
        "physical_parameters":{
            "m":MASS,"k":STIFFNESS,"c":DAMPING,
            "force":"0.2*sin(2*pi*t)",
            "original_q0":INITIAL_Q,"original_p0":INITIAL_P,
        },
        "time_stepper":"original_BDF2",
        "initializer":"exact_constitutive",
        "no_pseudo_time_initialization":True,
        "canonical_showcase":bool(args.showcase),
        "canonical_visuals":figures,
        "native_q0":float(q[0]),
        "native_p0":float(p[0]),
        "native_e_q0":float(eq[0]),
        "native_e_p0":float(ep[0]),
        "power_identity_scope":"BDF2_k_ge_2_no_unverified_startup",
        "residuals":residuals,
        "gates":{k:bool(v) for k,v in gates.items()},
        "H_initial":float(H[0]),"H_final":float(H[-1]),
        "energy_change":float(H[-1]-H[0]),
        "force_power_not_silently_omitted":bool(np.max(np.abs(P))>1e-10),
        "full_state_SHA256":hashlib.sha256(z.tobytes()).hexdigest(),
    }
    (args.out/f"rank{rank}.json").write_text(
        json.dumps(report,indent=2)+"\n",encoding="utf-8"
    )
    print("ORIGINAL LUMPED MASS-SPRING pH/DAE NATIVE REPORT "
          +json.dumps(report,sort_keys=True),flush=True)


if __name__=="__main__":
    main()
