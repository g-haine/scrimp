# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Issue #29: INDEPENDENT linear pH-DAE wave boundary-control cold restart.

Faithful physical bricks/FEMs from examples/wave.py (not Saint-Venant).
Zero q/p/controls at t0 give an exactly compatible stationary t0 state;
the real time-dependent boundary velocity input subsequently excites the
linear DAE. No sandbox init_step monkeypatch. The *real* production
DPHS.save_checkpoint/load_checkpoint APIs are exercised on new models
and new PETSc TS instances in 1 or 2 MPI ranks.

The experiment is not a proof of BDF/IMEX history restoration.
"""

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from mpi4py import MPI
from petsc4py import PETSc
import scrimp as S

from scrimp.utils.dae_checkpoint import model_provenance

DT = 0.002
TSTAR = 0.02
TFINAL = 0.024
MESH = 0.5


def digest(z):
    return hashlib.sha256(
        np.asarray(z, dtype="<f8").tobytes(order="C")
    ).hexdigest()


def build_linear_wave(t_final):
    """Original examples/wave.py linear 2-D wave DAE, smaller rectangle mesh."""
    s = S.DPHS("real")
    s.set_domain(S.Domain("Rectangle", {"L": 2.0, "l": 1.0, "h": MESH}))
    q = S.State("q", "Strain", "vector-field")
    p = S.State("p", "Linear momentum", "scalar-field")
    for item in (q, p):
        s.add_state(item)
    for item in (S.CoState("e_q", "Stress", q),
                 S.CoState("e_p", "Velocity", p)):
        s.add_costate(item)
    controls = [
        S.Control_Port(
            "Boundary control (bottom)", "U_B", "Normal force",
            "Y_B", "Velocity trace", "scalar-field", region=10,
        ),
        S.Control_Port(
            "Boundary control (right)", "U_R", "Normal force",
            "Y_R", "Velocity trace", "scalar-field", region=11,
        ),
        S.Control_Port(
            "Boundary control (top)", "U_T", "Normal force",
            "Y_T", "Velocity trace", "scalar-field", region=12,
        ),
        S.Control_Port(
            "Boundary control (left)", "U_L", "Velocity trace",
            "Y_L", "Normal force", "scalar-field", region=13,
            position="flow",
        ),
    ]
    for control in controls:
        s.add_control_port(control)
    for fem in [
        S.FEM("q", 1, "DG"),
        S.FEM("p", 2, "CG"),
        *[S.FEM(control.get_name(), 1, "DG")
          for control in controls],
    ]:
        s.add_FEM(fem)
    # Preserve original physically heterogeneous constitutive parameter.
    s.add_parameter(S.Parameter(
        "T", "Young's modulus", "tensor-field",
        "[[5+x,x*y],[x*y,2+y]]", "q",
    ))
    s.add_parameter(S.Parameter(
        "rho", "Mass density", "scalar-field", "3-x", "p",
    ))
    for brick in [
        S.Brick("M_q", "q.Test_q", [1], dt=True, position="flow"),
        S.Brick("M_p", "p*Test_p", [1], dt=True, position="flow"),
        S.Brick("M_Y_B", "Y_B*Test_Y_B", [10], position="flow"),
        S.Brick("M_Y_R", "Y_R*Test_Y_R", [11], position="flow"),
        S.Brick("M_Y_T", "Y_T*Test_Y_T", [12], position="flow"),
        S.Brick("M_Y_L", "U_L*Test_Y_L", [13], position="flow"),
        S.Brick("D", "Grad(e_p).Test_q", [1], position="effort"),
        S.Brick("-D^T", "-e_q.Grad(Test_p)", [1], position="effort"),
        S.Brick("B_B", "U_B*Test_p", [10], position="effort"),
        S.Brick("B_R", "U_R*Test_p", [11], position="effort"),
        S.Brick("B_T", "U_T*Test_p", [12], position="effort"),
        S.Brick("B_L", "Y_L*Test_p", [13], position="effort"),
        S.Brick("C_B", "-e_p*Test_Y_B", [10], position="effort"),
        S.Brick("C_R", "-e_p*Test_Y_R", [11], position="effort"),
        S.Brick("C_T", "-e_p*Test_Y_T", [12], position="effort"),
        S.Brick("C_L", "-e_p*Test_Y_L", [13], position="effort"),
        S.Brick("-M_e_q", "-e_q.Test_e_q", [1]),
        S.Brick("CR_q", "q.T.Test_e_q", [1]),
        S.Brick("-M_e_p", "-e_p*Test_e_p", [1]),
        S.Brick("CR_p", "p/rho*Test_e_p", [1]),
    ]:
        s.add_brick(brick)
    # At t0=0, all boundary inputs and differential states vanish.
    # An analytical compatible initial DAE state; no hidden projection.
    for ctrl,expr in zip(
        controls,
        ("0.", "0.", "0.", "0.2*sin(25*t)*sin(pi*y)"),
    ):
        s.set_control(ctrl.get_name(), expr)
    s.set_initial_value("q", "[0.,0.]")
    s.set_initial_value("p", "0.")
    s.set_time_scheme(
        ts_type="cn", ts_adapt_type="none",
        t_0=0., t_f=t_final,
        dt=DT, dt_save=DT*(1.-1e-9),
        init_step=False,
    )
    s.hamiltonian.set_name("Mechanical energy")
    s.hamiltonian.add_term(
        S.Term("Potential energy", "0.5*q.T.q", [1])
    )
    s.hamiltonian.add_term(
        S.Term("Kinetic energy", "0.5*p*p/rho", [1])
    )
    return s


def get_indices(system):
    n = int(system.gf_model.nbdof())
    dynamic = []
    for state in system.states.values():
        start,count = map(
            int, system.gf_model.interval_of_variable(state.get_name())
        )
        dynamic.extend(range(start,start+count))
    dyn = np.asarray(dynamic, dtype=int)
    if not dyn.size or np.unique(dyn).size != dyn.size:
        raise RuntimeError("Linear DAE dynamic FE state intervals invalid")
    alg = np.ones(n,dtype=bool)
    alg[dyn]=False
    if not np.any(alg):
        raise RuntimeError("Second test model is NOT a genuine DAE")
    return dyn, alg


def original_IFunction(system, mpi, t, state, state_dot):
    """Evaluate original real GetFEM IFunction on native global PETSc vectors."""
    if not np.all(np.isfinite(state)) or not np.all(np.isfinite(state_dot)):
        raise RuntimeError("Non-finite original linear DAE coefficients")
    vec=system.tangent_stiffness.createVecRight()
    vel=vec.duplicate()
    residual=vec.duplicate()
    try:
        lo,hi=vec.getOwnershipRange()
        if hi>lo:
            # PETSc can be compiled with 32-bit PetscInt even on
            # 64-bit Python/NumPy, where dtype=int is int64. Index
            # arrays passed to Vec.setValues MUST use PETSc.IntType.
            owned=np.arange(lo,hi,dtype=PETSc.IntType)
            vec.setValues(owned,state[lo:hi],addv=PETSc.InsertMode.INSERT_VALUES)
            vel.setValues(owned,state_dot[lo:hi],addv=PETSc.InsertMode.INSERT_VALUES)
        vec.assemble()
        vel.assemble()
        residual.set(0.)
        system.IFunction(None,float(t),vec,vel,residual)
        parts=mpi.allgather(
            (lo,hi,np.asarray(residual.getArray(readonly=True),dtype=float).copy())
        )
        output=np.zeros(int(system.gf_model.nbdof()))
        seen=np.zeros(output.size,dtype=int)
        for x,y,values in parts:
            output[x:y]=values
            seen[x:y]+=1
        if not np.all(seen==1):
            raise RuntimeError("Original linear DAE residual MPI ownership invalid")
        return output
    finally:
        residual.destroy()
        vel.destroy()
        vec.destroy()


def verify_trajectory(system, mpi, first, last, expected_steps):
    if (int(system.last_ts_converged_reason)<=0
            or int(system.last_ts_rejections)
            or int(system.last_ts_snes_failures)
            or int(system.last_ts_steps)!=expected_steps
            or abs(float(system.last_ts_time)-last)>1e-12):
        raise RuntimeError(
            "Original linear-wave PETSc TS stopped early/rejected steps"
        )
    times=np.asarray(system.solution["t"],dtype=float)
    states=np.asarray(system.solution["z"],dtype=float)
    n=int(system.gf_model.nbdof())
    expected=np.linspace(first,last,expected_steps+1)
    if (times.shape!=expected.shape
            or states.shape!=(expected.size,n)
            or not np.all(np.isfinite(states))
            or not np.allclose(times,expected,rtol=0.,atol=1e-12)):
        raise RuntimeError("Original linear-wave TS full trajectory unavailable")
    indices,alg = get_indices(system)
    zero=np.zeros(n)
    # Real original algebraic constraints at each accepted physical time.
    fa=[]
    cn=[]
    for k,t in enumerate(times):
        values=original_IFunction(system,mpi,t,states[k],zero)
        fa.append(float(np.linalg.norm(values[alg])))
        if k>0:
            deriv=2*(states[k]-states[k-1])/(times[k]-times[k-1])
            current=original_IFunction(system,mpi,t,states[k],deriv)
            previous=original_IFunction(system,mpi,times[k-1],states[k-1],zero)
            cn.append(float(np.linalg.norm(current+previous)))
    fa_max=max(fa)
    cn_max=max(cn)
    if fa_max>1e-8 or cn_max>1e-7:
        raise RuntimeError(
            f"Original real linear pH-DAE equations fail: Fa={fa_max:.9g}, "
            f"TSCN={cn_max:.9g}"
        )
    if first==0. and np.linalg.norm(states[-1,indices])<1e-8:
        raise RuntimeError(
            "Time-dependent control did not excite any differential wave DOFs"
        )
    return times,states,indices,{
        "maximum_original_TSCN_residual":cn_max,
        "maximum_original_Fa_residual":fa_max,
        "endpoint_Fa_residuals":fa,
        "full_step_TSCN_residuals":cn,
    }


def original_energy_and_boundary_power(system, times):
    """Use existing SCRIMP Hamiltonian + actual GetFEM control-port powers."""
    energies=np.asarray(system.get_Hamiltonian(),dtype=float)
    system.compute_powers()
    powers={}
    for name,port in system.ports.items():
        if port.get_algebraic() and not port.get_substituted():
            values=np.asarray(port.get_power(),dtype=float)
            if values.shape!=times.shape:
                raise RuntimeError("Original wave boundary-port trace invalid")
            powers[port.get_name()]=values
    if not powers:
        raise RuntimeError("Original linear DAE lacks real boundary powers")
    net=np.sum(np.stack(list(powers.values())),axis=0)
    # SCRIMP Port.compute evaluates net[k] (k>=1) from the
    # Crank–Nicolson midpoint z[k-1]/z[k], with t[k] passed to
    # the GetFEM algebraic port. Thus each net[k] is a STEP
    # power, NOT an instantaneous sample suitable for the
    # endpoint trapezoidal rule. net[0] is instantaneous
    # at the first accepted state and must NOT be used to
    # integrate the next CN step (especially after restart).
    # Original left-boundary port sign convention: supplied work
    # W=-sum(dt_k*net[k]) and ΔH=W for this conservative wave.
    work=np.r_[0.,np.cumsum(
        -np.diff(times)*net[1:]
    )]
    per_step_balance=np.diff(energies)+np.diff(times)*net[1:]
    max_balance=float(np.max(np.abs(per_step_balance)))
    if (not np.isfinite(max_balance) or max_balance>=1e-11):
        raise RuntimeError(
            "Original linear pH-DAE Crank–Nicolson power identity "
            f"ΔH+dt*P != 0; maximal per-step defect {max_balance:.9g}"
        )
    if (energies.shape!=times.shape or not np.all(np.isfinite(energies))
            or not np.all(np.isfinite(work))):
        raise RuntimeError("Original GetFEM Hamiltonian/control work invalid")
    return energies,work,net,max_balance


def run():
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--phase",choices=("seed","continuous","resume"),required=True)
    cli.add_argument("--checkpoint",type=Path,required=True)
    cli.add_argument("--out",type=Path,required=True)
    args=cli.parse_args()
    mpi=PETSc.COMM_WORLD.tompi4py()
    first=TSTAR if args.phase=="resume" else 0.
    last=TSTAR if args.phase=="seed" else TFINAL
    steps=int(round((last-first)/DT))
    system=build_linear_wave(last)
    if args.phase=="resume":
        system.load_checkpoint(args.checkpoint)
        if (abs(float(system.time_scheme["t_0"])-first)>1e-12
                or str(system.time_scheme["init_step"]).lower()!="false"):
            raise RuntimeError("Real production wave checkpoint failed to set t0")
    system.solve()
    t,z,dyn,diagnostics=verify_trajectory(
        system,mpi,first,last,steps
    )
    fingerprint=model_provenance(system)
    if args.phase=="resume":
        with args.checkpoint.open("rb") as handle:
            with np.load(handle,allow_pickle=False) as cp:
                original=np.asarray(cp["state"],dtype=float)
        if not np.array_equal(z[0],original):
            raise RuntimeError("Production DPHS linear DAE cold restore changed z*")
        status=system.last_checkpoint_report
        if (not isinstance(status,dict)
                or status.get("status")!="COMPLETE_DAE_STATE_COLD_RESTORED_IN_DPHS"
                or status.get("no_pseudo_BDF_or_projection") is not True
                or status.get(
                    "all_GetFEM_and_GlobalVec_coefficients_restored_bitwise"
                ) is not True
                or status.get("original_algebraic_residual_l2",1.)>1e-9):
            raise RuntimeError("Production DPHS linear DAE restore was projected")
    else:
        status={"status":"ORIGINAL_ANALYTICALLY_COMPATIBLE_ZERO_INITIAL_DAE",
                "init_step":False,"projected":False}
    if args.phase=="seed":
        system.save_checkpoint(args.checkpoint,extra={
            "source":"original_examples_wave_linear_DAE",
            "physical_time_origin":"absolute",
            "control_law":"0.2*sin(25*t)*sin(pi*y)",
            "source_MPI_ranks":mpi.size,
            "original_Fa_max":diagnostics["maximum_original_Fa_residual"],
        })
    h,w,power,discrete_power_defect=original_energy_and_boundary_power(system,t)
    report={
        "schema":"scrimp-v130-issue29-linear-DAE-public-API-v2",
        "status":"NATIVE_LINEAR_DAE_SEGMENT_NOT_YET_CROSS_MPI_COMPARED",
        "phase":args.phase,
        "MPI":mpi.size,"mesh_h":MESH,
        "t0":first,"tf":last,"dt":DT,
        "step_count":int(system.last_ts_steps),
        "PETSc_reason":int(system.last_ts_converged_reason),
        "PETSc_rejected":int(system.last_ts_rejections),
        "SNES_failures":int(system.last_ts_snes_failures),
        "ndof":int(system.gf_model.nbdof()),
        "dynamic_ndof":int(dyn.size),
        "algebraic_ndof":int(z.shape[1]-dyn.size),
        "fingerprint_sha256":fingerprint["fingerprint_sha256"],
        "initial_full_sha256":digest(z[0]),
        "final_full_sha256":digest(z[-1]),
        "original_equations":diagnostics,
        "Hamiltonian_initial":float(h[0]),
        "Hamiltonian_final":float(h[-1]),
        "boundary_work":float(w[-1]),
        "boundary_work_quadrature":"minus_CN_step_midpoint_port_power_times_dt",
        "maximum_discrete_energy_balance_defect":discrete_power_defect,
        "dH_minus_W":float(h[-1]-h[0]-w[-1]),
        "initialization_report":status,
        "production_API_used":args.phase in ("seed","resume"),
        "no_integrator_history_claimed":True,
        "not_SWE_model":True,
    }
    if mpi.rank==0:
        args.out.mkdir(parents=True,exist_ok=True)
        (args.out/"linear_report.json").write_text(
            json.dumps(report,indent=2)+"\n",encoding="utf-8"
        )
        np.savez_compressed(
            args.out/"linear_states.npz",
            times=t,states=z,H=h,work=w,port_power=power,
        )
        print(
            f"LINEAR DAE PHASE={args.phase} MPI={mpi.size} "
            f"steps={steps} total_FE={z.shape[1]} "
            f"Fa_max={diagnostics['maximum_original_Fa_residual']:.8g} "
            f"TSCN_max={diagnostics['maximum_original_TSCN_residual']:.8g} "
            f"dH-W={report['dH_minus_W']:.8g}",
            flush=True,
        )
    mpi.Barrier()


if __name__=="__main__":
    run()
