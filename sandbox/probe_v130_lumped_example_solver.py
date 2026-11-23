# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Actual original lumped example acceptance on PETSc LU backends.

No physics, weak forms, integrator, geometry or default times are
rewritten. A short *explicitly selected* time window is used for
the preliminary example/solver test, before full shipped defaults.
"""

import argparse
import json
from pathlib import Path

import numpy as np
from petsc4py import PETSc

from examples.mass_spring_damper import mass_spring_damper
from examples.nonlinear_pendulum import nonlinear_pendulum


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--example", choices=("mass_spring_damper",
                                         "nonlinear_pendulum"), required=True)
    cli.add_argument("--factor-solver",
                     choices=("default","mumps","superlu","superlu_dist"),
                     required=True)
    cli.add_argument("--out",type=Path,required=True)
    args=cli.parse_args()
    comm=PETSc.COMM_WORLD.tompi4py()
    if args.factor_solver=="superlu" and comm.size!=1:
        raise ValueError("Sequential SuperLU cannot be used at MPI>1")
    actual_solver = "mumps" if args.factor_solver=="default" else args.factor_solver
    if not PETSc.Sys.hasExternalPackage(actual_solver):
        raise RuntimeError("Missing PETSc backend "+actual_solver)

    if args.example=="mass_spring_damper":
        system=mass_spring_damper(
            final_time=.04,time_step=.01,plot=False,
            factor_solver=None if args.factor_solver=="default" else actual_solver,
        )
        names=("q","p")
    else:
        system=nonlinear_pendulum(
            final_time=.04,time_step=.01,initial_angle=.8,
            torque_amplitude=0.,plot=False,
            factor_solver=None if args.factor_solver=="default" else actual_solver,
        )
        names=("theta","p")

    requested=system.time_scheme.getString(
        "pc_factor_mat_solver_type",""
    ).lower()
    if requested!=actual_solver:
        raise RuntimeError("Actual SCRIMP solver configuration differs: "
                           +requested+" != "+actual_solver)
    t=np.asarray(system.solution["t"],dtype=float)
    states=np.asarray(system.solution["z"],dtype=float)
    H=np.asarray(system.get_Hamiltonian(),dtype=float)
    variables={
        name:np.asarray(system.get_solution(name),dtype=float).ravel()
        for name in names
    }
    if (
        system.last_ts_converged_reason<=0
        or system.last_ts_rejections
        or system.last_ts_snes_failures
        or system.last_ts_steps!=4
        or abs(float(system.last_ts_time)-.04)>1e-12
        or t.shape!=(5,) or states.shape[0]!=5
        or H.shape!=(5,)
        or any(v.shape!=(5,) for v in variables.values())
        or not np.allclose(t,np.linspace(0.,.04,5),rtol=0,atol=1e-11)
        or not np.all(np.isfinite(states))
        or not np.all(np.isfinite(H))
        or any(not np.all(np.isfinite(v)) for v in variables.values())
    ):
        raise RuntimeError("Original PETSc short example integration failed")

    if args.example=="mass_spring_damper":
        # H is the explicitly declared original spring + momentum energy.
        q,p=variables["q"],variables["p"]
        check=.5*5.*q*q+p*p/(2.*1.5)
        if (np.max(abs(H-check))>1e-8
                or abs(q[-1]-q[0])<1e-7):
            raise RuntimeError("Original driven oscillator energy/state incorrect")
    else:
        theta,p=variables["theta"],variables["p"]
        check=9.81*(1.-np.cos(theta))+.5*p*p
        if (np.max(abs(H-check))>1e-8
                or np.max(H-H[0])>1e-5
                or abs(theta[-1]-theta[0])<1e-7):
            raise RuntimeError("Original damped nonlinear pendulum energy incorrect")

    summaries=comm.allgather((t,states,H,variables))
    rank_diff=0.
    for other_t,other_z,other_h,other_vars in summaries:
        if other_z.shape!=states.shape:
            raise RuntimeError("Original full state differs in MPI shape")
        rank_diff=max(rank_diff,float(np.max(abs(other_z-states))),
                      float(np.max(abs(other_t-t))),
                      float(np.max(abs(other_h-H))))
        for name in names:
            rank_diff=max(
                rank_diff,
                float(np.max(abs(other_vars[name]-variables[name])))
            )
    if rank_diff>=1e-8:
        raise RuntimeError("Original full FE/DAE example differs by MPI rank")

    if comm.rank==0:
        args.out.mkdir(parents=True,exist_ok=True)
        out_npz=args.out/"original_example_trace.npz"
        out_json=args.out/"original_example_report.json"
        if out_npz.exists() or out_json.exists():
            raise FileExistsError("Do not overwrite prior native example outputs")
        np.savez_compressed(
            out_npz,times=t,states=states,H=H,
            **variables
        )
        report={
            "schema":"scrimp-v130-lumped-solver-pilot-native-v1",
            "status":"PASS_ORIGINAL_LUMPED_EXAMPLE_NATIVE",
            "example":args.example,"factor_solver":requested,
            "explicit_user_factor_solver":args.factor_solver != "default",
            "mpi_ranks":comm.size,"full_state_dofs":states.shape[1],
            "accepted_steps":system.last_ts_steps,
            "t_final":float(t[-1]),"dt":.01,
            "ts_converged_reason":system.last_ts_converged_reason,
            "ts_rejections":system.last_ts_rejections,
            "snes_failures":system.last_ts_snes_failures,
            "max_rank_local_full_state_discrepancy":rank_diff,
            "max_energy_formula_error":float(np.max(abs(check-H))),
            "physics_and_original_weak_forms_unchanged":True,
            "this_is_not_a_performance_benchmark":True,
        }
        out_json.write_text(json.dumps(report,indent=2)+"\n",
                            encoding="utf-8")
        print(
            "PASS ORIGINAL EXAMPLE "
            f"{args.example} LU={requested} MPI={comm.size} "
            f"DOFs={states.shape[1]} full_rank_diff={rank_diff:.9g}",
            flush=True,
        )
    comm.Barrier()


if __name__=="__main__":
    main()
