# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Actual short PETSc TS: original inviscid SWE weak forms, explicit G.

Unlike the already-PASSED full IFunction/IJacobian directional-FD gate,
this test executes REAL physical PETSc TS steps. It overrides only
mesh density and time horizon (and, for an independent nonvacuous
case, prescribed p0), never the original PDE/G weak forms.
A separate manufactured initial p=[y,-x] makes lagged G measurable.
The original published p=[0,0] is tested independently as baseline.
This records signed Hamiltonian drift without falsely asserting
second-order accuracy or exact structure-preserving split power.
"""

import argparse
import importlib.util
import json
from pathlib import Path
from unittest.mock import patch

import numpy as np
from mpi4py import MPI
import scrimp as S


def run_case(mesh, dt, steps, manufactured, mpi):
    path=Path(__file__).resolve().parents[1]/"examples/shallow_water.py"
    spec=importlib.util.spec_from_file_location(
        "scrimp_v130_shallow_ts_original_weak_forms", path
    )
    mod=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)

    original_domain=S.Domain
    original_set_time=S.DPHS.set_time_scheme
    original_set_initial=S.DPHS.set_initial_value
    original_prestep=S.DPHS._prestep_explicit_rhs
    state={"main_prestep":[],"max_main_cache_drift":0.,
           "main_residual_checks":0}
    live={"previous_cache":None}

    def domain(name,options):
        if name=="Rectangle":
            options={**options,"h":mesh}
        return original_domain(name,options)

    def initial(self,name,value):
        if manufactured and name=="p":
            # Smooth DG-compatible rotating momentum field; this
            # intentionally changes the published initial DATA
            # only in the separate manufactured test, not its PDE.
            value="[y,-x]"
        return original_set_initial(self,name,value)

    def scheme(self,**kwargs):
        # Retain original BDF2 and explicit pseudo relaxation;
        # substitute only bounded native runtime options.
        kwargs.update(
            init_step_ts_type="pseudo",
            ts_type="bdf",ts_bdf_order=2,
            pc_type="lu",pc_factor_mat_solver_type="mumps",
            ksp_type="gmres",
            t_f=dt*steps,dt=dt,dt_save=dt,
            ts_adapt_type="none",
        )
        return original_set_time(self,**kwargs)

    def prestep(self,ts):
        original_prestep(self,ts)
        if self._profile_phase!="main":
            return
        if not self._explicit_rhs_ready:
            raise RuntimeError("Physical TS has no frozen explicit RHS")
        if live["previous_cache"] is not None:
            live["previous_cache"].destroy()
        live["previous_cache"]=self._explicit_rhs.copy()
        state["main_prestep"].append({
            "time":float(ts.getTime()),
            "accepted_step_count":int(ts.getStepNumber()),
            "accepted_state_norm":float(ts.getSolution().norm()),
            "frozen_gyroscopic_rhs_norm":float(self._explicit_rhs.norm()),
            "frozen_rhs_time":float(self._explicit_rhs_time),
        })

    original_ifunction=S.DPHS.IFunction

    def ifunction(self,ts,t,z,zd,F):
        original_ifunction(self,ts,t,z,zd,F)
        if self._profile_phase!="main" or live["previous_cache"] is None:
            return
        delta=live["previous_cache"].copy()
        try:
            delta.axpy(-1,self._explicit_rhs)
            drift=float(delta.norm())
            state["max_main_cache_drift"]=max(
                state["max_main_cache_drift"],drift
            )
            state["main_residual_checks"]+=1
        finally:
            delta.destroy()

    def no_plot(self,**kwargs):
        return None

    def no_export(self,*args,**kwargs):
        return None

    with (
        patch.object(S,"Domain",side_effect=domain),
        patch.object(S.DPHS,"set_initial_value",initial),
        patch.object(S.DPHS,"set_time_scheme",scheme),
        patch.object(S.DPHS,"_prestep_explicit_rhs",prestep),
        patch.object(S.DPHS,"IFunction",ifunction),
        patch.object(S.DPHS,"plot_Hamiltonian",no_plot),
        patch.object(S.DPHS,"export_to_pv",no_export),
    ):
        try:
            system=mod.shallow_water_eq()
        finally:
            if live["previous_cache"] is not None:
                live["previous_cache"].destroy()
                live["previous_cache"]=None

    if not system.solve_done:
        raise RuntimeError("Production physical shallow-water solve failed")
    if (
        system.last_ts_converged_reason<=0
        or system.last_ts_steps!=steps
        or system.last_ts_rejections!=0
        or system.last_ts_snes_failures!=0
        or abs(system.last_ts_time-steps*dt)>1e-10
    ):
        raise RuntimeError(
            "Unsuccessful or unexpected production TS integration: "
            f"reason={system.last_ts_converged_reason}, "
            f"steps={system.last_ts_steps}, "
            f"rejections={system.last_ts_rejections}, "
            f"snes_failures={system.last_ts_snes_failures}, "
            f"t={system.last_ts_time}"
        )
    events=state["main_prestep"]
    if len(events)!=steps:
        raise RuntimeError(
            f"PreStep was invoked {len(events)} times for "
            f"{steps} accepted physical steps"
        )
    event_times=np.array([x["time"] for x in events])
    if (not np.all(np.diff(event_times)>0)
        or abs(event_times[0])>1e-12):
        raise RuntimeError("Explicit forcing did not follow accepted TS times")
    if state["max_main_cache_drift"]!=0:
        raise RuntimeError("Explicit RHS changed during physical Newton")
    if state["main_residual_checks"]<steps:
        raise RuntimeError("Physical Newton did not evaluate residuals")

    times=np.asarray(system.solution["t"],dtype=float)
    trajectory=np.asarray(system.solution["z"],dtype=float)
    if (len(times)!=len(trajectory)
        or len(times)<2 or times[-1]<steps*dt-1e-10
        or not np.all(np.diff(times)>0)
        or not np.all(np.isfinite(trajectory))):
        raise RuntimeError("Incomplete/nonfinite physical trajectory")
    H=np.asarray(system.get_Hamiltonian(),dtype=float)
    if (len(H)!=len(times) or not np.all(np.isfinite(H))
        or np.min(H)<=0):
        raise RuntimeError("Invalid original mechanical Hamiltonian")
    nonzero_gyroscopic=max(
        x["frozen_gyroscopic_rhs_norm"] for x in events
    )
    if manufactured and nonzero_gyroscopic<=1e-12:
        raise RuntimeError("Manufactured explicit G was vacuous")
    # Do not enforce H-H0 == 0 for a first-order lagged
    # explicit term. Record the signed defect for separate
    # time-refinement and physical-power analysis.
    return {
        "manufactured_nonzero_G_initial_data":manufactured,
        "p_initial_expression":"[y,-x]" if manufactured else "[0.,0.]",
        "time_scheme":"bdf",
        "time_steps":steps,"dt":dt,
        "time_end":float(system.last_ts_time),
        "TS_converged_reason":system.last_ts_converged_reason,
        "TS_accepted_steps":system.last_ts_steps,
        "TS_rejected_steps":system.last_ts_rejections,
        "TS_SN​​ES_failures":system.last_ts_snes_failures,
        "main_explicit_prestep":events,
        "main_IFunction_checks":state["main_residual_checks"],
        "max_frozen_explicit_RHS_drift":state["max_main_cache_drift"],
        "max_nonzero_explicit_G_norm":nonzero_gyroscopic,
        "times":times.tolist(),
        "full_state_trajectory":trajectory.tolist(),
        "hamiltonian":H.tolist(),
        "signed_H_final_minus_H_initial":float(H[-1]-H[0]),
        "relative_abs_H_change":float(abs(H[-1]-H[0])/H[0]),
        "exact_H_balance_claimed":False,
        "nonlinear_DAE_initialization_qualified":False,
    }


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--mesh",type=float,default=.2)
    p.add_argument("--dt",type=float,default=1e-4)
    p.add_argument("--steps",type=int,default=4)
    p.add_argument("--case",choices=("published","manufactured"),required=True)
    p.add_argument("--output-dir",type=Path,required=True)
    args=p.parse_args()
    if args.mesh<=0 or args.dt<=0 or args.steps<3:
        p.error("positive mesh/dt and at least 3 steps required")
    mpi=MPI.COMM_WORLD
    m=args.case=="manufactured"
    record=run_case(args.mesh,args.dt,args.steps,m,mpi)
    if mpi.Get_rank()==0:
        args.output_dir.mkdir(parents=True,exist_ok=True)
        filename=f"shallow_explicit_TS_{args.case}_{mpi.Get_size()}r.json"
        (args.output_dir/filename).write_text(
            json.dumps(record,indent=2)+"\n",encoding="utf-8"
        )
        print(
            f"SHALLOW PHYSICAL TS {args.case} MPI={mpi.Get_size()} "
            f"steps={record['TS_accepted_steps']} "
            f"reason={record['TS_converged_reason']} "
            f"G_norm={record['max_nonzero_explicit_G_norm']:.9g} "
            f"frozen_drift={record['max_frozen_explicit_RHS_drift']:.9g} "
            f"delta_H={record['signed_H_final_minus_H_initial']:.9g} "
            f"relative_delta_H={record['relative_abs_H_change']:.9g}",
            flush=True
        )


if __name__=="__main__":
    main()
