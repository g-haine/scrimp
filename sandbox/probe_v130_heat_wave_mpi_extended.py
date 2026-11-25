# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Strict, original Heat-Wave pH CN physics on real MPI1/MPI2 ranks.

Short stage: 4 x 0.001 CN, BOTH MPI1 and MPI2, full 7689-DOF trajectories.
Long stage: 3000 x 0.005 CN (15 s), MPI1, every accepted state saved.
The *canonical* heat_wave_eq runs, with only plots/PVD disabled. All
physical quadratures (including two independent oriented facets) run
on EVERY MPI rank. A snapshot every ten steps would NOT permit a CN
midpoint power identity; require dt_save=dt for THIS scientific stage.

Never replace GetFEM observables by synthetic/analytic fields, alter
DAE tolerances or weaken the strict physical acceptance conditions.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from mpi4py import MPI
from petsc4py import PETSc

from examples.heat_wave import heat_wave_eq
from measure_v130_hw_dual_facet_physics import (
    collect_heat_wave_dual_facet_physics,
)


def real(a, label):
    result=np.asarray(a)
    if np.iscomplexobj(result):
        if np.max(np.abs(result.imag))>1e-12:
            raise RuntimeError(f"Unexpected complex FE data: {label}")
        result=result.real
    result=np.asarray(result,dtype=np.float64)
    if not np.all(np.isfinite(result)):
        raise RuntimeError(f"Nonfinite GetFEM state or quantity: {label}")
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out",type=Path,required=True)
    parser.add_argument("--dt",type=float,required=True)
    parser.add_argument("--steps",type=int,required=True)
    parser.add_argument("--profile",choices=("gaussian","gallery_packet"),
                        default="gaussian")
    parser.add_argument("--store-full-states",action="store_true")
    args=parser.parse_args()
    comm=MPI.COMM_WORLD
    rank,size=comm.rank,comm.size
    if (size not in (1,2) or PETSc.COMM_WORLD.getSize()!=size
            or not np.isfinite(args.dt) or args.dt<=0
            or args.steps<4 or args.steps>10000
            or not args.out.is_dir()):
        raise RuntimeError("Invalid native MPI/time/output Heat-Wave contract")
    if not PETSc.Sys.hasExternalPackage("mumps"):
        raise RuntimeError("Pinned PETSc must have true MUMPS LU support")
    tf=float(args.steps*args.dt)

    # The exact same canonical Heat-Wave pH PDE and algebraic DAE
    # initialization as the previously qualified MPI1 example.
    # Only optional file-producing postprocessing is disabled.
    system=heat_wave_eq(
        heat_region=1,wave_region=2,factor_solver="mumps",
        final_time=tf,time_step=args.dt,save_step=args.dt,
        initial_wave_boundary="clamp",
        initial_wave_profile=args.profile,
        physics_only=True,
    )
    original=system.heat_wave_initialization_report
    if (original.get("method")!="clamp"
            or not original.get("boundary_compatibility_enforced")
            or original.get("boundary_region")!=20
            or original.get("after_boundary_l2_coefficients")!=0.
            or original.get("dae_tolerances_modified") is not False
            or original.get("boundary_DOFS",0)<=0):
        raise RuntimeError("Canonical true L2-projected FE wave clamp missing")

    compiler=getattr(system,"_internal_interface_compiler",None)
    if compiler is None or not compiler.completed:
        raise RuntimeError("Missing production conforming-interface compiler")
    heat_faces=int(system.ports["Interface Heat"].get_assembly_region())
    wave_faces=int(system.ports["Interface Wave"].get_assembly_region())
    if (heat_faces==wave_faces or 10 in (heat_faces,wave_faces)
            or heat_faces<=0 or wave_faces<=0):
        raise RuntimeError("True two-sided outward interface regions absent")

    times=real(system.solution["t"],"accepted CN times").reshape(-1)
    states=real(system.solution["z"],"original full FE history")
    H=real(system.get_Hamiltonian(),"original integrated Hamiltonian").reshape(-1)
    # Preserve exact native TS/monitor diagnostics BEFORE any strict
    # rejection. On source 42116b58, PETSc completed all 3000 steps,
    # but its time-threshold monitor dropped only accepted step 2999;
    # 3000 saved FE states cannot support 3000 separate CN identities.
    # Do not manufacture the missing state or relax the physical audit.
    measured_indices=np.rint(times/args.dt).astype(np.int64)
    present={int(i) for i in measured_indices
             if 0<=i<=args.steps}
    missing=sorted(set(range(args.steps+1))-present)
    preflight={
        "schema":"v130-heat-wave-original-CN-monitor-preflight-v1",
        "rank":rank,"MPI":size,
        "accepted_PETSc_TS_steps":int(system.last_ts_steps),
        "TS_converged_reason":int(system.last_ts_converged_reason),
        "TS_SNES_failures":int(system.last_ts_snes_failures),
        "TS_step_rejections":int(system.last_ts_rejections),
        "TS_final_time":float(system.last_ts_time),
        "expected_saved_FE_states":args.steps+1,
        "actual_saved_FE_states":int(len(times)),
        "H_values":int(len(H)),
        "FE_state_shape":list(states.shape),
        "missing_accepted_grid_indices":missing,
        "maximum_saved_clock_deviation":float(
            np.max(abs(times/args.dt-measured_indices))
        ) if len(times) else None,
        "first_saved_time":float(times[0]) if len(times) else None,
        "last_saved_time":float(times[-1]) if len(times) else None,
        "full_physics_gate_claimed":False,
    }
    (args.out/f"preflight-rank{rank}.json").write_text(
        json.dumps(preflight,indent=2)+"\n",encoding="utf-8"
    )
    if (times.shape!=(args.steps+1,)
            or H.shape!=times.shape or states.shape!=(len(times),7689)
            or not np.allclose(times,np.arange(args.steps+1)*args.dt,
                               atol=1e-10,rtol=0.)
            or system.last_ts_converged_reason<=0
            or system.last_ts_snes_failures!=0
            or system.last_ts_rejections!=0
            or system.last_ts_steps!=args.steps
            or abs(system.last_ts_time-tf)>1e-9
            or H[0]<=0.
            or np.max(abs(states[-1]-states[0]))<1e-12):
        raise RuntimeError(
            "Original GetFEM/PETSc fixed-CN full FE trace invalid: "
            f"TS_accepted={preflight['accepted_PETSc_TS_steps']} "
            f"saved={preflight['actual_saved_FE_states']}/"
            f"{preflight['expected_saved_FE_states']} "
            f"missing_grid_indices={preflight['missing_accepted_grid_indices'][:20]} "
            f"last_t={preflight['last_saved_time']} "
            f"final_TS={preflight['TS_final_time']}; "
            "see preflight-rankN.json (no synthetic FE states)"
        )

    traces,diagnostic=collect_heat_wave_dual_facet_physics(
        system,times,H,heat_faces=heat_faces,wave_faces=wave_faces,
        heat_bulk=1,wave_bulk=2,
        expected_dt=args.dt,expected_steps=args.steps,
    )
    scale=max(1.,float(np.max(np.abs(H))))
    step_scale=max(scale,float(np.max(
        abs(np.diff(times)*traces["CN_D_heat_flux"][1:])
    )))
    heat=real(traces["H_heat"],"heat energy")
    wave=real(traces["H_wave"],"wave energy")
    components=real(traces["H_components"],"independent components")
    d=real(traces["CN_D_heat_flux"][1:],"CN flux dissipation")
    pi=real(traces["CN_P_interface_physical"][1:],"signed interface")
    po=real(traces["CN_P_outer"][1:],"outer physical work")
    he=real(traces["CN_per_step_heat_balance_error"][1:],"heat step balance")
    wa=real(traces["CN_per_step_wave_balance_error"][1:],"wave step balance")
    total=real(traces["CN_per_step_H_balance_error"][1:],"total step balance")
    cum=real(traces["CN_cumulative_total_balance_error"],"cumulative balance")
    step_independent=np.diff(H)+np.diff(times)*(d-pi-po)
    # PASS requires independent complete histories, not just a few
    # end-point scalar values or PETSc's "converged" reason.
    gates={
        "genuine_MUMPS_Crank_Nicolson_all_steps":True,
        "strict_DAE_compatible_original_clamp":True,
        "separate_heat_and_wave_energies_nonnegative":
            bool(np.min(heat)>=-1e-9*scale and np.min(wave)>=-1e-9*scale),
        "original_H_equals_independent_components":
            bool(np.max(abs(H-components))<1e-8*scale),
        "positive_nontrivial_thermal_dissipation":
            bool(np.min(d)>=-1e-10*scale and np.max(d)>1e-10*scale),
        "oriented_two_sided_interface_power_cancels":
            bool(np.max(abs(pi))<1e-8*scale),
        "zero_outer_boundary_work":
            bool(np.max(abs(po))<1e-8*scale),
        "each_CN_heat_subsystem_balance":
            bool(np.max(abs(he))<1e-8*step_scale),
        "each_CN_wave_subsystem_balance":
            bool(np.max(abs(wa))<1e-8*step_scale),
        "each_CN_total_balance_independently":
            bool(max(np.max(abs(total)),np.max(abs(step_independent)))
                 <1e-8*step_scale),
        "complete_cumulative_energy_power_identity":
            bool(np.max(abs(cum))<1e-8*scale),
        "monotone_total_energy_at_every_accepted_step":
            bool(np.max(np.diff(H))<1e-8*step_scale),
    }
    if size==2:
        # Every MPI rank independently computes full original FE and
        # GetFEM integral traces. Consensus is tested here and again
        # from their individually archived original outputs on host.
        gathered=comm.allgather({
            "H_hash":hashlib.sha256(H.tobytes()).hexdigest(),
            "z_hash":hashlib.sha256(np.ascontiguousarray(states).tobytes()).hexdigest(),
            "physical_gates":gates,
            "sided_regions":(heat_faces,wave_faces),
        })
        rank_consensus=(gathered[0]==gathered[1])
    else:
        rank_consensus=True
    gates["all_MPI_ranks_agree_exactly"]=bool(rank_consensus)

    output={
        "schema":"scrimp-v130-heat-wave-MPI-and-long-CN-physical-v1",
        "status":("PASS_NATIVE_STRICT_DAE_CN_ORIENTED_POWER"
                  if all(gates.values()) else "FAIL_NATIVE_PHYSICAL_GATES"),
        "profile":args.profile,"MPI_ranks":size,"rank":rank,
        "dt":float(args.dt),"t_final":tf,"saved_all_accepted_steps":True,
        "steps":args.steps,"DOFs":int(states.shape[1]),
        "initial_wave_boundary":original,
        "initial_wave_profile":system.heat_wave_initialization_profile,
        "oriented_heat_faces":heat_faces,"oriented_wave_faces":wave_faces,
        "physical_gates":gates,
        "measurements":{
            "H0":float(H[0]),"Hf":float(H[-1]),
            "max_total_H_increase":float(np.max(np.diff(H))),
            "max_heat_power_defect":float(np.max(abs(he))),
            "max_wave_power_defect":float(np.max(abs(wa))),
            "max_total_power_defect":float(np.max(abs(step_independent))),
            "max_cumulative_power_defect":float(np.max(abs(cum))),
            "max_oriented_interface_power":float(np.max(abs(pi))),
            "max_outer_power":float(np.max(abs(po))),
            "max_thermal_dissipation":float(np.max(d)),
        },
        "state_history_sha256":hashlib.sha256(
            np.ascontiguousarray(states).tobytes()
        ).hexdigest(),
        "H_history_sha256":hashlib.sha256(
            np.ascontiguousarray(H).tobytes()
        ).hexdigest(),
        "all_rank_native_assembly_and_quadratures":True,
        "parallel_ParaView_export_qualified":False,
        "whispering_gallery_mode_qualified":False,
        "diagnostic_summary":diagnostic,
    }
    bundle={
        "t":times,"H":H,
        "z0":states[0].copy(),"zf":states[-1].copy(),
        **{"phys_"+k:real(v,k) for k,v in traces.items()
           if isinstance(v,np.ndarray) and v.shape==times.shape},
    }
    if args.store_full_states:
        bundle["z"]=states
    filename=args.out/f"rank{rank}.npz"
    record=args.out/f"rank{rank}.json"
    if filename.exists() or record.exists():
        raise FileExistsError("Original rank evidence may not be overwritten")
    np.savez_compressed(filename,**bundle)
    record.write_text(json.dumps(output,indent=2)+"\n",encoding="utf-8")
    comm.Barrier()
    if rank==0:
        print("HEAT_WAVE_CN_NATIVE",json.dumps({
            "status":output["status"],"MPI":size,"steps":args.steps,
            "profile":args.profile,"physical_gates":gates,
            "measurements":output["measurements"],
        }),flush=True)
    if not all(gates.values()):
        raise RuntimeError("Real MPI Heat-Wave physical balance failed, evidence kept")


if __name__=="__main__":
    main()
