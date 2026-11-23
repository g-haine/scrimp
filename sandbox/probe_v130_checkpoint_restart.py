# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See LICENSE for license information.

"""Issue #29: ORIGINAL physical SWE cold checkpoint/restart native pilot.

Phases: seed [0,0.02], continuous [0,0.024], resume [0.02,0.024].
All use original curved P2 viscous Saint-Venant weak bricks and true
C1 facet-tangent wall law. The only changed physical experiment vs
historical invalid rotating q0 is the separately validated rest ramp.

No production DPHS checkpoint is installed. Sandbox-only intercepted
init_step is the injection point AFTER PETSc/GetFEM assembly and BEFORE
the real original main TS. Resume uses a new real TS, all algebraic
and dynamic GetFEM coefficients, and absolute original control time.
No fake step, no BDF/pseudo initializer, no PETSc history claim.
"""

import argparse
import hashlib
import json
import math
from pathlib import Path
from unittest.mock import patch

import numpy as np
from mpi4py import MPI
from petsc4py import PETSc
import scrimp as S

from viscous_swe_2d_mpi import build_system, diagnostics, ROTATION
from probe_v130_issue28_structure import indices, layout
from probe_v130_issue28_one_cn_spinup import (
    native_residual, publish_sandbox_rest,
)
from probe_v130_issue28_weak_wall import ensure_real
from v130_checkpoint_state import (
    collective_read, collective_write, model_provenance, sha256,
)

MESHES=(.2, .1)
DT=.002
SPINUP_END=.02
FINAL_TIME=.024


def original_physical_model(mesh, final_time):
    """Exactly the same original physical weak bricks for all phases."""
    return build_system(
        gamma=0., mesh=mesh, geometry_order=2, height_order=2,
        final_time=final_time, save_dt=DT*(1.-1e-9), max_step=DT,
        ramp_start=.03, ramp_end=.06, init_step=True,
        ts_type="cn", adapt_type="none", init_dt=1e-9,
        rotation=ROTATION, initial_profile="rigid",
        wall_velocity_profile="facet_tangent",
        algebraic_initial_guess="zero",
        start_from_rest_duration=.02, init_method="bdf",
    )


def validate_actual_ts(system, initial_time, final_time, expected_steps):
    actual=(float(system.last_ts_time),
            int(system.last_ts_steps),
            int(system.last_ts_converged_reason),
            int(system.last_ts_rejections),
            int(system.last_ts_snes_failures))
    if (not math.isfinite(actual[0])
            or abs(actual[0]-final_time) > 1e-12
            or actual[1] != expected_steps or actual[2] <= 0
            or actual[3] or actual[4]):
        raise RuntimeError(
            f"Real PETSc TS did not accept exactly {expected_steps} "
            f"steps from {initial_time} to {final_time}: {actual!r}"
        )
    times=np.asarray(system.solution["t"],dtype=float)
    states=np.asarray([
        ensure_real(z,"accepted original physical GetFEM z").copy()
        for z in system.solution["z"]
    ])
    _,_,_,n=layout(system)
    expected=np.linspace(initial_time,final_time,expected_steps+1)
    if (times.shape != expected.shape
            or states.shape != (expected_steps+1,n)
            or not np.all(np.isfinite(states))
            or not np.allclose(times,expected,rtol=0.,atol=1e-12)):
        raise RuntimeError(
            "Physical PETSc TS accepted states not fully archived; "
            "cannot create an unqualified checkpoint"
        )
    return times,states


def validate_original_equations(system, mpi, times, states):
    differential, algebraic, iv, n = layout(system)
    indices_a=indices(algebraic,iv)
    indices_d=indices(differential,iv)
    zero=np.zeros(n)
    cn_norms=[]
    fa_norms=[]
    for k in range(len(times)-1):
        t0,t1=map(float,(times[k],times[k+1]))
        z0,z1=states[k],states[k+1]
        derivative=(z1-z0)/(t1-t0)
        _,old= native_residual(
            system,mpi,z0,zero,t0,iv,(differential,algebraic),
            with_vector=True
        )
        _,new= native_residual(
            system,mpi,z1,2.*derivative,t1,iv,
            (differential,algebraic),with_vector=True
        )
        total=old+new
        stage=float(np.linalg.norm(total))
        fa=float(np.linalg.norm(
            native_residual(
                system,mpi,z1,zero,t1,iv,
                (differential,algebraic),with_vector=True
            )[1][indices_a]
        ))
        if (not math.isfinite(stage) or stage > 1e-7
                or not math.isfinite(fa) or fa > 1e-8):
            raise RuntimeError(
                f"Original GetFEM/PETSc TSCN and DAE rejected step {k+1}: "
                f"||F_TSCN||={stage:.9g}, ||Fa||={fa:.9g}"
            )
        cn_norms.append(stage)
        fa_norms.append(fa)
    return differential, algebraic, iv, indices_d, {
        "maximum_actual_original_TSCN_residual":max(cn_norms),
        "maximum_actual_original_Fa_residual":max(fa_norms),
        "full_step_TSCN_residuals":cn_norms,
        "full_step_Fa_residuals":fa_norms,
    }


def publish_checkpoint_at_absolute_time(system, mpi, saved, path):
    """Sandbox ONLY. No pseudo/BDF or algebraic projection is called."""
    full_before=ensure_real(system.gf_model.from_variables(),
                            "newly reconstructed GetFEM model").copy()
    fingerprint=model_provenance(system)
    t, z, header=collective_read(
        path, mpi=mpi,current_provenance=fingerprint,ts_type="cn"
    )
    if abs(t-SPINUP_END)>1e-12 or len(z)!=len(full_before):
        raise RuntimeError("Checkpoint does not belong to t*=0.02 physical ramp")
    requested_t0=float(system.time_scheme["t_0"])
    if abs(t-requested_t0)>1e-12:
        raise RuntimeError("PETSc absolute t0 differs from verified checkpoint")
    if not np.array_equal(
        ensure_real(system.gf_model.from_variables(),"GetFEM untouched"),
        full_before
    ):
        raise RuntimeError("Checkpoint preflight modified GetFEM before validation")

    differential, algebraic, iv, n=layout(system)
    original_fa=native_residual(
        system,mpi,z,np.zeros(n),t,iv,(differential,algebraic)
    )
    if (not math.isfinite(original_fa["algebraic_l2"])
            or original_fa["algebraic_l2"]>1e-9):
        raise RuntimeError(
            "Original complete DAE checkpoint is incompatible with "
            "its exact physical boundary law: "
            f"Fa={original_fa['algebraic_l2']:.9g}"
        )
    # This publication is the only mutation of the initial FE state.
    # Absolute physical time is set BEFORE the original PETSc TS starts.
    system.gf_model.set_time(t)
    system.gf_model.to_variables(z.copy())
    system.GlobalVec.getArray()[:] = z
    accepted=np.array_equal(
        ensure_real(system.gf_model.from_variables(),"restored GetFEM"),z
    ) and np.array_equal(
        ensure_real(system.GlobalVec.getArray(readonly=True),
                    "restored replicated global vector"),z
    )
    if not mpi.allreduce(bool(accepted),op=MPI.LAND):
        raise RuntimeError("Bitwise full-coefficient GetFEM/PETSc restore failed")
    system.initialization_report={
        "status":"SANDBOX_CHECKPOINT_COLD_RESTORED_NOT_PRODUCTION",
        "method":"complete_accepted_state_at_absolute_time",
        "original_algebraic_residual_l2": original_fa["algebraic_l2"],
        "checkpoint_sha256": sha256(z.tobytes()),
        "getfem_globalvec_published_bitwise": True,
        "pseudo_or_bdf_initialization_performed":False,
        "getfem_initial_reconstruction_left_unchanged_until_validation":True,
    }
    saved["state"]=z.copy()
    saved["time"]=t
    saved["model_fingerprint_sha256"]=fingerprint["fingerprint_sha256"]
    saved["header"]=header
    if mpi.rank==0:
        print(
            f"CHECKPOINT RESTORE PASS ranks={mpi.size} "
            f"t*={t:.9g} Fa={original_fa['algebraic_l2']:.9g} "
            "all algebraic and physical variables restored bitwise",
            flush=True
        )


def run():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--phase",choices=("seed","continuous","resume"),required=True)
    p.add_argument("--mesh",type=float,choices=MESHES,required=True)
    p.add_argument("--checkpoint",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    mpi=PETSc.COMM_WORLD.tompi4py()
    t0=SPINUP_END if args.phase=="resume" else 0.
    tf=SPINUP_END if args.phase=="seed" else FINAL_TIME
    n_steps=int(round((tf-t0)/DT))
    system=original_physical_model(args.mesh,tf)
    saved={}
    if args.phase=="resume":
        system.set_time_scheme(
            t_0=SPINUP_END,t_f=tf,dt=DT,
            dt_save=DT*(1.-1e-9),init_step=True,
        )
        def injected(self):
            publish_checkpoint_at_absolute_time(
                self,mpi,saved,args.checkpoint
            )
    else:
        def injected(self):
            publish_sandbox_rest(self,mpi,saved)

    # Exact native solver with ONLY sandbox interception of the original
    # init_step. On resume no physical/algebraic initializer is called.
    with patch.object(S.DPHS,"init_step",injected):
        system.solve()
    times,states=validate_actual_ts(system,t0,tf,n_steps)
    dnames,anames,iv,dindices,audit=validate_original_equations(
        system,mpi,times,states
    )
    if args.phase=="resume":
        if not np.array_equal(states[0],saved["state"]):
            raise RuntimeError("Main TSCN first accepted state differs from checkpoint")
        init_report=system.initialization_report
        if (init_report.get("pseudo_or_bdf_initialization_performed") is not False
                or init_report.get("getfem_globalvec_published_bitwise") is not True):
            raise RuntimeError("Real cold restart silently reinitialized state")
    else:
        init_report=system.initialization_report
        if init_report.get("strategy")!="sandbox_exact_stationary_spinup_only":
            raise RuntimeError("Original physical initial state was not qualified")
    # Check model identity outside the solve and keep rank symmetry.
    fingerprint=model_provenance(system)
    original_end=np.asarray(states[-1],dtype="<f8")
    if args.phase=="seed":
        accepted=(
            int(system.last_ts_converged_reason)>0
            and int(system.last_ts_rejections)==0
            and int(system.last_ts_snes_failures)==0
            and abs(float(times[-1])-SPINUP_END)<=1e-12
            and audit["maximum_actual_original_Fa_residual"]<=1e-8
        )
        collective_write(
            args.checkpoint,mpi=mpi,z=original_end,
            t=float(times[-1]),provenance=fingerprint,ts_type="cn",
            accepted=accepted,
            extra={
                "source":"native_GetFEM_PETSc_real_TSCN_accepted_state",
                "nominal_mesh_size":args.mesh,
                "physical_spinup_duration":SPINUP_END,
                "physical_rotation_target":ROTATION,
                "physical_control_clock":"absolute",
                "original_TSCN_max":audit["maximum_actual_original_TSCN_residual"],
                "original_Fa_max":audit["maximum_actual_original_Fa_residual"],
                "source_mpi_ranks":mpi.size,
            },
        )
    power=diagnostics(system,gamma=0.)
    fields=("H","volume","volume_balance","supplied","dissipated",
            "p_normal","p_viscous","d_grad","d_div","balance")
    traces={key:np.asarray(power[key],dtype=float) for key in fields}
    for key,values in traces.items():
        if values.shape != times.shape or not np.all(np.isfinite(values)):
            raise RuntimeError(f"Original physical trace malformed: {key}")
    metrics={
        "schema":"scrimp-v130-issue29-native-cold-restart-pilot-v1",
        "status":"SANDBOX_NATIVE_SEGMENT_DIAGNOSTIC_NOT_PRODUCTION",
        "phase":args.phase,"mesh":args.mesh,"mpi_ranks":mpi.size,
        "source_time":t0,"actual_final_time":tf,"dt":DT,
        "actual_ts_steps":n_steps,"actual_ts_reason":int(system.last_ts_converged_reason),
        "actual_ts_rejections":int(system.last_ts_rejections),
        "actual_ts_snes_failures":int(system.last_ts_snes_failures),
        "model_fingerprint_sha256":fingerprint["fingerprint_sha256"],
        "original_h_p_state_change_l2":
            float(np.linalg.norm(states[-1,dindices]-states[0,dindices])),
        "complete_initial_state_sha256":sha256(states[0].tobytes()),
        "complete_final_state_sha256":sha256(states[-1].tobytes()),
        "exact_original_equations":audit,
        "initialization_report":init_report,
        "hamiltonian_initial":float(traces["H"][0]),
        "hamiltonian_final":float(traces["H"][-1]),
        "hamiltonian_change":float(traces["H"][-1]-traces["H"][0]),
        "original_boundary_work":float(traces["supplied"][-1]),
        "original_viscous_dissipation":float(traces["dissipated"][-1]),
        "original_H_plus_D_minus_W":float(traces["balance"][-1]),
        "original_mass_balance":float(traces["volume_balance"][-1]),
        "cold_restart_not_identical_TS_history":True,
        "general_high_index_initializer_qualified":False,
        "production_checkpoint_qualified":False,
    }
    if mpi.rank==0:
        args.out.mkdir(parents=True,exist_ok=True)
        (args.out/"segment_report.json").write_text(
            json.dumps(metrics,indent=2)+"\n",encoding="utf-8"
        )
        np.savez_compressed(
            args.out/"segment_states.npz",
            times=times,states=states,**traces
        )
        print(
            f"CHECKPOINT PHASE={args.phase} mesh={args.mesh:g} "
            f"MPI={mpi.size} steps={n_steps} "
            f"Fa_max={audit['maximum_actual_original_Fa_residual']:.8g} "
            f"TSCN_max={audit['maximum_actual_original_TSCN_residual']:.8g} "
            f"power={metrics['original_H_plus_D_minus_W']:.8g}",
            flush=True,
        )
    mpi.Barrier()


if __name__=="__main__":
    run()
