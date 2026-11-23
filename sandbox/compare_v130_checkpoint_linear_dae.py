# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See LICENSE for license information.

"""Independent #29 production DPHS API linear pH-DAE wave restart audit.

Same original accepted mesh/FEM and absolute control clock; 8 native
real GetFEM/PETSc runs: two seeds, two uninterrupted, four cross-MPI
cold restarts. Scientific proof for this second model is PENDING until
native evidence is supplied. No BDF/IMEX history restart is claimed.
"""

import argparse
import importlib.util
import json
import math
from pathlib import Path

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location(
    "_scrimp_ckpt_linear_host",ROOT/"scrimp"/"utils"/"dae_checkpoint.py",
)
ckpt=importlib.util.module_from_spec(spec)
spec.loader.exec_module(ckpt)
RANKS=(1,2)
DT=.002
TSTAR=.02
TFINAL=.024
SCHEMA="scrimp-v130-issue29-linear-DAE-public-API-v1"


def require(condition,reason):
    if not condition:
        raise ValueError(reason)


def case_path(phase,ranks,source=None):
    if phase=="resume":
        require(source in RANKS,"Missing source rank for linear DAE restart")
        return f"wave_restart_from{source}_to{ranks}"
    require(phase in ("seed","continuous"),"Unknown DAE phase")
    return f"wave_{phase}_mpi{ranks}"


def checkpoint(root,source):
    path=Path(root)/f"wave_source_mpi{source}"/"accepted_state.npz"
    require(path.is_file(),"Missing linear DAE source checkpoint")
    with path.open("rb") as stream:
        with np.load(stream,allow_pickle=False) as x:
            raw=x["manifest"]
            require(raw.ndim==1 and raw.dtype==np.dtype("uint8"),
                    "Invalid linear DAE checkpoint manifest")
            meta=json.loads(raw.tobytes().decode("utf-8"))
    t,z,meta=ckpt.read_checkpoint(
        path,current_provenance=meta["model"],ts_type="cn",
    )
    require(abs(t-TSTAR)<1e-12,"Wrong linear DAE checkpoint absolute t*")
    extra=meta.get("extra",{})
    require(extra.get("source_MPI_ranks")==source
            and extra.get("source")=="original_examples_wave_linear_DAE"
            and extra.get("physical_time_origin")=="absolute"
            and extra.get("control_law")=="0.2*sin(25*t)*sin(pi*y)",
            "Wrong physical linear DAE checkpoint control/time provenance")
    return z,meta


def segment(root,phase,ranks,source=None):
    root=Path(root)
    name=case_path(phase,ranks,source)
    folder=root/name
    rp=folder/"linear_report.json"
    zp=folder/"linear_states.npz"
    log=root/f"{name}.log"
    require(rp.is_file() and zp.is_file() and log.is_file(),
            "Missing native linear DAE segment or raw log: "+name)
    r=json.loads(rp.read_text(encoding="utf-8"))
    begin=TSTAR if phase=="resume" else 0.
    end=TSTAR if phase=="seed" else TFINAL
    n=int(round((end-begin)/DT))
    require(
        r.get("schema")==SCHEMA
        and r.get("status")==
            "NATIVE_LINEAR_DAE_SEGMENT_NOT_YET_CROSS_MPI_COMPARED"
        and r.get("phase")==phase
        and r.get("MPI")==ranks
        and r.get("mesh_h")==.5
        and abs(r.get("t0",-1)-begin)<1e-12
        and abs(r.get("tf",-1)-end)<1e-12
        and r.get("step_count")==n
        and r.get("PETSc_reason",0)>0
        and r.get("PETSc_rejected")==0
        and r.get("SNES_failures")==0
        and r.get("not_SWE_model") is True
        and r.get("no_integrator_history_claimed") is True
        and r.get("production_API_used") is (phase in ("seed","resume")),
        "Actual linear DAE PETSc phase, provenance or public API mismatch"
    )
    require(
        f"LINEAR DAE PHASE={phase} MPI={ranks} " in log.read_text(
            encoding="utf-8",errors="replace"
        ),"Original native linear DAE completion marker absent"
    )
    audit=r.get("original_equations",{})
    require(len(audit.get("endpoint_Fa_residuals",[]))==n+1
            and len(audit.get("full_step_TSCN_residuals",[]))==n,
            "Original linear DAE all-steps residual evidence missing")
    for val in audit["endpoint_Fa_residuals"]:
        require(math.isfinite(val) and val<=1e-8,
                "Original linear DAE algebraic constraint failure")
    for val in audit["full_step_TSCN_residuals"]:
        require(math.isfinite(val) and val<=1e-7,
                "Original linear DAE TSCN failure")
    with zp.open("rb") as handle:
        with np.load(handle,allow_pickle=False) as npz:
            require(set(npz.files)=={
                "times","states","H","work","port_power",
            },"Missing linear DAE full FE trajectory or physical powers")
            times=np.asarray(npz["times"],dtype=float)
            states=np.asarray(npz["states"],dtype=float)
            h=np.asarray(npz["H"],dtype=float)
            w=np.asarray(npz["work"],dtype=float)
            powers=np.asarray(npz["port_power"],dtype=float)
    require(times.shape==(n+1,)
            and states.ndim==2
            and states.shape==(n+1,r.get("ndof"))
            and r["dynamic_ndof"]>0
            and r["algebraic_ndof"]>0
            and r["dynamic_ndof"]+r["algebraic_ndof"]==r["ndof"]
            and all(x.shape==times.shape for x in (h,w,powers))
            and all(np.all(np.isfinite(x)) for x in (
                times,states,h,w,powers
            ))
            and np.allclose(
                times,np.linspace(begin,end,n+1),atol=1e-12,rtol=0.
            ),"Linear DAE original full vector/physical timeline incomplete")
    require(np.all(np.abs(w-np.r_[0,np.cumsum(
        .5*np.diff(times)*(powers[:-1]+powers[1:])
    )])<1e-12),"Linear DAE independent boundary-work integration mismatch")
    require(abs(w[-1])>1e-12,
            "Linear DAE physical boundary drive and response vanished")
    require(abs(float(h[0])-r["Hamiltonian_initial"])<1e-12
            and abs(float(h[-1])-r["Hamiltonian_final"])<1e-12
            and abs(float(w[-1])-r["boundary_work"])<1e-12
            and abs(float(h[-1]-h[0]-w[-1])-r["dH_minus_W"])<1e-12,
            "Native linear DAE Hamiltonian/work report differs from traces")
    require(ckpt.sha256(states[0].tobytes())==r["initial_full_sha256"]
            and ckpt.sha256(states[-1].tobytes())==r["final_full_sha256"],
            "Linear DAE accepted full FE state checksum mismatch")
    if phase=="seed":
        state,meta=checkpoint(root,ranks)
        require(np.array_equal(state,states[-1])
                and meta["model"]["fingerprint_sha256"]==
                    r["fingerprint_sha256"],
                "Production DPHS saved checkpoint differs from seed native FE")
    if phase=="resume":
        state,meta=checkpoint(root,source)
        report=r["initialization_report"]
        require(np.array_equal(state,states[0])
                and meta["model"]["fingerprint_sha256"]==
                    r["fingerprint_sha256"]
                and report.get("status")==
                    "COMPLETE_DAE_STATE_COLD_RESTORED_IN_DPHS"
                and report.get(
                    "all_GetFEM_and_GlobalVec_coefficients_restored_bitwise"
                ) is True
                and report.get("no_pseudo_BDF_or_projection") is True
                and report.get("PETSc_TS_history_restored") is False
                and report.get("original_algebraic_residual_l2",float("inf"))<1e-9,
                "Actual DPHS linear DAE resume projected/changed q0 or violated Fa")
    return dict(report=r,times=times,states=states,H=h,
                work=w,power=powers)


def compare(root):
    root=Path(root)
    manifest=json.loads(
        (root/"manifest.json").read_text(encoding="utf-8")
    )
    require(manifest.get("schema")==
            "scrimp-v130-issue29-linear-DAE-runner-v1"
            and manifest.get("original_example")=="wave_2d_linear_DAE"
            and manifest.get("cold_state_restart_only") is True
            and manifest.get("physical_final_time")==TFINAL,
            "Native second-model checkpoint runner manifest invalid")
    baselines={r:segment(root,"continuous",r) for r in RANKS}
    seeds={r:segment(root,"seed",r) for r in RANKS}
    fingerprints={item["report"]["fingerprint_sha256"]
                  for item in [*baselines.values(),*seeds.values()]}
    require(len(fingerprints)==1,
            "Independent linear DAE original mesh, weak form or FE differ")
    checks={}
    for source in RANKS:
        for target in RANKS:
            restart=segment(root,"resume",target,source)
            baseline=baselines[target]
            seed=seeds[source]
            require(restart["report"]["fingerprint_sha256"] in fingerprints,
                    "Cross MPI linear DAE uses physically different model")
            full_error=float(np.linalg.norm(
                baseline["states"][-1]-restart["states"][-1]
            ))
            require(full_error<1e-8,
                    f"Linear DAE continuous/restarted full FE mismatch "
                    f"{source}->{target}: {full_error:.9g}")
            # Correct provenance for the full original conserved quantity
            # and physical boundary work, including the restart junction.
            hjump=abs(float(seed["H"][-1]-restart["H"][0]))
            dH_split=float(
                seed["H"][-1]-seed["H"][0]+
                restart["H"][-1]-restart["H"][0]
            )
            wsplit=float(seed["work"][-1]+restart["work"][-1])
            hdifference=abs(
                dH_split-(baseline["H"][-1]-baseline["H"][0])
            )
            wdifference=abs(wsplit-baseline["work"][-1])
            balance_difference=abs(
                (dH_split-wsplit)-
                (baseline["H"][-1]-baseline["H"][0]-
                 baseline["work"][-1])
            )
            require(hjump<1e-9
                    and hdifference<1e-8
                    and wdifference<1e-8
                    and balance_difference<1e-8,
                    "Linear DAE Hamiltonian/boundary-work not continuous "
                    f"on MPI {source}->{target}")
            checks[f"{source}_to_{target}"]={
                "full_original_FE_vector_l2":full_error,
                "original_Hamiltonian_restart_jump":hjump,
                "continuous_vs_split_Hamiltonian_change_l1":hdifference,
                "continuous_vs_split_boundary_work_l1":wdifference,
                "continuous_vs_split_signed_energy_defect_l1":balance_difference,
                "DAE_constraints_verified_each_real_TS_step":True,
                "cold_restart_exact_history_claimed":False,
            }
            print(
                f"LINEAR DAE CHECKPOINT PASS MPI={source}->{target} "
                f"full_FE_error={full_error:.8g} "
                f"power_diff={balance_difference:.8g}",
                flush=True,
            )
    result={
        "schema":"scrimp-v130-issue29-linear-DAE-comparison-v1",
        "status":"NATIVE_LINEAR_PHDae_DPHS_COLD_RESTART_PASSED",
        "cases":checks,
        "actual_original_model":"examples/wave.py linear mixed boundary DAE",
        "production_DPHS_save_load":True,
        "SWE_or_sloshing_model":False,
        "no_PETSc_TS_history_claim":True,
        "generic_checkpoint_supported_CN_only":True,
    }
    (root/"linear_dae_comparison.json").write_text(
        json.dumps(result,indent=2)+"\n",encoding="utf-8",
    )
    return result


def main():
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--evidence",type=Path,required=True)
    cli.add_argument("--verify-case",default=None)
    cli.add_argument("--phase",choices=("seed","continuous","resume"))
    cli.add_argument("--ranks",type=int,choices=RANKS)
    cli.add_argument("--source-ranks",type=int,choices=RANKS)
    args=cli.parse_args()
    if args.verify_case:
        require(args.verify_case==case_path(
            args.phase,args.ranks,args.source_ranks
        ),"Wrong requested linear DAE native phase name")
        segment(args.evidence,args.phase,args.ranks,args.source_ranks)
        print("PASS: original linear DAE native phase "+args.verify_case)
    else:
        compare(args.evidence)


if __name__=="__main__":
    main()
