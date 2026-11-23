# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Issue #29: independent strict native GetFEM/PETSc checkpoint restart gate.

Compares physical 0->.024 full TSCN with 0->.02 + NEW TS .02->.024,
including MPI1->1,1->2,2->1,2->2, all original FE coefficients,
algebraic constraints and true original h-weighted port/viscous power.

Pilot tests cold restart, not restoring PETSc TS history or qualifying
the production DPHS API or a general high-index initializer.
"""

import argparse
import json
import math
from pathlib import Path

import numpy as np

from v130_checkpoint_state import read_checkpoint

MESHES = (.2, .1)
RANKS = (1, 2)
DT = .002
TSTAR = .02
TEND = .024
SCHEMA = "scrimp-v130-issue29-native-cold-restart-pilot-v1"


def require(test, message):
    if not test:
        raise ValueError(message)


def finite(v, name):
    require(type(v) in (int,float) and math.isfinite(v),
            f"Invalid {name} in physical checkpoint evidence")
    return float(v)


def segment(root, mesh, name, rank, source_rank=None):
    loc = f"h{str(mesh).replace('.', 'p')}"
    if name == "resume":
        path = root/f"{loc}_restart_from{source_rank}_to{rank}"
    else:
        path = root/f"{loc}_{name}_mpi{rank}"
    rp = path/"segment_report.json"
    zp = path/"segment_states.npz"
    require(rp.is_file() and zp.is_file(),
            f"Missing real native original checkpoint segment {path}")
    report=json.loads(rp.read_text(encoding="utf-8"))
    require(report.get("schema")==SCHEMA
            and report.get("status")==
            "SANDBOX_NATIVE_SEGMENT_DIAGNOSTIC_NOT_PRODUCTION"
            and report.get("phase")==name
            and report.get("mpi_ranks")==rank
            and report.get("mesh")==mesh,
            "Unexpected original physical checkpoint model/provenance")
    expected_begin=TSTAR if name=="resume" else 0.
    expected_end=TSTAR if name=="seed" else TEND
    requested=int(round((expected_end-expected_begin)/DT))
    require(report.get("actual_ts_reason",0)>0
            and report.get("actual_ts_steps")==requested
            and report.get("actual_ts_rejections")==0
            and report.get("actual_ts_snes_failures")==0
            and abs(finite(report.get("source_time"),"segment start")-
                    expected_begin)<1e-12
            and abs(finite(report.get("actual_final_time"),"segment end")-
                    expected_end)<1e-12,
            "PETSc original physical segment not genuinely accepted")
    constraints=report.get("exact_original_equations",{})
    require(len(constraints.get("full_step_TSCN_residuals",[]))==requested
            and len(constraints.get("full_step_Fa_residuals",[]))==requested,
            "Original GetFEM per-step physical DAE diagnostics missing")
    for k in range(requested):
        require(finite(
            constraints["full_step_TSCN_residuals"][k],"TSCN step"
        )<1e-7 and finite(
            constraints["full_step_Fa_residuals"][k],"algebraic step"
        )<1e-8, "Real original TSCN/algebraic equations fail")
    require(finite(
        constraints.get("maximum_actual_original_TSCN_residual"),"max TSCN"
    )<1e-7 and finite(
        constraints.get("maximum_actual_original_Fa_residual"),"max Fa"
    )<1e-8,"Max original DAE/TSCN constraint violated")
    require(report.get("cold_restart_not_identical_TS_history") is True
            and report.get("general_high_index_initializer_qualified") is False
            and report.get("production_checkpoint_qualified") is False,
            "Unjustified PETSc history/production initialization qualification")
    if name=="resume":
        init=report.get("initialization_report",{})
        require(init.get("status")==
                "SANDBOX_CHECKPOINT_COLD_RESTORED_NOT_PRODUCTION"
                and init.get("getfem_globalvec_published_bitwise") is True
                and init.get("pseudo_or_bdf_initialization_performed") is False
                and finite(init.get("original_algebraic_residual_l2"),
                           "initial cold checkpoint Fa")<1e-9,
                "Checkpoint was reprojected or incompatible at restoration")
    else:
        require(report.get("initialization_report",{}).get("strategy")==
                "sandbox_exact_stationary_spinup_only",
                "Original t0 physical startup not separately native-qualified")
    with np.load(zp,allow_pickle=False) as stored:
        required={"times","states","H","volume","volume_balance",
                  "supplied","dissipated","p_normal","p_viscous",
                  "d_grad","d_div","balance"}
        require(required<=set(stored.files),
                "Original full FE coefficient and power arrays incomplete")
        times=np.asarray(stored["times"],dtype=float)
        states=np.asarray(stored["states"],dtype=float)
        arrays={key:np.asarray(stored[key],dtype=float)
                for key in required-{"times","states"}}
    require(times.shape==(requested+1,)
            and states.shape[0]==requested+1
            and states.ndim==2
            and states.shape[1]>0
            and np.all(np.isfinite(states))
            and np.allclose(
                times,np.linspace(expected_begin,expected_end,requested+1),
                rtol=0.,atol=1e-12
            ),"Wrong or incomplete native absolute physical timeline")
    for name_, v in arrays.items():
        require(v.shape==times.shape and np.all(np.isfinite(v)),
                "Original h-weighted physical power trace invalid: "+name_)
    require(np.allclose(
        arrays["balance"], arrays["H"]-arrays["H"][0]
        +arrays["dissipated"]-arrays["supplied"],
        rtol=0,atol=1e-11
    ),"Original Hamiltonian/work/dissipation identity not reproducible")
    require(np.array_equal(
        states[0].astype("<f8"), states[0]
    ),"Invalid original FE coefficient scalar format")
    scalars={
        "dH":float(arrays["H"][-1]-arrays["H"][0]),
        "W":float(arrays["supplied"][-1]),
        "D":float(arrays["dissipated"][-1]),
        "volume_defect":float(arrays["volume_balance"][-1]),
    }
    require(abs(scalars["volume_defect"])<1e-9,
            "Physical mass balance failed")
    scalars["balance"]=scalars["dH"]+scalars["D"]-scalars["W"]
    for report_key,expected in (
        ("hamiltonian_change",scalars["dH"]),
        ("original_boundary_work",scalars["W"]),
        ("original_viscous_dissipation",scalars["D"]),
        ("original_H_plus_D_minus_W",scalars["balance"]),
    ):
        require(math.isclose(
            finite(report.get(report_key),report_key),
            expected,rel_tol=1e-8,abs_tol=1e-11
        ),"Original physical segment power report mismatch: "+report_key)
    return {"report":report,"times":times,"states":states,
            "physics":scalars,"traces":arrays}


def compare(root):
    root=Path(root)
    manifest=json.loads((root/"manifest.json").read_text(encoding="utf-8"))
    require(manifest.get("schema")=="scrimp-v130-issue29-runner-v1"
            and manifest.get("original_physics_untouched") is True
            and manifest.get("cold_state_restart_only") is True
            and manifest.get("source_spinup_time")==TSTAR
            and manifest.get("physical_final_time")==TEND,
            "Native checkpoint pilot manifest/provenance incorrect")
    final={}
    for mesh in MESHES:
        baseline={r:segment(root,mesh,"continuous",r) for r in RANKS}
        seed={r:segment(root,mesh,"seed",r) for r in RANKS}
        rows={}
        for source in RANKS:
            checkpoint=root/(
                f"h{str(mesh).replace('.', 'p')}_source_mpi{source}"
            )/"accepted_state.npz"
            require(checkpoint.is_file(),"Native original seed checkpoint missing")
            # Validate manifest, accepted full FE byte hash and
            # version before comparing actual native seed coefficients.
            with np.load(checkpoint,allow_pickle=False) as x:
                head=json.loads(x["manifest"].tobytes().decode("utf-8"))
            time,z,head=read_checkpoint(
                checkpoint,current_provenance=head["model"]
            )
            require(abs(time-TSTAR)<1e-12
                    and z.size==seed[source]["states"].shape[1]
                    and np.array_equal(z,seed[source]["states"][-1]),
                    "Saved accepted FE checkpoint differs from source native TS")
            for target in RANKS:
                replay=segment(root,mesh,"resume",target,source_rank=source)
                require(np.array_equal(
                    replay["states"][0],z
                ),"Cold GetFEM/PETSc restart silently changed any FE coefficient")
                original=baseline[target]
                end_norm=float(np.linalg.norm(
                    original["states"][-1]-replay["states"][-1]
                ))
                base_change=float(np.linalg.norm(
                    original["states"][-1]-original["states"][0]
                ))
                tolerance=5e-8+2e-6*base_change
                require(end_norm<=tolerance,
                        f"Physical continuous/split full FE state mismatch "
                        f"mesh={mesh} MPI {source}->{target}: "
                        f"{end_norm:.9g}>{tolerance:.9g}")
                split={
                    name:seed[source]["physics"][name]+replay["physics"][name]
                    for name in ("dH","W","D")
                }
                split["balance"]=split["dH"]+split["D"]-split["W"]
                diffs={
                    name:abs(split[name]-original["physics"][name])
                    for name in ("dH","W","D","balance")
                }
                for name,v in diffs.items():
                    require(v<2e-6,
                            f"Physical continuous/split h-weighted {name} "
                            f"mismatch after MPI {source}->{target}")
                # H must be continuous at the exact stored restart time.
                h_jump=abs(
                    seed[source]["traces"]["H"][-1]-
                    replay["traces"]["H"][0]
                )
                require(h_jump<1e-9,
                        "Physical Hamiltonian jumped at cold restart")
                rows[f"{source}_to_{target}"]={
                    "full_GetFEM_endpoint_l2":end_norm,
                    "allowed_full_state_norm":tolerance,
                    "Hamiltonian_restart_jump":h_jump,
                    "continuous_vs_split_signed_power_differences":diffs,
                    "all_algebraic_variables_restored_bitwise":True,
                    "original_control_absolute_time_preserved":True,
                    "PETSc_TS_history_restored":False,
                }
                print(
                    f"CHECKPOINT PASS mesh={mesh:g} MPI={source}->{target} "
                    f"full_FE_error={end_norm:.8g} "
                    f"power_diff={diffs['balance']:.8g}",
                    flush=True,
                )
        final[f"h{mesh:g}"]=rows
    output={
        "schema":"scrimp-v130-issue29-native-checkpoint-comparison-v1",
        "status":"NATIVE_COLD_CHECKPOINT_SPLIT_VS_CONTINUOUS_PILOT_PASSED",
        "cases":final,
        "production_DPHS_checkpoint_API_qualified":False,
        "history_exact_BDF_IMEX_restoration_qualified":False,
        "general_high_index_DAE_initializer_qualified":False,
        "second_independent_linear_DAE_model_qualified":False,
        "interpretation":(
            "Sandbox accepted-state GetFEM/PETSc cold restart and 1/2 cross "
            "MPI replays match the true original physical continuous "
            "Saint-Venant trajectory in whole FE vector, original DAE "
            "and weighted H+viscosity-boundary-work; "
            "not a qualified production API or PETSc integrator-history reset."
        ),
    }
    (root/"checkpoint_comparison.json").write_text(
        json.dumps(output,indent=2)+"\n",encoding="utf-8"
    )
    print("COMPLETE: native original physical cold restart pilot; NO production checkpoint API qualification.")
    return output


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--evidence",type=Path,required=True)
    args=p.parse_args()
    compare(args.evidence)


if __name__=="__main__":
    main()
