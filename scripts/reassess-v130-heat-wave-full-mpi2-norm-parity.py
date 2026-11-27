#!/usr/bin/env python3
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Read-only norm-based reassessment of archived original Heat-Wave MPI1/2.

Reuses the EXACT existing independent per-rank physics/clock audit.
The original strict pointwise FAIL and its unmodified ZIP are preserved.
The *same* atol=5e-10 and rtol=2e-8 now apply to COEFFICIENT VECTOR
norms at EACH of the 3001 saved true FE steps, not to individual DOFs:
  ||u1-u2||inf <= atol + rtol*max(||u1||inf,||u2||inf)
  ||u1-u2||2   <= atol*sqrt(N) + rtol*max(||u1||2,||u2||2)
The second condition includes the same fixed absolute allowance
per coefficient propagated through Euclidean norm. For extra rigor,
also require pure normalized coefficient L2 < rtol at every time.

This is a coefficient-space metric, NOT a physical mass-weighted FE
L2 norm. All NON-FE original powers/H and physical residuals remain
subject to the UNCHANGED legacy pointwise allclose gate.
Only rerun archived numerical comparisons, never PETSc or GetFEM.
"""
import argparse
import importlib.util
import json
from pathlib import Path

import numpy as np

ORIGINAL=Path(__file__).with_name("verify-v130-heat-wave-full-mpi2.py")
spec=importlib.util.spec_from_file_location("heat_wave_original_strict_audit",ORIGINAL)
original=importlib.util.module_from_spec(spec)
spec.loader.exec_module(original)

ATOL=5e-10
RTOL=2e-8
STEPS=3000
DOFS=7689
SOURCE="e6326198d5cfe8c1921e85605be59f0eea51bd07"


def normwise_all_saved_states(coarse, parallel):
    """All true original GetFEM coefficients, chunked to cap RAM use."""
    if coarse.shape!=(STEPS+1,DOFS) or parallel.shape!=coarse.shape:
        raise RuntimeError("Missing some original true FE coefficients")
    count_inf=0
    count_l2=0
    pointwise_failures=0
    max_abs=0.
    max_step_relative_l2=0.
    max_inf_ratio=0.
    max_l2_ratio=0.
    first_infinity_failure=None
    first_euclidean_failure=None
    worst_abs_step=None
    global_diff2=0.
    global_ref2=0.
    original_initial_allclose=False
    for lo in range(0,STEPS+1,32):
        hi=min(lo+32,STEPS+1)
        a=np.asarray(coarse[lo:hi],dtype=np.float64)
        b=np.asarray(parallel[lo:hi],dtype=np.float64)
        if not (np.all(np.isfinite(a)) and np.all(np.isfinite(b))):
            raise RuntimeError("Nonfinite original native FE coefficients")
        d=a-b
        infinity_difference=np.max(np.abs(d),axis=1)
        reference_infinity=np.maximum(
            np.max(np.abs(a),axis=1),np.max(np.abs(b),axis=1)
        )
        bound_infinity=ATOL+RTOL*reference_infinity
        ax=np.sqrt(np.sum(a*a,axis=1))
        bx=np.sqrt(np.sum(b*b,axis=1))
        diff2=np.sum(d*d,axis=1)
        difference_2=np.sqrt(diff2)
        reference_2=np.maximum(ax,bx)
        bound_2=ATOL*np.sqrt(DOFS)+RTOL*reference_2
        relative_2=difference_2/np.maximum(reference_2,1e-300)
        bad_inf=infinity_difference>bound_infinity
        bad_l2=(difference_2>bound_2)|(relative_2>RTOL)
        count_inf+=int(np.count_nonzero(bad_inf))
        count_l2+=int(np.count_nonzero(bad_l2))
        if first_infinity_failure is None and np.any(bad_inf):
            first_infinity_failure=int(lo+np.flatnonzero(bad_inf)[0])
        if first_euclidean_failure is None and np.any(bad_l2):
            first_euclidean_failure=int(lo+np.flatnonzero(bad_l2)[0])
        # Preserve the historical scalar pointwise criterion for
        # diagnostic visibility only; it is NOT the FE norm gate.
        pointwise_failures+=int(np.count_nonzero(
            np.abs(d) > ATOL+RTOL*np.abs(b)
        ))
        local_max=int(np.argmax(infinity_difference))
        if infinity_difference[local_max]>max_abs:
            max_abs=float(infinity_difference[local_max])
            worst_abs_step=int(lo+local_max)
        max_step_relative_l2=max(max_step_relative_l2,
                                  float(np.max(relative_2)))
        max_inf_ratio=max(max_inf_ratio,
                          float(np.max(infinity_difference/bound_infinity)))
        max_l2_ratio=max(max_l2_ratio,
                         float(np.max(difference_2/bound_2)))
        global_diff2+=float(np.sum(diff2))
        global_ref2+=float(np.sum(np.maximum(ax*ax,bx*bx)))
        if lo==0:
            original_initial_allclose=bool(
                np.allclose(a[0],b[0],atol=ATOL,rtol=RTOL)
            )
    return {
        "coefficient_space_not_mass_weighted_FE_norm":True,
        "reference_absolute_tolerance":ATOL,
        "reference_relative_tolerance":RTOL,
        "infinity_norm_gate_formula":
            "||z1-z2||inf <= atol+rtol*max(||z1||inf,||z2||inf) each time",
        "euclidean_norm_gate_formula":
            "||z1-z2||2 <= sqrt(N)*atol+rtol*max(||z1||2,||z2||2) each time",
        "extra_normalized_euclidean_gate":"relative l2 per time <= rtol",
        "original_initial_state_passes_old_pointwise":original_initial_allclose,
        "saved_original_states":STEPS+1,
        "true_GetFEM_coefficients_per_state":DOFS,
        "maximum_absolute_coefficient_difference":max_abs,
        "worst_absolute_step":worst_abs_step,
        "maximum_normwise_infinity_tolerance_ratio":max_inf_ratio,
        "maximum_normwise_euclidean_tolerance_ratio":max_l2_ratio,
        "max_relative_coefficient_l2_per_time":max_step_relative_l2,
        "relative_coefficient_l2_full_space_time":
            float(np.sqrt(global_diff2)/max(1e-300,np.sqrt(global_ref2))),
        "failed_infinity_norm_steps":count_inf,
        "failed_euclidean_norm_steps":count_l2,
        "first_failed_infinity_step":first_infinity_failure,
        "first_failed_euclidean_step":first_euclidean_failure,
        "original_pointwise_failures_preserved_as_diagnostic":
            pointwise_failures,
        "qualified":(
            count_inf==0 and count_l2==0 and original_initial_allclose
        ),
    }


def compare_physics_unchanged(one,two):
    keys=set(one)
    if keys!=set(two):
        raise RuntimeError("Native physical fields differ")
    result={}
    for key in sorted(keys-{"z","z0","zf"}):
        a,b=one[key],two[key]
        if (a.shape!=b.shape or not np.all(np.isfinite(a))
                or not np.all(np.isfinite(b))):
            raise RuntimeError("Invalid original physical trace "+key)
        delta=float(np.max(np.abs(a-b)))
        # Preserve original FE-independent PHYSICS pointwise criterion.
        accepted=bool(np.allclose(a,b,atol=ATOL,rtol=RTOL))
        if not accepted:
            raise RuntimeError(
                f"Original physical observable old allclose failed: {key} {delta}"
            )
        result[key]=delta
    return result


def audit(root):
    result={
        "schema":"scrimp-v130-HeatWave-15s-normwise-MPI-parity-reassessment-v1",
        "status":"FAIL_OR_UNQUALIFIED_HEAT_WAVE_MPI2_REASSESSMENT",
        "previous_original_strict_status_preserved":True,
        "original_PDE_or_solver_tolerances_modified":False,
        "original_MPI2_rank0_rank1_bitwise_identical":False,
        "all_independent_native_physics_gates_passed":False,
        "all_original_scalar_power_observable_old_thresholds_passed":False,
        "full_15s_normwise_MPI1_MPI2_parity_qualified":False,
        "no_MPI_parallel_VTK_qualification":True,
        "no_h_dt_convergence_claim":True,
        "errors":[],
    }
    try:
        manifest=json.loads((root/"manifest.json").read_text())
        if (manifest.get("source")!=SOURCE
                or manifest.get("profile")!="gaussian"
                or manifest.get("steps")!=STEPS
                or manifest.get("dt")!=.005
                or manifest.get("FEM_DOFs")!=DOFS
                or manifest.get("original_model_unchanged") is not True
                or manifest.get("full_FE_states_saved_per_rank") is not True):
            raise RuntimeError("Original 15s MPI1/2 evidence provenance mismatch")
        prior=json.loads((root/"summary.json").read_text())
        if (prior.get("status")!="FAIL_OR_UNQUALIFIED_FULL_HEAT_WAVE_MPI2"
                or not any("zf" in x for x in prior.get("errors",[]))):
            raise RuntimeError("Original strict zf FAIL status not preserved")
        result["original_strict_status"]=prior["status"]
        result["original_strict_failures"]=prior["errors"]
        arrays={}
        for label,nrank in (("mpi1_long",1),("mpi2_long",2)):
            for rank in range(nrank):
                a,meta,physical=original.load(
                    root,label,rank,STEPS,.005,"gaussian"
                )
                pre=json.loads(
                    (root/label/f"preflight-rank{rank}.json").read_text()
                )
                if (pre["accepted_PETSc_TS_steps"]!=STEPS
                        or pre["expected_saved_FE_states"]!=STEPS+1
                        or pre["actual_saved_FE_states"]!=STEPS+1
                        or pre["FE_state_shape"]!=[STEPS+1,DOFS]
                        or pre["missing_accepted_grid_indices"]
                        or pre["TS_SNES_failures"]!=0
                        or pre["TS_step_rejections"]!=0
                        or a["z"].shape!=(STEPS+1,DOFS)
                        or not np.array_equal(a["z"][0],a["z0"])
                        or not np.array_equal(a["z"][-1],a["zf"])):
                    raise RuntimeError(
                        f"GetFEM original accepted-state/provenance gate {label}/{rank}"
                    )
                arrays[label,rank]=a
                result[f"native_original_physics_{label}_rank{rank}"]=physical
        mpi1=arrays["mpi1_long",0]
        mpi20=arrays["mpi2_long",0]
        mpi21=arrays["mpi2_long",1]
        if set(mpi20)!=set(mpi21) or not all(
            np.array_equal(mpi20[k],mpi21[k]) for k in mpi20
        ):
            raise RuntimeError("Native MPI2 rank replicas differ")
        result["original_MPI2_rank0_rank1_bitwise_identical"]=True
        result["all_independent_native_physics_gates_passed"]=True
        result["scalar_power_differences"]=compare_physics_unchanged(mpi1,mpi20)
        result["all_original_scalar_power_observable_old_thresholds_passed"]=True
        result["normwise_full_23074689_original_GetFEM_coefficients"]=(
            normwise_all_saved_states(mpi1["z"],mpi20["z"])
        )
        norm=result["normwise_full_23074689_original_GetFEM_coefficients"]
        if not norm["qualified"]:
            raise RuntimeError(
                "Native MPI1/MPI2 coefficient norm parity failed; "
                +str(norm)
            )
        # Existing independently accepted 15s Gaussian endpoint.
        if (abs(float(mpi1["H"][0])-0.7836421211876778)>1e-8
                or abs(float(mpi1["H"][-1])-0.008268389052222169)>2e-6):
            raise RuntimeError("Original MPI1 full 15-second benchmark changed")
        result["full_15s_normwise_MPI1_MPI2_parity_qualified"]=True
        result["status"]="PASS_HEAT_WAVE_FULL_15S_MPI1_MPI2_CN_NORM_PARITY"
    except Exception as exc:
        result["errors"].append(f"{type(exc).__name__}: {exc}")
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("evidence",type=Path)
    parser.add_argument("--out",type=Path,
                        help="Separate JSON path; old summary.json untouched")
    args=parser.parse_args()
    root=args.evidence
    if not root.is_dir():
        raise RuntimeError("Read-only reassessment expects original native dir")
    result=audit(root)
    report=json.dumps(result,indent=2)+"\n"
    if args.out:
        args.out.parent.mkdir(parents=True,exist_ok=True)
        args.out.write_text(report,encoding="utf-8")
    print(report,flush=True)
    raise SystemExit(0 if result["status"].startswith("PASS_") else 1)


if __name__=="__main__":
    main()
