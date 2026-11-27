"""Post-mortem of a completed, but not yet accepted, 15s Heat-Wave run.

Read only ORIGINAL archived native GetFEM vectors and power traces.
Do not solve another PDE, modify a tolerance, or label MPI2 as PASS.

Usage:
    python sandbox/diagnose_v130_heat_wave_full_mpi_parity.py \
        sandbox/outputs/v130-heat-wave-full-15s-MPI1-MPI2-... \
        --out /tmp/heat-wave-MPI-parity-diagnostic.json

The input is the original evidence *directory* (preferred to re-reading
a potentially huge double-compressed ZIP). This program never changes
the original archive, recorded results or their failed acceptance status.
"""
import argparse
import json
from pathlib import Path

import numpy as np


ABS_REFERENCE=5e-10
REL_REFERENCE=2e-8
STATE_SHAPE=(3001,7689)


def as_finite(array,name):
    x=np.asarray(array,dtype=np.float64)
    if not np.all(np.isfinite(x)):
        raise RuntimeError("Nonfinite native evidence: "+name)
    return x


def read_npz(folder,name):
    f=folder/name
    if not f.is_file():
        raise FileNotFoundError(str(f))
    with np.load(f,allow_pickle=False) as data:
        return {key:as_finite(data[key],key) for key in data.files}


def state_mismatch_summary(coarse,parallel,time):
    """Chunked full space-time audit, no invented FE transformations."""
    if coarse.shape!=STATE_SHAPE or parallel.shape!=STATE_SHAPE:
        raise RuntimeError("Expected all 3001 genuine 7689-DOF snapshots")
    worst=(-1.,None,None)
    max_scaled=(-1.,None,None)
    violation_count=0
    near_zero_violations=0
    cumulative_diff_sq=0.
    cumulative_ref_sq=0.
    cumulative_other_sq=0.
    steps_max_abs=np.zeros(coarse.shape[0])
    steps_rel_l2=np.zeros(coarse.shape[0])
    steps_violations=np.zeros(coarse.shape[0],dtype=np.int64)
    for lo in range(0,len(time),64):
        hi=min(lo+64,len(time))
        x=coarse[lo:hi]
        y=parallel[lo:hi]
        difference=np.abs(x-y)
        threshold=ABS_REFERENCE+REL_REFERENCE*np.abs(y)
        threshold=np.maximum(
            threshold,
            ABS_REFERENCE+REL_REFERENCE*np.abs(x),
        )
        exceeds=difference>threshold
        norms_x=np.sum(x*x,axis=1)
        norms_y=np.sum(y*y,axis=1)
        norms_d=np.sum((x-y)*(x-y),axis=1)
        cumulative_diff_sq+=float(np.sum(norms_d))
        cumulative_ref_sq+=float(np.sum(norms_x))
        cumulative_other_sq+=float(np.sum(norms_y))
        steps_max_abs[lo:hi]=np.max(difference,axis=1)
        steps_rel_l2[lo:hi]=(
            np.sqrt(norms_d) /
            np.maximum(1e-300,np.maximum(np.sqrt(norms_x),np.sqrt(norms_y)))
        )
        steps_violations[lo:hi]=np.count_nonzero(exceeds,axis=1)
        violation_count+=int(np.count_nonzero(exceeds))
        near_zero_violations+=int(np.count_nonzero(
            exceeds & (np.maximum(np.abs(x),np.abs(y))<1e-3)
        ))
        position=np.unravel_index(int(np.argmax(difference)),difference.shape)
        val=float(difference[position])
        if val>worst[0]:worst=(val,int(lo+position[0]),int(position[1]))
        scaled=difference/threshold
        position_scaled=np.unravel_index(
            int(np.argmax(scaled)),scaled.shape
        )
        factor=float(scaled[position_scaled])
        if factor>max_scaled[0]:
            max_scaled=(factor,int(lo+position_scaled[0]),int(position_scaled[1]))

    v,i,j=worst
    factor,k,m=max_scaled
    final_difference=coarse[-1]-parallel[-1]
    top=np.argsort(np.abs(final_difference))[-12:][::-1]
    return {
        "maximum_absolute_state_difference":v,
        "worst_step":i,"worst_time":float(time[i]),"worst_DOF":j,
        "original_MPI1_value_at_worst":float(coarse[i,j]),
        "original_MPI2_value_at_worst":float(parallel[i,j]),
        "tolerance_at_worst":float(
            ABS_REFERENCE+REL_REFERENCE*max(abs(coarse[i,j]),abs(parallel[i,j]))
        ),
        "maximum_pointwise_tolerance_ratio":factor,
        "worst_scaled_step":k,"worst_scaled_time":float(time[k]),
        "worst_scaled_DOF":m,
        "absolute_difference_at_scaled_worst":
            float(abs(coarse[k,m]-parallel[k,m])),
        "original_MPI1_value_at_scaled_worst":float(coarse[k,m]),
        "original_MPI2_value_at_scaled_worst":float(parallel[k,m]),
        "pointwise_failed_comparisons":violation_count,
        "pointwise_compared_coefficients":int(coarse.size),
        "failed_comparisons_with_both_values_below_1e-3":
            near_zero_violations,
        "relative_l2_over_all_saved_FE_states":float(
            np.sqrt(cumulative_diff_sq) /
            max(1e-300,np.sqrt(max(cumulative_ref_sq,cumulative_other_sq)))
        ),
        "maximum_relative_l2_per_saved_state":
            float(np.max(steps_rel_l2)),
        "time_of_max_relative_l2":
            float(time[int(np.argmax(steps_rel_l2))]),
        "first_failed_step":(
            int(np.flatnonzero(steps_violations)[0])
            if violation_count else None
        ),
        "last_failed_step":(
            int(np.flatnonzero(steps_violations)[-1])
            if violation_count else None
        ),
        "sampled_time_evolution":[
            {"step":n,
             "t":float(time[n]),
             "max_abs":float(steps_max_abs[n]),
             "relative_l2":float(steps_rel_l2[n]),
             "pointwise_failures":int(steps_violations[n])}
            for n in (0,1,4,20,100,200,500,1000,1500,2000,2500,3000)
            if n<len(time)
        ],
        "t_final_largest_original_DOF_differences":[
            {"DOF":int(j),
             "MPI1":float(coarse[-1,j]),
             "MPI2":float(parallel[-1,j]),
             "difference":float(final_difference[j])}
            for j in top
        ],
    }


def scalar_trace_comparison(a,b):
    out={}
    for key in sorted(a):
        if key not in b:
            raise RuntimeError("MPI1/MPI2 field absent: "+key)
        if key=="z":continue
        x=a[key];y=b[key]
        if x.shape!=y.shape:
            raise RuntimeError("Different physical trace shapes: "+key)
        diff=np.abs(x-y)
        index=int(np.argmax(diff))
        out[key]={
            "maximum_absolute_difference":float(diff.flat[index]),
            "maximum_relative_l2":float(
                np.linalg.norm((x-y).ravel())/
                max(1e-300,
                    np.linalg.norm(x.ravel()),np.linalg.norm(y.ravel()))
            ),
            "pointwise_old_allclose":bool(np.allclose(
                x,y,atol=ABS_REFERENCE,rtol=REL_REFERENCE
            )),
            "max_MPI1_magnitude":float(np.max(np.abs(x))),
            "max_MPI2_magnitude":float(np.max(np.abs(y))),
        }
    return out


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("evidence",type=Path,
                        help="Existing original evidence directory, not a new run")
    parser.add_argument("--out",type=Path,
                        help="Optional separate JSON path; original archive unchanged")
    args=parser.parse_args()
    root=args.evidence
    if not root.is_dir():
        raise RuntimeError("Specify existing original evidence directory (not ZIP)")
    manifest=json.loads((root/"manifest.json").read_text())
    if (manifest.get("profile")!="gaussian"
            or manifest.get("dt")!=.005
            or manifest.get("steps")!=3000
            or manifest.get("FEM_DOFs")!=7689):
        raise RuntimeError("Not the original full 15-second Heat-Wave comparison")
    paths=["mpi1_long","mpi2_long"]
    for name in paths:
        if (root/(name+".exitcode")).read_text().strip()!="0":
            raise RuntimeError("Not a completed native run: "+name)
    one=read_npz(root/"mpi1_long","rank0.npz")
    two=read_npz(root/"mpi2_long","rank0.npz")
    replica=read_npz(root/"mpi2_long","rank1.npz")
    if set(one)!=set(two) or set(two)!=set(replica):
        raise RuntimeError("Missing original native observable arrays")
    if one["z"].shape!=STATE_SHAPE or two["z"].shape!=STATE_SHAPE:
        raise RuntimeError("Incomplete native FE history")
    if not all(np.array_equal(two[key],replica[key]) for key in two):
        raise RuntimeError("MPI2 ranks are NOT bitwise identical")
    t=one["t"]
    if not (t.shape==(3001,)
            and np.allclose(t,np.arange(3001)*.005,atol=1e-10,rtol=0.)):
        raise RuntimeError("Original 3001-step clock changed")
    if not np.allclose(t,two["t"],atol=1e-10,rtol=0.):
        raise RuntimeError("MPI1/MPI2 native time clocks differ")
    z_report=state_mismatch_summary(one["z"],two["z"],t)
    scalars=scalar_trace_comparison(one,two)
    summary={
        "schema":"scrimp-v130-heat-wave-15s-original-MPI-parity-DIAGNOSTIC-v1",
        "status":"POST_MORTEM_ONLY_NOT_A_RELEASE_ACCEPTANCE",
        "original_native_source":manifest["source"],
        "original_evidence_directory":str(root),
        "original_report_unchanged":True,
        "physical_model_or_tolerances_modified":False,
        "MPI2_rank0_rank1_exact_identical":True,
        "absolute_comparison_atol":ABS_REFERENCE,
        "relative_comparison_rtol":REL_REFERENCE,
        "state":z_report,
        "physical_and_scalar_traces":scalars,
        "all_original_mpi_physical_observable_traces_pass_old_allclose":
            all(v["pointwise_old_allclose"] for v in scalars.values()),
    }
    text=json.dumps(summary,indent=2)+"\n"
    if args.out:
        args.out.parent.mkdir(parents=True,exist_ok=True)
        args.out.write_text(text,encoding="utf-8")
    print(text)


if __name__=="__main__":
    main()
