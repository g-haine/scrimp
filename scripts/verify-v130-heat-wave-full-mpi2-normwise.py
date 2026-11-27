#!/usr/bin/env python3
"""Reusable full 15-second original Heat-Wave CN MPI1/2 release checker.

This is the accepted normwise coefficient-space criterion validated
offline against the complete native source-e6326198 dataset. Unlike
the historical *strict pointwise* verifier, it does NOT require a
previous failed summary.json, and can accept future genuine runs.

Uses the SAME original native metadata, per-rank physical checks and
true GetFEM quadratures from verify-v130-heat-wave-full-mpi2.py,
plus the tested coefficient-vector normwise comparator from
reassess-v130-heat-wave-full-mpi2-norm-parity.py.

Never edits any physical simulation, PETSc solver tolerances or
historical files unless an explicit --out path is provided.
"""
import argparse
import importlib.util
import json
from pathlib import Path

import numpy as np


def _import_file(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    result=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


DIR=Path(__file__).resolve().parent
previous=_import_file(
    "scrimp_v130_heat_wave_orig_physics_verify",
    DIR/"verify-v130-heat-wave-full-mpi2.py",
)
norm=_import_file(
    "scrimp_v130_heat_wave_normwise_metrics",
    DIR/"reassess-v130-heat-wave-full-mpi2-norm-parity.py",
)


def audit(root):
    """Fail closed on any absent/failing native physical or normwise data."""
    result={
        "schema":"scrimp-v130-HeatWave-15s-reusable-normwise-MPI-acceptance-v1",
        "status":"FAIL_OR_UNQUALIFIED_HEAT_WAVE_FULL_MPI_NORM_PARITY",
        "errors":[],
        "original_FE_model_unchanged":False,
        "original_MPI2_rank0_rank1_bitwise_identical":False,
        "all_original_MPI1_and_MPI2_native_CN_physics_gates_passed":False,
        "original_scalar_and_power_pointwise_parity_qualified":False,
        "full_15s_normwise_MPI1_MPI2_parity_qualified":False,
        "per_DOF_pointwise_parity_claimed":False,
        "MPI_parallel_VTK_writer_qualified":False,
        "mass_weighted_FE_L2_convergence_claimed":False,
    }
    try:
        manifest=json.loads((root/"manifest.json").read_text())
        if (
            not manifest.get("source")
            or manifest.get("image")!="scrimp-mpi:v1.3.0-rc3"
            or manifest.get("image_revision") !=
                "1d6d0c90115213298047b28237a7f38087563a35"
            or manifest.get("canonical_example")!="examples/heat_wave.py"
            or manifest.get("original_model_unchanged") is not True
            or manifest.get("full_FE_states_saved_per_rank") is not True
            or manifest.get("profile")!="gaussian"
            or manifest.get("dt")!=.005
            or manifest.get("steps")!=3000
            or manifest.get("FEM_DOFs")!=7689
        ):
            raise RuntimeError("Unqualified native Heat-Wave provenance")
        result["source"]=manifest["source"]
        result["native_image"]=manifest["image"]
        result["original_FE_model_unchanged"]=True
        arrays={}
        for label,ranks in (("mpi2_long",2),("mpi1_long",1)):
            for rank in range(ranks):
                data,meta,physics=previous.load(
                    root,label,rank,3000,.005,"gaussian"
                )
                pre=json.loads(
                    (root/label/f"preflight-rank{rank}.json").read_text()
                )
                if (
                    pre["accepted_PETSc_TS_steps"]!=3000
                    or pre["expected_saved_FE_states"]!=3001
                    or pre["actual_saved_FE_states"]!=3001
                    or pre["FE_state_shape"]!=[3001,7689]
                    or pre["missing_accepted_grid_indices"]
                    or pre["TS_SNES_failures"]!=0
                    or pre["TS_step_rejections"]!=0
                    or data["z"].shape!=(3001,7689)
                    or not np.array_equal(data["z0"],data["z"][0])
                    or not np.array_equal(data["zf"],data["z"][-1])
                ):
                    raise RuntimeError(
                        f"Missing original accepted GetFEM FE states {label}/{rank}"
                    )
                arrays[label,rank]=data
                result[f"native_{label}_rank{rank}_physics"]=physics
        x,y,z=(arrays["mpi1_long",0],arrays["mpi2_long",0],
               arrays["mpi2_long",1])
        if set(y)!=set(z) or not all(
            np.array_equal(y[field],z[field]) for field in y
        ):
            raise RuntimeError("Original MPI2 rank0/rank1 numerical records differ")
        result["original_MPI2_rank0_rank1_bitwise_identical"]=True
        result["all_original_MPI1_and_MPI2_native_CN_physics_gates_passed"]=True
        # Keep strict original *per scalar* physics power comparisons.
        result["original_scalar_and_power_differences"]=(
            norm.compare_physics_unchanged(x,y)
        )
        result["original_scalar_and_power_pointwise_parity_qualified"]=True
        result["all_saved_GetFEM_coefficient_vector_norms"]=(
            norm.normwise_all_saved_states(x["z"],y["z"])
        )
        errors=result["all_saved_GetFEM_coefficient_vector_norms"]
        if not errors["qualified"]:
            raise RuntimeError("Genuine full GetFEM normwise MPI parity failed")
        if (abs(float(x["H"][0])-0.7836421211876778)>1e-8
                or abs(float(x["H"][-1])-0.008268389052222169)>2e-6):
            raise RuntimeError("Unexpected native 15s Gaussian reference")
        result["full_15s_normwise_MPI1_MPI2_parity_qualified"]=True
        result["status"]="PASS_HEAT_WAVE_FULL_15S_MPI1_MPI2_CN_NORM_PARITY"
    except Exception as exc:
        result["errors"].append(f"{type(exc).__name__}: {exc}")
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("evidence",type=Path)
    parser.add_argument("--out",type=Path,
                        help="Write the result here; no default source mutation")
    args=parser.parse_args()
    report=audit(args.evidence)
    serialized=json.dumps(report,indent=2)+"\n"
    if args.out is not None:
        args.out.parent.mkdir(parents=True,exist_ok=True)
        args.out.write_text(serialized)
    print(serialized,flush=True)
    raise SystemExit(0 if report["status"].startswith("PASS_") else 1)


if __name__=="__main__":
    main()
