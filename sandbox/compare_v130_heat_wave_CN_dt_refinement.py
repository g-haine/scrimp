# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Full-FE physical time-refinement diagnostic for canonical Heat-Wave CN.

Four genuine GetFEM/PETSc original-pH model simulations at the EXACT
SAME conforming spatial mesh h=.05 and t_f=.2, with:
  dt=0.01, 0.005, 0.0025, 0.00125 (20, 40, 80, 160 steps).
All rank0 true FE vectors are retained; no spatial interpolation,
ODE replacement, exact-solution fabrication or energy-only proxy.

For each nested triple, compare time-aligned original full FE coefficient
histories and REPORT the Richardson/self-convergence ratio without
asserting a given temporal order. In particular, this diagnostic is
NOT a validation of asymptotic second-order accuracy, spatial error,
analytic PDE convergence or an independent exact solution.

Each native physical gate and all original GetFEM subsystem power
balances must pass independently. No change to canonical physics.
"""
import argparse
import json
from pathlib import Path
import zipfile

import numpy as np


CASES=(
    ("dt_0p01",.01,20),
    ("dt_0p005",.005,40),
    ("dt_0p0025",.0025,80),
    ("dt_0p00125",.00125,160),
)
FE_DOF=7689
T_FINAL=.2


def finite(name,value):
    a=np.asarray(value,dtype=np.float64)
    if not np.all(np.isfinite(a)):
        raise ValueError("Nonfinite genuine GetFEM array: "+name)
    return a


def load_case(root,label,dt,steps):
    if (root/f"{label}.exitcode").read_text().strip()!="0":
        raise RuntimeError(f"Original native CN calculation failed {label}")
    folder=root/label
    meta=json.loads((folder/"rank0.json").read_text())
    pre=json.loads((folder/"preflight-rank0.json").read_text())
    if not (
        meta.get("schema")=="scrimp-v130-heat-wave-MPI-and-long-CN-physical-v1"
        and meta.get("status")=="PASS_NATIVE_STRICT_DAE_CN_ORIENTED_POWER"
        and meta.get("rank")==0
        and meta.get("MPI_ranks")==1
        and meta.get("DOFs")==FE_DOF
        and meta.get("steps")==steps
        and meta.get("profile")=="gaussian"
        and abs(meta.get("dt",99)-dt)<1e-14
        and abs(meta.get("t_final",99)-T_FINAL)<1e-10
        and meta.get("saved_all_accepted_steps") is True
        and all(meta.get("physical_gates",{}).values())
        and meta.get("all_rank_native_assembly_and_quadratures") is True
        and meta.get("initial_wave_boundary",{}).get("method")=="clamp"
        and meta["initial_wave_boundary"].get("boundary_compatibility_enforced")
        and meta["initial_wave_boundary"].get("after_boundary_l2_coefficients")==0.
        and pre.get("accepted_PETSc_TS_steps")==steps
        and pre.get("TS_SNES_failures")==0
        and pre.get("TS_step_rejections")==0
        and pre.get("actual_saved_FE_states")==steps+1
        and pre.get("expected_saved_FE_states")==steps+1
        and pre.get("FE_state_shape")==[steps+1,FE_DOF]
        and pre.get("missing_accepted_grid_indices")==[]
    ):
        raise RuntimeError(f"Original native GetFEM CN / physical gate failed: {label}")
    with np.load(folder/"rank0.npz",allow_pickle=False) as f:
        a={name:finite(name,f[name]) for name in f.files}
    if (a.get("z",np.array([])).shape!=(steps+1,FE_DOF)
            or a.get("t",np.array([])).shape!=(steps+1,)
            or a.get("H",np.array([])).shape!=(steps+1,)
            or not np.allclose(a["t"],np.arange(steps+1)*dt,
                               rtol=0,atol=1e-10)
            or not np.array_equal(a["z"][0],a["z0"])
            or not np.array_equal(a["z"][-1],a["zf"])):
        raise RuntimeError(f"Missing original consecutive full FE states: {label}")
    if (abs(float(a["H"][0])-0.7836421211876778)>1e-8
            or np.min(a["H"])<=0
            or np.max(np.diff(a["H"]))>1e-8):
        raise RuntimeError(f"Original coupled wave-heat initial energy invalid: {label}")
    for field in (
        "phys_CN_per_step_H_balance_error",
        "phys_CN_per_step_heat_balance_error",
        "phys_CN_per_step_wave_balance_error",
        "phys_CN_cumulative_total_balance_error",
        "phys_CN_P_interface_physical",
        "phys_CN_P_outer",
    ):
        if field not in a or a[field].shape!=(steps+1,):
            raise RuntimeError(f"Missing canonical independent physical power: {label}/{field}")
        if np.max(np.abs(a[field]))>1e-8:
            raise RuntimeError(f"Original CN energy/port power failure: {label}/{field}")
    return a,meta


def nested_triple(coarse,medium,fine,label):
    """Compute both true original full-FE differences on SAME time grid."""
    z0=coarse["z"]
    z1=medium["z"][::2]
    z2=fine["z"][::4]
    if z0.shape!=z1.shape or z0.shape!=z2.shape:
        raise RuntimeError("Real FE nested time-grid dimensions mismatch "+label)
    if not (np.allclose(coarse["t"],medium["t"][::2],rtol=0.,atol=1e-10)
            and np.allclose(coarse["t"],fine["t"][::4],rtol=0.,atol=1e-10)):
        raise RuntimeError("Original PETSc CN accepted times not nested "+label)
    if (np.max(np.abs(z0[0]-z1[0]))>1e-9
            or np.max(np.abs(z0[0]-z2[0]))>1e-9):
        raise RuntimeError("Original physical initial GetFEM state changed with dt "+label)
    e01=float(np.linalg.norm(z0-z1))
    e12=float(np.linalg.norm(z1-z2))
    if not np.isfinite(e01+e12) or min(e01,e12)<=1e-14:
        raise RuntimeError("Degenerate native CN refinement difference "+label)
    ratio=e01/e12
    result={
        "dt_coarse":float(coarse["t"][1]-coarse["t"][0]),
        "compared_true_saved_original_GetFEM_states":int(len(z0)),
        "real_FE_DOF_per_state":FE_DOF,
        "coarse_minus_half_full_coefficient_Frobenius":e01,
        "half_minus_quarter_full_coefficient_Frobenius":e12,
        "Richardson_ratio":ratio,
        "observed_order_log2_ratio":float(np.log2(ratio)),
        "coarse_half_max_abs_coefficient":float(np.max(np.abs(z0-z1))),
        "half_quarter_max_abs_coefficient":float(np.max(np.abs(z1-z2))),
        "spatial_mesh_kept_fixed":True,
        "analytic_PDE_exact_solution_not_claimed":True,
    }
    # Distinct from a solution error: this is self-convergence
    # of the original time integrator for the same FE semidiscretization.
    # An observed ratio in [3,5] NEVER proves order two. Stiff CN/DAE
    # transients may not lie in an asymptotic regime; this is descriptive.
    result["time_refinement_differences_decrease"] = bool(e12<e01)
    result["temporal_order_claimed_or_required"] = False
    return result


def audit(root):
    report={
        "schema":"scrimp-v130-original-heatwave-CN-temporal-refinement-v1",
        "status":"FAIL_OR_UNQUALIFIED_HEATWAVE_CN_PHYSICAL_REFINEMENT",
        "errors":[],
        "same_original_conforming_two_sided_GetFEM_problem":True,
        "same_original_mesh_h":.05,
        "fixed_physical_final_time":T_FINAL,
        "only_MUMPS_MPI1":True,
        "no_spatial_convergence_claim":True,
        "no_asymptotic_temporal_order_claim":True,
        "all_true_accepted_CN_FE_states_required":True,
        "no_synthetic_or_interpolated_FE_solution":True,
        "cases":{},
        "Richardson_self_convergence":{},
    }
    try:
        manifest=json.loads((root/"manifest.json").read_text())
        if not (
            manifest.get("original_model_unchanged") is True
            and manifest.get("profile")=="gaussian"
            and manifest.get("space_mesh_size_h")==.05
            and manifest.get("original_FE_DOF")==FE_DOF
            and manifest.get("final_time")==T_FINAL
            and manifest.get("solver")=="mumps"
            and manifest.get("MPI")==1
        ):
            raise RuntimeError("Native refinement manifest / provenance invalid")
        output={}
        for label,dt,steps in CASES:
            a,meta=load_case(root,label,dt,steps)
            output[label]=a
            report["cases"][label]={
                "dt":dt,"steps":steps,"saved_states":steps+1,
                "H0":float(a["H"][0]),"H_final":float(a["H"][-1]),
                "all_original_physical_gates":meta["physical_gates"],
                "max_CN_total_power_residual":float(np.max(
                    np.abs(a["phys_CN_per_step_H_balance_error"])
                )),
            }
        for i in (0,1):
            n0,n1,n2=(CASES[j][0] for j in (i,i+1,i+2))
            report["Richardson_self_convergence"][f"{n0}_over_{n1}_over_{n2}"]=(
                nested_triple(output[n0],output[n1],output[n2],n0)
            )
        report["status"]="PASS_HEATWAVE_CN_FIXED_MESH_PHYSICAL_REFINEMENT_DIAGNOSTIC"
    except Exception as e:
        report["errors"].append(f"{type(e).__name__}: {e}")
    return report


def main():
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument("evidence",type=Path)
    args=cli.parse_args()
    root=args.evidence
    result=audit(root)
    (root/"summary.json").write_text(json.dumps(result,indent=2)+"\n")
    archive=Path(str(root)+".zip")
    with zipfile.ZipFile(archive,"x",compression=zipfile.ZIP_DEFLATED,
                         compresslevel=4,allowZip64=True) as z:
        for file in sorted(root.rglob("*")):
            if file.is_file():
                z.write(file,str(file.relative_to(root)))
    with zipfile.ZipFile(archive) as z:
        bad=z.testzip()
        if bad is not None:
            raise RuntimeError("Corrupt original CN temporal-refinement evidence "+bad)
    print(json.dumps({**result,"evidence":str(root),"archive":str(archive),
                      "ZIP_CRC_pass":True},indent=2),flush=True)
    raise SystemExit(0 if result["status"].startswith("PASS_") else 1)


if __name__=="__main__":
    main()
