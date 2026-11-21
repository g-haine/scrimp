# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Recheck archived #28 original rotating P2 physical-RHS/Y0 diagnostics.

Pure Python host-only, no Docker, PETSc or GetFEM; the source native
pilots are deliberately read-only, and a Y0 candidate is NOT an accepted
consistent initial traction without independently satisfied Fa.
"""

import argparse
import json
import math
from pathlib import Path

EXPECTED_IMAGE = "1d6d0c90115213298047b28237a7f38087563a35"


def require(value, reason):
    if not value:
        raise ValueError(reason)


def finite(value, name):
    require(isinstance(value, (float, int)) and not isinstance(value, bool)
            and math.isfinite(value), f"Nonfinite {name}: {value!r}")
    return float(value)


def compare(evidence):
    folder = Path(evidence)
    require(folder.is_dir(), "Missing saved host evidence directory")
    manifest = json.loads((folder/"manifest.json").read_text(encoding="utf-8"))
    require(manifest.get("schema") == "v130-issue28-physical-rhs-Y0-pilot-v1",
            "Wrong native campaign schema")
    require(manifest.get("image_revision") == EXPECTED_IMAGE,
            "Wrong immutable Docker revision")
    require(manifest.get("physical_Y0_accepted") is False
            and manifest.get("no_physical_time_steps") is True,
            "This gate may not claim physical Y0 initialization")
    summary = {}
    for case in ("rest", "rotating"):
        runs = []
        for ranks in (1, 2):
            path = folder/f"result_{case}_{ranks}"/"physical_rhs_y0_pilot.json"
            require(path.is_file(), f"Missing native report {path}")
            r = json.loads(path.read_text(encoding="utf-8"))
            require(r.get("schema") == "scrimp-v130-issue28-physical-RHS-Y0-pilot-v1",
                    f"Native schema wrong {case}/{ranks}")
            require(r.get("status") == "DIAGNOSTIC_TRUE_RHS_SOLVED_NOT_INITIALIZED",
                    f"Physical Y0 improperly accepted {case}/{ranks}")
            require(r.get("case") == case and r.get("mpi_ranks") == ranks,
                    f"Case/rank mismatch {case}/{ranks}")
            require(r.get("dofs") == 3300
                    and r.get("dynamic_dofs") == 1452
                    and r.get("regular_algebraic_dofs") == 1716
                    and r.get("traction_Y_dofs") == 132,
                    f"Original FE block size changed {case}/{ranks}")
            require(r.get("original_P2_geometry") is True
                    and r.get("original_q0_preserved_bitwise") is True,
                    f"Original P2 geometry or physical initial q0 changed {case}/{ranks}")
            require(r.get("bordered_KSP_reason", 0) > 0,
                    f"Physical-RHS bordered solver failed {case}/{ranks}")
            require(finite(r.get("bordered_relative_true_rhs_residual"), "linear residual")
                    < 1e-8,
                    f"Original bordered physical RHS not numerically solved {case}/{ranks}")
            require(finite(r.get(
                "original_native_IFunction_Fd_with_candidate_derivative_Y0_l2"
            ), "original dynamic residual") < 1e-8,
                    f"Actual original GetFEM dynamic residual too high {case}/{ranks}")
            require(finite(r.get(
                "original_native_IFunction_Fa_with_candidate_Y0_l2"
            ), "original algebraic residual") >= 0.,
                    "Invalid original algebraic residual")
            require(finite(r.get("candidate_traction_Y0_l2"), "Y0 norm") >= 0.,
                    "Invalid candidate traction norm")
            require(r.get("independent_rotating_algebraic_compatibility_qualified") is False
                    and r.get("physical_rotating_Y0_accepted") is False
                    and r.get("mesh_independent_physical_constraint_check_qualified") is False,
                    "Unproven physical data compatibility/initialization advertised")
            require(r.get("true_initialization_TS_executed") is False
                    and r.get("true_physical_TS_executed") is False,
                    "Native time stepping occurred in matrix pilot")
            p=r.get("initial_projection",{})
            require(finite(p.get("native_final_Fa_l2"), "post-projection algebraic Fa")>=0.,
                    "Projection norm invalid")
            require(p.get("exact_discrete_rotating_compatibility_certified") is False,
                    "Incompatible physical projection certification")
            if case=="rest":
                require(p["native_final_Fa_l2"] < 1e-11,
                        "Known exact stationary algebraic solution not reproduced")
                require(r["candidate_traction_Y0_l2"] < 1e-7,
                        "Zero-rest physical traction did not reproduce zero")
            runs.append(r)
        a,b=runs
        fa1=a["initial_projection"]["native_final_Fa_l2"]
        fa2=b["initial_projection"]["native_final_Fa_l2"]
        require(math.isclose(fa1,fa2,rel_tol=1e-3,abs_tol=1e-8),
                f"MPI algebraic projection mismatch {case}: {fa1},{fa2}")
        Y1=a["candidate_traction_Y0_l2"]
        Y2=b["candidate_traction_Y0_l2"]
        # Record Y0 rank sensitivity, but DO NOT accept a physical
        # rotating Y0 based on rank parity at an incompatible state.
        parity=abs(Y1-Y2)/max(1.,abs(Y1),abs(Y2))
        summary[case]={
            "MPI1_projected_native_Fa_l2":fa1,
            "MPI2_projected_native_Fa_l2":fa2,
            "MPI1_candidate_Y0_l2":Y1,
            "MPI2_candidate_Y0_l2":Y2,
            "candidate_Y0_norm_MPI_relative_difference":parity,
            "candidate_Y0_norm_MPI_parity_below_1e-3":bool(parity <= 1e-3),
            "physical_rotating_compatibility_established":False,
        }
        print(
            f"{case}: native physical RHS diagnostic MPI1/2 complete; "
            f"Fa=({fa1:.9g},{fa2:.9g}), "
            f"candidate |Y0|=({Y1:.9g},{Y2:.9g}), "
            f"Y0 norm relative MPI difference={parity:.9g}",
            flush=True,
        )
    result={
        "schema":"scrimp-v130-issue28-physical-rhs-Y0-comparator-v1",
        "status":"COMPLETE_DIAGNOSTIC_NOT_INITIALIZER",
        "cases":summary,
        "rotating_physical_initial_Y0_qualified":False,
        "actual_rotating_constraint_compatibility_qualified":False,
        "no_main_or_pseudo_TS":True,
    }
    (folder/"physical_rhs_y0_comparison.json").write_text(
        json.dumps(result,indent=2)+"\n",encoding="utf-8"
    )
    print("COMPLETE: true physical RHS pilot is diagnostic only; rotating Y0 NOT initialized.")
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--evidence",required=True,type=Path)
    args=parser.parse_args()
    compare(args.evidence)


if __name__=="__main__":
    main()
