# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""#28 negative control: reject a computed algebraic trial below numeric accuracy.

Build exactly the shipped rectangular examples/shallow_water.py initial
DAE (h0,p0 and physical boundary conditions UNMODIFIED) and opt into
production nonlinear_affine. Set *only* an unrealistically strict
algebraic numerical accuracy 1e-16, below the known ~2e-12 achievable
floating-point residual. The correction must be computed, rejected,
and rolled back in BOTH GetFEM AND the PETSc GlobalVec plus cached
explicit RHS. No pseudo-time or physical TS.

This deliberately tests numerical acceptance/rollback, NOT physical
incompatibility of initial differential data.
"""

import argparse
import json
from pathlib import Path

import numpy as np
from mpi4py import MPI
from petsc4py import PETSc

from scrimp.utils.dae_initialization import InconsistentDAEInitialData
from check_v130_shallow_water_mass_jacobian import _build_original


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    mpi = PETSc.COMM_WORLD.tompi4py()

    system = _build_original(0.1)
    before = np.asarray(system.gf_model.from_variables(), dtype=float).copy()
    global_before = np.asarray(system.GlobalVec.getArray(readonly=True)).copy()
    old_ready = bool(getattr(system, "_explicit_rhs_ready", False))
    old_time = getattr(system, "_explicit_rhs_time", None)
    old_rhs = getattr(system, "_explicit_rhs", None)
    old_rhs_contents = (
        np.asarray(old_rhs.getArray(readonly=True)).copy()
        if old_rhs is not None else None
    )

    system.time_scheme["init_step_ts_type"] = "nonlinear_affine"
    system.time_scheme["init_algebraic_accuracy"] = "1e-16"
    caught = None
    try:
        system.init_step()
    except InconsistentDAEInitialData as exc:
        caught = str(exc)
    else:
        raise AssertionError(
            "FAIL: precision 1e-16 should reject this original "
            "shallow-water finite-precision algebraic trial"
        )

    report = system.initialization_report
    if report.get("status") != "REJECTED":
        raise RuntimeError("Rejected nonlinear affine trial not recorded")
    if report.get("strategy") != "nonlinear_affine":
        raise RuntimeError("Wrong production path")
    if "did not reach strict accuracy" not in caught:
        raise RuntimeError(f"Wrong failure reason: {caught}")
    residual = float(report["native_final_algebraic_residual_l2"])
    target = float(report["accuracy_target"])
    if not (0 < target <= 1e-16 and residual > target):
        raise RuntimeError(
            f"Negative precision control not exercised: residual={residual}, "
            f"target={target}"
        )
    if int(report.get("linear_solver_converged_reason", 0)) <= 0:
        raise RuntimeError("Algebraic correction was never computed")
    if int(report.get("algebraic_null_columns", -1)) != 0:
        raise RuntimeError("Incorrect structural rejection")
    after = np.asarray(system.gf_model.from_variables(), dtype=float)
    global_after = np.asarray(system.GlobalVec.getArray(readonly=True))
    getfem_unchanged = bool(np.array_equal(before, after))
    petsc_unchanged = bool(np.array_equal(global_before, global_after))
    cache_ready_unchanged = (
        bool(getattr(system, "_explicit_rhs_ready", False)) == old_ready
    )
    cache_time_unchanged = getattr(system, "_explicit_rhs_time", None) == old_time
    new_rhs = getattr(system, "_explicit_rhs", None)
    rhs_contents_unchanged = (
        True if old_rhs_contents is None else
        new_rhs is not None and np.array_equal(
            np.asarray(new_rhs.getArray(readonly=True)),
            old_rhs_contents,
        )
    )
    checks = (getfem_unchanged, petsc_unchanged, cache_ready_unchanged,
              cache_time_unchanged, rhs_contents_unchanged)
    if not mpi.allreduce(all(checks), op=MPI.LAND):
        raise RuntimeError(
            "FAIL rollback after *computed* rejected algebraic candidate: "
            f"getfem={getfem_unchanged}, GlobalVec={petsc_unchanged}, "
            f"ready={cache_ready_unchanged}, time={cache_time_unchanged}, "
            f"explicit_contents={rhs_contents_unchanged}"
        )
    result = {
        "schema": "scrimp-v130-issue28-negative-accuracy-v1",
        "status": "PASS_EXPECTED_NUMERICAL_REJECTION",
        "mpi_ranks": int(mpi.size),
        "model": "original examples/shallow_water.py",
        "physical_initial_h0_p0_modified": False,
        "physical_boundary_laws_modified": False,
        "numerical_accuracy_target": target,
        "computed_candidate_Fa_l2": residual,
        "algebraic_correction_computed": True,
        "getfem_bitwise_rolled_back": True,
        "petsc_GlobalVec_bitwise_rolled_back": True,
        "explicit_cache_ready_restored": True,
        "explicit_cache_time_restored": True,
        "explicit_cache_previous_contents_preserved": True,
        "expected_exception": "InconsistentDAEInitialData",
        "exception_message": caught,
        "production_initialization_report": report,
        "pseudo_TS_executed": False,
        "main_TS_executed": False,
        "physically_incompatible_initial_data_demonstrated": False,
        "rotating_high_index_initialization_qualified": False,
    }
    if mpi.rank == 0:
        args.out.mkdir(parents=True, exist_ok=True)
        (args.out / "negative_accuracy.json").write_text(
            json.dumps(result, indent=2) + "\n", encoding="utf-8"
        )
        print(
            f"ISSUE28 NEGATIVE ACCURACY MPI={mpi.size}: "
            f"PASS_EXPECTED_NUMERICAL_REJECTION "
            f"computed |Fa|={residual:.9g} > target={target:.9g}; "
            "GetFEM+PETSc+cache rolled back",
            flush=True,
        )
    mpi.Barrier()


if __name__ == "__main__":
    main()
