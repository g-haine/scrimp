# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""#28 native NEGATIVE control: reject original high-index rotating traction Y.

Calls the REAL production opt-in nonlinear_affine DPHS.init_step()
through the original rotating viscous Saint-Venant system.solve().
This DAE has dF_a/dY = 0 but nonzero dF_d/dY, as established by
the independent MPI1/2 issue-28 structural audit. Hence naive square
algebraic Newton MUST reject the singular block before factorization.

Run on original nonzero rotation, true P2 GEOMETRIC disk, original
h0,p0 and the corrected clockwise GWFL wall control. No pseudo-BDF,
accepted physical TS, relaxed threshold or altered differential state.
Verify rollback including GetFEM algebraic variables and explicit cache.
"""

import argparse
import json
from pathlib import Path

import numpy as np
from mpi4py import MPI
from petsc4py import PETSc

from scrimp.utils.nonlinear_affine_dae import (
    UnsupportedNonlinearAlgebraicDAE,
)
from viscous_swe_2d_mpi import build_system, ROTATION


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    mpi = PETSc.COMM_WORLD.tompi4py()
    system = build_system(
        gamma=0.0,
        mesh=0.2,
        geometry_order=2,
        height_order=2,
        final_time=0.02,
        save_dt=0.002,
        max_step=0.002,
        ramp_start=0.003,
        ramp_end=0.009,
        init_step=True,
        ts_type="bdf",
        adapt_type="none",
        init_dt=1e-9,
        rotation=ROTATION,
        initial_profile="tangent",
        wall_velocity_profile="facet_tangent",
        algebraic_initial_guess="zero",
        init_method="bdf",
    )
    n = int(system.gf_model.nbdof())
    original = np.asarray(
        system.gf_model.from_variables(), dtype=np.float64
    ).copy()
    ys, yn = map(int, system.gf_model.interval_of_variable("Y"))
    assert yn > 0 and np.isfinite(original).all()

    # Override only the initial strategy, NOT the model or wall law.
    system.time_scheme["init_step_ts_type"] = "nonlinear_affine"
    if system.time_scheme.getString("init_step_ts_type") != "nonlinear_affine":
        raise RuntimeError("Production opt-in not applied")

    caught = None
    try:
        system.solve()
    except UnsupportedNonlinearAlgebraicDAE as exc:
        caught = str(exc)
    else:
        raise RuntimeError(
            "FAIL: singular rotating Y was accepted by nonlinear_affine!"
        )

    report = getattr(system, "initialization_report", None)
    if not isinstance(report, dict):
        raise RuntimeError("No native initializer rejection report")
    if report.get("strategy") != "nonlinear_affine":
        raise RuntimeError("Wrong production initializer was exercised")
    if report.get("status") != "REJECTED":
        raise RuntimeError("Production initializer was not marked REJECTED")
    if "null columns" not in (caught or ""):
        raise RuntimeError(f"Expected algebraic null-column rejection, got: {caught}")
    null_columns = int(report.get("algebraic_null_columns", -1))
    if null_columns < yn:
        raise RuntimeError(
            f"Expected at least {yn} Y null columns; found {null_columns}"
        )
    after = np.asarray(
        system.gf_model.from_variables(), dtype=np.float64
    ).copy()
    model_rolled_back = bool(np.array_equal(after, original))
    if not mpi.allreduce(model_rolled_back, op=MPI.LAND):
        raise RuntimeError("Rejected initializer changed original GetFEM state")
    if getattr(system, "_explicit_rhs_ready", False):
        raise RuntimeError("Rejected initial guess left frozen explicit forcing")
    if report.get("no_pseudo_ts_executed") is not True:
        raise RuntimeError("Initializer incorrectly entered historical pseudo-TS")
    if report.get("no_physical_ts_executed") is not True:
        raise RuntimeError("Initializer incorrectly entered physical TS")
    if n != len(original):
        raise RuntimeError("Wrong native model size")
    result = {
        "schema": "scrimp-v130-issue28-negative-high-index-v1",
        "status": "PASS_EXPECTED_REJECTION",
        "mpi_ranks": int(mpi.size),
        "geometry_order": 2,
        "omega": float(ROTATION),
        "original_fixed_h0_p0_preserved": True,
        "full_original_GetFEM_model_rolled_back": True,
        "algebraic_Y_dofs": int(yn),
        "algebraic_Y_start": int(ys),
        "detected_algebraic_null_columns": null_columns,
        "expected_exception": "UnsupportedNonlinearAlgebraicDAE",
        "exception_message": caught,
        "production_initialization_report": report,
        "physical_TS_executed": False,
        "pseudo_TS_executed": False,
        "nonlinear_affine_success_claimed": False,
        "rotating_high_index_Y0_qualified": False,
    }
    if mpi.rank == 0:
        args.out.mkdir(parents=True, exist_ok=True)
        (args.out / "negative_high_index.json").write_text(
            json.dumps(result, indent=2) + "\n", encoding="utf-8"
        )
        print(
            f"ISSUE28 NEGATIVE HIGH INDEX MPI={mpi.size}: "
            f"PASS_EXPECTED_REJECTION null_columns={null_columns}, "
            f"Y_dofs={yn}, original h0,p0 and algebraic state rolled back",
            flush=True,
        )
    mpi.Barrier()


if __name__ == "__main__":
    main()
