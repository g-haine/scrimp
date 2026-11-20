# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Issue #28: production-path nonlinear_affine DAE initialization, NO TS.

Load the shipped original rectangular shallow-water model and native
operators exactly as in the independently passed read-only diagnostic.
Opt into SCRIMP's new production DPHS.init_step type, without running
any physical time step. Check the committed algebraic result, original
differential state preservation, MPI parity, and explicit G cache.
This is NOT a general nonlinear algebraic Newton or high-index gate.
"""

import argparse
import json
from pathlib import Path

import numpy as np
from mpi4py import MPI
from petsc4py import PETSc

from check_v130_shallow_water_mass_jacobian import _build_original


def _indices(names, sys):
    return np.sort(np.concatenate([
        np.arange(*(
            sys.gf_model.interval_of_variable(name)[0],
            sum(sys.gf_model.interval_of_variable(name))
        ), dtype=np.int64)
        for name in names
    ]))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    mpi = PETSc.COMM_WORLD.tompi4py()
    sys = _build_original(.1)
    if sys.linear_mass or sys.linear_stiffness:
        raise RuntimeError("The original shallow-water example is no longer nonlinear")
    before = np.asarray(sys.gf_model.from_variables(), dtype=np.float64).copy()
    dnames = tuple(state.get_name() for state in sys.states.values())
    if dnames != ("h", "p"):
        raise RuntimeError("Unexpected original dynamic field labels")
    d = _indices(dnames, sys)
    original_g = sys.bricks["G"].get_explicit()
    if not original_g:
        raise RuntimeError("The original gyroscopic G must remain explicit=True")

    # Use the *actual production initializer dispatch* and publish
    # accepted algebraic coefficients in the real GetFEM model.
    sys.time_scheme["init_step_ts_type"] = "nonlinear_affine"
    sys.init_step()
    report = dict(sys.initialization_report)
    after = np.asarray(sys.gf_model.from_variables(), dtype=np.float64).copy()
    dynamic_unchanged = bool(np.array_equal(before[d], after[d]))
    if not mpi.allreduce(dynamic_unchanged, op=MPI.LAND):
        raise RuntimeError("Production initializer modified prescribed initial h0,p0")

    if not sys._explicit_rhs_ready:
        raise RuntimeError("The initial explicit gyro was not frozen")

    # Re-evaluate actual production original IFunction at the
    # PUBLISHED initial state without a physical time step.
    z = sys.tangent_stiffness.createVecRight()
    zdot = z.duplicate()
    residual = z.duplicate()
    lo, hi = z.getOwnershipRange()
    if hi > lo:
        z.setValues(np.arange(lo, hi, dtype=PETSc.IntType), after[lo:hi])
    z.assemble()
    zdot.set(0.)
    try:
        sys.IFunction(None, 0., z, zdot, residual)
        names = list(dnames)
        for port in sys.ports.values():
            if port.get_algebraic():
                names.extend((port.get_flow(), port.get_effort()))
            else:
                names.append(port.get_effort())
        all_indices = np.concatenate([
            np.arange(sys.gf_model.interval_of_variable(name)[0],
                      sum(sys.gf_model.interval_of_variable(name)), dtype=np.int64)
            for name in names if name
        ])
        algebraic_indices = np.sort(np.setdiff1d(
            all_indices, d, assume_unique=True,
        ))
        owned_algebraic = algebraic_indices[
            (algebraic_indices >= lo) & (algebraic_indices < hi)
        ]
        ia = PETSc.IS().createGeneral(
            np.asarray(owned_algebraic, dtype=PETSc.IntType),
            comm=z.getComm(),
        )
        try:
            sub = residual.getSubVector(ia)
            try:
                verify_norm = float(sub.norm())
            finally:
                residual.restoreSubVector(ia, sub)
        finally:
            ia.destroy()
    finally:
        residual.destroy()
        zdot.destroy()
        z.destroy()

    if abs(verify_norm-report["native_final_algebraic_residual_l2"]) > 1e-11:
        raise RuntimeError("Native published algebraic state no longer matches report")
    if verify_norm > 1e-11:
        raise RuntimeError(
            f"Native production initializer did not meet strict target: {verify_norm:.9g}"
        )
    if report["status"] != "PASS":
        raise RuntimeError("Production initializer reported non-PASS status")
    if report["algebraic_null_columns"] != 0:
        raise RuntimeError("Unexpected null column in original shallow-water Kaa")

    result = {
        "schema": "scrimp-v130-issue28-production-nonlinear-affine-shallow-v1",
        "status": "PASS_NATIVE_PRODUCTION_INITIALIZER",
        "mpi_ranks": int(mpi.size),
        "model": "original shipped examples/shallow_water.py",
        "mesh": .1,
        "dofs": int(len(after)),
        "original_dynamic_h0_p0_bitwise": True,
        "published_algebraic_dofs_modified": bool(
            np.any(before != after)
        ),
        "actual_published_IFunction_Fa_l2": verify_norm,
        "real_production_init_step_report": report,
        "explicit_gyro_from_original_shallow_water": original_g,
        "explicit_gyro_cache_ready": bool(sys._explicit_rhs_ready),
        "main_physical_TS_executed": False,
        "pseudo_initial_TS_executed": False,
        "generic_nonlinear_algebraic_Newton_qualified": False,
        "high_index_rotating_Y_initialization_qualified": False,
    }
    if not result["published_algebraic_dofs_modified"]:
        raise RuntimeError("Production initializer failed to publish an algebraic correction")
    if mpi.rank == 0:
        args.out.mkdir(parents=True, exist_ok=True)
        (args.out/"production_nonlinear_affine.json").write_text(
            json.dumps(result, indent=2)+"\n", encoding="utf-8"
        )
        print(
            f"ISSUE28 PRODUCTION NONLINEAR_AFFINE MPI={mpi.size}: "
            f"PASS |Fa(initial)|={report['native_initial_algebraic_residual_l2']:.9g}, "
            f"|Fa(published)|={verify_norm:.9g}, "
            f"affine={report['affine_algebraic_identity_error_l2']:.9g}, "
            f"h0,p0 untouched, no physical TS",
            flush=True,
        )
    mpi.Barrier()


if __name__ == "__main__":
    main()
