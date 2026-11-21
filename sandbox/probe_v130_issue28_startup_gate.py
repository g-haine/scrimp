# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See LICENSE for license information.

"""#28/#6: fail-closed opt-in initializer vs compatible rest-to-spinup t0.

Two DISTINCT diagnostic cases, both using the unchanged original
curved P2 viscous Saint-Venant bricks and native GetFEM IFunction:

- incompatible_rotation: run the ACTUAL DPHS.init_step dispatch with
  opt-in nonlinear_affine on original h0,p0 and already rotating wall.
  High-index null algebraic columns must cause a typed, collective,
  rollback-safe UNSUPPORTED rejection BEFORE any physical/pseudo TS.
  This structural rejection is not by itself a proof of physical
  inconsistency; independent exact block elimination already provided it.

- compatible_rest_spinup: physical omega=0.2, initial p0=0, original
  physical facet tangent wall ramped from rest with a C1 smoothstep.
  Build the analytically known t0 algebraic state e_h=rho*g*h,
  e_p=0, Y_n=-e_h, U_n=U=Y=0; verify native FULL IFunction
  at z_dot0=0, all original algebraic + dynamic blocks. Check
  actuation appears O(t^2) near t=0 rather than impulsively.
  This is a compatibility positive-control CANDIDATE, NOT a
  qualified production nonlinear/high-index initializer or TS run.

Never alter the originally prescribed rotating q0 in the negative
case. Never accept an unqualified Y0 or perform a physical TS step.
"""

import argparse
import json
from pathlib import Path
from unittest.mock import patch

import numpy as np
from mpi4py import MPI
from petsc4py import PETSc
import scrimp as S
from scrimp.utils.nonlinear_affine_dae import UnsupportedNonlinearAlgebraicDAE

from viscous_swe_2d_mpi import (
    build_system, H_INITIAL, RHO, GRAVITY, ROTATION,
)
from probe_v130_issue28_structure import (
    AuditComplete, indices, indexset, layout, norm_sub,
)
from probe_v130_issue28_weak_wall import ensure_real

SPINUP_DURATION = 0.02
RAMP_START = 0.03
RAMP_END = 0.06
DELTA = 0.0001


def _global_values(vec, mpi):
    return np.concatenate(mpi.allgather(
        ensure_real(vec.getArray(readonly=True), "distributed PETSc Vec")
        .copy()
    ))


def _all_reduced_exact(mpi, ok, description):
    if not mpi.allreduce(bool(ok), op=MPI.LAND):
        raise RuntimeError(f"Native startup gate changed {description}")


def _positive_rest_gate(system, mpi, mesh):
    names_d, names_a, iv, n = layout(system)
    if names_d != ("h", "p") or not all(
        k in iv for k in ("e_h", "e_p", "Y_n", "U_n", "U", "Y")
    ):
        raise RuntimeError("Unexpected original SWE algebraic layout")

    original = ensure_real(system.gf_model.from_variables(), "original q0")
    dynamic = indices(names_d, iv)
    if not np.all(original[iv["p"][0]:sum(iv["p"])] == 0.):
        raise RuntimeError("Spinup must start from zero original momentum")
    if not np.array_equal(
        original[iv["h"][0]:sum(iv["h"])],
        np.full(iv["h"][1], H_INITIAL),
    ):
        raise RuntimeError("Spinup must preserve prescribed constant h0")

    prepared = original.copy()
    eh, ch = iv["e_h"]
    yn, cy = iv["Y_n"]
    if ch == 0 or cy == 0:
        raise RuntimeError("Missing interior pressure or normal traction FEM")
    prepared[eh:eh+ch] = RHO*GRAVITY*H_INITIAL
    prepared[yn:yn+cy] = -RHO*GRAVITY*H_INITIAL
    _all_reduced_exact(mpi, np.array_equal(
        prepared[dynamic], original[dynamic]
    ), "original h0 or p0")

    lo, hi = system.tangent_stiffness.getOwnershipRange()
    created = []

    def keep(obj):
        created.append(obj)
        return obj

    snap_global = system.GlobalVec.copy()
    cache_ready = getattr(system, "_explicit_rhs_ready", False)
    cache_time = getattr(system, "_explicit_rhs_time", None)
    cache_rhs = getattr(system, "_explicit_rhs", None)
    cache_saved = cache_rhs.copy() if cache_rhs is not None else None
    try:
        z = keep(system.tangent_stiffness.createVecRight())
        if hi > lo:
            z.setValues(
                np.arange(lo, hi, dtype=PETSc.IntType), prepared[lo:hi]
            )
        z.assemble()
        zd = keep(z.duplicate())
        zd.set(0.)
        residual = keep(z.duplicate())

        system.IFunction(None, 0., z, zd, residual)
        full = _global_values(residual, mpi)
        block_norms = {
            name: float(np.linalg.norm(full[start:start+width]))
            for name, (start, width) in iv.items()
        }
        algebraic = indices(names_a, iv)
        fd = float(np.linalg.norm(full[dynamic]))
        fa = float(np.linalg.norm(full[algebraic]))
        total = float(np.linalg.norm(full))
        if any(not np.isfinite(x) for x in (fd, fa, total)):
            raise FloatingPointError("Native rest IFunction is not finite")
        if max(fd, fa, total) > 1e-10:
            raise RuntimeError(
                f"Original positive spinup rest not DAE-compatible at t0: "
                f"Fd={fd:.9g}, Fa={fa:.9g}, full={total:.9g}"
            )
        # At 0<t<ramp_start, the original decreasing envelope is one.
        # The rising C1 smoothstep is 3s^2-2s^3 (s=t/tau);
        # its derivative at zero is zero. Verify the ORIGINAL
        # time-dependent GetFEM port-U form responds accordingly.
        u_index = keep(indexset(
            indices(("U",), iv), lo, hi, z.getComm()
        ))
        def sample(t):
            system.IFunction(None, t, z, zd, residual)
            return norm_sub(residual, u_index)
        u_1 = sample(DELTA)
        u_2 = sample(2*DELTA)
        a = DELTA/SPINUP_DURATION
        smooth1 = 3*a*a-2*a*a*a
        smooth2 = 3*(2*a)**2-2*(2*a)**3
        expected = smooth2/smooth1
        observed = u_2/u_1 if u_1 else float("inf")
        if (u_1 <= 1e-15 or not np.isfinite(observed) or
                abs(observed/expected-1.) > 2e-6):
            raise RuntimeError(
                "Native wall actuation is not O(t²) with smoothstep "
                f"near t0: U(t)={u_1:.9g}, U(2t)={u_2:.9g}, "
                f"ratio={observed:.9g} expected={expected:.9g}"
            )

        # No physical state or initial algebraic variables are PUBLISHED:
        # this report qualifies only the t0 feasibility candidate.
        report = {
            "schema": "scrimp-v130-issue28-startup-gate-v1",
            "status": "COMPATIBLE_REST_SPINUP_T0_DIAGNOSTIC_ONLY",
            "nominal_mesh": mesh,
            "mpi_ranks": mpi.size, "rotation_after_spinup": ROTATION,
            "spinup_duration": SPINUP_DURATION,
            "geometry_order": 2, "height_order": 2,
            "original_h0": H_INITIAL, "original_p0_zero": True,
            "original_physical_wall_facet_tangent": True,
            "initial_wall_rotation_exact_zero": True,
            "initial_wall_rotation_derivative_zero": True,
            "initial_algebraic_coenergy_e_h": RHO*GRAVITY*H_INITIAL,
            "initial_algebraic_normal_traction_Y_n":
                -RHO*GRAVITY*H_INITIAL,
            "native_IFunction_t0_Fd_l2": fd,
            "native_IFunction_t0_Fa_l2": fa,
            "native_IFunction_t0_full_l2": total,
            "native_IFunction_t0_blocks_l2": block_norms,
            "native_F_U_at_delta_l2": float(u_1),
            "native_F_U_at_two_delta_l2": float(u_2),
            "native_wall_smoothstep_ratio": float(observed),
            "analytic_wall_smoothstep_ratio": float(expected),
            "original_differential_h0_p0_preserved_bitwise": True,
            "entire_GetFEM_and_GlobalVec_rollback_bitwise": True,
            "no_pseudo_or_physical_TS": True,
            "physical_time_evolution_qualified": False,
            "production_nonlinear_initializer_qualified": False,
            "physical_rotating_Y0_accepted": False,
            "interpretation": (
                "The deliberately distinct physical rest-to-rotation "
                "experiment has a complete native zero-derivative "
                "compatible t0 candidate; the original physical wall "
                "rotates only after t0, with zero initial derivative. "
                "No differential time step, differentiated constraint "
                "solve or production DAE initial state was accepted."
            ),
        }
        if mpi.rank == 0:
            print(
                f"ISSUE28 STARTUP POSITIVE mesh={mesh:g} MPI={mpi.size} "
                f"native_Fa={fa:.9g} native_Fd={fd:.9g} "
                f"U_ratio={observed:.9g} expected={expected:.9g}",
                flush=True,
            )
        return report
    finally:
        system.gf_model.set_time(0.)
        system.gf_model.to_variables(original.copy())
        snap_global.copy(system.GlobalVec)
        _all_reduced_exact(
            mpi, np.array_equal(
                _global_values(system.GlobalVec, mpi),
                _global_values(snap_global, mpi),
            ), "PETSc GlobalVec restored coefficients",
        )
        system._explicit_rhs_ready = cache_ready
        system._explicit_rhs_time = cache_time
        if cache_rhs is not None and cache_saved is not None:
            cache_saved.copy(cache_rhs)
        _all_reduced_exact(mpi, np.array_equal(
            ensure_real(system.gf_model.from_variables(), "restored GetFEM"),
            original,
        ), "GetFEM original coefficients")
        for obj in reversed(created):
            obj.destroy()
        snap_global.destroy()
        if cache_saved is not None:
            cache_saved.destroy()


def _negative_initializer_gate(system, mpi, mesh, original_init_step):
    original = ensure_real(system.gf_model.from_variables(), "negative q0")
    global_saved = system.GlobalVec.copy()
    old_method = str(system.time_scheme.getString(
        "init_step_ts_type", "algebraic"
    ))
    initial_ready = getattr(system, "_explicit_rhs_ready", False)
    initial_time = getattr(system, "_explicit_rhs_time", None)
    initial_rhs = getattr(system, "_explicit_rhs", None)
    initial_rhs_snapshot = (
        initial_rhs.copy() if initial_rhs is not None else None
    )

    try:
        # Exercise REAL production dispatch, but explicit OPT-IN only.
        # Default linear initializer and historical pseudo/BDF
        # are intentionally unchanged.
        system.time_scheme["init_step_ts_type"] = "nonlinear_affine"
        try:
            original_init_step(system)
        except UnsupportedNonlinearAlgebraicDAE as exc:
            reason = str(exc)
        else:
            raise RuntimeError(
                "High-index original rotating model unexpectedly accepted "
                "by production nonlinear_affine DAE initializer"
            )
        report = system.initialization_report
        if not isinstance(report, dict):
            raise RuntimeError("Native initializer left no rejection report")
        if (report.get("status") != "REJECTED" or
                report.get("high_index_boundary_multipliers_supported") is not False or
                report.get("no_pseudo_ts_executed") is not True or
                report.get("no_physical_ts_executed") is not True or
                not isinstance(report.get("algebraic_null_columns"), int) or
                report["algebraic_null_columns"] < 1 or
                "null columns" not in reason):
            raise RuntimeError(
                f"Production opt-in initializer rejected for an "
                f"unexpected reason: {reason}; report={report}"
            )
        _all_reduced_exact(mpi, np.array_equal(
            ensure_real(system.gf_model.from_variables(), "rejected q0"),
            original,
        ), "original rejected GetFEM state")
        before = _global_values(global_saved, mpi)
        after = _global_values(system.GlobalVec, mpi)
        _all_reduced_exact(
            mpi, np.array_equal(before, after), "rejected PETSc GlobalVec"
        )
        result = {
            "schema": "scrimp-v130-issue28-startup-gate-v1",
            "status": "HIGH_INDEX_OPT_IN_REJECTED_BEFORE_TS",
            "nominal_mesh": mesh,
            "mpi_ranks": mpi.size,
            "rotation": ROTATION,
            "geometry_order": 2, "height_order": 2,
            "original_rotating_q0_fixed": True,
            "original_wall_already_rotating_at_t0": True,
            "actual_production_init_step_dispatch": "nonlinear_affine",
            "expected_exception_type": "UnsupportedNonlinearAlgebraicDAE",
            "rejection_reason": reason,
            "native_algebraic_null_columns":
                report["algebraic_null_columns"],
            "production_initialization_report": report,
            "original_differential_h0_p0_preserved_bitwise": True,
            "entire_GetFEM_and_GlobalVec_rollback_bitwise": True,
            "no_pseudo_or_physical_TS": True,
            "physical_time_evolution_qualified": False,
            "production_nonlinear_initializer_qualified": False,
            "physical_rotating_Y0_accepted": False,
            "interpretation": (
                "The existing production opt-in nonlinear_affine path "
                "refuses a structurally singular higher-index rotating "
                "model with the expected typed exception, full initial "
                "state rollback, and no pseudo/BDF/physical TS. The "
                "separately qualified exact weak block gate, not "
                "this structural check, established the fixed-q0 "
                "physical incompatibility."
            ),
        }
        if mpi.rank == 0:
            print(
                f"ISSUE28 STARTUP NEGATIVE mesh={mesh:g} MPI={mpi.size} "
                f"null_columns={report['algebraic_null_columns']} "
                "typed_rejection=PASS rollback=PASS",
                flush=True,
            )
        return result
    finally:
        system.time_scheme["init_step_ts_type"] = old_method
        system.gf_model.set_time(0.)
        system.gf_model.to_variables(original.copy())
        # If production rejection was incorrect, preserve original
        # state regardless; do not let the diagnostic contaminate it.
        global_saved.copy(system.GlobalVec)
        system._explicit_rhs_ready = initial_ready
        system._explicit_rhs_time = initial_time
        if initial_rhs is not None and initial_rhs_snapshot is not None:
            initial_rhs_snapshot.copy(initial_rhs)
        global_saved.destroy()
        if initial_rhs_snapshot is not None:
            initial_rhs_snapshot.destroy()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--case", required=True,
        choices=("compatible_rest_spinup", "incompatible_rotation"),
    )
    parser.add_argument("--mesh", type=float, required=True, choices=(.2, .1))
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    mpi = PETSc.COMM_WORLD.tompi4py()
    positive = args.case == "compatible_rest_spinup"
    system = build_system(
        gamma=0., mesh=args.mesh, geometry_order=2, height_order=2,
        final_time=.08, save_dt=.002, max_step=.002,
        ramp_start=RAMP_START if positive else .003,
        ramp_end=RAMP_END if positive else .009,
        init_step=True, ts_type="bdf", adapt_type="none",
        init_dt=1e-9, rotation=ROTATION,
        initial_profile="rigid" if positive else "tangent",
        wall_velocity_profile="facet_tangent",
        algebraic_initial_guess="zero",
        start_from_rest_duration=SPINUP_DURATION if positive else 0.,
        init_method="bdf",
    )

    result = {}
    native_init = S.DPHS.init_step

    def intercept(self):
        if positive:
            result["report"] = _positive_rest_gate(self, mpi, args.mesh)
        else:
            result["report"] = _negative_initializer_gate(
                self, mpi, args.mesh, native_init
            )
        raise AuditComplete

    with patch.object(S.DPHS, "init_step", intercept):
        try:
            system.solve()
        except AuditComplete:
            pass
        else:
            raise RuntimeError("Native startup audit ran a physical time step")

    if "report" not in result:
        raise RuntimeError("Native startup diagnostic did not execute")
    if mpi.rank == 0:
        args.out.mkdir(parents=True, exist_ok=True)
        (args.out/"startup_gate.json").write_text(
            json.dumps(result["report"], indent=2)+"\n",
            encoding="utf-8",
        )
    mpi.Barrier()


if __name__ == "__main__":
    main()
