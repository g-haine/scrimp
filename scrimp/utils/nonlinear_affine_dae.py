# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Opt-in nonlinear-state / *affine-algebraic* consistent initialization.

The nonlinear full DPHS may have state-dependent mass and nonlinear
dynamic bricks. For fixed original differential q0 we accept only an
algebraic block whose Jacobian is regular and whose true original
IFunction satisfies the affine identity to numerical precision.

One distributed PETSc direct Kaa solve, no pseudo-time TS. The original
differential coefficients are preserved BITWISE. An incompatible or
singular algebraic block is rejected with *complete state rollback*.
This does NOT handle high-index boundary multipliers (such as Y in the
rotating viscous disk) requiring differentiated constraints.

No user data or physical constraints are silently projected.
"""

import numpy as np
from mpi4py import MPI
from petsc4py import PETSc

from scrimp.utils.dae_initialization import (
    _options, _intervals, InconsistentDAEInitialData,
)


class UnsupportedNonlinearAlgebraicDAE(ValueError):
    """This DAE requires more than a regular affine algebraic solve."""


def initialize_nonlinear_affine_algebraic(system):
    """Initialize an opt-in regular affine algebraic block at fixed q0.

    Called only after the ordinary DPHS native mass/stiffness/rhs
    assembly. Publishes a candidate iff native IFunction, affine
    identity, fixed prescribed states and direct KSP acceptance all
    succeed. Leaves an initialization report on success or rejection.
    """

    mpi = system.stiffness.getComm().tompi4py()
    atol, rtol, accuracy, _, _, _ = _options(system)
    n, dynamic, a = _intervals(system, allow_nonlinear=True)
    if system.linear_mass and system.linear_stiffness:
        raise UnsupportedNonlinearAlgebraicDAE(
            "Use init_step_ts_type=algebraic for fully linear DPHS"
        )

    # The reference is only the user's initially prescribed GetFEM
    # coefficients, before any nonlinear solves or BDF pseudo stepping.
    initial_raw = np.asarray(system.gf_model.from_variables())
    if np.iscomplexobj(initial_raw):
        raise ValueError("Real-valued nonlinear affine initializer only")
    original = np.asarray(initial_raw, dtype=np.float64).copy()
    if original.shape != (n,) or not np.all(np.isfinite(original)):
        raise ValueError("Nonfinite initial GetFEM coefficients")

    d = np.concatenate([
        np.arange(start, start+width, dtype=np.int64)
        for start, width in dynamic.values()
    ])
    d.sort()
    if len(a) + len(d) != n:
        raise ValueError("Nonlinear DAE unknown partition is incomplete")
    lo, hi = system.tangent_stiffness.getOwnershipRange()
    if system.tangent_mass.getOwnershipRange() != (lo, hi):
        raise ValueError("Nonlinear DAE mass/stiffness ownership mismatch")
    t0 = float(system.time_scheme["t_0"])
    objects = []
    accepted = False
    # A rejected candidate must not remain in either the GetFEM model
    # OR the PETSc sequential GlobalVec used by later SCRIMP routines.
    # The candidate is scattered to GlobalVec *before* the final native
    # residual/accuracy acceptance checks; save the pretrial vector.
    global_snapshot = system.GlobalVec.copy()
    report = {
        "strategy": "nonlinear_affine",
        "status": "REJECTED",
        "dofs": n,
        "differential_dofs": int(len(d)),
        "algebraic_dofs": int(len(a)),
        "differential_initial_exact": False,
        "no_pseudo_ts_executed": True,
        "no_physical_ts_executed": True,
        "full_nonlinear_algebraic_newton_qualified": False,
        "high_index_boundary_multipliers_supported": False,
    }
    system.initialization_report = report

    # Preserve and restore any previous cached accepted-state explicit
    # forcing on failure; otherwise freeze it at the ORIGINAL initial
    # state, never at an algebraic Newton trial.
    explicit_old_ready = getattr(system, "_explicit_rhs_ready", False)
    explicit_old_time = getattr(system, "_explicit_rhs_time", None)
    explicit_old = getattr(system, "_explicit_rhs", None)
    explicit_snapshot = None
    if explicit_old is not None:
        explicit_snapshot = explicit_old.copy()

    def keep(obj):
        objects.append(obj)
        return obj

    try:
        z = keep(system.tangent_stiffness.createVecRight())
        if hi > lo:
            z.setValues(
                np.arange(lo, hi, dtype=PETSc.IntType),
                original[lo:hi],
                addv=PETSc.InsertMode.INSERT_VALUES,
            )
        z.assemble()
        zd = keep(z.duplicate())
        zd.set(0.)
        residual = keep(z.duplicate())
        ai = a[(a >= lo) & (a < hi)]
        ia = keep(PETSc.IS().createGeneral(
            np.asarray(ai, dtype=PETSc.IntType),
            comm=z.getComm(),
        ))
        if len(ai):
            local_ptr, _, entries = system.tangent_mass.getValuesCSR()
            local_max = max(
                (float(np.max(np.abs(entries[local_ptr[i-lo]:local_ptr[i-lo+1]])))
                 for i in ai if local_ptr[i-lo+1] > local_ptr[i-lo]),
                default=0.,
            )
        else:
            local_max = 0.
        mass_max = float(mpi.allreduce(local_max, op=MPI.MAX))
        report["max_abs_algebraic_mass_rows"] = mass_max
        if mass_max > 1e-12:
            raise UnsupportedNonlinearAlgebraicDAE(
                "Algebraic mass rows are nonzero: "
                f"max={mass_max:.9g}"
            )

        if any(b.get_explicit() for b in system.bricks.values()):
            system._freeze_explicit_rhs(t0, z)
            if not system._explicit_rhs_ready:
                raise RuntimeError("Explicit force not frozen at original q0")
        system.IFunction(None, t0, z, zd, residual)
        r0 = residual.getSubVector(ia)
        try:
            norm0 = float(r0.norm())
        finally:
            residual.restoreSubVector(ia, r0)
        tolerance = float(atol + rtol * norm0)
        target = float(min(tolerance, accuracy))
        report.update({
            "native_initial_algebraic_residual_l2": norm0,
            "feasibility_tolerance": tolerance,
            "accuracy_target": target,
        })

        # Assemble the actual nonlinear tangent at the ORIGINAL q0.
        # Only its algebraic-algebraic part is used: no state update.
        system._sync_lumped_trial_state(t0, z)
        if not system.linear_stiffness:
            system.assemble_nl_stiffness()
        source = (
            system.tangent_stiffness if not system.linear_stiffness
            else system.stiffness
        )
        A = keep(source.createSubMatrix(ia, ia))
        if A.getSize() != (len(a), len(a)):
            raise UnsupportedNonlinearAlgebraicDAE(
                "Algebraic Jacobian block is not square"
            )

        # Detect the simplest high-index obstruction *before* MUMPS.
        # A zero algebraic column represents a multiplier only coupled
        # to dynamic equations. The rotating viscous disk has this
        # exact structure; a full square Jaa Newton is forbidden.
        _, col, val = A.getValuesCSR()
        local_col_max = np.zeros(len(a), dtype=np.float64)
        if len(val):
            np.maximum.at(local_col_max, np.asarray(col), np.abs(val))
        col_max = np.zeros(len(a), dtype=np.float64)
        mpi.Allreduce(local_col_max, col_max, op=MPI.MAX)
        reference = float(col_max.max(initial=0.))
        null_cols = np.flatnonzero(col_max <= max(1e-14, 1e-12*reference))
        report["algebraic_null_columns"] = int(len(null_cols))
        if len(null_cols):
            raise UnsupportedNonlinearAlgebraicDAE(
                f"Algebraic Jacobian has {len(null_cols)} numerically "
                "null columns; a differentiated-constraint / higher-index "
                "boundary initializer may be required"
            )

        rhs = keep(A.createVecLeft())
        sub = residual.getSubVector(ia)
        try:
            rhs.getArray()[:] = -np.asarray(sub.getArray(readonly=True))
        finally:
            residual.restoreSubVector(ia, sub)
        rhs.assemble()
        delta = keep(A.createVecRight())
        delta.set(0.)

        ksp = keep(PETSc.KSP().create(comm=A.getComm()))
        ksp.setOptionsPrefix("scrimp_dae_nonlinear_affine_")
        ksp.setOperators(A)
        ksp.setType(PETSc.KSP.Type.PREONLY)
        ksp.getPC().setType(PETSc.PC.Type.LU)
        ksp.getPC().setFactorSolverType("mumps")
        ksp.solve(rhs, delta)
        ksp_reason = int(ksp.getConvergedReason())
        report.update({
            "linear_solver": "PETSc PREONLY + MUMPS",
            "linear_solver_converged_reason": ksp_reason,
            "algebraic_correction_l2": float(delta.norm()),
        })
        if ksp_reason <= 0:
            raise InconsistentDAEInitialData(
                "Algebraic Kaa factorization/solve did not converge; "
                f"PETSc reason={ksp_reason}"
            )
        da = np.asarray(delta.getArray(readonly=True)).copy()
        if not bool(mpi.allreduce(
            np.all(np.isfinite(da)), op=MPI.LAND
        )):
            raise FloatingPointError("Nonfinite algebraic correction")
        if len(ai) != len(da):
            raise RuntimeError("PETSc local algebraic block ownership mismatch")
        if len(ai):
            z.setValues(
                np.asarray(ai, dtype=PETSc.IntType), da,
                addv=PETSc.InsertMode.ADD_VALUES,
            )
        z.assemble()
        system.Scatter[0].scatter(z, system.GlobalVec)
        all_candidate = np.asarray(
            system.GlobalVec.getArray(readonly=True).real,
            dtype=np.float64,
        ).copy()
        exact = bool(np.array_equal(all_candidate[d], original[d]))
        if not mpi.allreduce(exact, op=MPI.LAND):
            raise RuntimeError("Nonlinear algebraic solve changed prescribed q0")

        # The real nonlinear production IFunction, not a proxy or the
        # linearized residual, decides acceptance.
        system.IFunction(None, t0, z, zd, residual)
        rr = residual.getSubVector(ia)
        try:
            norm_final = float(rr.norm())
            predicted = keep(A.createVecLeft())
            A.mult(delta, predicted)
            predicted.axpy(-1., rhs)
            predicted.axpy(-1., rr)
            affine_error = float(predicted.norm())
        finally:
            residual.restoreSubVector(ia, rr)
        report.update({
            "native_final_algebraic_residual_l2": norm_final,
            "affine_algebraic_identity_error_l2": affine_error,
            "differential_initial_exact": exact,
            "accuracy_target_met": bool(norm_final <= target),
        })
        affine_limit = 1e-10 * max(1., norm0)
        report["affine_identity_tolerance"] = affine_limit
        if not np.isfinite(norm_final) or not np.isfinite(affine_error):
            raise FloatingPointError("Nonfinite native nonlinear DAE residual")
        if affine_error > affine_limit:
            raise UnsupportedNonlinearAlgebraicDAE(
                "Algebraic residual not affine at fixed original "
                f"differential data: error={affine_error:.9g}"
            )
        # Accuracy target is intentionally stricter than feasibility.
        # This mirrors the existing linear initializer's separate
        # feasibility / accuracy policy but here both are required.
        if norm_final > target:
            raise InconsistentDAEInitialData(
                "Nonlinear affine algebraic solve did not reach strict "
                f"accuracy: native |Fa|={norm_final:.9g}, "
                f"target={target:.9g}"
            )

        # Publish only accepted algebraic coefficients, never h0,p0
        # changes or an intermediate Newton iterate.
        system.gf_model.set_time(t0)
        system.gf_model.to_variables(all_candidate.copy())
        report["status"] = "PASS"
        accepted = True
        return report
    except Exception as exc:
        report["status"] = "REJECTED"
        report["rejection_reason"] = f"{type(exc).__name__}: {exc}"
        raise
    finally:
        if not accepted:
            # Restore both representations, even when final validation
            # rejects a *computed* (but unaccepted) algebraic candidate.
            global_snapshot.copy(system.GlobalVec)
            system.gf_model.set_time(t0)
            system.gf_model.to_variables(original.copy())
            system._explicit_rhs_ready = explicit_old_ready
            system._explicit_rhs_time = explicit_old_time
            if explicit_old is not None and explicit_snapshot is not None:
                explicit_snapshot.copy(explicit_old)
            elif getattr(system, "_explicit_rhs", None) is not None:
                system._explicit_rhs_ready = False
        for obj in reversed(objects):
            obj.destroy()
        global_snapshot.destroy()
        if explicit_snapshot is not None:
            explicit_snapshot.destroy()
