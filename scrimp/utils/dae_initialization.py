# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Opt-in, linear-only, fixed-differential-data DAE initialization.

Mathematics: hold all initial differential coefficients exactly fixed,
find only algebraic corrections using distributed PETSc K_aa and KSP
LSQR, and accept the result only if the *original* IFunction's
algebraic rows satisfy the requested tolerance. No pseudo-time TS.

This deliberately refuses lumped, nonlinear, and incomplete variable
layouts. Substituted FE ports are classified by their actual GetFEM
unknown variables. See docs/V130_DAE_ALGEBRAIC_INITIALIZER_DESIGN.md.
"""

import numpy as np
from mpi4py import MPI
from petsc4py import PETSc


class InconsistentDAEInitialData(RuntimeError):
    """No algebraically consistent initial state was found within tolerance."""


def _options(system):
    """Validate separate DAE feasibility and accuracy settings.

    An accuracy target tight enough for reproducible MPI trajectories
    must not silently change the actual feasibility acceptance criterion.
    """
    opts = system.time_scheme
    atol = float(opts.getString("init_algebraic_atol", "1e-10"))
    rtol = float(opts.getString("init_algebraic_rtol", "1e-8"))
    accuracy = float(opts.getString("init_algebraic_accuracy", "1e-11"))
    passes = int(opts.getString("init_algebraic_max_passes", "3"))
    ksp_rtol = float(opts.getString("init_algebraic_ksp_rtol", "1e-13"))
    ksp_atol = float(opts.getString("init_algebraic_ksp_atol", "1e-13"))
    for name, value in (
        ("init_algebraic_atol", atol),
        ("init_algebraic_rtol", rtol),
        ("init_algebraic_accuracy", accuracy),
        ("init_algebraic_ksp_rtol", ksp_rtol),
        ("init_algebraic_ksp_atol", ksp_atol),
    ):
        if not np.isfinite(value) or value < 0:
            raise ValueError(f"{name} must be finite and nonnegative")
    if passes < 1 or passes > 20:
        raise ValueError("init_algebraic_max_passes must be between 1 and 20")
    if not (ksp_rtol or ksp_atol):
        raise ValueError("At least one algebraic KSP tolerance must be > 0")
    return atol, rtol, accuracy, passes, ksp_rtol, ksp_atol


def _intervals(system, *, allow_nonlinear=False):
    """Classify every GetFEM unknown once and reject all unknown layouts.

    Nonlinear mass/stiffness is permitted ONLY when explicitly requested
    by the separate opt-in nonlinear-affine initializer. The existing
    nonlinear-affine initializer retains its substituted-port guard.
    """
    if not allow_nonlinear and (not system.linear_mass or not system.linear_stiffness):
        raise ValueError(
            "Algebraic initializer supports only linear mass and stiffness"
        )
    if system.lumped_variables:
        raise ValueError("Algebraic initializer does not support lumped ports")
    if allow_nonlinear and any(p.get_substituted() for p in system.ports.values()):
        raise ValueError(
            "Nonlinear-affine initializer does not yet support substituted ports"
        )

    n = int(system.gf_model.nbdof())
    dynamic, algebraic = {}, {}
    for state in system.states.values():
        name = state.get_name()
        dynamic[name] = tuple(map(
            int, system.gf_model.interval_of_variable(name)
        ))
    for port in system.ports.values():
        if port.get_algebraic():
            names = [port.get_flow()]
            if not port.get_substituted():
                names.append(port.get_effort())
        elif port.get_substituted():
            names = []
        else:
            names = [port.get_effort()]
        for name in names:
            if not name:
                continue
            if name in algebraic or name in dynamic:
                raise ValueError(
                    f"Repeated differential/algebraic interval {name}"
                )
            algebraic[name] = tuple(map(
                int, system.gf_model.interval_of_variable(name)
            ))
    if not dynamic or not algebraic:
        raise ValueError("Algebraic initializer requires both variable types")
    all_intervals = (*dynamic.values(), *algebraic.values())
    total = sum(width for _, width in all_intervals)
    if total != n:
        raise ValueError("Unclassified GetFEM DOFs in DAE initializer")
    coverage = np.zeros(n, dtype=np.int8)
    for start, width in all_intervals:
        if start < 0 or width <= 0 or start + width > n:
            raise ValueError("GetFEM variable interval out of bounds")
        coverage[start:start + width] += 1
    if np.any(coverage != 1):
        raise ValueError("DAE variable intervals overlap or leave gaps")
    aidx = np.sort(np.concatenate([
        np.arange(start, start + width, dtype=np.int64)
        for start, width in algebraic.values()
    ]))
    return n, dynamic, aidx


def _collective_error(mpi, local_error):
    """Prevent one-rank-only validation failures before PETSc collectives."""
    errors = mpi.allgather(local_error)
    if any(errors):
        raise ValueError(
            "DAE initialization structural error on MPI rank(s): "
            + "; ".join(
                f"{r}: {message}" for r, message in enumerate(errors)
                if message
            )
        )


def initialize_linear_algebraic(system):
    """Initialize compatible algebraic variables while preserving states.

    Must be called from DPHS.init_step *after* mass/stiffness/rhs/J/F
    assembly and *before* the main TS. Stores a serializable report in
    system.initialization_report, including on feasibility rejection.

    Returns the report on success. Errors are never silently converted
    into a pseudo-TS fallback.
    """
    mpi = system.stiffness.getComm().tompi4py()
    local_error = None
    try:
        atol, rtol, accuracy, max_passes, kr, ka = _options(system)
        n, dynamic, algebraic_indices = _intervals(system)
        initial_raw = np.asarray(system.gf_model.from_variables())
        if np.iscomplexobj(initial_raw):
            raise ValueError(
                "Algebraic initializer currently supports real states only"
            )
        initial = np.asarray(
            initial_raw, dtype=np.float64
        ).copy()
        if initial.shape != (n,) or not np.all(np.isfinite(initial)):
            raise ValueError("Invalid GetFEM real initial coefficient vector")
        lo, hi = map(int, system.stiffness.getOwnershipRange())
        if tuple(system.mass.getOwnershipRange()) != (lo, hi):
            raise ValueError("PETSc mass/stiffness ownership mismatch")
        owned = algebraic_indices[
            (algebraic_indices >= lo) & (algebraic_indices < hi)
        ]
        if len(owned) and (
            owned.min() < np.iinfo(PETSc.IntType).min
            or owned.max() > np.iinfo(PETSc.IntType).max
        ):
            raise ValueError("Algebraic index exceeds PetscInt width")
        ptr, cols, vals = system.mass.getValuesCSR()
        if len(ptr) != hi - lo + 1:
            raise ValueError("Wrong PETSc mass CSR row layout")
        local_mass_max = max((
            float(np.max(np.abs(
                vals[ptr[i-lo]:ptr[i-lo+1]]
            )))
            for i in owned if ptr[i-lo+1] > ptr[i-lo]
        ), default=0.0)
    except (ValueError, TypeError, AttributeError, IndexError) as exc:
        local_error = f"{type(exc).__name__}: {exc}"
    _collective_error(mpi, local_error)
    mass_max = float(mpi.allreduce(local_mass_max, op=MPI.MAX))
    if mass_max > 1e-12:
        raise ValueError(
            "Classified algebraic mass rows are not zero: "
            f"maximum abs entry {mass_max:.9g}"
        )

    t0 = float(system.time_scheme["t_0"])
    is_alg = None
    kaa = None
    ksp = None
    z = None
    zdot = None
    residual = None
    rhs = None
    delta = None
    accepted = False
    try:
        owned_idx = np.asarray(owned, dtype=PETSc.IntType)
        is_alg = PETSc.IS().createGeneral(
            owned_idx, comm=system.stiffness.getComm()
        )
        z = system.stiffness.createVecRight()
        if hi > lo:
            z.setValues(
                range(lo, hi), initial[lo:hi],
                addv=PETSc.InsertMode.INSERT_VALUES,
            )
        z.assemble()
        zdot = z.duplicate()
        zdot.set(0)
        residual = z.duplicate()
        residual.set(0)

        system.gf_model.set_time(t0)
        system.gf_model.to_variables(initial.copy())

        def original_residual():
            system.IFunction(None, t0, z, zdot, residual)
            sub = residual.getSubVector(is_alg)
            try:
                return float(sub.norm(PETSc.NormType.NORM_2))
            finally:
                residual.restoreSubVector(is_alg, sub)

        initial_norm = original_residual()
        tolerance = atol + rtol * initial_norm
        target = min(tolerance, accuracy)
        kaa = system.stiffness.createSubMatrix(is_alg, is_alg)
        rhs = kaa.createVecRight()
        delta = rhs.duplicate()
        a_lo, a_hi = rhs.getOwnershipRange()
        if a_hi - a_lo != len(owned):
            raise ValueError(
                "Distributed algebraic IS/submatrix ownership mismatch"
            )

        ksp = PETSc.KSP().create(comm=system.stiffness.getComm())
        ksp.setOptionsPrefix("scrimp_dae_init_")
        ksp.setOperators(kaa)
        ksp.setType(PETSc.KSP.Type.LSQR)
        ksp.getPC().setType(PETSc.PC.Type.NONE)
        ksp.setTolerances(
            rtol=kr, atol=ka,
            max_it=max(1000, 5 * len(algebraic_indices)),
        )
        ksp.setFromOptions()
        history = []
        best_state = initial.copy()
        best_norm = initial_norm
        prev_norm = initial_norm
        for iteration in range(1, max_passes + 1):
            entries = np.asarray(residual.getArray(readonly=True))
            rhs.getArray()[:] = -entries[owned - lo]
            rhs.assemble()
            delta.set(0)
            ksp.solve(rhs, delta)
            ksp_reason = int(ksp.getConvergedReason())
            ksp_iters = int(ksp.getIterationNumber())
            local_finite = bool(np.all(np.isfinite(
                delta.getArray(readonly=True)
            )))
            if not mpi.allreduce(local_finite, op=MPI.LAND):
                raise FloatingPointError("Non-finite PETSc LSQR correction")
            if len(owned):
                z.setValues(
                    owned_idx,
                    np.asarray(delta.getArray(readonly=True)).copy(),
                    addv=PETSc.InsertMode.ADD_VALUES,
                )
            z.assemble()
            system.Scatter[0].scatter(z, system.GlobalVec)
            candidate = np.asarray(
                system.GlobalVec.getArray(readonly=True).real,
                dtype=np.float64,
            ).copy()
            for name, (start, width) in dynamic.items():
                if not np.array_equal(
                    candidate[start:start+width],
                    initial[start:start+width],
                ):
                    raise RuntimeError(
                        f"Algebraic initializer modified prescribed {name}"
                    )
            system.gf_model.set_time(t0)
            system.gf_model.to_variables(candidate.copy())
            actual = original_residual()
            history.append({
                "pass": iteration,
                "ksp_reason": ksp_reason,
                "ksp_iterations": ksp_iters,
                "ksp_reported_residual_norm": float(
                    ksp.getResidualNorm()
                ),
                "algebraic_correction_l2": float(
                    delta.norm(PETSc.NormType.NORM_2)
                ),
                "native_algebraic_residual_l2": actual,
            })
            if not np.isfinite(actual):
                raise FloatingPointError("Non-finite native algebraic residual")
            if actual < best_norm:
                best_state = candidate.copy()
                best_norm = actual
            if actual <= target or ksp_reason <= 0:
                break
            if actual >= prev_norm * (1.0 - 1e-6):
                break
            prev_norm = actual

        # An unsuccessful final pass must not overwrite a better earlier
        # state, nor may a rejected result contaminate the model.
        if not np.array_equal(candidate, best_state):
            z.getArray()[:] = best_state[lo:hi]
            z.assemble()
            system.Scatter[0].scatter(z, system.GlobalVec)
            candidate = np.asarray(
                system.GlobalVec.getArray(readonly=True).real,
                dtype=np.float64,
            ).copy()
            system.gf_model.set_time(t0)
            system.gf_model.to_variables(candidate.copy())
        final_norm = original_residual()
        if abs(final_norm - best_norm) > max(1e-13, 1e-6*best_norm):
            raise RuntimeError(
                "Restored PETSc initial state no longer matches best residual"
            )

        report = {
            "strategy": "algebraic",
            "status": "PASS" if final_norm <= tolerance else "UNRESOLVED",
            "dofs": n,
            "algebraic_dofs": int(len(algebraic_indices)),
            "native_initial_residual_l2": initial_norm,
            "native_final_residual_l2": final_norm,
            "feasibility_tolerance": tolerance,
            "accuracy_target": target,
            "accuracy_target_met": bool(final_norm <= target),
            "mass_algebraic_max_abs": mass_max,
            "differential_initial_exact": True,
            "pseudo_ts_executed": False,
            "ksp_type": str(ksp.getType()),
            "ksp_options_prefix": ksp.getOptionsPrefix(),
            "passes": history,
        }
        system.initialization_report = report
        if final_norm > tolerance:
            raise InconsistentDAEInitialData(
                "Fixed prescribed dynamic initial data do not admit "
                "an algebraic state at the requested numerical tolerance "
                f"(||F_a||={final_norm:.9g}, tolerance={tolerance:.9g}); "
                "consider revising boundary/initial data or explicitly "
                "requesting historical pseudo/BDF initial relaxation"
            )
        accepted = True
        return report
    finally:
        if not accepted:
            # This path must NEVER publish rejected algebraic variables.
            system.gf_model.set_time(t0)
            system.gf_model.to_variables(initial.copy())
        for obj in (ksp, delta, rhs, kaa, is_alg, residual, zdot, z):
            if obj is not None:
                obj.destroy()
