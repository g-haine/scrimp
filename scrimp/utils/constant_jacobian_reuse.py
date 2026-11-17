# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Opt-in reuse of LU factors for a CERTIFIED constant implicit Jacobian.

For a time-invariant linear system with fixed-step Crank--Nicolson,
J = K + (2/dt) M is constant throughout the physical TS. The
consistency/initialization TS remains entirely unchanged.

This module deliberately uses no PETSc import: its admission gate can
be unit-tested separately. It is NOT sufficient for arbitrary linear
systems with time-varying coefficients, event-driven changes or varying
time step. The user certifies these conditions when opting in.
"""
import math


def enable_constant_jacobian_lu_reuse(
    ts, *, all_bricks_linear, linear_mass, linear_stiffness,
    lumped_variables, mumps_policy, initial_time, final_time,
    configured_dt,
):
    """Ask PETSc KSP to retain the FIRST full LU factorization.

    Must be called after TS.setFromOptions() and before TS.solve(),
    on all MPI ranks. The independent initialization TS was already
    destroyed and is never affected by this optimization.

    Return descriptive configuration details for diagnostics.
    """
    if (not all_bricks_linear or not linear_mass
            or not linear_stiffness or bool(lumped_variables)):
        raise ValueError(
            "Constant-Jacobian LU reuse requires a fully linear "
            "PDE without lumped variables. Time-dependent matrix "
            "coefficients are NEVER permitted in this opt-in mode."
        )
    if mumps_policy != "never":
        raise ValueError(
            "Constant-Jacobian LU reuse requires "
            "SCRIMP_MUMPS_RESET_POLICY=never"
        )
    ts_type = ts.getType()
    if ts_type != "cn":
        raise ValueError(
            "Constant-Jacobian LU reuse currently supports only "
            "fixed-step Crank-Nicolson (ts_type=cn)"
        )
    adapt_type = ts.getAdapt().getType()
    if adapt_type != "none":
        raise ValueError(
            "Constant-Jacobian LU reuse requires "
            "ts_adapt_type=none"
        )
    t0 = float(initial_time)
    tf = float(final_time)
    dt = float(configured_dt)
    if not (math.isfinite(t0) and math.isfinite(tf)
            and math.isfinite(dt) and tf > t0 and dt > 0):
        raise ValueError("Require finite positive fixed CN time interval")
    steps = round((tf - t0) / dt)
    if (steps < 1 or not math.isclose(
        tf - t0, steps * dt, rel_tol=0., abs_tol=1e-11 * max(1., abs(tf))
    )):
        raise ValueError(
            "A partial final CN step changes the TS shift and "
            "invalidates constant-Jacobian LU reuse"
        )
    if not math.isclose(
        float(ts.getTimeStep()), dt, rel_tol=1e-12, abs_tol=1e-14
    ):
        raise ValueError(
            "Effective TS time step differs from configured CN dt"
        )
    ksp = ts.getSNES().getKSP()
    pc = ksp.getPC()
    if (
        ksp.getType() != "preonly"
        or pc.getType() != "lu"
        or pc.getFactorSolverType() != "mumps"
    ):
        raise ValueError(
            "Constant-Jacobian reuse is initially supported only "
            "for KSPPREONLY + PCLU + MUMPS"
        )
    # This preserves the first PC construction, then omits further
    # PC setup even if SNES reissues the SAME numerical Jacobian.
    # Crucially it does NOT suppress the first factorization.
    ksp.setReusePreconditioner(True)
    return {
        "ts_type": ts_type,
        "adapt_type": adapt_type,
        "steps": int(steps),
        "dt": dt,
        "ksp_type": ksp.getType(),
        "pc_type": pc.getType(),
        "factor_solver": pc.getFactorSolverType(),
        "method": "KSPSetReusePreconditioner",
        "certification": (
            "User expressly certifies a time-invariant linear "
            "operator, constant coefficients, no TS events that "
            "modify matrices, and no adaptive time stepping."
        ),
    }
