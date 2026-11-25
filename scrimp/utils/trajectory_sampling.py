# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Roundoff-safe PETSc TS output scheduling (not time integration).

PETSc's accepted step at t=1.49 can be represented very slightly
below the algebraic threshold last_saved_t + dt_save. An exact
comparison then silently drops a genuine accepted state. If the
postprocessor subsequently interprets each pair of SAVED states
as a consecutive CN step, a physically false power balance follows.

Only forgive a small ULP-scale ROUNDING error, additionally capped
to one-millionth of dt_save, so the policy cannot fill a genuinely
missing sampling time or interpolate an uncomputed state.
"""

import math


def saved_time_is_due(time, scheduled_time, dt_save):
    """Whether an actual accepted TS sample reaches its next save target.

    Never manufactures/interpolates a time or FE state.
    Accepts only a limited 16-ULP float comparison error, capped
    to 1e-6 of the positive requested save spacing. This remains
    safe both when dt_save=dt and when deliberate subsampling is used.
    """
    t, target, spacing = (
        float(time), float(scheduled_time), float(dt_save)
    )
    if (not all(map(math.isfinite, (t, target, spacing)))
            or spacing <= 0.0):
        raise ValueError(
            "Accepted TS sampling needs finite times and positive dt_save"
        )
    magnitude = max(1.0, abs(t), abs(target))
    margin = min(16.0 * math.ulp(magnitude), 1.e-6 * spacing)
    return target - t <= margin


def fixed_cn_accepted_step_is_due(
    step, *, dt, dt_save, ts_type, ts_adapt_type,
):
    """Save EACH genuine accepted TS step in an explicitly fixed CN run.

    When dt_save == dt and adaptive stepping is disabled, PETSc's
    monitor is called only after accepted TS steps. Rechecking the
    accumulated (previous saved t + dt_save) threshold can skip a
    perfectly legitimate accepted CN state due to long-horizon
    floating-point accumulation (seen at step 2999, t≈14.995,
    in the 3000-step Heat–Wave MUMPS pilot).

    This routine changes ONLY trajectory sampling, not physical time,
    the TS step size, FE solution, CN algebra or DAE tolerances.
    Keep saved_time_is_due() for adaptive schemes, subsampling, and
    any other time scheme. The caller's physical audit must STILL
    verify strict real step count and every recorded timestamp.
    """
    delta=float(dt)
    save=float(dt_save)
    if not (math.isfinite(delta) and math.isfinite(save)
            and delta>0. and save>0.):
        raise ValueError("Fixed-CN sample spacing must be positive and finite")
    return bool(
        int(step)>0
        and str(ts_type).lower()=="cn"
        and str(ts_adapt_type).lower()=="none"
        and delta==save
    )
