# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Time-quadrature conventions for already-assembled SCRIMP port powers.

Port.compute() now declares its sampling:
- cn_step: P[0] instantaneous, P[k>=1] at the original GetFEM
  CN mean state 0.5*(z[k-1]+z[k]) and model time t[k].
- instantaneous: all P[k] at (t[k], z[k]) for any PETSc scheme,
  suitable for approximate nodal quadrature.

This module integrates the raw, **signed** GetFEM flow-effort products.
It NEVER chooses the physical sign of a port, which depends on the
model's boundary orientation. It does not enforce energy conservation.
"""

import numpy as np


def integrate_port_power(times, samples, *, convention):
    """Return cumulative signed port-power integrals.

    Args:
        times: Strictly increasing accepted physical times (including t0).
        samples: GetFEM port powers at those stored indices, finite.
        convention: One of:
            'cn_step' -- P[k] is already the value for interval
                [t[k-1], t[k]] and is integrated as dt[k]*P[k].
                P[0] is INITIAL INSTANTANEOUS power, never used as
                the preceding step-power (critical for checkpoints).
            'instantaneous_trapezoid' -- integrate truly instantaneous
                nodal powers by an explicitly approximate trapezoid,
                with no claim to recover internal BDF/IMEX stages.
            'legacy_trapezoid' -- retained numerical alias for
                historical callers of this helper. The sampling
                convention must still be known by the caller.

    Returns:
        NumPy float array with W[0]=0 and one cumulative value per
        accepted physical time, with original port sign unchanged.
    """
    t = np.asarray(times, dtype=float)
    p = np.asarray(samples, dtype=float)
    if (t.ndim != 1 or p.ndim != 1 or len(t) != len(p)
            or len(t) < 1 or not np.all(np.isfinite(t))
            or not np.all(np.isfinite(p))):
        raise ValueError("Invalid physical time grid or port power samples")
    dt = np.diff(t)
    if not np.all(dt > 0.):
        raise ValueError("Expected strictly increasing accepted physical times")
    if convention == "cn_step":
        increments = dt * p[1:]
    elif convention in ("instantaneous_trapezoid", "legacy_trapezoid"):
        increments = .5 * dt * (p[:-1]+p[1:])
    else:
        raise ValueError("Unsupported port-power integration convention")
    if not np.all(np.isfinite(increments)):
        raise ValueError("Nonfinite integrated port powers")
    return np.r_[0., np.cumsum(increments)]
