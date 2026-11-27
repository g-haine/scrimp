# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

# SPDX-License-Identifier: GPL-3.0-or-later
"""Resolve numerically honest port-power sampling from the actual TS history.

CN step powers are evaluated using a mean of TWO CONSECUTIVE accepted
states by Port.compute. They must never be reconstructed from two
subsampled snapshots: dt_save > dt is a common user configuration.
Non-CN schemes use instantaneous nodal powers; their trapezoidal
integration is only approximate, not a claim of PETSc BDF/ARKIMEX
stage quadrature or discrete passivity.
"""
import numpy as np


def select_port_power_sampling(ts_type, saved_times, *, accepted_steps=None,
                               dt=None, adapt_type=None, requested=None):
    """Return 'cn_step' or 'instantaneous' or reject an unsafe request.

    If the real PETSc accepted-step count is unknown, CN is treated
    conservatively: instantaneous nodal samples rather than invented
    CN-step quantities. For fixed-step CN additionally check all saved
    time spacings equal the declared physical dt to roundoff.
    """
    if requested not in (None, "auto", "instantaneous", "cn_step"):
        raise ValueError("Invalid port power sampling override")
    t=np.asarray(saved_times,dtype=float)
    if (t.ndim!=1 or not len(t) or not np.all(np.isfinite(t))
            or (len(t)>1 and not np.all(np.diff(t)>0))):
        raise ValueError("Invalid accepted PETSc time history")
    scheme=str(ts_type).lower()
    valid_cn=False
    if scheme=="cn" and accepted_steps is not None:
        n=int(accepted_steps)
        valid_cn=(n>=1 and len(t)==n+1)
        if valid_cn and str(adapt_type).lower()=="none":
            try:
                requested_dt=float(dt)
            except (TypeError,ValueError):
                valid_cn=False
            else:
                valid_cn=(np.isfinite(requested_dt) and requested_dt>0
                          and np.allclose(np.diff(t),requested_dt,
                                          rtol=2e-10,atol=2e-12))
    if requested=="cn_step" and not valid_cn:
        raise ValueError(
            "CN step power requires every consecutive accepted TS state "
            "(and matching fixed dt when TS adaptation is disabled); "
            "use instantaneous sampling for subsampled results"
        )
    if requested=="instantaneous":
        return "instantaneous"
    return "cn_step" if valid_cn else "instantaneous"
