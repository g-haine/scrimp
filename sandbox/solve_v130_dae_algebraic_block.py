# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Numerical least-squares feasibility probe for a fixed-dynamic DAE block.

Solve K_aa * delta_a = -F_a with dynamic q,p frozen. A small *actual*
residual of the original IFunction is still required downstream;
convergence of LSMR alone is not DAE consistency. For an ill-conditioned
or rank-deficient matrix, failure to meet tolerance is INCONCLUSIVE
about mathematical incompatibility unless separately certified.
This module intentionally has no GetFEM/PETSc import dependencies.
"""

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import lsmr


def solve_algebraic_block(matrix, residual, *, atol=1e-10, rtol=1e-8):
    """Return candidate algebraic correction and transparent diagnostics.

    Parameters are a square, assembled algebraic Jacobian and the actual
    algebraic residual evaluated with *fixed* differential unknowns.
    The upper bound on the algebraic residual comes from this iterative
    candidate, not necessarily from a unique or true minimum.
    """
    A = csr_matrix(matrix)
    r = np.asarray(residual).reshape(-1)
    if A.shape != (r.size, r.size) or r.size == 0:
        raise ValueError("Expected a non-empty square algebraic block")
    if (not np.all(np.isfinite(A.data))
            or not np.all(np.isfinite(r))):
        raise ValueError("Non-finite algebraic matrix or residual")
    norm_before = float(np.linalg.norm(r))
    tolerance = float(atol + rtol * norm_before)
    if not np.isfinite(tolerance) or tolerance <= 0:
        raise ValueError("Positive finite algebraic tolerance is required")

    # A Jacobi-style diagonal scaling is NOT generally valid for this
    # saddle-point block: diagonal entries can be structurally zero.
    # LSMR handles rectangular/singular-like structure without changing
    # the actual constraint equations. Report condition estimate and
    # termination flags instead of asserting success from the API call.
    it_limit = max(1000, 5 * r.size)
    result = lsmr(
        A, -r, atol=1e-13, btol=1e-13,
        conlim=1e12, maxiter=it_limit,
    )
    correction, istop, iterations, normr, normar, normA, condA, normx = result
    corrected_linear = A @ correction + r
    norm_after = float(np.linalg.norm(corrected_linear))
    if (not np.all(np.isfinite(correction))
            or not np.isfinite(norm_after)):
        raise FloatingPointError("Non-finite algebraic LSMR candidate")

    qualified = bool(norm_after <= tolerance)
    if qualified:
        decision = "LINEAR_COMPATIBLE"
    elif istop in (1, 2) and condA < 1e12:
        decision = "NONZERO_LEAST_SQUARES_RESIDUAL"
    else:
        decision = "UNRESOLVED_RANK_OR_CONDITIONING"
    return correction, {
        "algebraic_dimension": int(r.size),
        "algebraic_jacobian_nnz": int(A.nnz),
        "raw_constraint_residual_l2": norm_before,
        "linearized_corrected_residual_l2": norm_after,
        "requested_tolerance": tolerance,
        "linearized_compatible_within_tolerance": qualified,
        "lsmr_istop": int(istop),
        "lsmr_iterations": int(iterations),
        "lsmr_reported_residual": float(normr),
        "lsmr_normal_equation_residual": float(normar),
        "lsmr_operator_norm_estimate": float(normA),
        "lsmr_condition_estimate": float(condA),
        "algebraic_correction_l2": float(normx),
        "linear_feasibility_decision": decision,
        "interpretation": (
            "Iterative feasibility evidence; a nonzero residual "
            "under poor conditioning is not proof of incompatibility. "
            "Recheck the actual production IFunction on the candidate."
        ),
    }
