# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Compare distributed sparse operators modulo explicitly stored zeros.

PETSc COO assembly can retain additional *stored zero* coefficients
when the global matrix is split between MPI ranks. These are not
different numerical matrix entries and must not cause a false negative
in a genuine MPI1/MPI2 assembled-operator parity test.

Use canonical CSR, preserving and reporting raw structural information
rather than silently treating different stored-nnz values as identical.
"""

import numpy as np


def compare_assembled_csr(reference, parallel, expected_n):
    """Return numerical parity and auditable raw/canonical CSR evidence.

    Requires identical post-canonicalization sparsity and compares actual
    finite coefficient values with a strict absolute/relative tolerance.
    The bitwise flag has its literal meaning only *after eliminating*
    stored zeros and combining duplicated sparse entries.
    """
    A = reference.tocsr(copy=True)
    B = parallel.tocsr(copy=True)
    for K in (A, B):
        K.sum_duplicates()
        K.sort_indices()
        if not np.all(np.isfinite(K.data)):
            raise ValueError("Nonfinite GetFEM/PETSc stiffness coefficients")

    raw_nnz_a, raw_nnz_b = int(A.nnz), int(B.nnz)
    zeros_a = int(np.count_nonzero(A.data == 0))
    zeros_b = int(np.count_nonzero(B.data == 0))
    A.eliminate_zeros()
    B.eliminate_zeros()
    A.sort_indices()
    B.sort_indices()

    delta = (A - B).tocsr()
    delta.eliminate_zeros()
    max_difference = float(np.max(np.abs(delta.data))) if delta.nnz else 0.
    max_reference = float(np.max(np.abs(A.data))) if A.nnz else 0.
    threshold = 1e-12 + 1e-11 * max_reference
    same_canonical_pattern = bool(
        A.shape == B.shape
        and np.array_equal(A.indptr, B.indptr)
        and np.array_equal(A.indices, B.indices)
    )
    same_canonical_bytes = bool(
        same_canonical_pattern and np.array_equal(A.data, B.data)
    )
    passed = bool(
        A.shape == B.shape == (int(expected_n), int(expected_n))
        and same_canonical_pattern
        and max_difference <= threshold
    )
    return {
        "passed": passed,
        "global_shape": list(A.shape),
        "MPI1_stored_nnz": raw_nnz_a,
        "MPI2_stored_nnz": raw_nnz_b,
        "MPI1_explicit_zeros": zeros_a,
        "MPI2_explicit_zeros": zeros_b,
        "MPI1_numerical_nnz": int(A.nnz),
        "MPI2_numerical_nnz": int(B.nnz),
        "raw_stored_nnz_identical": raw_nnz_a == raw_nnz_b,
        "canonical_sparsity_identical": same_canonical_pattern,
        "canonical_CSR_bitwise_identical": same_canonical_bytes,
        "difference_nnz": int(delta.nnz),
        "max_abs_difference": max_difference,
        "threshold": threshold,
        "comparison_convention": (
            "CSR after sum_duplicates and eliminate_zeros, "
            "numerical values compared with explicit tolerance"
        ),
    }
