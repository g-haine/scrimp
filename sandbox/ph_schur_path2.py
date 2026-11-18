# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See LICENSE for license information.
#
# github: https://github.com/g-haine/scrimp

"""Experimental sparse length-two-path approximation of a constraint Schur.

Given original A in (principal z, constraint lambda) blocks [H B; C D],
let J = diag(H), E = H-J. The truncated diagonal expansion

    H^-1 ~ J^-1 - J^-1 E J^-1

induces the sparse path-preserving approximation

    S_hat = D - 2*C*J^-1*B + C*J^-1*H*J^-1*B.

Unlike PETSc SELFP (D-C*J^-1*B), this retains paths through one
off-diagonal coupling of H. It is NOT a convergent Neumann expansion
without additional assumptions, nor a guaranteed nonsingular PC.

Only sparse distributed PETSc MatMatMult and diagonal scaling. A
hard safety cap on constraint DoFs prevents accidental quadratic fill
on large models. No explicit exact Schur, dense gather, or system edit.
"""

from time import perf_counter

import numpy as np
from mpi4py import MPI
from petsc4py import PETSc


def build_path2_schur(A, iz, il, *, max_constraint_dofs=512):
    """Return sparse PETSc USER Schur Pmat and diagnostics (caller owns Mat)."""
    comm = MPI.COMM_WORLD
    nlam = int(il.getSize())
    if nlam <= 0 or nlam > max_constraint_dofs:
        raise ValueError(
            f"Path2 diagnostic limited to 1..{max_constraint_dofs} "
            f"constraint DoFs (found {nlam}); no fallback to dense assembly"
        )
    start = perf_counter()
    H = A.createSubMatrix(iz, iz)
    B = A.createSubMatrix(iz, il)
    C = A.createSubMatrix(il, iz)
    D = A.createSubMatrix(il, il)
    resources = [H, B, C, D]
    output = None
    try:
        diag = H.getDiagonal()
        resources.append(diag)
        values = diag.getArray(readonly=True)
        local_zero = int(np.count_nonzero(np.abs(values) <= 1e-14))
        zero = comm.allreduce(local_zero, op=MPI.SUM)
        if zero:
            raise ValueError(
                f"Principal H has {zero} near-zero diagonal entries; "
                "the graph-selected path2 approximation is not admissible"
            )
        inv = diag.copy()
        resources.append(inv)
        inv.reciprocal()

        # All intermediate products remain PETSc sparse and distributed.
        Bj = B.copy()
        resources.append(Bj)
        Bj.diagonalScale(L=inv)
        Cj = C.copy()
        resources.append(Cj)
        Cj.diagonalScale(R=inv)
        term0 = Cj.matMult(B)       # C J^-1 B
        resources.append(term0)
        HBj = H.matMult(Bj)        # H J^-1 B
        resources.append(HBj)
        term1 = Cj.matMult(HBj)    # C J^-1 H J^-1 B
        resources.append(term1)
        output = D.copy()
        output.axpy(
            -2.0, term0,
            structure=PETSc.Mat.Structure.DIFFERENT_NONZERO_PATTERN,
        )
        output.axpy(
            1.0, term1,
            structure=PETSc.Mat.Structure.DIFFERENT_NONZERO_PATTERN,
        )
        output.assemble()
        norm = float(output.norm(PETSc.NormType.FROBENIUS))
        if not np.isfinite(norm) or norm <= 0:
            raise ValueError("Path2 Schur approximation is zero/nonfinite")
        d = output.getDiagonal()
        try:
            vals = d.getArray(readonly=True)
            local_bad = int(np.count_nonzero(np.abs(vals) < 1e-14))
            bad = comm.allreduce(local_bad, op=MPI.SUM)
        finally:
            d.destroy()
        elapsed = comm.allreduce(perf_counter()-start, op=MPI.MAX)
        diag_report = {
            "path2_build_seconds_max": float(elapsed),
            "path2_frobenius": norm,
            "path2_near_zero_diagonal": int(bad),
            "path2_constraint_dofs": nlam,
            "approximation": "D-2*C*J^-1*B+C*J^-1*H*J^-1*B",
            "path2_nnz_used": int(output.getInfo()["nz_used"]),
        }
        return output, diag_report
    except BaseException:
        if output is not None:
            output.destroy()
        raise
    finally:
        for obj in reversed(resources):
            obj.destroy()
