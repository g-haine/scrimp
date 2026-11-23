#!/usr/bin/env python3
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp
"""PETSc LU backend preflight on a real complex, distributed AIJ matrix.

This is not an SCRIMP example, a scaling benchmark, or a GetFEM
qualification. It only checks that installed native PETSc backends
can factor and solve the *same* nontrivial parallel sparse system.
"""

import argparse
import json
import math
from pathlib import Path
import time

import numpy as np
from petsc4py import PETSc


def run():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--backend", choices=("mumps", "superlu", "superlu_dist"),
        required=True,
    )
    parser.add_argument("--size", type=int, default=384)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()

    comm = PETSc.COMM_WORLD
    ranks = comm.getSize()
    rank = comm.getRank()
    n = args.size
    if n < 16 or n < 4 * ranks:
        raise ValueError("Matrix too small for a meaningful MPI test")
    if args.backend == "superlu" and ranks != 1:
        raise ValueError("Sequential SuperLU is only tested at MPI=1")
    if not PETSc.Sys.hasExternalPackage(args.backend):
        raise RuntimeError("PETSc external package unavailable: " + args.backend)
    if not args.out.parent.is_dir():
        raise FileNotFoundError("Host-owned output directory not mounted")

    # Original native PETSc scalar type is complex. Deliberately use a
    # nonsymmetric, complex, nonsingular tridiagonal operator so this
    # cannot pass as an accidentally real/symmetric/dense LU test.
    if not np.issubdtype(np.dtype(PETSc.ScalarType), np.complexfloating):
        raise RuntimeError("Expected the qualified complex PETSc image")

    matrix = PETSc.Mat().createAIJ(size=(n, n), nnz=3, comm=comm)
    matrix.setUp()
    lo, hi = matrix.getOwnershipRange()
    for i in range(lo, hi):
        cols = [i]
        vals = [4.0 + (0.08 + 0.001 * i) * 1j]
        if i > 0:
            cols.append(i - 1)
            vals.append(-0.8 + 0.03j)
        if i + 1 < n:
            cols.append(i + 1)
            vals.append(-1.1 - 0.02j)
        matrix.setValues(i, cols, vals)
    matrix.assemblyBegin()
    matrix.assemblyEnd()
    matrix_type = matrix.getType()
    if ranks > 1 and matrix_type != "mpiaij":
        raise RuntimeError(f"Expected a truly distributed MPIAIJ: {matrix_type}")
    if ranks == 1 and matrix_type not in ("seqaij", "mpiaij"):
        raise RuntimeError(f"Unexpected sparse matrix type: {matrix_type}")

    exact = matrix.createVecRight()
    for i in range(lo, hi):
        exact.setValue(i, 0.75 + (i % 11) / 19.0 + (i % 5) * 0.031j)
    exact.assemblyBegin()
    exact.assemblyEnd()
    rhs = matrix.createVecLeft()
    matrix.mult(exact, rhs)

    ksp = PETSc.KSP().create(comm=comm)
    ksp.setType("preonly")
    pc = ksp.getPC()
    pc.setType("lu")
    pc.setFactorSolverType(args.backend)
    ksp.setOperators(matrix)

    started = time.perf_counter()
    ksp.setUp()
    setup_seconds_local = time.perf_counter() - started

    solution = matrix.createVecRight()
    started = time.perf_counter()
    ksp.solve(rhs, solution)
    solve_seconds_local = time.perf_counter() - started
    reason = int(ksp.getConvergedReason())
    chosen = pc.getFactorSolverType()
    if reason <= 0:
        raise RuntimeError(f"{args.backend}: KSP did not converge: {reason}")
    if chosen != args.backend:
        raise RuntimeError(f"Unexpected solver: expected {args.backend}, got {chosen}")

    difference = solution.copy()
    difference.axpy(-1.0, exact)
    relative_solution_error = float(difference.norm() / exact.norm())
    residual = matrix.createVecLeft()
    matrix.mult(solution, residual)
    residual.axpy(-1.0, rhs)
    relative_residual = float(residual.norm() / rhs.norm())
    if (not math.isfinite(relative_solution_error)
            or not math.isfinite(relative_residual)
            or relative_solution_error >= 1e-11
            or relative_residual >= 1e-11):
        raise RuntimeError(
            f"{args.backend} MPI={ranks}: actual LU accuracy failure "
            f"relative solution={relative_solution_error:.9g}, "
            f"relative residual={relative_residual:.9g}"
        )

    # Distinct wall-clock samples are recorded only for diagnosis.
    # A single 384x384 factorization is emphatically NOT a benchmark.
    result = {
        "schema": "scrimp-v130-installed-petsc-lu-preflight-v1",
        "status": "PASS_ACTUAL_COMPLEX_SPARSE_LU",
        "backend": args.backend,
        "ranks": ranks,
        "matrix_size": n,
        "matrix_type": matrix_type,
        "nonzero_complex_matrix": True,
        "getfem_or_scrimp_example_qualified": False,
        "performance_comparison_qualified": False,
        "selected_factor_solver": chosen,
        "KSP_reason": reason,
        "relative_solution_l2_error": relative_solution_error,
        "relative_residual_l2": relative_residual,
        "local_setup_seconds_rank0": setup_seconds_local,
        "local_solve_seconds_rank0": solve_seconds_local,
    }
    comm.barrier()
    if rank == 0:
        if args.out.exists():
            raise FileExistsError("Never overwrite an earlier solver observation")
        args.out.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
        print(
            f"PASS PETSc LU backend={args.backend} MPI={ranks} "
            f"mat={matrix_type} n={n} "
            f"relative_error={relative_solution_error:.4g} "
            f"relative_residual={relative_residual:.4g}",
            flush=True,
        )
    comm.barrier()


if __name__ == "__main__":
    run()
