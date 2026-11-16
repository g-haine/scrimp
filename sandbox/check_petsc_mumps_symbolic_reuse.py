# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Isolate PETSc/MUMPS symbolic reuse with the same MatCopy/MatAXPY sequence
as DPHS.IJacobian, without GetFEM, a time integrator, or a physical model.

The default 'native' mode deliberately does NOT call PC.reset(): PETSc's
PCSetUp should notice matrix nonzero-state changes and rebuild symbolic LU
only when needed. 'always' reproduces the old conservative workaround.

Six numerical Jacobians are tested: three with fixed sparsity, two with
an added second-neighbour coupling, then one without that coupling.
The full solution and relative algebraic residual are checked every time.
The local stored CSR pattern is compared exactly (no probabilistic hash),
and MPI.LOR reports changes on any rank.

Use PETSc's -log_view output (enabled by --petsc-log) to inspect the counts
of MatLUFactorSym, MatLUFactorNum and MatSolve; the first should be much
smaller in native mode for fixed-pattern updates. The test does not prove
GetFEM's nonlinear Jacobian respects PETSc's structure-state bookkeeping.

Example:
  mpiexec -n 4 python -u sandbox/check_petsc_mumps_symbolic_reuse.py \
      --mode native --petsc-log
"""

import argparse
import json

from mpi4py import MPI
import numpy as np
from petsc4py import PETSc


def stiffness_matrix(n, cycle):
    """A well-conditioned distributed operator with a controlled pattern."""
    matrix = PETSc.Mat().createAIJ(
        size=(n, n), nnz=5, comm=PETSc.COMM_WORLD
    )
    matrix.setUp()
    first, last = matrix.getOwnershipRange()
    extra = cycle in (3, 4)
    for row in range(first, last):
        cols = [row]
        vals = [8.0 + 0.2 * cycle]
        for distance, value in ((1, -1.0), (2, -0.15)):
            if distance == 2 and not extra:
                continue
            for col in (row - distance, row + distance):
                if 0 <= col < n:
                    cols.append(col)
                    vals.append(value)
        matrix.setValues([row], cols, [vals])
    matrix.assemble()
    return matrix


def mass_matrix(n, cycle):
    matrix = PETSc.Mat().createAIJ(
        size=(n, n), nnz=1, comm=PETSc.COMM_WORLD
    )
    matrix.setUp()
    first, last = matrix.getOwnershipRange()
    if last > first:
        diagonal = 1.0 + 0.01 * cycle
        for row in range(first, last):
            matrix.setValue(row, row, diagonal)
    matrix.assemble()
    return matrix


def _pattern(matrix):
    """Exact, rank-local stored CSR structure; coefficients are excluded."""
    indptr, indices, _values = matrix.getValuesCSR()
    return np.array(indptr, copy=True), np.array(indices, copy=True)


def _changed(before, after):
    return (
        not np.array_equal(before[0], after[0])
        or not np.array_equal(before[1], after[1])
    )


def run(mode, n, tolerance, petsc_log):
    mpi = MPI.COMM_WORLD
    rank, size = mpi.Get_rank(), mpi.Get_size()
    n = max(n, 24 * size)
    if petsc_log:
        PETSc.Log.begin()

    # The SAME PETSc Pmat, PC, and KSP persist across six solves.
    operator = PETSc.Mat().createAIJ(
        size=(n, n), nnz=5, comm=PETSc.COMM_WORLD
    )
    operator.setUp()
    operator.assemble()
    ksp = PETSc.KSP().create(comm=PETSc.COMM_WORLD)
    ksp.setType("preonly")
    pc = ksp.getPC()
    pc.setType("lu")
    pc.setFactorSolverType("mumps")
    ksp.setOperators(operator)

    exact = operator.createVecRight()
    rhs = operator.createVecLeft()
    numerical = operator.createVecRight()
    residual = operator.createVecLeft()
    exact.set(1.0)

    previous_pattern = None
    events = []
    try:
        for cycle in range(6):
            stiffness = stiffness_matrix(n, cycle)
            mass = mass_matrix(n, cycle)
            try:
                if mode == "always":
                    pc.reset()

                # Mirrors DPHS.IJacobian: never claim equal sparsity.
                structure = PETSc.Mat.Structure.DIFFERENT_NONZERO_PATTERN
                stiffness.copy(operator, structure=structure)
                operator.axpy(0.5 + 0.025 * cycle, mass,
                              structure=structure)
                operator.assemble()

                current_pattern = _pattern(operator)
                changed_local = (
                    previous_pattern is not None
                    and _changed(previous_pattern, current_pattern)
                )
                changed_any = bool(mpi.allreduce(changed_local, op=MPI.LOR))
                previous_pattern = current_pattern

                operator.mult(exact, rhs)
                numerical.set(0.0)
                ksp.solve(rhs, numerical)
                operator.mult(numerical, residual)
                residual.axpy(-1.0, rhs)
                residual_relative = residual.norm() / max(1.0, rhs.norm())
                numerical.axpy(-1.0, exact)
                state_relative = numerical.norm() / max(1.0, exact.norm())
                reason = int(ksp.getConvergedReason())
                if reason <= 0 or max(residual_relative, state_relative) > tolerance:
                    raise RuntimeError(
                        "MUMPS mismatch at cycle {}: reason={}, residual={}, "
                        "state={}".format(
                            cycle, reason, residual_relative, state_relative
                        )
                    )
                if cycle == 3 and not changed_any:
                    raise RuntimeError(
                        "The test failed to introduce a new stored nonzero "
                        "at cycle 3; the pattern-change test is inconclusive"
                    )
                if cycle in (1, 2, 4) and changed_any:
                    raise RuntimeError(
                        "Unexpected stored nonzero pattern change at "
                        f"cycle {cycle}"
                    )

                events.append({
                    "cycle": cycle,
                    "pattern_changed": changed_any,
                    "relative_residual": float(residual_relative),
                    "relative_state_error": float(state_relative),
                    "ksp_reason": reason,
                })
            finally:
                mass.destroy()
                stiffness.destroy()
    finally:
        residual.destroy()
        numerical.destroy()
        rhs.destroy()
        exact.destroy()
        ksp.destroy()
        operator.destroy()

    if rank == 0:
        print(json.dumps({
            "benchmark": "PETSc-only LU/MUMPS symbolic reuse",
            "mpi_size": size, "matrix_size": n, "mode": mode,
            "cycles": events,
        }, indent=2), flush=True)
    if petsc_log:
        PETSc.Log.view()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=("native", "always"),
                        default="native")
    parser.add_argument("--size", type=int, default=240)
    parser.add_argument("--tol", type=float, default=1e-10)
    parser.add_argument("--petsc-log", action="store_true")
    args = parser.parse_args()
    if args.size <= 0 or args.tol <= 0 or not np.isfinite(args.tol):
        parser.error("size and tol must be finite and positive")
    run(args.mode, args.size, args.tol, args.petsc_log)
