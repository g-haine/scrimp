# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Minimal native MPI PETSc/MUMPS factorization-and-reset diagnosis.

This program deliberately does NOT import SCRIMP, GetFEM or Gmsh.
It constructs a small, well-conditioned MPIAIJ tridiagonal matrix,
solves it with LU/MUMPS, then explicitly resets KSP (which destroys
its MUMPS factor). The goal is to isolate the native failure seen in
SCRIMP's pseudo-TS initializer at MatDestroy_MUMPS -> TSReset.

Examples:
  mpiexec -n 1 python -u sandbox/test_petsc_mumps_reset.py --solver mumps
  mpiexec -n 2 python -u sandbox/test_petsc_mumps_reset.py --solver mumps
  mpiexec -n 2 python -u sandbox/test_petsc_mumps_reset.py --solver bjacobi

A successful minimal run does not rule out a MUMPS defect for the
real, 3300x3300 indefinite Saint-Venant operator.
"""

import argparse
import sys

from petsc4py import PETSc


def trace(rank, size, message):
    print(f"[rank {rank}/{size}] {message}", flush=True)


def build_matrix(n, comm):
    A = PETSc.Mat().createAIJ(size=(n, n), nnz=3, comm=comm)
    A.setUp()
    lo, hi = A.getOwnershipRange()
    for i in range(lo, hi):
        cols = [i]
        vals = [4.0]
        if i > 0:
            cols.append(i - 1)
            vals.append(-1.0)
        if i + 1 < n:
            cols.append(i + 1)
            vals.append(-1.0)
        A.setValues([i], cols, [vals], addv=PETSc.InsertMode.INSERT_VALUES)
    A.assemble()
    return A


def run(solver, size, cycles):
    comm = PETSc.COMM_WORLD
    rank, nranks = comm.getRank(), comm.getSize()
    n = max(size, 16 * nranks)
    trace(rank, nranks, f"starting: solver={solver} N={n}, cycles={cycles}")

    for cycle in range(cycles):
        A = build_matrix(n, comm)
        b, x = A.createVecRight(), A.createVecLeft()
        ksp = PETSc.KSP().create(comm=comm)
        try:
            b.set(1.0)
            x.set(0.0)
            ksp.setOperators(A)
            pc = ksp.getPC()
            if solver == "mumps":
                ksp.setType("preonly")
                pc.setType("lu")
                pc.setFactorSolverType("mumps")
            else:
                ksp.setType("gmres")
                pc.setType("bjacobi")
                ksp.setTolerances(rtol=1e-10, max_it=400)

            trace(rank, nranks, f"cycle={cycle} before KSP.solve")
            ksp.solve(b, x)
            trace(
                rank, nranks,
                f"cycle={cycle} KSP.solve completed, "
                f"reason={int(ksp.getConvergedReason())}, "
                f"its={ksp.getIterationNumber()}",
            )
            if ksp.getConvergedReason() <= 0:
                raise RuntimeError("Minimal PETSc KSP did not converge")

            trace(rank, nranks, f"cycle={cycle} before KSP.reset")
            ksp.reset()
            trace(rank, nranks, f"cycle={cycle} after KSP.reset")
        finally:
            trace(rank, nranks, f"cycle={cycle} before KSP.destroy")
            ksp.destroy()
            trace(rank, nranks, f"cycle={cycle} after KSP.destroy")
            b.destroy()
            x.destroy()
            A.destroy()

    trace(rank, nranks, "all cycles finished")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--solver", choices=("mumps", "bjacobi"),
                        default="mumps")
    parser.add_argument("--size", type=int, default=160)
    parser.add_argument("--cycles", type=int, default=3)
    args = parser.parse_args()
    if args.size <= 0 or args.cycles <= 0:
        parser.error("size and cycles must be positive")
    try:
        run(args.solver, args.size, args.cycles)
    except Exception as exc:
        print(f"minimal PETSc test failed: {exc!r}", file=sys.stderr, flush=True)
        raise
