# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Minimal real GetFEM MPI-2 RHS-vs-tangent scale experiment.

Unlike the full SWE, each 1D example uses ONE GetFEM brick. The
mathematical RHS sign convention is checked against its own assembled
tangent via a symmetric state perturbation. No TS/SNES/KSP/MUMPS.

Run identically on 1, 2, and 4 ranks:
  mpiexec -n 1 python -u sandbox/check_getfem_mpi_rhs.py
  mpiexec -n 2 python -u sandbox/check_getfem_mpi_rhs.py
  mpiexec -n 4 python -u sandbox/check_getfem_mpi_rhs.py

For N(u) = integral(u**p Test_u) on [0,1] at u=c=1:
  1^T DN(u) 1 = p,  1^T (-D rhs(u)) 1 = p.
If GetFEM's full RHS has an unintended P-fold factor after MPI
assembly, only the latter will scale with P.
"""

import hashlib
import json

import numpy as np
from mpi4py import MPI
from petsc4py import PETSc

from scrimp.utils.linalg import convert_gmm_to_petsc
import getfem as gf


def evaluate_case(mf, mim, expression, c=1.0, epsilon=1e-5):
    model = gf.Model("real")
    model.add_fem_variable("u", mf)
    model.add_nonlinear_term(mim, expression)
    n = model.nbdof()

    def rhs_at(u):
        model.set_variable("u", np.full(n, u))
        model.assembly(option="build_rhs")
        return np.asarray(model.rhs(), dtype=float).ravel().copy()

    plus, minus = rhs_at(c + epsilon), rhs_at(c - epsilon)
    derivative = -(plus - minus) / (2 * epsilon)

    model.set_variable("u", np.full(n, c))
    model.assembly(option="build_matrix")

    # Look at the GetFEM tangent BEFORE SCRIMP's 'replicated' transfer.
    # If GetFEM truly assembles a partitioned matrix, each rank's
    # raw contribution should generally differ, and the sum over MPI
    # ranks should give the theoretical integral (1,2,3).
    # If each rank already has the entire matrix, every raw action is
    # 1/2/3, and summing matrices would produce P times too much.
    raw = gf.Spmat("copy", model.tangent_matrix())
    raw.transpose()  # GetFEM CSC of transpose -> row-wise CSR
    row_ptr, col_idx = (np.asarray(x) for x in raw.csc_ind())
    raw_values = np.asarray(raw.csc_val())
    row_action = np.asarray([
        np.sum(raw_values[row_ptr[i]:row_ptr[i + 1]])
        for i in range(n)
    ])
    # Include the entire sparse CSR structure and values in the
    # fingerprint. A row-sum fingerprint alone is NOT proof that two
    # tangent matrices are identical: differing sparse matrices can
    # have the same action on the constant vector.
    full_digest = hashlib.blake2b(digest_size=16)
    for piece in (row_ptr, col_idx, raw_values):
        full_digest.update(str(piece.dtype).encode("ascii"))
        full_digest.update(
            np.ascontiguousarray(piece).view(np.uint8).tobytes()
        )
    rank_matrix_summaries = MPI.COMM_WORLD.allgather({
        "rank": MPI.COMM_WORLD.rank,
        "raw_tangent_integral": float(np.real(row_action.sum())),
        "raw_tangent_abs_sum": float(np.sum(np.abs(raw_values))),
        "raw_tangent_nnz": int(len(raw_values)),
        "raw_tangent_digest": hashlib.blake2b(
            np.ascontiguousarray(row_action).view(np.uint8),
            digest_size=8,
        ).hexdigest(),
        "raw_tangent_full_csr_digest": full_digest.hexdigest(),
    })

    matrix = PETSc.Mat().createAIJ(size=(n, n), comm=PETSc.COMM_WORLD)
    matrix.setUp()
    matrix.setOption(PETSc.Mat.Option.NEW_NONZERO_ALLOCATION_ERR, False)
    try:
        convert_gmm_to_petsc(model.tangent_matrix(), matrix)
        x, y = matrix.createVecRight(), matrix.createVecLeft()
        try:
            x.set(1.0)
            matrix.mult(x, y)
            jacobian_integral = float(np.real(x.dot(y)))
            lo, hi = y.getOwnershipRange()
            local_difference = (
                derivative[lo:hi] - np.asarray(y.getArray(readonly=True)).real
            )
            own_difference_sq = float(np.dot(local_difference, local_difference))
            difference_norm = float(
                MPI.COMM_WORLD.allreduce(own_difference_sq, op=MPI.SUM) ** 0.5
            )
        finally:
            x.destroy()
            y.destroy()
    finally:
        matrix.destroy()

    all_rhs = MPI.COMM_WORLD.allgather(derivative)
    rank_delta = max(
        float(np.linalg.norm(array - all_rhs[0])) for array in all_rhs
    )
    return {
        "expression": expression,
        "dofs": n,
        "getfem_negative_rhs_derivative_integral":
            float(np.sum(derivative)),
        "petsc_tangent_integral": jacobian_integral,
        "rhs_to_tangent_integral_ratio":
            float(np.sum(derivative)) / max(abs(jacobian_integral), 1e-14),
        "rhs_vs_tangent_vector_difference_norm": difference_norm,
        "rhs_rank_copy_difference": rank_delta,
        "raw_getfem_tangent_by_rank": rank_matrix_summaries,
        "sum_raw_getfem_tangent_integrals": sum(
            item["raw_tangent_integral"]
            for item in rank_matrix_summaries
        ),
        "raw_rhs_minus_tangent_ratio_at_each_rank": [
            float(np.sum(derivative)) /
            max(abs(item["raw_tangent_integral"]), 1e-14)
            for item in rank_matrix_summaries
        ],
    }


def main():
    world = MPI.COMM_WORLD
    mesh = gf.Mesh("cartesian", np.linspace(0.0, 1.0, 17))
    mf = gf.MeshFem(mesh, 1)
    mf.set_fem(gf.Fem("FEM_PK(1,1)"))
    mim = gf.MeshIm(mesh, gf.Integ("IM_GAUSS1D(5)"))

    cases = [
        evaluate_case(mf, mim, "u*Test_u"),
        evaluate_case(mf, mim, "u*u*Test_u"),
        evaluate_case(mf, mim, "u*u*u*Test_u"),
    ]
    if world.rank == 0:
        print(json.dumps({
            "mpi_ranks": world.size,
            "getfem_mpi_level": gf.util_mpi_parallelism_level(),
            "cases": cases,
        }, indent=2), flush=True)


if __name__ == "__main__":
    main()
