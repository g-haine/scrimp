# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Opt-in native PETSc/MUMPS constitutive-block and Schur-identity probe.

This does NOT install or benchmark a new SCRIMP solver. It leaves
the full GetFEM pH-DAE and the original MUMPS temporal solver unchanged.

Using the actual assembled K and M, form A(sigma)=K+sigma*M at two
different positive PETSc TS shifts. Partition actual GetFEM unknowns
into independent coenergies e and remaining x (dynamic + ALL ports).

A(sigma) in (e,x) row/column order:
    [-Q, C]     where -Q = K[e,e]
    [ D, G]

For every sigma, test that A[e,e]==-Q and check the *uncondensed*
Schur action Sx = (G + D*Q^{-1}*C)*x against A*[Q^{-1}Cx,x].
The full boundary multiplier variables stay in x. An actual MPI
MUMPS KSP solves Q*w=C*x; its operator is never reset across shifts
or right-hand sides. No global matrix gather, dense inversion or
changes to time discretization. All matrix norms are collective.

A valid result establishes an exact algebraic Schur identity for
the sampled vectors, NOT a scalable Schur preconditioner and NOT
a full spectral proof. Its purpose is to qualify the next prototype.
"""

import argparse
import json
import os
from pathlib import Path

from mpi4py import MPI
import numpy as np
from petsc4py import PETSc

from bench_acoustic_2d import _build_wave
from scrimp.utils.ph_block_layout import discover_layout


def _owned_variable_indices(fields, lo, hi, *, coenergy):
    """Sorted original GLOBAL PETSc indices, locally owned by a rank."""
    indices = []
    for field in fields:
        if (field["group"] == "coenergy") != coenergy:
            continue
        first = max(int(lo), int(field["start"]))
        last = min(int(hi), int(field["stop"]))
        if first < last:
            indices.extend(range(first, last))
    return np.asarray(indices, dtype=PETSc.IntType)


def _filled_right(mat, *, seed):
    """Deterministic manufactured distributed vector in submatrix order."""
    x = mat.createVecRight()
    lo, hi = x.getOwnershipRange()
    rows = np.arange(lo, hi, dtype=PETSc.IntType)
    values = (
        np.sin(0.00131 * (rows + 1) * seed)
        + 0.41 * np.cos(0.0041 * (rows + 3) * seed)
    )
    if len(rows):
        x.setValues(rows, values)
    x.assemble()
    return x


def _difference_rel(lhs, rhs, *, mpi):
    lhs, rhs = np.asarray(lhs), np.asarray(rhs)
    if lhs.shape != rhs.shape:
        raise ValueError("Schur and full-PETSc ownership distributions differ")
    delta = lhs-rhs
    squared_diff = mpi.allreduce(
        float(np.vdot(delta, delta).real), op=MPI.SUM
    )
    squared_ref = mpi.allreduce(
        float(np.vdot(rhs, rhs).real), op=MPI.SUM
    )
    return float(np.sqrt(squared_diff / max(squared_ref, 1e-300)))


def probe(args):
    mpi = MPI.COMM_WORLD
    rank, ranks = mpi.Get_rank(), mpi.Get_size()
    if not (args.mesh > 0 and args.dt > 0):
        raise ValueError("Mesh and dt must be positive")
    if not (
        np.isfinite(args.shift1) and np.isfinite(args.shift2)
        and args.shift1 > 0 and args.shift2 > 0
        and args.shift1 != args.shift2
    ):
        raise ValueError("Require distinct finite positive PETSc shifts")
    if os.environ.get("SCRIMP_GETFEM_MATRIX_MODE") != "distributed":
        raise ValueError("Require SCRIMP_GETFEM_MATRIX_MODE=distributed")
    if os.environ.get("SCRIMP_GETFEM_BRIDGE_DISTRIBUTED") != "coo":
        raise ValueError("Require SCRIMP_GETFEM_BRIDGE_DISTRIBUTED=coo")

    if args.petsc_log:
        PETSc.Log.begin()

    system = _build_wave(args.mesh)
    system.set_time_scheme(
        ts_type="cn", ts_adapt_type="none",
        ts_exact_final_time="matchstep",
        dt=args.dt, dt_save=args.dt,
        t_f=2.0 * args.dt,
        ksp_type="preonly", pc_type="lu",
        pc_factor_mat_solver_type="mumps",
    )
    system.solve()
    if (
        system.last_ts_converged_reason <= 0
        or system.last_ts_steps != 2
        or system.last_ts_rejections != 0
        or system.last_ts_snes_failures != 0
    ):
        raise RuntimeError("Reference MUMPS/CN PDE integration failed")
    if not system.linear_mass or not system.linear_stiffness:
        raise RuntimeError("Require fully linear time-invariant FE blocks")

    fields = discover_layout(
        system.gf_model,
        states=system.states,
        costates=system.costates,
        ports=system.ports,
        controls=system.controls,
    )
    K, M = system.stiffness, system.mass
    n = int(system.gf_model.nbdof())
    first, last = K.getOwnershipRange()
    local_e = _owned_variable_indices(
        fields, first, last, coenergy=True
    )
    local_x = _owned_variable_indices(
        fields, first, last, coenergy=False
    )
    if local_e.size + local_x.size != last-first:
        raise RuntimeError("A PETSc row is missing from the pH partition")
    n_e = mpi.allreduce(int(len(local_e)), op=MPI.SUM)
    n_x = mpi.allreduce(int(len(local_x)), op=MPI.SUM)
    if n_e + n_x != n or min(n_e, n_x) <= 0:
        raise RuntimeError("Incomplete coenergy/complement global split")

    is_e = PETSc.IS().createGeneral(local_e, comm=PETSc.COMM_WORLD)
    is_x = PETSc.IS().createGeneral(local_x, comm=PETSc.COMM_WORLD)

    neg_Q = K.createSubMatrix(is_e, is_e)
    Q = neg_Q.copy()
    Q.scale(-1.0)
    Q.assemble()
    qnorm = float(Q.norm(PETSc.NormType.FROBENIUS))
    if qnorm <= 0:
        raise ValueError("Empty constitutive -Q block")
    symmetry = bool(Q.isSymmetric(tol=1e-10))
    if not symmetry:
        raise ValueError(
            "Sampled coenergy block is not symmetric: cannot yet "
            "treat it as the constitutive mass matrix Q"
        )
    diagonal = Q.getDiagonal()
    diag = diagonal.getArray(readonly=True)
    local_diag_min = float(np.min(diag)) if len(diag) else np.inf
    min_diag = float(mpi.allreduce(local_diag_min, op=MPI.MIN))
    diagonal.destroy()
    if min_diag <= 0:
        raise ValueError(
            "Constitutive Q has nonpositive diagonal; do not "
            "assume a coercive constitutive block"
        )

    # One MPI-capable direct KSP for the time-INDEPENDENT Q block.
    # Crucial: no setOperators/reset of this KSP when sigma changes.
    ksp = PETSc.KSP().create(comm=PETSc.COMM_WORLD)
    ksp.setType("preonly")
    ksp.setTolerances(rtol=1e-12, atol=1e-14)
    pc = ksp.getPC()
    pc.setType("lu")
    pc.setFactorSolverType("mumps")
    ksp.setOperators(Q)

    series = []
    for j, sigma in enumerate((args.shift1, args.shift2)):
        # Reconstruct the exact ORIGINAL PETSc implicit linear Jacobian
        # (not the hypothetical reduced physical operator).
        A = K.copy()
        A.axpy(
            float(sigma), M,
            structure=PETSc.Mat.Structure.DIFFERENT_NONZERO_PATTERN,
        )
        A.assemble()

        Aee = A.createSubMatrix(is_e, is_e)
        diff = Aee.copy()
        diff.axpy(
            1.0, Q,
            structure=PETSc.Mat.Structure.DIFFERENT_NONZERO_PATTERN,
        )
        relative_constitutive_change = (
            float(diff.norm(PETSc.NormType.FROBENIUS)) / qnorm
        )
        if relative_constitutive_change > 2e-11:
            raise ValueError(
                "Constitutive -Q block depends on TS shift "
                f"{sigma}: relative change={relative_constitutive_change:g}"
            )

        C = A.createSubMatrix(is_e, is_x)
        D = A.createSubMatrix(is_x, is_e)
        G = A.createSubMatrix(is_x, is_x)

        # The Schur apply has the sign +D Q^{-1} C because
        # A_ee=-Q. The full uncondensed state [e;x] satisfies
        # algebraic e-rows to solver tolerance.
        x = _filled_right(G, seed=j+1)
        q_rhs = C.createVecLeft()
        C.mult(x, q_rhs)
        e = Q.createVecRight()
        ksp.solve(q_rhs, e)
        if ksp.getConvergedReason() <= 0:
            raise RuntimeError("Constitutive MUMPS/KSP did not converge")
        q_check = Q.createVecLeft()
        Q.mult(e, q_check)
        q_check.axpy(-1.0, q_rhs)
        q_rhs_norm = float(q_rhs.norm())
        q_residual = float(q_check.norm()) / max(q_rhs_norm, 1e-300)

        reduced = G.createVecLeft()
        G.mult(x, reduced)
        added = D.createVecLeft()
        D.mult(e, added)
        reduced.axpy(1.0, added)

        if e.getLocalSize() != local_e.size or x.getLocalSize() != local_x.size:
            raise RuntimeError(
                "PETSc MatCreateSubMatrix changed the expected "
                "rank-local IS distribution"
            )
        original = A.createVecRight()
        original.set(0.0)
        if local_e.size:
            original.setValues(
                local_e, e.getArray(readonly=True)
            )
        if local_x.size:
            original.setValues(
                local_x, x.getArray(readonly=True)
            )
        original.assemble()
        full_y = A.createVecLeft()
        A.mult(original, full_y)

        full_x_result = full_y.getValues(local_x)
        full_e_residual = full_y.getValues(local_e)
        schur_action_error = _difference_rel(
            reduced.getArray(readonly=True),
            full_x_result, mpi=mpi,
        )
        # For the homogeneous algebraic equation normalize the
        # global residual by ||Cx||, not by zero.
        alg_local_sq = float(np.vdot(
            full_e_residual, full_e_residual
        ).real)
        algebraic_residual_rel = float(
            np.sqrt(mpi.allreduce(alg_local_sq, op=MPI.SUM))
            / max(q_rhs_norm, 1e-300)
        )

        if max(
            q_residual, schur_action_error, algebraic_residual_rel
        ) > 5e-9:
            raise RuntimeError(
                "Exact constitutive Schur identity failed: "
                f"Q residual={q_residual}, "
                f"Schur action={schur_action_error}, "
                f"DAE algebraic row={algebraic_residual_rel}"
            )

        series.append({
            "sigma": float(sigma),
            "relative_Aee_plus_Q": relative_constitutive_change,
            "Q_MUMPS_manufactured_relative_residual": q_residual,
            "exact_schur_action_relative_error": schur_action_error,
            "full_DAE_algebraic_relative_residual": algebraic_residual_rel,
            "KSP_iterations": int(ksp.getIterationNumber()),
            "KSP_converged_reason": int(ksp.getConvergedReason()),
            "MUMPS_Q_operator_reused": True,
        })
        if rank == 0:
            print(
                f"PASS sigma={sigma:.9g}: "
                f"Aee+Q={relative_constitutive_change:.2e}, "
                f"Qsolve={q_residual:.2e}, "
                f"Schur={schur_action_error:.2e}, "
                f"algebraic={algebraic_residual_rel:.2e}",
                flush=True,
            )

        for obj in (
            original, full_y, reduced, added, x, q_rhs, e, q_check,
            C, D, G, Aee, diff, A,
        ):
            obj.destroy()

    # Confirm Q is a genuine invertible block on manufactured vectors
    # independent of the choice of sigma (different RHS, same factors).
    truth = _filled_right(Q, seed=7)
    rhs = Q.createVecLeft()
    Q.mult(truth, rhs)
    sol = Q.createVecRight()
    ksp.solve(rhs, sol)
    if ksp.getConvergedReason() <= 0:
        raise RuntimeError("Manufactured MUMPS solution did not converge")
    error = sol.copy()
    error.axpy(-1.0, truth)
    manufactured_error = (
        float(error.norm()) / max(float(truth.norm()), 1e-300)
    )
    if manufactured_error > 5e-9:
        raise RuntimeError(
            f"Constitutive MUMPS manufactured solution differs: "
            f"{manufactured_error:g}"
        )

    report = {
        "schema": "scrimp-constitutive-schur-probe-v1",
        "source_commit": os.environ.get("SCRIMP_BENCH_GIT_SHA", ""),
        "mpi_ranks": ranks,
        "mesh_h": args.mesh,
        "dofs": n,
        "coenergy_dofs": n_e,
        "complement_dofs": n_x,
        "coenergy_fields": [
            f["name"] for f in fields if f["group"] == "coenergy"
        ],
        "complement_fields": [
            f["name"] for f in fields if f["group"] != "coenergy"
        ],
        "constitutive_Q_frobenius_norm": qnorm,
        "constitutive_Q_symmetric": symmetry,
        "constitutive_Q_min_diagonal": min_diag,
        "constitutive_Q_manufactured_solution_relative_error":
            manufactured_error,
        "Q_solver": "PETSc KSPPREONLY + PCLU + MUMPS",
        "Q_factor_operator_configuration_count": 1,
        "Q_reused_across_shifts_and_rhs": True,
        "shift_results": series,
        "reference_pde": {
            "system": "wave.py (full mixed DAE)",
            "integrator": "CN 2 fixed steps only to assemble matrices",
            "mumps_full_system_unchanged": True,
            "converged_reason": int(system.last_ts_converged_reason),
            "accepted_steps": int(system.last_ts_steps),
        },
        "caveats": [
            "Schur action is verified WITHOUT reducing production DAE.",
            "Symmetry + positive diagonal do NOT prove global SPD.",
            "MUMPS is retained for full-solver reference and Q solves.",
            "No new iterative preconditioner installed or timed yet.",
            "The two shift values emulate adaptive TS Jacobian shifts "
            "on the same time-invariant FE operators; they do NOT "
            "claim a complete adaptive-integration trajectory test.",
        ],
    }
    if rank == 0:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(
            json.dumps(report, indent=2) + "\n", encoding="utf-8"
        )
        print(
            f"PASS: full pH-DAE Schur block identity at two TS shifts, "
            f"1 MUMPS Q operator, {n_e}/{n} coenergy DoFs",
            flush=True,
        )
        print(f"Saved: {args.output}", flush=True)

    for obj in (error, truth, rhs, sol, ksp, Q, neg_Q, is_e, is_x):
        obj.destroy()
    if args.petsc_log:
        PETSc.Log.view()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mesh", type=float, default=0.08)
    parser.add_argument("--dt", type=float, default=0.002)
    parser.add_argument(
        "--shift1", type=float, default=1000.,
        help="First positive PETSc IJacobian shift",
    )
    parser.add_argument(
        "--shift2", type=float, default=2000.,
        help="Second positive PETSc IJacobian shift",
    )
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--petsc-log", action="store_true",
        help="Print PETSc LU symbolic/numeric/MatSolve counters",
    )
    probe(parser.parse_args())


if __name__ == "__main__":
    main()
