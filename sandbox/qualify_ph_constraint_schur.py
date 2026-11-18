# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See LICENSE for license information.
#
# github: https://github.com/g-haine/scrimp

"""Qualify graph-discovered zero-self blocks using actual distributed PETSc solves.

For a candidate lambda selected purely by the graph audit, partition the
UNMODIFIED A(sigma) in the order (z, lambda):
    [ H  B ]
    [ C  D ].
Test MUMPS invertibility of H and a matrix-free exact Schur identity:
    S_lambda v = D v - C H^{-1} B v.
This is a read-only algebraic experiment; no condensed production system,
no global dense gather, and no assumption about boundary-variable names.
A zero self-block is necessary for this particular probe, NOT proof of a
Lagrange multiplier. A failed H solve is a rejected partition.
"""
import argparse
import json
import os
from pathlib import Path

import numpy as np
from mpi4py import MPI
from petsc4py import PETSc

from bench_acoustic_2d import _build_wave
from audit_ph_block_graph import audit_blocks
from probe_constitutive_mumps_schur import _filled_right
from scrimp.utils.ph_block_layout import discover_layout


def locally_owned_indices(fields, lo, hi, selected, *, invert=False):
    result = []
    for field in fields:
        keep = (field["name"] in selected)
        if invert:
            keep = not keep
        if not keep:
            continue
        start = max(lo, field["start"])
        stop = min(hi, field["stop"])
        if start < stop:
            result.extend(range(start, stop))
    return np.asarray(sorted(result), dtype=PETSc.IntType)


def rel_norm(diff, reference):
    return float(diff.norm()) / max(float(reference.norm()), 1e-300)


def qualify(A, fields, candidate, *, tol):
    mpi = MPI.COMM_WORLD
    lo, hi = A.getOwnershipRange()
    selected = {candidate}
    iz_local = locally_owned_indices(fields, lo, hi, selected, invert=True)
    il_local = locally_owned_indices(fields, lo, hi, selected)
    if len(iz_local) + len(il_local) != hi - lo:
        raise ValueError("Candidate split does not cover owned matrix rows")
    iz = PETSc.IS().createGeneral(iz_local, comm=PETSc.COMM_WORLD)
    il = PETSc.IS().createGeneral(il_local, comm=PETSc.COMM_WORLD)
    H = A.createSubMatrix(iz, iz)
    B = A.createSubMatrix(iz, il)
    C = A.createSubMatrix(il, iz)
    D = A.createSubMatrix(il, il)
    ks = PETSc.KSP().create(comm=PETSc.COMM_WORLD)
    ks.setType("preonly")
    ks.setErrorIfNotConverged(False)
    ks.setOperators(H)
    pc = ks.getPC()
    pc.setType("lu")
    pc.setFactorSolverType("mumps")
    report = {
        "candidate": candidate,
        "principal_dofs": int(mpi.allreduce(len(iz_local), op=MPI.SUM)),
        "candidate_dofs": int(mpi.allreduce(len(il_local), op=MPI.SUM)),
        "diagonal_lambda_norm": float(D.norm(PETSc.NormType.FROBENIUS)),
        "H_solve_reason": None,
        "H_manufactured_relative_error": None,
        "H_manufactured_relative_residual": None,
        "schur_action_relative_error": None,
        "algebraic_principal_row_residual": None,
        "qualified": False,
    }
    if report["principal_dofs"] == 0 or report["candidate_dofs"] == 0:
        raise ValueError("Empty Schur partition")
    try:
        # The H solve checks invertibility numerically on a manufactured
        # RHS. It does not by itself prove full rank or spectral properties.
        truth = _filled_right(H, seed=13)
        rhs = H.createVecLeft()
        H.mult(truth, rhs)
        found = H.createVecRight()
        try:
            ks.solve(rhs, found)
            reason = int(ks.getConvergedReason())
            report["H_solve_reason"] = reason
            if reason <= 0:
                report["rejected_reason"] = "MUMPS failed on principal block"
                return report
            err = found.copy()
            err.axpy(-1.0, truth)
            report["H_manufactured_relative_error"] = rel_norm(err, truth)
            err.destroy()
            check = H.createVecLeft()
            H.mult(found, check)
            check.axpy(-1.0, rhs)
            report["H_manufactured_relative_residual"] = rel_norm(check, rhs)
            check.destroy()
        finally:
            for obj in (truth, rhs, found):
                obj.destroy()

        v = _filled_right(D, seed=17)
        bv = B.createVecLeft()
        B.mult(v, bv)
        w = H.createVecRight()
        ks.solve(bv, w)
        if ks.getConvergedReason() <= 0:
            report["rejected_reason"] = "H inverse failed on coupling RHS"
            for obj in (v, bv, w):
                obj.destroy()
            return report

        reduced = D.createVecLeft()
        D.mult(v, reduced)
        correction = C.createVecLeft()
        C.mult(w, correction)
        reduced.axpy(-1.0, correction)

        # Form original distributed [z;lambda] with z=-H^-1 B lambda,
        # preserving native GetFEM ordering on the complete A.
        if w.getLocalSize() != len(iz_local) or v.getLocalSize() != len(il_local):
            raise RuntimeError("PETSc submatrix ownership differs from selected IS")
        full = A.createVecRight()
        full.set(0.0)
        if len(iz_local):
            full.setValues(iz_local, -w.getArray(readonly=True))
        if len(il_local):
            full.setValues(il_local, v.getArray(readonly=True))
        full.assemble()
        afull = A.createVecLeft()
        A.mult(full, afull)
        # Ownership-compatible vectors in submatrix order.
        observed_lambda = D.createVecLeft()
        observed_z = H.createVecLeft()
        if len(il_local):
            observed_lambda.setValues(
                np.arange(*observed_lambda.getOwnershipRange(), dtype=PETSc.IntType),
                afull.getValues(il_local),
            )
        if len(iz_local):
            observed_z.setValues(
                np.arange(*observed_z.getOwnershipRange(), dtype=PETSc.IntType),
                afull.getValues(iz_local),
            )
        observed_lambda.assemble()
        observed_z.assemble()
        delta = observed_lambda.copy()
        delta.axpy(-1.0, reduced)
        report["schur_action_relative_error"] = rel_norm(delta, reduced)
        report["algebraic_principal_row_residual"] = (
            float(observed_z.norm()) / max(float(bv.norm()), 1e-300)
        )
        report["qualified"] = (
            report["H_manufactured_relative_error"] < tol
            and report["H_manufactured_relative_residual"] < tol
            and report["schur_action_relative_error"] < tol
            and report["algebraic_principal_row_residual"] < tol
        )
        if not report["qualified"]:
            report["rejected_reason"] = "Schur algebraic accuracy gate"
        for obj in (v, bv, w, reduced, correction, full, afull,
                    observed_lambda, observed_z, delta):
            obj.destroy()
        return report
    except PETSc.Error as exc:
        report["rejected_reason"] = f"PETSc error during H factor/solve: {exc}"
        return report
    finally:
        for obj in (ks, H, B, C, D, iz, il):
            obj.destroy()


def run(args):
    mpi = MPI.COMM_WORLD
    if not (args.mesh > 0 and args.dt > 0 and
            args.shift1 > 0 and args.shift2 > 0 and args.shift1 != args.shift2):
        raise ValueError("Invalid mesh, time step or Jacobian shifts")
    if (os.getenv("SCRIMP_GETFEM_MATRIX_MODE") != "distributed"
            or os.getenv("SCRIMP_GETFEM_BRIDGE_DISTRIBUTED") != "coo"):
        raise RuntimeError("Require distributed GetFEM-to-PETSc COO bridge")
    system = _build_wave(args.mesh)
    system.set_time_scheme(
        ts_type="cn", ts_adapt_type="none",
        ts_exact_final_time="matchstep",
        dt=args.dt, dt_save=args.dt, t_f=2*args.dt,
        ksp_type="preonly", pc_type="lu",
        pc_factor_mat_solver_type="mumps",
    )
    system.solve()
    if (system.last_ts_converged_reason <= 0
            or system.last_ts_steps != 2
            or system.last_ts_rejections or system.last_ts_snes_failures):
        raise RuntimeError("MUMPS/CN reference integration failed")
    fields = discover_layout(
        system.gf_model, states=system.states,
        costates=system.costates, ports=system.ports,
        controls=system.controls,
    )
    reports = []
    all_pass = True
    for sigma in (args.shift1, args.shift2):
        A = system.stiffness.copy()
        A.axpy(sigma, system.mass,
               structure=PETSc.Mat.Structure.DIFFERENT_NONZERO_PATTERN)
        A.assemble()
        graph = audit_blocks(A, fields, atol=args.atol)
        candidates = [x["name"] for x in graph["zero_self_block_candidates"]]
        results = [qualify(A, fields, name, tol=args.tol) for name in candidates]
        for item in results:
            if mpi.rank == 0:
                print(
                    f"{'PASS' if item['qualified'] else 'REJECT'} "
                    f"sigma={sigma:g} candidate={item['candidate']}: "
                    f"H_reason={item['H_solve_reason']} "
                    f"H_error={item['H_manufactured_relative_error']} "
                    f"Schur_error={item['schur_action_relative_error']}",
                    flush=True,
                )
        all_pass = all_pass and all(x["qualified"] for x in results)
        reports.append({"sigma": sigma, "candidates": results})
        A.destroy()
    result = {
        "schema": "scrimp-ph-constraint-schur-qualification-v1",
        "source_commit": os.getenv("SCRIMP_BENCH_GIT_SHA", ""),
        "mpi_ranks": mpi.size,
        "dofs": int(system.gf_model.nbdof()),
        "mesh": args.mesh,
        "shifts": reports,
        "all_candidates_qualified": all_pass,
        "limitations": [
            "Candidate selection is algebraic, independent of field name.",
            "Only manufactured H solves and exact Schur action tested.",
            "Numerical success does not prove H is nonsingular for all models.",
            "No preconditioner enabled, no rank assertion or physics trajectory gate.",
            "Submatrices remain distributed; no dense global gather.",
        ],
    }
    if mpi.rank == 0:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
        print("Saved:", args.output, flush=True)
    if not all_pass:
        raise RuntimeError("One or more algebraic candidate splits not qualified")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--mesh", type=float, default=0.08)
    p.add_argument("--dt", type=float, default=0.002)
    p.add_argument("--shift1", type=float, default=1000.0)
    p.add_argument("--shift2", type=float, default=2000.0)
    p.add_argument("--atol", type=float, default=0.0)
    p.add_argument("--tol", type=float, default=5e-9)
    p.add_argument("--output", type=Path, required=True)
    run(p.parse_args())


if __name__ == "__main__":
    main()
