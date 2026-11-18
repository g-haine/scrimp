# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Opt-in PETSc FGMRES/PCFIELDSPLIT trial on the *unmodified* wave DAE.

This is a baseline block triangular PC, NOT yet a Schur/AMG solver.
The complementary dynamic+boundary saddle block is approximated by
rank-local block Jacobi/ILU. Failure to converge is reported explicitly.

The full PETSc operator is A(sigma)=K+sigma*M, at two distinct shifts;
the exact distributed full-MUMPS solution of a manufactured RHS provides
the reference. No dense gather, substituted coenergy or production change.
"""
import argparse
import json
import os
from pathlib import Path
from time import perf_counter

import numpy as np
from mpi4py import MPI
from petsc4py import PETSc

from bench_acoustic_2d import _build_wave
from probe_constitutive_mumps_schur import _owned_variable_indices, _filled_right
from scrimp.utils.ph_block_layout import discover_layout, classify_ph_formulation


def solve_reference(A, rhs):
    ksp = PETSc.KSP().create(comm=PETSc.COMM_WORLD)
    ksp.setType("preonly")
    pc = ksp.getPC()
    pc.setType("lu")
    pc.setFactorSolverType("mumps")
    ksp.setOperators(A)
    result = A.createVecRight()
    tic = perf_counter()
    ksp.solve(rhs, result)
    elapsed = MPI.COMM_WORLD.allreduce(perf_counter() - tic, op=MPI.MAX)
    reason = int(ksp.getConvergedReason())
    ksp.destroy()
    if reason <= 0:
        result.destroy()
        raise RuntimeError(f"Full-system MUMPS reference failed: {reason}")
    return result, elapsed


def block_trial(A, rhs, is_e, is_x, rtol, max_it, x_pc, e_pc):
    ksp = PETSc.KSP().create(comm=PETSc.COMM_WORLD)
    ksp.setOptionsPrefix("ph_trial_")
    ksp.setType("fgmres")
    ksp.setTolerances(rtol=rtol, atol=1e-13, max_it=max_it)
    ksp.setGMRESRestart(60)
    ksp.setErrorIfNotConverged(False)
    pc = ksp.getPC()
    pc.setType("fieldsplit")
    # Keep the failed triangular block PC as a reproducible control.
    # Schur mode uses A00=-Q and a shift-dependent SELF-P Schur PC:
    # S_p = G - D diag(-Q)^(-1) C. The true Schur action retains
    # the complete constitutive inverse via the e-subsolver.
    if x_pc == "schur-selfp":
        pc.setFieldSplitType(PETSc.PC.CompositeType.SCHUR)
        pc.setFieldSplitSchurFactType(PETSc.PC.FieldSplitSchurFactType.FULL)
        pc.setFieldSplitSchurPreType(PETSc.PC.FieldSplitSchurPreType.SELFP)
    else:
        pc.setFieldSplitType(PETSc.PC.CompositeType.MULTIPLICATIVE)
    pc.setFieldSplitIS(("e", is_e), ("x", is_x))
    opts = PETSc.Options()
    settings = {
        "ph_trial_fieldsplit_e_ksp_type": "preonly",
        "ph_trial_fieldsplit_e_pc_type": e_pc,
        "ph_trial_fieldsplit_x_ksp_type": "preonly",
        "ph_trial_fieldsplit_x_pc_type": x_pc,
    }
    if e_pc == "lu":
        settings["ph_trial_fieldsplit_e_pc_factor_mat_solver_type"] = "mumps"
    elif e_pc == "bjacobi":
        settings.update({
            "ph_trial_fieldsplit_e_sub_ksp_type": "preonly",
            "ph_trial_fieldsplit_e_sub_pc_type": "ilu",
        })
    else:
        raise ValueError("Unsupported coenergy preconditioner")
    if x_pc == "bjacobi":
        settings.update({
            "ph_trial_fieldsplit_x_sub_ksp_type": "preonly",
            "ph_trial_fieldsplit_x_sub_pc_type": "ilu",
        })
    elif x_pc in ("lu", "schur-selfp"):
        settings["ph_trial_fieldsplit_x_pc_type"] = "lu"
        settings["ph_trial_fieldsplit_x_pc_factor_mat_solver_type"] = "mumps"
    else:
        raise ValueError("Unsupported complement block preconditioner")
    try:
        for key, val in settings.items():
            opts[key] = val
        ksp.setOperators(A)
        ksp.setFromOptions()
        result = A.createVecRight()
        tic = perf_counter()
        ksp.solve(rhs, result)
        elapsed = MPI.COMM_WORLD.allreduce(
            perf_counter() - tic, op=MPI.MAX
        )
        reason = int(ksp.getConvergedReason())
        its = int(ksp.getIterationNumber())

        diagnostics = {"outer_pc_failed_reason": None, "sub_ksp": []}
        try:
            diagnostics["outer_pc_failed_reason"] = str(pc.getFailedReason())
        except (AttributeError, RuntimeError, ValueError) as exc:
            diagnostics["outer_pc_failed_reason"] = f"unavailable: {exc}"
        try:
            for name, sub in zip(("e", "x"), pc.getFieldSplitSubKSP()):
                entry = {
                    "field": name,
                    "ksp_type": sub.getType(),
                    "pc_type": sub.getPC().getType(),
                    "converged_reason": int(sub.getConvergedReason()),
                }
                try:
                    entry["pc_failed_reason"] = str(sub.getPC().getFailedReason())
                except (AttributeError, RuntimeError, ValueError) as exc:
                    entry["pc_failed_reason"] = f"unavailable: {exc}"
                diagnostics["sub_ksp"].append(entry)
        except (AttributeError, RuntimeError, ValueError) as exc:
            diagnostics["sub_ksp_error"] = str(exc)
        return result, elapsed, reason, its, diagnostics
    finally:
        ksp.destroy()
        for key in settings:
            del opts[key]



def block_diagonal_diagnostics(A, index_set, mpi):
    block = A.createSubMatrix(index_set, index_set)
    try:
        diagonal = block.getDiagonal()
        values = np.asarray(diagonal.getArray(readonly=True))
        local_zero = int(np.count_nonzero(np.abs(values) < 1e-14))
        local_size = int(values.size)
        diagonal.destroy()
        return {
            "size": int(mpi.allreduce(local_size, op=MPI.SUM)),
            "near_zero_diagonal": int(mpi.allreduce(local_zero, op=MPI.SUM)),
            "frobenius_norm": float(block.norm(PETSc.NormType.FROBENIUS)),
        }
    finally:
        block.destroy()


def relative_residual(A, rhs, x):
    y = rhs.duplicate()
    A.mult(x, y)
    y.axpy(-1, rhs)
    rel = float(y.norm()) / max(float(rhs.norm()), 1e-300)
    y.destroy()
    return rel


def run(args):
    mpi = MPI.COMM_WORLD
    rank = mpi.Get_rank()
    if args.mesh <= 0 or args.dt <= 0 or args.rtol <= 0:
        raise ValueError("mesh, dt, rtol must be positive")
    if args.max_it <= 0 or args.shift1 <= 0 or args.shift2 <= 0 or args.shift1 == args.shift2:
        raise ValueError("Need positive iteration limit and distinct positive shifts")
    if os.environ.get("SCRIMP_GETFEM_MATRIX_MODE") != "distributed":
        raise ValueError("Distributed GetFEM matrix mode is mandatory")
    if os.environ.get("SCRIMP_GETFEM_BRIDGE_DISTRIBUTED") != "coo":
        raise ValueError("Distributed COO bridge is mandatory")
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
        raise RuntimeError("Reference two-step MUMPS/CN integration failed")
    if not system.linear_mass or not system.linear_stiffness:
        raise RuntimeError("Expected time-independent linear mass and stiffness")
    fields = discover_layout(
        system.gf_model, states=system.states,
        costates=system.costates, ports=system.ports,
        controls=system.controls,
    )
    classification = classify_ph_formulation(fields, costates=system.costates)
    if classification["formulation"] != "implicit_constitutive":
        raise RuntimeError("Independent constitutive coenergy DoFs required")
    K, M = system.stiffness, system.mass
    first, last = K.getOwnershipRange()
    local_e = _owned_variable_indices(fields, first, last, coenergy=True)
    local_x = _owned_variable_indices(fields, first, last, coenergy=False)
    if local_e.size + local_x.size != last-first:
        raise RuntimeError("Split does not cover every locally owned row")
    is_e = PETSc.IS().createGeneral(local_e, comm=PETSc.COMM_WORLD)
    is_x = PETSc.IS().createGeneral(local_x, comm=PETSc.COMM_WORLD)
    results = []
    for j, sigma in enumerate((args.shift1, args.shift2)):
        A = K.copy()
        A.axpy(sigma, M,
               structure=PETSc.Mat.Structure.DIFFERENT_NONZERO_PATTERN)
        A.assemble()
        true_x = _filled_right(A, seed=j+11)
        rhs = A.createVecLeft()
        A.mult(true_x, rhs)
        ref, ref_seconds = solve_reference(A, rhs)
        ref_residual = relative_residual(A, rhs, ref)
        trial, trial_seconds, reason, its, diagnostics = block_trial(
            A, rhs, is_e, is_x, args.rtol, args.max_it, args.x_pc, args.e_pc
        )
        trial_residual = relative_residual(A, rhs, trial)
        difference = trial.copy()
        difference.axpy(-1.0, ref)
        rel_error = float(difference.norm()) / max(float(ref.norm()), 1e-300)
        success = (
            reason > 0 and trial_residual <= args.accept_residual
            and rel_error <= args.accept_error
        )
        item = {
            "sigma": sigma, "converged": bool(success),
            "fgmres_reason": reason, "fgmres_iterations": its,
            "complement_pc": args.x_pc,
            "coenergy_pc": args.e_pc,
            "diagnostics": diagnostics,
            "coenergy_diagonal": block_diagonal_diagnostics(A, is_e, mpi),
            "complement_diagonal": block_diagonal_diagnostics(A, is_x, mpi),
            "mumps_reference_residual": ref_residual,
            "fgmres_true_relative_residual": trial_residual,
            "relative_state_error_vs_mumps": rel_error,
            "mumps_seconds_max": ref_seconds,
            "fieldsplit_seconds_max": trial_seconds,
        }
        results.append(item)
        if rank == 0 and not success:
            print('DIAGNOSTIC:', json.dumps({'sigma': sigma, 'split': diagnostics, 'e': item['coenergy_diagonal'], 'x': item['complement_diagonal']}, sort_keys=True), flush=True)
        if rank == 0:
            print(
                f"{'PASS' if success else 'FAIL'} sigma={sigma:g}: "
                f"FGMRES its={its}, reason={reason}, "
                f"residual={trial_residual:.3e}, "
                f"state diff={rel_error:.3e}, "
                f"ref={ref_seconds:.3g}s, trial={trial_seconds:.3g}s",
                flush=True,
            )
        for obj in (difference, trial, ref, rhs, true_x, A):
            obj.destroy()
    report = {
        "schema": "scrimp-ph-fieldsplit-baseline-v1",
        "source_commit": os.environ.get("SCRIMP_BENCH_GIT_SHA", ""),
        "mpi_ranks": mpi.Get_size(), "mesh": args.mesh,
        "dofs": int(system.gf_model.nbdof()),
        "fields": fields,
        "formulation": classification["formulation"],
        "solver": (f"FGMRES/PCFIELDSPLIT ({'schur full, SELF-P' if args.x_pc == 'schur-selfp' else 'multiplicative'}); e: {args.e_pc}; x: {args.x_pc}"),
        "reference": "full PETSc MUMPS on identical A(sigma) and RHS",
        "requested_fgmres_rtol": args.rtol,
        "accept_residual": args.accept_residual,
        "accept_error": args.accept_error,
        "trials": results,
        "all_converged": all(x["converged"] for x in results),
        "limitations": [
            "Algebraic manufactured linear-solve comparison, not adaptive TS",
            "Triangular BJACOBI/ILU and LU-on-G are diagnostic baselines only",
            "Schur SELF-P uses MUMPS on approximate Schur, with configurable coenergy PC",
            "Timings include PC setup; no scaling or memory claims",
            "No Schur approximation/AMG configured; no production solver modified",
        ],
    }
    if rank == 0:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2)+"\n", encoding="utf-8")
        print("Saved:", args.output, flush=True)
    is_e.destroy()
    is_x.destroy()
    if not report["all_converged"]:
        raise RuntimeError("Block baseline did not meet the reference accuracy gate")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--mesh", type=float, default=.08)
    p.add_argument("--dt", type=float, default=.002)
    p.add_argument("--shift1", type=float, default=1000)
    p.add_argument("--shift2", type=float, default=2000)
    p.add_argument("--rtol", type=float, default=1e-12)
    p.add_argument("--max-it", type=int, default=400)
    p.add_argument("--e-pc", choices=("lu", "bjacobi"), default="lu",
                   help="Coenergy block PC: MUMPS LU or BJACOBI/ILU")
    p.add_argument("--x-pc", choices=("bjacobi", "lu", "schur-selfp"), default="bjacobi",
                   help="Complement PC: BJACOBI/ILU, LU on G, or Schur SELF-P + LU")
    p.add_argument("--accept-residual", type=float, default=5e-8)
    p.add_argument("--accept-error", type=float, default=5e-8)
    p.add_argument("--output", type=Path, required=True)
    run(p.parse_args())


if __name__ == "__main__":
    main()
