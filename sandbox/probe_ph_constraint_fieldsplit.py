# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See LICENSE for license information.
#
# github: https://github.com/g-haine/scrimp

"""Experimental graph-selected constraint Schur solve of the full pH-DAE.

Select every zero-self-block candidate from a distributed variable graph.
For each candidate, run FGMRES with PETSc full Schur factorization in the
(original) (principal, candidate) field partition. MUMPS on the principal
block and on PETSc FULL or SELF-P are controls, not a scalable implementation.
No production solver change, no name-specific rule and no dense gather.
Compare the true residual and complete state with full-system MUMPS.
"""
import argparse
import json
import os
from pathlib import Path
from time import perf_counter

from mpi4py import MPI
from petsc4py import PETSc

from bench_acoustic_2d import _build_wave
from audit_ph_block_graph import audit_blocks
from probe_constitutive_mumps_schur import _filled_right
from probe_ph_fieldsplit_baseline import solve_reference, relative_residual
from qualify_ph_constraint_schur import locally_owned_indices
from ph_schur_path2 import build_path2_schur
from scrimp.utils.ph_block_layout import discover_layout


def trial(A, rhs, iz, il, rtol, max_it, schur_pre, inner_restart_cap, repeat_rhs):
    ksp = PETSc.KSP().create(comm=PETSc.COMM_WORLD)
    ksp.setOptionsPrefix("ph_constraint_")
    ksp.setType("fgmres")
    ksp.setTolerances(rtol=rtol, atol=1e-13, max_it=max_it)
    ksp.setGMRESRestart(60)
    ksp.setErrorIfNotConverged(False)
    pc = ksp.getPC()
    pc.setType("fieldsplit")
    pc.setFieldSplitType(PETSc.PC.CompositeType.SCHUR)
    pc.setFieldSplitSchurFactType(PETSc.PC.FieldSplitSchurFactType.FULL)
    path2 = None
    path2_report = None
    if schur_pre == "path2":
        path2, path2_report = build_path2_schur(A, iz, il)
        if MPI.COMM_WORLD.rank == 0:
            print(
                f"PATH2 built: {json.dumps(path2_report, sort_keys=True)}",
                flush=True,
            )
        pc.setFieldSplitSchurPreType(
            PETSc.PC.FieldSplitSchurPreType.USER, path2
        )
    else:
        pc.setFieldSplitSchurPreType({
            "full": PETSc.PC.FieldSplitSchurPreType.FULL,
            "selfp": PETSc.PC.FieldSplitSchurPreType.SELFP,
            "self": PETSc.PC.FieldSplitSchurPreType.SELF,
        }[schur_pre])
    pc.setFieldSplitIS(("principal", iz), ("constraint", il))
    opts = PETSc.Options()
    options = {
        "ph_constraint_fieldsplit_principal_ksp_type": "preonly",
        "ph_constraint_fieldsplit_principal_pc_type": "lu",
        "ph_constraint_fieldsplit_principal_pc_factor_mat_solver_type": "mumps",
        "ph_constraint_fieldsplit_constraint_ksp_type": (
            "gmres" if schur_pre in ("self", "path2") else "preonly"
        ),
        "ph_constraint_fieldsplit_constraint_pc_type": (
            "none" if schur_pre == "self" else "lu"
        ),
    }
    inner_restart = min(int(il.getSize()) + 5, inner_restart_cap)
    if schur_pre in ("self", "path2"):
        options["ph_constraint_fieldsplit_constraint_ksp_gmres_restart"] = str(inner_restart)
        options["ph_constraint_fieldsplit_constraint_ksp_rtol"] = "1e-12"
        options["ph_constraint_fieldsplit_constraint_ksp_atol"] = "1e-14"
        options["ph_constraint_fieldsplit_constraint_ksp_max_it"] = str(max(100, 2 * inner_restart))
    else:
        options["ph_constraint_fieldsplit_constraint_pc_factor_mat_solver_type"] = "mumps"
    if schur_pre == "path2":
        options["ph_constraint_fieldsplit_constraint_pc_factor_mat_solver_type"] = "mumps"
    result = A.createVecRight()
    try:
        for key, val in options.items():
            opts[key] = val
        ksp.setOperators(A)
        ksp.setFromOptions()
        # Instrument the nested KSPs *before* solving. Count inner
        # iterations across all Schur applications, not only the last
        # inner solve, and report the number of subsolver invocations.
        setup_start = perf_counter()
        ksp.setUp()
        setup_seconds = MPI.COMM_WORLD.allreduce(perf_counter()-setup_start, op=MPI.MAX)
        sub_stats = {
            name: {"solve_calls": 0, "iterations_total": 0,
                   "max_iterations_per_solve": 0}
            for name in ("principal", "constraint")
        }
        for name, sub in zip(
            ("principal", "constraint"), pc.getFieldSplitSubKSP()
        ):
            def monitor(subksp, it, rnorm, name=name):
                stats = sub_stats[name]
                if it == 0:
                    stats["solve_calls"] += 1
                else:
                    stats["iterations_total"] += 1
                    stats["max_iterations_per_solve"] = max(
                        stats["max_iterations_per_solve"], int(it)
                    )
            sub.setMonitor(monitor)
        # Same assembled Jacobian, same KSP and PC, independent
        # manufactured RHS on subsequent solves. Reuse is intentional.
        # PETSc may defer numerical factorization to the first solve:
        # report warm and cold timings separately, not as pure solve cost.
        repetitions = []
        elapsed = None
        for repetition in range(repeat_rhs):
            manufactured = None
            this_rhs = rhs
            if repetition:
                manufactured = _filled_right(A, seed=101+repetition)
                this_rhs = A.createVecLeft()
                A.mult(manufactured, this_rhs)
            before = {
                name: dict(stat) for name, stat in sub_stats.items()
            }
            result.set(0.0)
            start = perf_counter()
            ksp.solve(this_rhs, result)
            seconds = MPI.COMM_WORLD.allreduce(
                perf_counter()-start, op=MPI.MAX
            )
            if elapsed is None:
                elapsed = seconds
                first_result = result.copy()
            sub_solvers = list(pc.getFieldSplitSubKSP())
            nested_ok = (
                len(sub_solvers) == 2
                and all(
                    int(sub.getConvergedReason()) > 0
                    and str(sub.getPC().getFailedReason()) == "0"
                    for sub in sub_solvers
                )
            )
            true_res = relative_residual(A, this_rhs, result)
            state_error = None
            if manufactured is not None:
                diff = result.copy()
                diff.axpy(-1.0, manufactured)
                state_error = (
                    float(diff.norm())
                    / max(float(manufactured.norm()), 1e-300)
                )
                diff.destroy()
            repetitions.append({
                "index": repetition,
                "seconds_max": seconds,
                "outer_reason": int(ksp.getConvergedReason()),
                "outer_iterations": int(ksp.getIterationNumber()),
                "nested_converged": nested_ok,
                "true_relative_residual": true_res,
                "manufactured_relative_state_error": state_error,
                "subsolver_calls_delta": {
                    name: sub_stats[name]["solve_calls"]
                          - before[name]["solve_calls"]
                    for name in sub_stats
                },
                "subsolver_iterations_delta": {
                    name: sub_stats[name]["iterations_total"]
                          - before[name]["iterations_total"]
                    for name in sub_stats
                },
            })
            if manufactured is not None:
                manufactured.destroy()
                this_rhs.destroy()
        result.destroy()
        result = first_result
        diagnostics = {
            "explicit_ksp_setup_seconds_max": setup_seconds,
            "repeated_rhs": repetitions,
            "outer_converged_reason": int(ksp.getConvergedReason()),
            "iterations": int(ksp.getIterationNumber()),
            "pc_failed_reason": str(pc.getFailedReason()),
            "subsolver_cumulative_monitor": sub_stats,
            "inner_gmres_restart_configured": inner_restart if schur_pre in ("self", "path2") else None,
            "sparse_path2_preconditioner": path2_report,
        }
        try:
            diagnostics["sub_solvers"] = [
                {
                    "name": name,
                    "pc_type": sub.getPC().getType(),
                    "reason": int(sub.getConvergedReason()),
                    "pc_failed_reason": str(sub.getPC().getFailedReason()),
                }
                for name, sub in zip(
                    ("principal", "constraint"), pc.getFieldSplitSubKSP()
                )
            ]
        except (AttributeError, RuntimeError, ValueError) as exc:
            diagnostics["sub_solver_inspection_error"] = str(exc)
        return result, elapsed, diagnostics
    except Exception:
        result.destroy()
        raise
    finally:
        ksp.destroy()
        if path2 is not None:
            path2.destroy()
        for key in options:
            del opts[key]


def run(args):
    mpi = MPI.COMM_WORLD
    if args.repeat_rhs < 1 or args.repeat_rhs > 20:
        raise ValueError("Require 1..20 manufactured RHS repetitions")
    if args.inner_restart_cap < 2:
        raise ValueError("Inner restart cap must be at least 2")
    if not (args.mesh > 0 and args.dt > 0 and
            args.shift1 > 0 and args.shift2 > 0 and
            args.shift1 != args.shift2 and args.rtol > 0):
        raise ValueError("Invalid mesh, time step, shifts or tolerance")
    if (os.getenv("SCRIMP_GETFEM_MATRIX_MODE") != "distributed" or
            os.getenv("SCRIMP_GETFEM_BRIDGE_DISTRIBUTED") != "coo"):
        raise RuntimeError("Distributed GetFEM COO assembly required")
    sys = _build_wave(args.mesh)
    sys.set_time_scheme(
        ts_type="cn", ts_adapt_type="none",
        ts_exact_final_time="matchstep",
        dt=args.dt, dt_save=args.dt, t_f=2*args.dt,
        ksp_type="preonly", pc_type="lu",
        pc_factor_mat_solver_type="mumps",
    )
    sys.solve()
    if (sys.last_ts_converged_reason <= 0 or sys.last_ts_steps != 2 or
            sys.last_ts_rejections or sys.last_ts_snes_failures):
        raise RuntimeError("Reference CN trajectory failed")
    if not sys.linear_mass or not sys.linear_stiffness:
        raise RuntimeError("This probe requires time-independent FE matrices")
    fields = discover_layout(
        sys.gf_model, states=sys.states, costates=sys.costates,
        ports=sys.ports, controls=sys.controls,
    )
    results = []
    all_pass = True
    for j, sigma in enumerate((args.shift1, args.shift2)):
        A = sys.stiffness.copy()
        A.axpy(sigma, sys.mass,
               structure=PETSc.Mat.Structure.DIFFERENT_NONZERO_PATTERN)
        A.assemble()
        graph = audit_blocks(A, fields)
        candidates = [v["name"] for v in graph["zero_self_block_candidates"]]
        if not candidates:
            raise RuntimeError("No qualified graph candidate; do not report vacuous PASS")
        sigma_runs = []
        for candidate in candidates:
            lo, hi = A.getOwnershipRange()
            local_z = locally_owned_indices(
                fields, lo, hi, {candidate}, invert=True
            )
            local_l = locally_owned_indices(
                fields, lo, hi, {candidate}, invert=False
            )
            if len(local_z)+len(local_l) != hi-lo:
                raise RuntimeError("Field split does not cover owned rows")
            iz = PETSc.IS().createGeneral(local_z, comm=PETSc.COMM_WORLD)
            il = PETSc.IS().createGeneral(local_l, comm=PETSc.COMM_WORLD)
            truth = _filled_right(A, seed=j+21)
            rhs = A.createVecLeft()
            A.mult(truth, rhs)
            ref, ref_t = solve_reference(A, rhs)
            computed, trial_t, details = trial(
                A, rhs, iz, il, args.rtol, args.max_it, args.schur_pre,
                args.inner_restart_cap, args.repeat_rhs
            )
            error = computed.copy()
            error.axpy(-1.0, ref)
            delta = float(error.norm())/max(float(ref.norm()), 1e-300)
            residual = relative_residual(A, rhs, computed)
            sub_solvers = details.get("sub_solvers", [])
            inner_converged = (
                len(sub_solvers) == 2
                and all(sub["reason"] > 0 and sub["pc_failed_reason"] == "0"
                        for sub in sub_solvers)
            )
            all_repeated_ok = all(
                entry["outer_reason"] > 0
                and entry["nested_converged"]
                and entry["true_relative_residual"] < args.accept_residual
                and (entry["manufactured_relative_state_error"] is None
                     or entry["manufactured_relative_state_error"] < args.accept_error)
                for entry in details["repeated_rhs"]
            )
            passed = (all_repeated_ok
                      and details["outer_converged_reason"] > 0
                      and inner_converged
                      and delta < args.accept_error
                      and residual < args.accept_residual)
            item = {
                "candidate": candidate,
                "schur_preconditioner": args.schur_pre,
                "constraint_dofs": int(
                    mpi.allreduce(len(local_l), op=MPI.SUM)),
                "principal_dofs": int(
                    mpi.allreduce(len(local_z), op=MPI.SUM)),
                "pass": bool(passed),
                "all_nested_ksp_converged": bool(inner_converged),
                "all_repeated_rhs_pass": bool(all_repeated_ok),
                "relative_state_error_vs_full_mumps": delta,
                "true_relative_residual": residual,
                "trial_seconds_max": trial_t,
                "full_mumps_seconds_max": ref_t,
                "ksp": details,
            }
            sigma_runs.append(item)
            if mpi.rank == 0:
                print(
                    f"NESTED sigma={sigma:g} candidate={candidate}: "
                    f"{json.dumps(details['subsolver_cumulative_monitor'], sort_keys=True)}",
                    flush=True,
                )
            all_pass = all_pass and passed
            if mpi.rank == 0:
                print(
                    f"{'PASS' if passed else 'FAIL'} sigma={sigma:g} "
                    f"candidate={candidate}: its={details['iterations']} "
                    f"reason={details['outer_converged_reason']} "
                    f"residual={residual:.3e} diff={delta:.3e} "
                    f"trial={trial_t:.3g}s ref={ref_t:.3g}s "
                    f"inner_restart={details['inner_gmres_restart_configured']} "
                    f"inner_reasons={[v['reason'] for v in sub_solvers]}",
                    flush=True,
                )
            for obj in (iz, il, truth, rhs, ref, computed, error):
                obj.destroy()
        results.append({"sigma": sigma, "candidates": sigma_runs})
        A.destroy()
    report = {
        "schema": "scrimp-ph-constraint-fieldsplit-v1",
        "source_commit": os.getenv("SCRIMP_BENCH_GIT_SHA", ""),
        "mpi_ranks": mpi.size,
        "mesh": args.mesh,
        "dofs": int(sys.gf_model.nbdof()),
        "fgmres_rtol": args.rtol,
        "repeat_rhs": args.repeat_rhs,
        "schur_preconditioner": args.schur_pre,
        "results": results,
        "all_pass": bool(all_pass),
        "limitations": [
            "Full original Jacobian solved on manufactured RHS only.",
            "All candidate blocks detected from graph; no name-specific rule.",
            "MUMPS on principal H; in SELF mode Schur is implicit with GMRES.",
            "FULL explicit Schur construction is diagnostic and previously errored in MUMPS.",
            "Nested monitor counts exclude setup and are per rank; report compares solve wall times.",
            "No production or time-trajectory FGMRES changes.",
            "No performance/scalability or physics claims from this probe.",
        ],
    }
    if mpi.rank == 0:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2)+"\n", encoding="utf-8")
        print("Saved:", args.output, flush=True)
    if not all_pass:
        raise RuntimeError("Graph-selected constraint Schur trial failed accuracy gate")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--mesh", type=float, default=0.08)
    p.add_argument("--dt", type=float, default=0.002)
    p.add_argument("--shift1", type=float, default=1000)
    p.add_argument("--shift2", type=float, default=2000)
    p.add_argument("--rtol", type=float, default=1e-12)
    p.add_argument("--max-it", type=int, default=400)
    p.add_argument("--repeat-rhs", type=int, default=1,
                   help="Number of independent manufactured solves reusing the same KSP/PC")
    p.add_argument("--inner-restart-cap", type=int, default=256,
                   help="Cap for dimension-aware GMRES Schur restart; diagnostic only")
    p.add_argument("--schur-pre", choices=("self", "path2", "full", "selfp"), default="self",
                   help="SELF baseline, PATH2 sparse Schur USER Pmat, FULL/SELFP controls")
    p.add_argument("--accept-error", type=float, default=5e-8)
    p.add_argument("--accept-residual", type=float, default=5e-8)
    p.add_argument("--output", type=Path, required=True)
    run(p.parse_args())


if __name__ == "__main__":
    main()
