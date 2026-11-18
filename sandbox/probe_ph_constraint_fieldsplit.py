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
block and on PETSc SELF-P are controls, not a scalable implementation.
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
from scrimp.utils.ph_block_layout import discover_layout


def trial(A, rhs, iz, il, rtol, max_it):
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
    pc.setFieldSplitSchurPreType(PETSc.PC.FieldSplitSchurPreType.SELFP)
    pc.setFieldSplitIS(("principal", iz), ("constraint", il))
    opts = PETSc.Options()
    options = {
        "ph_constraint_fieldsplit_principal_ksp_type": "preonly",
        "ph_constraint_fieldsplit_principal_pc_type": "lu",
        "ph_constraint_fieldsplit_principal_pc_factor_mat_solver_type": "mumps",
        "ph_constraint_fieldsplit_constraint_ksp_type": "preonly",
        "ph_constraint_fieldsplit_constraint_pc_type": "lu",
        "ph_constraint_fieldsplit_constraint_pc_factor_mat_solver_type": "mumps",
    }
    result = A.createVecRight()
    try:
        for key, val in options.items():
            opts[key] = val
        ksp.setOperators(A)
        ksp.setFromOptions()
        start = perf_counter()
        ksp.solve(rhs, result)
        elapsed = MPI.COMM_WORLD.allreduce(
            perf_counter()-start, op=MPI.MAX
        )
        diagnostics = {
            "outer_converged_reason": int(ksp.getConvergedReason()),
            "iterations": int(ksp.getIterationNumber()),
            "pc_failed_reason": str(pc.getFailedReason()),
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
        for key in options:
            del opts[key]


def run(args):
    mpi = MPI.COMM_WORLD
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
                A, rhs, iz, il, args.rtol, args.max_it
            )
            error = computed.copy()
            error.axpy(-1.0, ref)
            delta = float(error.norm())/max(float(ref.norm()), 1e-300)
            residual = relative_residual(A, rhs, computed)
            passed = (details["outer_converged_reason"] > 0
                      and delta < args.accept_error
                      and residual < args.accept_residual)
            item = {
                "candidate": candidate,
                "constraint_dofs": int(
                    mpi.allreduce(len(local_l), op=MPI.SUM)),
                "principal_dofs": int(
                    mpi.allreduce(len(local_z), op=MPI.SUM)),
                "pass": bool(passed),
                "relative_state_error_vs_full_mumps": delta,
                "true_relative_residual": residual,
                "trial_seconds_max": trial_t,
                "full_mumps_seconds_max": ref_t,
                "ksp": details,
            }
            sigma_runs.append(item)
            all_pass = all_pass and passed
            if mpi.rank == 0:
                print(
                    f"{'PASS' if passed else 'FAIL'} sigma={sigma:g} "
                    f"candidate={candidate}: its={details['iterations']} "
                    f"reason={details['outer_converged_reason']} "
                    f"residual={residual:.3e} diff={delta:.3e} "
                    f"trial={trial_t:.3g}s ref={ref_t:.3g}s",
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
        "results": results,
        "all_pass": bool(all_pass),
        "limitations": [
            "Full original Jacobian solved on manufactured RHS only.",
            "All candidate blocks detected from graph; no name-specific rule.",
            "MUMPS used on principal H and on approximate Schur SELF-P.",
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
    p.add_argument("--accept-error", type=float, default=5e-8)
    p.add_argument("--accept-residual", type=float, default=5e-8)
    p.add_argument("--output", type=Path, required=True)
    run(p.parse_args())


if __name__ == "__main__":
    main()
