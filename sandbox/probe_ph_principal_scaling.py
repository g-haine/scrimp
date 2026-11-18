# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See LICENSE for license information.
#
# github: https://github.com/g-haine/scrimp

"""MPI strong-scaling probe for the graph-selected principal DAE block H.

Compare MUMPS global LU with MPI-local BJACOBI/ILU and ASM/ILU, on
IDENTICAL manufactured RHS, preserving a single KSP and PC across
repeated solves. Report setup/cold/warm times separately, true relative
residual and manufactured solution error, plus process peak RSS.

This is deliberately a standalone diagnostic of H, not a substitute
for the complete pH-DAE solver. No general claim of admissibility of
local ILU for indefinite H is made. Failed methods remain in the JSON.
"""
import argparse
import json
import os
import resource
from pathlib import Path
from time import perf_counter

from mpi4py import MPI
from petsc4py import PETSc

from audit_ph_block_graph import audit_blocks
from bench_acoustic_2d import _build_wave
from probe_constitutive_mumps_schur import _filled_right
from qualify_ph_constraint_schur import locally_owned_indices
from scrimp.utils.ph_block_layout import discover_layout


def max_seconds(start, comm):
    return float(comm.allreduce(perf_counter() - start, op=MPI.MAX))


def memory_kib(comm):
    local = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return {
        "process_peak_rss_kib_max": comm.allreduce(local, op=MPI.MAX),
        "process_peak_rss_kib_sum": comm.allreduce(local, op=MPI.SUM),
        "warning": "Linux ru_maxrss high-water includes GetFEM build and reference CN solve",
    }


def run_method(H, method, *, repeat, rtol, max_it, tol):
    comm = MPI.COMM_WORLD
    report = {
        "method": method, "setup_seconds_max": None,
        "repetitions": [], "all_pass": False,
        "memory_after_method": None,
    }
    ksp = PETSc.KSP().create(comm=PETSc.COMM_WORLD)
    ksp.setOptionsPrefix("h_scaling_")
    ksp.setErrorIfNotConverged(False)
    pc = ksp.getPC()
    if method == "mumps":
        ksp.setType("preonly")
        pc.setType("lu")
        pc.setFactorSolverType("mumps")
    else:
        ksp.setType("fgmres")
        ksp.setTolerances(rtol=rtol, atol=1e-14, max_it=max_it)
        ksp.setGMRESRestart(min(100, max_it))
        pc.setType("bjacobi" if method == "bjacobi" else "asm")
        if method == "asm":
            pc.setASMOverlap(1)
        # PETSc on MPIAIJ uses local sequential submatrices for each
        # subdomain; explicit ILU only on those local submatrices.
        opts = PETSc.Options()
        opts["h_scaling_sub_ksp_type"] = "preonly"
        opts["h_scaling_sub_pc_type"] = "ilu"
    try:
        ksp.setOperators(H)
        start = perf_counter()
        ksp.setUp()
        report["setup_seconds_max"] = max_seconds(start, comm)
        for j in range(repeat):
            truth = _filled_right(H, seed=31+j)
            rhs = H.createVecLeft()
            H.mult(truth, rhs)
            found = H.createVecRight()
            try:
                start = perf_counter()
                error_message = None
                try:
                    ksp.solve(rhs, found)
                except PETSc.Error as exc:
                    error_message = str(exc)
                seconds = max_seconds(start, comm)
                reason = int(ksp.getConvergedReason())
                if error_message is None:
                    diff = found.copy()
                    diff.axpy(-1.0, truth)
                    state_error = float(diff.norm()) / max(float(truth.norm()), 1e-300)
                    diff.destroy()
                    check = H.createVecLeft()
                    H.mult(found, check)
                    check.axpy(-1.0, rhs)
                    residual = float(check.norm()) / max(float(rhs.norm()), 1e-300)
                    check.destroy()
                else:
                    state_error = None
                    residual = None
                passed = (
                    error_message is None and reason > 0
                    and state_error < tol and residual < tol
                )
                report["repetitions"].append({
                    "index": j, "seconds_max": seconds,
                    "ksp_iterations": int(ksp.getIterationNumber()),
                    "converged_reason": reason,
                    "pc_failed_reason": str(pc.getFailedReason()),
                    "true_relative_residual": residual,
                    "manufactured_relative_error": state_error,
                    "pass": bool(passed), "exception": error_message,
                })
                # Avoid repeated failures with an unusable local PC.
                if error_message is not None or reason < 0:
                    break
            finally:
                truth.destroy()
                rhs.destroy()
                found.destroy()
        report["all_pass"] = (
            len(report["repetitions"]) == repeat
            and all(item["pass"] for item in report["repetitions"])
        )
        report["memory_after_method"] = memory_kib(comm)
        return report
    finally:
        ksp.destroy()
        if method != "mumps":
            for key in ("h_scaling_sub_ksp_type", "h_scaling_sub_pc_type"):
                del PETSc.Options()[key]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--mesh", type=float, default=0.05)
    ap.add_argument("--shift", type=float, default=1000)
    ap.add_argument("--repeat", type=int, default=4)
    ap.add_argument("--rtol", type=float, default=1e-10)
    ap.add_argument("--tol", type=float, default=5e-8)
    ap.add_argument("--max-it", type=int, default=400)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    comm = MPI.COMM_WORLD
    if args.mesh <= 0 or args.shift <= 0 or not 1 <= args.repeat <= 20:
        raise ValueError("Invalid mesh/shift/repetitions")
    if (os.getenv("SCRIMP_GETFEM_MATRIX_MODE") != "distributed"
            or os.getenv("SCRIMP_GETFEM_BRIDGE_DISTRIBUTED") != "coo"):
        raise RuntimeError("Distributed GetFEM-PETSc COO bridge required")

    system = _build_wave(args.mesh)
    system.set_time_scheme(
        ts_type="cn", ts_adapt_type="none",
        ts_exact_final_time="matchstep",
        dt=0.002, dt_save=0.002, t_f=0.004,
        ksp_type="preonly", pc_type="lu",
        pc_factor_mat_solver_type="mumps",
    )
    system.solve()
    if system.last_ts_converged_reason <= 0 or system.last_ts_steps != 2:
        raise RuntimeError("CN/MUMPS reference failed")
    fields = discover_layout(
        system.gf_model, states=system.states,
        costates=system.costates, ports=system.ports,
        controls=system.controls,
    )
    A = system.stiffness.copy()
    A.axpy(
        args.shift, system.mass,
        structure=PETSc.Mat.Structure.DIFFERENT_NONZERO_PATTERN,
    )
    A.assemble()
    graph = audit_blocks(A, fields)
    candidates = [c["name"] for c in graph["zero_self_block_candidates"]]
    if not candidates:
        raise RuntimeError("No candidate from graph; cannot define H")
    if len(candidates) != 1:
        raise RuntimeError(
            "Multiple candidates require graph aggregation; refuse arbitrary H"
        )
    lo, hi = A.getOwnershipRange()
    local_z = locally_owned_indices(
        fields, lo, hi, {candidates[0]}, invert=True
    )
    iz = PETSc.IS().createGeneral(local_z, comm=PETSc.COMM_WORLD)
    H = A.createSubMatrix(iz, iz)
    results = []
    for method in ("mumps", "bjacobi", "asm"):
        result = run_method(
            H, method, repeat=args.repeat, rtol=args.rtol,
            max_it=args.max_it, tol=args.tol,
        )
        results.append(result)
        if comm.rank == 0:
            reps = result["repetitions"]
            print(
                f"METHOD {method}: pass={result['all_pass']} "
                f"setup={result['setup_seconds_max']:.5g}s "
                f"times={[round(t['seconds_max'], 6) for t in reps]} "
                f"iterations={[t['ksp_iterations'] for t in reps]} "
                f"reasons={[t['converged_reason'] for t in reps]}",
                flush=True,
            )
    report = {
        "schema": "scrimp-ph-principal-H-scaling-v1",
        "source_commit": os.getenv("SCRIMP_BENCH_GIT_SHA", ""),
        "mpi_ranks": comm.size, "mesh_h": args.mesh, "sigma": args.shift,
        "total_dofs": int(A.getSize()[0]), "principal_dofs": int(H.getSize()[0]),
        "constraint_candidate": candidates[0],
        "repeat_rhs": args.repeat,
        "methods": results,
        "warning": [
            "The benchmark includes no full-system iterative pH time integration.",
            "MPI strong scaling must compare the same mesh and identical RHS.",
            "ILU(BJACOBI/ASM) may fail for an indefinite DAE principal block.",
            "RSS is process high-water, not isolated solver allocations.",
            "Production reference remains unchanged.",
        ],
    }
    if comm.rank == 0:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2)+"\n", encoding="utf-8")
        print("Saved:", args.output, flush=True)
    H.destroy()
    iz.destroy()
    A.destroy()
    # Do not report a globally successful scalable solver just because
    # the independent direct MUMPS control passed.
    if not results[0]["all_pass"]:
        raise RuntimeError("H MUMPS reference did not pass")


if __name__ == "__main__":
    main()
