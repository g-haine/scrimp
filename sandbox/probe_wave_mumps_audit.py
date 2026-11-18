# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# See LICENSE for license information.
"""Audit PETSc/MUMPS on the complete wave CN Jacobian (not principal H).

Separate KSPSetUp (which includes LU symbolic/numeric in this case)
from the first and repeated triangular solves. Capture effective MUMPS
ICNTL/INFOG/RINFOG and PETSc log events in text logs. Ordering variants
run in SEPARATE MPI jobs so no factor/option state leaks across tests.

No production solver changes. Reference state is manufactured x*, with
b=A*x*, and both residual and state error are checked.
"""
import argparse
import json
import os
import resource
from pathlib import Path
from time import perf_counter

from mpi4py import MPI
from petsc4py import PETSc

from bench_acoustic_2d import _build_wave
from probe_constitutive_mumps_schur import _filled_right
from replay_petsc_mumps_snapshot import load_matrix


def interval(start):
    return float(MPI.COMM_WORLD.allreduce(perf_counter()-start, op=MPI.MAX))


def rss():
    local = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return {
        "peak_per_process_kib_max": int(MPI.COMM_WORLD.allreduce(local, op=MPI.MAX)),
        "sum_process_peak_kib": int(MPI.COMM_WORLD.allreduce(local, op=MPI.SUM)),
    }


def build_operator(args):
    if args.matrix:
        matrix = load_matrix(args.matrix, PETSc.COMM_WORLD)
        if matrix.getSize()[0] != matrix.getSize()[1]:
            raise ValueError("Snapshot matrix must be square")
        return matrix, "petsc_binary_snapshot"
    if (os.getenv("SCRIMP_GETFEM_MATRIX_MODE") != "distributed"
            or os.getenv("SCRIMP_GETFEM_BRIDGE_DISTRIBUTED") != "coo"):
        raise RuntimeError("This audit requires distributed COO GetFEM bridge")
    system = _build_wave(args.mesh)
    system.set_time_scheme(
        ts_type="cn", ts_adapt_type="none",
        ts_exact_final_time="matchstep",
        dt=0.002, dt_save=0.002, t_f=0.004,
        ksp_type="preonly", pc_type="lu",
        pc_factor_mat_solver_type="mumps",
    )
    system.solve()
    if (system.last_ts_converged_reason <= 0
            or system.last_ts_steps != 2
            or system.last_ts_rejections or system.last_ts_snes_failures):
        raise RuntimeError("Reference CN preparation failed")
    matrix = system.stiffness.copy()
    matrix.axpy(
        args.shift, system.mass,
        structure=PETSc.Mat.Structure.DIFFERENT_NONZERO_PATTERN,
    )
    matrix.assemble()
    return matrix, "assembled_wave_cn_operator"


def stats(factor):
    output = {}
    for name, indices, method in (
        ("icntl", (7, 20, 21, 28, 29), "getMumpsIcntl"),
        ("infog", (1, 2, 3, 4, 5, 6, 7, 8, 9, 12, 16, 17, 18, 19, 20,
                   21, 22, 23, 24, 25, 27, 28, 29), "getMumpsInfog"),
        ("rinfog", (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12), "getMumpsRinfog"),
    ):
        values = {}
        getter = getattr(factor, method, None)
        if getter is None:
            output[name+"_unsupported"] = method
            continue
        for index in indices:
            try:
                values[str(index)] = getter(index)
            except (RuntimeError, PETSc.Error, ValueError) as exc:
                values[str(index)] = "unavailable: " + str(exc)[:120]
        output[name] = values
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mesh", type=float, default=0.025)
    parser.add_argument("--shift", type=float, default=1000)
    parser.add_argument("--matrix", type=str,
                        help="Optional PETSc binary snapshot instead of wave build")
    parser.add_argument("--ordering", choices=("default", "sequential", "parallel"),
                        default="default")
    parser.add_argument("--repeat", type=int, default=4)
    parser.add_argument("--tol", type=float, default=5e-8)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    PETSc.Log.begin()
    comm = MPI.COMM_WORLD
    if args.mesh <= 0 or args.repeat < 2 or args.repeat > 12:
        raise ValueError("mesh > 0 and repeat in [2,12] required")
    # PETSc MUMPS options control the actual factor Mat. Crucially these
    # must be applied BEFORE factor creation and NOT recycled between jobs.
    opts = PETSc.Options()
    if args.ordering == "sequential":
        opts["mat_mumps_icntl_28"] = "1"
    elif args.ordering == "parallel":
        opts["mat_mumps_icntl_28"] = "2"
        opts["mat_mumps_icntl_29"] = "0"  # request MUMPS automatic parallel ordering

    start = perf_counter()
    matrix, origin = build_operator(args)
    build_s = interval(start)
    if matrix.getSize()[0] != matrix.getSize()[1]:
        raise RuntimeError("Non-square operator")
    before_factor = rss()
    ksp = PETSc.KSP().create(comm=PETSc.COMM_WORLD)
    ksp.setType("preonly")
    ksp.setErrorIfNotConverged(False)
    pc = ksp.getPC()
    pc.setType("lu")
    pc.setFactorSolverType("mumps")
    ksp.setOperators(matrix)
    ksp.setFromOptions()
    report = {
        "schema": "scrimp-wave-mumps-diagnostics-v1",
        "commit": os.getenv("SCRIMP_BENCH_GIT_SHA", ""),
        "mpi_ranks": comm.size,
        "mesh": args.mesh, "matrix_source": origin,
        "ordering_requested": args.ordering,
        "matrix_size": list(matrix.getSize()),
        "matrix_nnz_used_sum": int(comm.allreduce(int(matrix.getInfo()["nz_used"]), op=MPI.SUM)),
        "wave_preparation_seconds_max": build_s,
        "peak_rss_before_factor": before_factor,
        "factor_setup_seconds_max": None,
        "mumps_after_setup": None,
        "solves": [], "all_pass": False,
        "limitations": [
            "KSPSetUp includes MUMPS symbolic/numeric phases, not phase-resolved.",
            "For phase-resolved PETSc times use -log_view from the emitted PETSc log.",
            "Process RSS high-water includes wave setup, not isolated factor allocations.",
            "Full-system manufactured solves, not a time-integration speedup claim.",
            "Parallel ordering availability depends on container MUMPS build.",
            "Unhonored ICNTL(28) request is NON-PASS despite accurate solves.",
        ],
    }
    try:
        stage_factor = PETSc.Log.Stage("MUMPS_AUDIT_FACTOR")
        stage_factor.push()
        start = perf_counter()
        try:
            ksp.setUp()
            report["factor_setup_seconds_max"] = interval(start)
        finally:
            stage_factor.pop()
        factor = pc.getFactorMatrix()
        report["mumps_after_setup"] = stats(factor)
        actual28 = report["mumps_after_setup"].get("icntl", {}).get("28")
        expected28 = {"sequential": 1, "parallel": 2}.get(args.ordering)
        report["requested_icntl_28"] = expected28
        report["effective_icntl_28"] = actual28
        report["ordering_request_honored"] = (
            expected28 is None or actual28 == expected28
        )
        if comm.rank == 0:
            print(
                f"ORDERING-CHECK requested={args.ordering} "
                f"ICNTL28={actual28} expected={expected28} "
                f"honored={report['ordering_request_honored']}",
                flush=True,
            )
        stage_solve = PETSc.Log.Stage("MUMPS_AUDIT_SOLVES")
        stage_solve.push()
        for rep in range(args.repeat):
            truth = _filled_right(matrix, seed=501+rep)
            rhs = matrix.createVecLeft()
            solution = matrix.createVecRight()
            matrix.mult(truth, rhs)
            start = perf_counter()
            ksp.solve(rhs, solution)
            elapsed = interval(start)
            residual_vec = matrix.createVecLeft()
            matrix.mult(solution, residual_vec)
            residual_vec.axpy(-1.0, rhs)
            delta = solution.copy()
            delta.axpy(-1.0, truth)
            error = float(delta.norm()) / max(float(truth.norm()), 1e-300)
            residual = float(residual_vec.norm()) / max(float(rhs.norm()), 1e-300)
            reason = int(ksp.getConvergedReason())
            passed = reason > 0 and error < args.tol and residual < args.tol
            report["solves"].append({
                "index": rep, "seconds_max": elapsed,
                "converged_reason": reason, "state_relative_error": error,
                "true_relative_residual": residual, "pass": bool(passed),
            })
            for obj in (truth, rhs, solution, residual_vec, delta):
                obj.destroy()
        stage_solve.pop()
        report["peak_rss_after_solves"] = rss()
        report["all_pass"] = (
            all(row["pass"] for row in report["solves"])
            and report["ordering_request_honored"]
        )
    finally:
        ksp.destroy()
        matrix.destroy()
    if comm.rank == 0:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2)+"\n", encoding="utf-8")
        print(
            "MUMPS-AUDIT ordering={ordering} ranks={ranks} "
            "dofs={dofs} setup={setup:.6g}s solves={solves} "
            "pass={passed}".format(
                ordering=args.ordering, ranks=comm.size,
                dofs=report["matrix_size"][0],
                setup=report["factor_setup_seconds_max"],
                solves=[round(x["seconds_max"], 6) for x in report["solves"]],
                passed=report["all_pass"],
            ), flush=True,
        )
        print("Saved:", args.output, flush=True)
    # Built-in PETSc log separates symbolic/numeric LU and MatSolve events
    # within the tagged audit stages. Captured by the Docker run .log.
    PETSc.Log.view()
    if not report["all_pass"]:
        raise RuntimeError("MUMPS accuracy gate failed or requested ordering was not honored")


if __name__ == "__main__":
    main()
