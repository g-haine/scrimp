# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Version-neutral acoustic wave benchmark on one fixed physical problem.

This standalone script lives in the CURRENT sandbox, but loads SCRIMP
and examples/wave.py from the selected /workspace/scrimp mount. For
historical v1.2.0 it never changes the solver, Jacobian, time
integration, or PETSc/MUMPS reset behavior.

Full run wall times (from just before set_time_scheme/solve until
after solve) and individual accepted CN-step times are measured.
For a step n>0, we time from AFTER monitor(n-1) to BEFORE monitor(n),
excluding result-writing/monitor overhead. This separates first-step
factorization from steady iterations. All rank-local times are reduced
AFTER solve; no MPI timing collectives occur inside PETSc callbacks.
The true t_final and complete distributed PETSc Vec are collected
once at the beginning and end for numerical parity.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import resource
import time

from mpi4py import MPI
from petsc4py import PETSc
import numpy as np
import scrimp as S

from bench_acoustic_2d import _build_wave


def _gather_global(local, size, mpi):
    """Collect once after the timed solve, NOT for every timestep."""
    partitions = mpi.gather(local, root=0)
    if mpi.Get_rank() != 0:
        return None
    result = np.empty(size, dtype=np.float64)
    end = 0
    for lo, hi, values in sorted(partitions, key=lambda p: p[0]):
        if lo != end or len(values) != hi - lo:
            raise ValueError("Invalid PETSc final-state ownership")
        result[lo:hi] = values
        end = hi
    if end != size:
        raise ValueError("Incomplete PETSc final-state ownership")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--variant", choices=("release_120", "current"), required=True
    )
    parser.add_argument(
        "--policy", choices=("release_native", "never", "on_change"),
        required=True,
    )
    parser.add_argument("--mesh", type=float, default=0.08)
    parser.add_argument("--dt", type=float, default=0.002)
    parser.add_argument("--final-time", type=float, default=0.04)
    parser.add_argument("--save-dt", type=float, default=0.01)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--source-sha", required=True)
    parser.add_argument("--source-tag", default="")
    parser.add_argument("--repetition", type=int, default=1)
    parser.add_argument("--profile", action="store_true",
                        help="Current revision: SCRIMP phase timings")
    parser.add_argument("--petsc-log", action="store_true",
                        help="Current revision: native PETSc staged log")
    args = parser.parse_args()
    mpi = MPI.COMM_WORLD
    rank, ranks = mpi.Get_rank(), mpi.Get_size()
    if args.variant == "release_120" and (args.profile or args.petsc_log):
        parser.error("legacy release v1.2.0 uses only neutral per-step timings")
    if args.profile:
        os.environ["SCRIMP_PROFILE"] = "1"
    if args.petsc_log:
        PETSc.Log.begin()

    if not all(x > 0 for x in (
        args.mesh, args.dt, args.final_time, args.save_dt
    )):
        parser.error("Mesh, dt, final time and save cadence must be positive")
    target_steps = int(round(args.final_time / args.dt))
    if abs(target_steps * args.dt - args.final_time) > 1e-12:
        parser.error("T must be an integer number of fixed CN steps")
    if target_steps < 5:
        parser.error(
            "At least five fixed CN steps are required to separate "
            "startup from meaningful early/late steady-step windows"
        )
    if args.variant == "release_120":
        if ranks != 1 or args.policy != "release_native":
            parser.error("v1.2.0 is ONLY benchmarked on 1 MPI rank")
    else:
        if args.policy not in ("never", "on_change"):
            parser.error("Current COO requires never or on_change")
        if os.environ.get("SCRIMP_GETFEM_MATRIX_MODE") != "distributed":
            raise RuntimeError("Current COO must use distributed GetFEM")
        if os.environ.get(
            "SCRIMP_GETFEM_BRIDGE_DISTRIBUTED"
        ) != "coo":
            raise RuntimeError("Current benchmark must use COO bridge")
        if os.environ.get(
            "SCRIMP_MUMPS_RESET_POLICY"
        ) != args.policy:
            raise RuntimeError("Unexpected current MUMPS policy")
    if not Path(S.__file__).resolve().as_posix().startswith(
        "/workspace/scrimp/scrimp/"
    ):
        raise RuntimeError(
            f"SCRIMP not imported from the selected Docker source: {S.__file__}"
        )
    mpi.Barrier()
    build_start = time.perf_counter()
    system = _build_wave(args.mesh)
    mpi.Barrier()
    build_stop = time.perf_counter()
    n_global = int(system.gf_model.nbdof())

    observations = {"initial": None, "final": None, "steps": None,
                    "rejected": None, "snes_failures": None,
                    "last_callback_time": None}
    original_monitor = system.monitor
    # Each rank records intervals locally. The previous monitor exit
    # and the next monitor entry bracket only PETSc step work, not
    # archival, output or the user's monitor callback itself.
    step_clock = {
        "previous_monitor_exit_ns": None,
        "indices": [], "times": [], "durations_s": [],
        "monitor_durations_s": [],
    }

    def measured_monitor(TS, i, t, z, *a, **kw):
        begin_ns = time.perf_counter_ns()
        is_main = not kw.get("initial_step", False) and i >= 0
        if is_main and i > 0:
            previous_end = step_clock["previous_monitor_exit_ns"]
            if previous_end is None:
                raise RuntimeError(
                    "Main TS step monitor lacks its step-zero boundary"
                )
            elapsed = (begin_ns - previous_end) * 1e-9
            step_clock["indices"].append(int(i))
            step_clock["times"].append(float(t))
            step_clock["durations_s"].append(float(elapsed))
        # The initial-consistency TS uses `initial_step=True` in both
        # versions. Do not count that pseudostep as a physical CN step.
        if not kw.get("initial_step", False) and i >= 0:
            observations["steps"] = int(TS.getStepNumber())
            observations["rejected"] = int(TS.getStepRejections())
            observations["snes_failures"] = int(TS.getSNESFailures())
            observations["last_callback_time"] = float(t)
            if i == 0 or (
                abs(float(t) - args.final_time) <= 1e-10
            ):
                lo, hi = z.getOwnershipRange()
                local = (
                    int(lo), int(hi),
                    np.asarray(z.getArray(readonly=True).real,
                               dtype=np.float64).copy(),
                )
                if i == 0:
                    observations["initial"] = local
                if abs(float(t) - args.final_time) <= 1e-10:
                    observations["final"] = local
        try:
            return original_monitor(TS, i, t, z, *a, **kw)
        finally:
            if is_main:
                end_ns = time.perf_counter_ns()
                step_clock["monitor_durations_s"].append(
                    (end_ns - begin_ns) * 1e-9
                )
                step_clock["previous_monitor_exit_ns"] = end_ns

    system.monitor = measured_monitor
    if args.petsc_log:
        # PETSc events separated by initial consistency and CN main TS.
        initial_stage = PETSc.Log.Stage("WAVE_INITIAL_CONSISTENCY")
        main_stage = PETSc.Log.Stage("WAVE_MAIN_CN")
        original_init_step = system.init_step
        original_memory_phase = system._memory_phase
        in_main = [False]

        def timed_init(*a, **kw):
            initial_stage.push()
            try:
                return original_init_step(*a, **kw)
            finally:
                initial_stage.pop()

        def timed_memory_phase(name, *a, **kw):
            if name == "before_main_ts_solve":
                main_stage.push()
                in_main[0] = True
            try:
                return original_memory_phase(name, *a, **kw)
            finally:
                if name == "after_main_ts_solve" and in_main[0]:
                    main_stage.pop()
                    in_main[0] = False

        system.init_step = timed_init
        system._memory_phase = timed_memory_phase

    # PETSc TS sees precisely the same fixed CN grid in both releases.
    mpi.Barrier()
    start = time.perf_counter()
    system.set_time_scheme(
        ts_type="cn", ts_adapt_type="none",
        ts_exact_final_time="matchstep",
        dt=args.dt, dt_save=args.save_dt,
        ts_adapt_dt_max=args.dt,
        t_f=args.final_time,
        ksp_type="preonly", pc_type="lu",
        pc_factor_mat_solver_type="mumps",
    )
    system.solve()
    mpi.Barrier()
    end = time.perf_counter()
    phase_profile = system.get_profile() if args.profile else None

    # Legacy SCRIMP v1.2.0 has no `last_ts_*` metadata; use the actual
    # standard PETSc TS.monitor as evidence, not a fictitious reason code.
    complete = (
        observations["final"] is not None
        and observations["steps"] == target_steps
        and observations["rejected"] == 0
        and observations["snes_failures"] == 0
        and abs(
            observations["last_callback_time"] - args.final_time
        ) <= 1e-10
    )
    if not all(mpi.allgather(bool(complete))):
        raise RuntimeError(
            "Physical time integration did not complete the required "
            f"{target_steps} fixed CN steps at t={args.final_time}; "
            f"last monitor metadata={observations}"
        )
    if args.variant == "current":
        if (getattr(system, "last_ts_converged_reason", 0) <= 0
                or getattr(system, "last_ts_steps", -1) != target_steps
                or getattr(system, "_mumps_reset_policy", None)
                != args.policy):
            raise RuntimeError("Current PETSc convergence/policy gate failed")

    # All accepted steps, across ALL ranks, must have a measured time
    # interval. Validate before gathering. Never use Tsolve/N as an
    # approximation for true post-startup per-step solve time.
    local_ok = (
        step_clock["indices"] == list(range(1, target_steps + 1))
        and len(step_clock["durations_s"]) == target_steps
        and len(step_clock["monitor_durations_s"]) == target_steps + 1
        and all(
            np.isfinite(step) and step > 0
            for step in step_clock["durations_s"]
        )
        and all(
            abs(float(t) - i * args.dt) <= 1e-9
            for i, t in zip(
                step_clock["indices"], step_clock["times"]
            )
        )
    )
    if not all(mpi.allgather(bool(local_ok))):
        raise RuntimeError(
            "Incomplete or inconsistent main-TS per-step timing: "
            + repr(step_clock)
        )
    rank_step_times = mpi.gather(
        step_clock["durations_s"], root=0
    )
    rank_monitor_times = mpi.gather(
        step_clock["monitor_durations_s"], root=0
    )

    initial = _gather_global(observations["initial"], n_global, mpi)
    final = _gather_global(observations["final"], n_global, mpi)
    build_wall = mpi.allreduce(
        build_stop - build_start, op=MPI.MAX
    )
    solve_wall = mpi.allreduce(end - start, op=MPI.MAX)
    max_rss = mpi.allreduce(
        int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss),
        op=MPI.MAX,
    )
    if rank == 0:
        rank_step_matrix = np.asarray(
            rank_step_times, dtype=np.float64
        )
        if rank_step_matrix.shape != (ranks, target_steps):
            raise RuntimeError("Incomplete per-rank accepted-step timings")
        # A distributed timestep completes only when the slowest rank
        # is ready; each step's wall is therefore max over ranks.
        step_wall = np.max(rank_step_matrix, axis=0)
        if not np.all(np.isfinite(step_wall)) or np.any(step_wall <= 0):
            raise RuntimeError("Invalid global PETSc step time")
        steady = step_wall[1:]  # exclude first factorization/warmup
        half = len(steady) // 2
        early = steady[:half]
        late = steady[-half:]
        args.out.parent.mkdir(parents=True, exist_ok=True)
        if not np.all(np.isfinite(initial)) or not np.all(np.isfinite(final)):
            raise RuntimeError("Nonfinite initial/final state")
        state_file = args.out.with_name(args.out.stem + "_state.npz")
        np.savez_compressed(
            state_file,
            initial=initial, final=final,
        )
        report = {
            "benchmark": "examples/wave.py:wave_eq()",
            "variant": args.variant,
            "source_sha": args.source_sha,
            "source_tag": args.source_tag,
            "repetition": args.repetition,
            "mpi_ranks": ranks, "mesh_h": args.mesh,
            "dofs": n_global, "scheme": "cn",
            "ts_adapt_type": "none",
            "dt": args.dt, "time_final_requested": args.final_time,
            "time_final_monitor": observations["last_callback_time"],
            "save_dt": args.save_dt,
            "steps": observations["steps"],
            "rejections": observations["rejected"],
            "snes_failures": observations["snes_failures"],
            "ts_convergence_reason":
                int(getattr(system, "last_ts_converged_reason"))
                if hasattr(system, "last_ts_converged_reason") else None,
            "matrix_backend": (
                "release_v1.2.0_native"
                if args.variant == "release_120" else "distributed_coo"
            ),
            "mumps_policy": args.policy,
            "matrix_memory_mode": (
                os.environ.get("SCRIMP_MATRIX_MEMORY_MODE", "legacy")
                if args.variant == "current" else "release_native"
            ),
            "trajectory_backend": "memory",
            "profile": phase_profile,
            "wall_build_s": float(build_wall),
            "wall_solve_s": float(solve_wall),
            "wall_total_s": float(build_wall + solve_wall),
            "step_timing_method": (
                "pre-monitor(n) minus post-monitor(n-1), in monotonic "
                "wall time per rank; maximum over MPI ranks per step; "
                "excludes both adjacent monitor callbacks"
            ),
            "step_warmup_excluded": 1,
            "step_indices": list(range(1, target_steps + 1)),
            "step_physical_times": step_clock["times"],
            "step_times_max_rank_s": step_wall.tolist(),
            "step_times_by_rank_s": rank_step_times,
            "monitor_times_by_rank_s": rank_monitor_times,
            "first_step_wall_s": float(step_wall[0]),
            "steady_step_count": int(steady.size),
            "steady_step_mean_s": float(np.mean(steady)),
            "steady_step_median_s": float(np.median(steady)),
            "steady_step_early_mean_s": float(np.mean(early)),
            "steady_step_late_mean_s": float(np.mean(late)),
            "steady_step_window_size": half,
            "max_rank_hwm_kib": int(max_rss),
            "initial_sha256": hashlib.sha256(
                initial.tobytes()
            ).hexdigest(),
            "final_sha256": hashlib.sha256(
                final.tobytes()
            ).hexdigest(),
            "state_file": state_file.name,
            "source_path": str(Path(S.__file__).resolve()),
            "notes": (
                "End-state captured from original PETSc TS.monitor at "
                "t_f; same single snapshot copy on all variants. "
                "v1.2.0 does not expose a post-TS convergence-reason "
                "attribute, so its native reason is recorded as null."
            ),
        }
        scratch = args.out.with_suffix(".json.partial")
        scratch.write_text(
            json.dumps(report, indent=2) + "\n", encoding="utf-8"
        )
        os.replace(scratch, args.out)
        print(
            f"PASS {args.variant}/{args.policy}/{ranks}r "
            f"DOFs={n_global}, steps={target_steps}, "
            f"solve={solve_wall:.3f}s, "
            f"first_step={step_wall[0]:.4f}s, "
            f"steady_step={np.mean(steady):.4f}s, "
            f"RSS={max_rss/1024:.1f}MiB",
            flush=True,
        )

    if args.petsc_log:
        PETSc.Log.view()


if __name__ == "__main__":
    main()
