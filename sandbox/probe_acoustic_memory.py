# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Measure actual per-rank process memory during a long acoustic TS run.

The same examples/wave.py physics, fixed-step CN and MUMPS policy are
used for 'dense' and 'sparse' saving. The only intended difference is
dt_save, not the integration step. The DPHS monitor is WRAPPED on each
instance, not replaced, so state storage, GetFEM updates and PETSc
callbacks preserve their original behavior.

No extra MPI collectives are introduced inside PETSc callbacks; each
rank writes its own JSON samples. /proc/self/status measures process
RSS, NOT PETSc-only or MUMPS-only memory. Saved-state payload bytes and
PETSc Mat.getInfo() storage are reported independently and are NOT
assumed to partition all resident memory.

Run only inside scrimp-mpi:local with GetFEM 5.5.
"""
import argparse
import json
import os
from pathlib import Path
import resource
import time

import numpy as np
from mpi4py import MPI
from petsc4py import PETSc

from bench_acoustic_2d import _build_wave


_PROC_KEYS = (
    "VmRSS", "VmHWM", "VmSize", "VmData", "VmSwap",
    "RssAnon", "RssFile", "RssShmem",
)


def _proc_memory_kib():
    """Read local Linux /proc memory without psutil or MPI collectives."""
    result = {}
    for line in Path("/proc/self/status").read_text(
        encoding="ascii"
    ).splitlines():
        name, _, value = line.partition(":")
        if name in _PROC_KEYS:
            words = value.strip().split()
            if len(words) < 2 or words[1] != "kB":
                raise RuntimeError(f"Unexpected /proc memory line: {line}")
            result[name + "_kib"] = int(words[0])
    if "VmRSS_kib" not in result or "VmHWM_kib" not in result:
        raise RuntimeError("Linux /proc memory counters unavailable")
    result["ru_maxrss_kib"] = resource.getrusage(
        resource.RUSAGE_SELF
    ).ru_maxrss
    return result


def _stored_state_info(system):
    """Count NumPy array payload, not Python/PETSc/MUMPS allocations."""
    saved = system.solution["z"]
    if isinstance(saved, np.ndarray):
        if not saved.size:
            return 0, 0
        return int(saved.shape[0]), int(saved.nbytes)
    return len(saved), sum(int(np.asarray(x).nbytes) for x in saved)


def _mat_info(system):
    """Diagnostic of persistent PETSc matrix allocation, excluding LU."""
    result = {}
    for name in (
        "mass", "stiffness", "nl_mass", "nl_stiffness",
        "tangent_mass", "tangent_stiffness", "J",
    ):
        mat = getattr(system, name, None)
        if mat is None:
            continue
        try:
            data = mat.getInfo()  # Local PETSc MatInfo, no MPI reduction.
            result[name] = {
                key: float(data[key])
                for key in ("memory", "nz_used", "nz_allocated", "mallocs")
                if key in data
            }
        except (PETSc.Error, ValueError, TypeError, AttributeError) as exc:
            result[name] = {"unavailable": type(exc).__name__}
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mesh", type=float, default=0.04)
    parser.add_argument("--dt", type=float, default=0.002)
    parser.add_argument("--final-time", type=float, default=0.5)
    parser.add_argument("--mode", choices=("dense", "sparse"), required=True)
    parser.add_argument("--sample-every", type=int, default=10)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    if min(args.mesh, args.dt, args.final_time) <= 0:
        parser.error("mesh, dt, final-time must be positive")
    if args.sample_every < 1:
        parser.error("--sample-every must be >= 1")
    expected_steps = round(args.final_time / args.dt)
    if not np.isclose(
        expected_steps * args.dt, args.final_time,
        rtol=0, atol=1e-12
    ):
        parser.error("Final time must be an integer multiple of dt")

    method = os.environ.get(
        "SCRIMP_GETFEM_BRIDGE_DISTRIBUTED", "coo"
    ).strip().lower()
    matrix_mode = os.environ.get(
        "SCRIMP_GETFEM_MATRIX_MODE", "replicated"
    ).strip().lower()
    if matrix_mode != "distributed" or method != "coo":
        raise RuntimeError(
            "Memory audit requires distributed GetFEM and default/explicit COO"
        )
    policy = os.environ.get("SCRIMP_MUMPS_RESET_POLICY", "never")
    if policy != "never":
        raise RuntimeError(
            "Linear acoustic memory audit requires MUMPS never"
        )

    os.environ["SCRIMP_PROFILE"] = "1"
    os.environ.setdefault("SCRIMP_INIT_FACTOR_SOLVER_TYPE", "mumps")
    mpi = MPI.COMM_WORLD
    rank = mpi.Get_rank()
    size = mpi.Get_size()
    begun = time.perf_counter()
    samples = []
    system = None

    def snap(phase, step, t):
        nstates, state_bytes = (
            _stored_state_info(system)
            if system is not None else (0, 0)
        )
        samples.append({
            "phase": phase,
            "step": int(step),
            "simulation_time": float(t),
            "wall_elapsed_seconds": time.perf_counter() - begun,
            "stored_snapshots": nstates,
            "stored_state_bytes": state_bytes,
            **_proc_memory_kib(),
        })

    snap("before_build", -1, 0.)
    mpi.Barrier()
    system = _build_wave(args.mesh)
    snap("after_build", -1, 0.)

    # 'sparse' saves only initial/final state, because its dt_save exceeds
    # the full integration interval. It still executes every TS step.
    save_dt = (
        args.dt if args.mode == "dense"
        else 2.0 * args.final_time
    )
    system.set_time_scheme(
        ts_type="cn",
        ts_adapt_type="none",
        ts_exact_final_time="matchstep",
        dt=args.dt,
        dt_save=save_dt,
        t_f=args.final_time,
        ksp_type="preonly",
        pc_type="lu",
        pc_factor_mat_solver_type="mumps",
    )

    original_monitor = system.monitor

    def instrumented_monitor(TS, i, t, z, **kwargs):
        # Keep SCRIMP's exact monitor behavior. Never inject MPI
        # collectives or modify saved arrays, GetFEM state or PETSc Vec.
        value = original_monitor(TS, i, t, z, **kwargs)
        if (not kwargs.get("initial_step", False)
                and (i <= 0 or i % args.sample_every == 0)):
            snap("main_monitor", i, t)
        return value

    system.monitor = instrumented_monitor
    snap("before_solve", -1, 0.)
    mpi.Barrier()
    t0 = time.perf_counter()
    system.solve()
    mpi.Barrier()
    duration = time.perf_counter() - t0
    snap("after_solve", -1, float(system.last_ts_time))

    times = np.asarray(system.solution["t"], dtype=np.float64)
    states = np.asarray(system.solution["z"])
    dofs = int(system.gf_model.nbdof())
    final_state = np.array(states[-1], copy=True)
    convergence = {
        "reason": int(system.last_ts_converged_reason),
        "steps": int(system.last_ts_steps),
        "rejections": int(system.last_ts_rejections),
        "snes_failures": int(system.last_ts_snes_failures),
        "last_ts_time": float(system.last_ts_time),
    }
    if (convergence["reason"] <= 0
            or convergence["steps"] != expected_steps
            or convergence["rejections"] != 0
            or convergence["snes_failures"] != 0
            or not np.isclose(
                convergence["last_ts_time"], args.final_time,
                atol=1e-12, rtol=0,
            )
            or not np.isclose(
                times[-1], args.final_time, atol=1e-12, rtol=0
            )
            or states.ndim != 2 or states.shape[1] != dofs
            or not np.all(np.isfinite(final_state))):
        raise RuntimeError(
            "Invalid acoustic solve/endpoint/trajectory: "
            + json.dumps(convergence)
        )

    count, payload = _stored_state_info(system)
    sample_stats = _mat_info(system)
    profile = system.get_profile()  # Collective, only AFTER solve.
    report = {
        "benchmark": "memory audit of examples/wave.py",
        "source_commit": os.environ.get("SCRIMP_BENCH_GIT_SHA", ""),
        "bridge_method": method,
        "getfem_matrix_mode": matrix_mode,
        "mumps_reset_policy": policy,
        "mpi_size": size,
        "mpi_rank": rank,
        "mesh": args.mesh,
        "dofs": dofs,
        "dt": args.dt,
        "dt_save": save_dt,
        "final_time": args.final_time,
        "save_mode": args.mode,
        "expected_steps": expected_steps,
        "sample_every": args.sample_every,
        "ts_convergence": convergence,
        "solve_wall_seconds_rank": duration,
        "snapshots": count,
        "stored_state_bytes": payload,
        "last_saved_time": float(times[-1]),
        "final_state_l2": float(np.linalg.norm(final_state)),
        "final_state_dtype": str(final_state.dtype),
        "mat_info_after_solve_local": sample_stats,
        "profile": profile,
        "samples": samples,
        "notes": [
            "VmRSS, RssAnon, VmHWM are Linux per-process counters.",
            "VmHWM and ru_maxrss include native PETSc/MUMPS and Python.",
            "Mat.getInfo().memory is matrix-local, not factorization RAM.",
            "Stored states are NumPy payload only; no double counting.",
            "The 'sparse' run changes only dt_save, not the PDE or dt.",
            "Samples after monitor reflect completed state saves.",
        ],
    }
    args.output_dir.mkdir(parents=True, exist_ok=True)
    stem = f"memory_{args.mode}_rank{rank}"
    (args.output_dir / (stem + ".json")).write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    np.savez_compressed(
        args.output_dir / (stem + "_state.npz"),
        initial=np.asarray(states[0]),
        final=final_state,
        first_saved_time=float(times[0]),
        final_time=float(times[-1]),
    )
    mpi.Barrier()
    if rank == 0:
        print(json.dumps({
            "mode": args.mode,
            "mpi_size": size,
            "dofs": dofs,
            "steps": convergence["steps"],
            "snapshots_rank0": count,
            "payload_mib_rank0": payload / (1024.**2),
            "wall_seconds_rank0": duration,
            "report_pattern": f"memory_{args.mode}_rank*.json",
        }, indent=2), flush=True)


if __name__ == "__main__":
    main()
