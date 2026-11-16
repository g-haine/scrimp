# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Reproducible 2D acoustic-wave strong-scaling benchmark.

Builds the *actual* examples/wave.py model, intercepting its solve()
immediately before the long published run. Overrides only Rectangle
mesh spacing, PETSc time/sparse solver options and output policy.
No plots, no edits to the existing physical formulation.

Run the SAME --mesh and --final-time at ranks 1, 2, 4, 8 to obtain
strong-scaling results. A mesh refined with P is a DIFFERENT experiment.
"""
import argparse
import importlib.util
import json
import os
from pathlib import Path
import resource
import time
from unittest.mock import patch

import numpy as np
from mpi4py import MPI
from petsc4py import PETSc
import scrimp as S


class _StopBeforeOriginalSolve(Exception):
    def __init__(self, system):
        self.system = system


def _build_wave(mesh):
    source = Path(__file__).resolve().parents[1] / "examples" / "wave.py"
    spec = importlib.util.spec_from_file_location("_scrimp_bench_acoustic_2d", source)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    original_domain = S.Domain

    def domain_with_mesh(kind, parameters, *args, **kwargs):
        params = dict(parameters)
        if kind == "Rectangle":
            params["h"] = mesh
        return original_domain(kind, params, *args, **kwargs)

    def stop_before_solve(system):
        raise _StopBeforeOriginalSolve(system)

    with patch.object(S, "Domain", domain_with_mesh), patch.object(
        S.DPHS, "solve", stop_before_solve
    ):
        try:
            module.wave_eq()
        except _StopBeforeOriginalSolve as captured:
            return captured.system
    raise RuntimeError("The existing acoustic wave example did not call DPHS.solve")


def _norms(system):
    if not system.solution["z"]:
        return {}
    result = {}
    last = np.asarray(system.solution["z"][-1], dtype=float)
    for name in ("q", "p", "e_q", "e_p", "U_B", "U_R", "U_T", "U_L"):
        try:
            start, n = system.gf_model.interval_of_variable(name)
        except (RuntimeError, ValueError, KeyError):
            continue
        result[name] = float(np.linalg.norm(last[int(start):int(start + n)]))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mesh", type=float, default=0.08,
                        help="Fixed 2D Rectangle mesh spacing for strong scaling")
    parser.add_argument("--final-time", type=float, default=0.02)
    parser.add_argument("--dt", type=float, default=0.001)
    parser.add_argument("--save-dt", type=float, default=0.01)
    parser.add_argument("--output", type=Path, default=None,
                        help="Write JSON report from MPI rank 0")
    # PETSc options, including -log_view, have already been parsed by
    # petsc4py before this parser. Leave their argv entries intact.
    args, _petsc_options = parser.parse_known_args()
    if min(args.mesh, args.final_time, args.dt, args.save_dt) <= 0:
        parser.error("mesh/final-time/dt/save-dt must all be positive")

    mpi = MPI.COMM_WORLD
    p = mpi.Get_size()
    rank = mpi.Get_rank()
    os.environ.setdefault("SCRIMP_INIT_FACTOR_SOLVER_TYPE", "mumps")

    mpi.Barrier()
    start = time.perf_counter()
    system = _build_wave(args.mesh)
    mpi.Barrier()
    built = time.perf_counter()
    dofs = system.gf_model.nbdof()

    # Fixed-step Crank-Nicolson. Petsc LU/MUMPS settings are the
    # same for all ranks, so speedups are not time-adapt artifacts.
    system.set_time_scheme(
        ts_type="cn",
        ts_adapt_type="none",
        ts_exact_final_time="matchstep",
        dt=args.dt,
        dt_save=args.save_dt,
        t_f=args.final_time,
        ksp_type="preonly",
        pc_type="lu",
        pc_factor_mat_solver_type="mumps",
    )
    system.solve()
    mpi.Barrier()
    finished = time.perf_counter()

    build_wall = mpi.allreduce(built-start, op=MPI.MAX)
    solve_wall = mpi.allreduce(finished-built, op=MPI.MAX)
    total_wall = mpi.allreduce(finished-start, op=MPI.MAX)
    # Linux ru_maxrss is KiB, *per rank*. Max does not measure total
    # memory across ranks but exposes per-process bottlenecks.
    rss_rank_kib = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    rss_peak_rank_kib = mpi.allreduce(rss_rank_kib, op=MPI.MAX)
    t_petsc = float(system.last_ts_time)
    t_saved = float(system.solution["t"][-1])
    report = {
        "benchmark": "examples/wave.py:wave_eq() (2D acoustic waves)",
        "mpi_size": p,
        "degrees_of_freedom": int(dofs),
        "rectangle_mesh_h": args.mesh,
        "ts_type": "cn",
        "ts_adapt_type": "none",
        "requested_final_time": args.final_time,
        "actual_final_time": t_petsc,
        "last_saved_time": t_saved,
        "ts_converged_reason": int(system.last_ts_converged_reason),
        "ts_steps": int(system.last_ts_steps),
        "ts_rejected_steps": int(system.last_ts_rejections),
        "ts_snes_failures": int(system.last_ts_snes_failures),
        "wall_build_seconds_max": build_wall,
        "wall_solve_seconds_max": solve_wall,
        "wall_total_seconds_max": total_wall,
        "peak_rank_rss_kib": int(rss_peak_rank_kib),
        "final_variable_norms": _norms(system),
        "numpy_version": np.__version__,
    }
    if rank == 0:
        text = json.dumps(report, indent=2)
        print(text, flush=True)
        if args.output is not None:
            args.output.parent.mkdir(parents=True, exist_ok=True)
            args.output.write_text(text + "\n", encoding="utf-8")
    tol = 32.0 * abs(np.spacing(args.final_time))
    if (report["ts_converged_reason"] <= 0
            or abs(t_petsc - args.final_time) > tol
            or abs(t_saved - args.final_time) > tol):
        raise RuntimeError(
            "Acoustic benchmark failed MATCHSTEP endpoint gate: "
            f"PETSc={t_petsc}, saved={t_saved}, requested={args.final_time}"
        )


if __name__ == "__main__":
    main()
