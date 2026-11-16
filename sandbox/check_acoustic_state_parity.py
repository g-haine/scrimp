# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Compare complete acoustic wave trajectories across MPI/reset policies.

Uses exactly the published examples/wave.py variational formulation via
bench_acoustic_2d._build_wave. Files contain every saved state coefficient
(including complex values) and every saved timestamp, not just norms.

Use --output for the reference run and --reference to compare against it.
Rank-local trajectories are checked against rank 0 after TS has finished;
no additional collectives are inserted in PETSc callbacks.
"""

import argparse
import json
import os
from pathlib import Path

from mpi4py import MPI
import numpy as np

from bench_acoustic_2d import _build_wave


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mesh", type=float, default=0.08)
    parser.add_argument("--dt", type=float, default=0.002)
    parser.add_argument("--final-time", type=float, default=0.04)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--reference", type=Path)
    parser.add_argument("--rtol", type=float, default=1e-10)
    parser.add_argument("--atol", type=float, default=1e-11)
    args = parser.parse_args()
    if min(args.mesh, args.dt, args.final_time) <= 0:
        parser.error("mesh, dt and final-time must be positive")
    if min(args.atol, args.rtol) < 0:
        parser.error("atol/rtol must be non-negative")

    mpi = MPI.COMM_WORLD
    rank = mpi.Get_rank()
    os.environ.setdefault("SCRIMP_INIT_FACTOR_SOLVER_TYPE", "mumps")
    system = _build_wave(args.mesh)
    system.set_time_scheme(
        ts_type="cn", ts_adapt_type="none",
        ts_exact_final_time="matchstep",
        dt=args.dt, dt_save=args.dt, t_f=args.final_time,
        ksp_type="preonly", pc_type="lu",
        pc_factor_mat_solver_type="mumps",
    )
    system.solve()
    times = np.asarray(system.solution["t"], dtype=np.float64)
    states = np.asarray(system.solution["z"])
    dofs = int(system.gf_model.nbdof())

    if (times.ndim != 1 or states.ndim != 2
            or len(times) != len(states) or states.shape[1] != dofs):
        raise RuntimeError(
            "Expected complete saved (time, global state) arrays, got "
            f"{times.shape}, {states.shape}, dofs={dofs}"
        )
    if (int(system.last_ts_converged_reason) <= 0
            or int(system.last_ts_rejections) != 0
            or int(system.last_ts_snes_failures) != 0
            or abs(times[-1] - args.final_time)
            > 32 * abs(np.spacing(args.final_time))):
        raise RuntimeError("Acoustic solve failed convergence/endpoint gate")

    # The complete, ordered global solution must agree across ranks.
    root_times = mpi.bcast(times if rank == 0 else None, root=0)
    root_states = mpi.bcast(states if rank == 0 else None, root=0)
    if not np.array_equal(times, root_times):
        raise RuntimeError("Saved timesteps differ between MPI ranks")
    rank_diff = float(np.max(np.abs(states - root_states)))
    rank_diff_max = mpi.allreduce(rank_diff, op=MPI.MAX)
    if not np.allclose(states, root_states, rtol=args.rtol, atol=args.atol):
        raise RuntimeError(
            f"Complete state differs between MPI ranks: {rank_diff_max:g}"
        )

    comparison = None
    if args.reference is not None:
        with np.load(args.reference, allow_pickle=False) as saved:
            expected_times = np.array(saved["times"], copy=True)
            expected_states = np.array(saved["states"], copy=True)
        if (not np.array_equal(expected_times, times)
                or expected_states.shape != states.shape):
            raise RuntimeError(
                "Reference trajectory differs in time grid or dimensions"
            )
        difference = states - expected_states
        comparison = {
            "reference": str(args.reference),
            "max_absolute_coefficient_difference":
                float(np.max(np.abs(difference))),
            "relative_trajectory_l2_difference": float(
                np.linalg.norm(difference.ravel())
                / max(1.0, np.linalg.norm(expected_states.ravel()))
            ),
        }
        if not np.allclose(states, expected_states, rtol=args.rtol,
                           atol=args.atol):
            raise RuntimeError(
                "Complete acoustic trajectory differs from reference: "
                + json.dumps(comparison)
            )

    if rank == 0:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(
            args.output, times=times, states=states,
            mesh=np.array(args.mesh), dt=np.array(args.dt),
            mpi_size=np.array(mpi.Get_size()),
        )
        print(json.dumps({
            "mpi_size": mpi.Get_size(),
            "dofs": dofs,
            "snapshots": len(times),
            "final_time": float(times[-1]),
            "max_rank_to_rank_state_difference": rank_diff_max,
            "reset_pc_on_jacobian":
                os.environ.get("SCRIMP_MUMPS_RESET_PC_ON_JACOBIAN", "0"),
            "output": str(args.output),
            "comparison": comparison,
        }, indent=2), flush=True)


if __name__ == "__main__":
    main()
