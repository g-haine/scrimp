# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""SHORT physical 2D acoustic parity gate for streamed trajectories.

Run both 'memory' and 'async_binary' on the same h, CN time grid,
MUMPS policy, PETSc ranks and source code. Compare complete
time-indexed states and Hamiltonian/power postprocessing later.
This smoke intentionally reconstructs the entire tiny trajectory
on rank 0 to verify it; production asynchronous storage does not.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path

import numpy as np
from mpi4py import MPI

from bench_acoustic_2d import _build_wave
from scrimp.utils.async_binary_trajectory import LazyTrajectory


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=("memory", "async_binary"), required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--mesh", type=float, default=0.08)
    parser.add_argument("--dt", type=float, default=0.001)
    parser.add_argument("--final-time", type=float, default=0.006)
    args = parser.parse_args()

    mpi = MPI.COMM_WORLD
    rank = mpi.Get_rank()
    size = mpi.Get_size()
    os.environ["SCRIMP_TRAJECTORY_BACKEND"] = args.mode
    os.environ["SCRIMP_GETFEM_MATRIX_MODE"] = "distributed"
    os.environ["SCRIMP_MUMPS_RESET_POLICY"] = "never"
    os.environ.setdefault("SCRIMP_INIT_FACTOR_SOLVER_TYPE", "mumps")

    system = _build_wave(args.mesh)
    system.set_time_scheme(
        ts_type="cn", ts_adapt_type="none",
        ts_exact_final_time="matchstep",
        dt=args.dt, dt_save=args.dt,
        ts_adapt_dt_max=args.dt,
        t_f=args.final_time,
        ksp_type="preonly", pc_type="lu",
        pc_factor_mat_solver_type="mumps",
    )
    system.solve()
    if bool(isinstance(system.solution["z"], LazyTrajectory)) != (
        args.mode == "async_binary"
    ):
        raise RuntimeError("Unexpected in-memory/disk trajectory backend")
    if (system.last_ts_converged_reason <= 0
            or system.last_ts_rejections != 0
            or system.last_ts_snes_failures != 0
            or abs(system.last_ts_time - args.final_time) > 1e-12):
        raise RuntimeError("PETSc CN integration invalid")

    times = np.asarray(system.solution["t"])
    if not len(times) or abs(times[-1] - args.final_time) > 1e-12:
        raise RuntimeError("Missing precise final saved time")
    if len(times) != len(system.solution["z"]):
        raise RuntimeError("Incorrect lazy trajectory length")

    # Test only: materialize every saved state for direct coefficient parity.
    states = np.stack(
        [system.solution["z"][k] for k in range(len(times))], axis=0
    )
    if states.shape != (len(times), system.gf_model.nbdof()):
        raise RuntimeError("Full archived state was not reconstructed")
    fingerprint = hashlib.sha256(states.tobytes()).hexdigest()
    if len(set(mpi.allgather(fingerprint))) != 1:
        raise RuntimeError("Reconstructed states differ between MPI ranks")

    system.compute_Hamiltonian()
    hamiltonian = np.asarray(system.get_Hamiltonian())
    system.compute_powers()
    power_names = []
    power_arrays = []
    for name, port in system.ports.items():
        if port.get_is_computed():
            power_names.append(name)
            power_arrays.append(np.asarray(port.get_power()))
    if not np.all(np.isfinite(hamiltonian)):
        raise RuntimeError("Nonfinite Hamiltonian")
    if any(len(power) != len(times) for power in power_arrays):
        raise RuntimeError("Power postprocessing lost saved times")

    if rank == 0:
        args.output_dir.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(
            args.output_dir / f"{args.mode}_physical.npz",
            times=times, states=states, hamiltonian=hamiltonian,
            **{f"power_{i}": item for i, item in enumerate(power_arrays)},
        )
        report = {
            "mode": args.mode,
            "source_commit": os.environ.get("SCRIMP_BENCH_GIT_SHA"),
            "mpi_size": size,
            "dofs": system.gf_model.nbdof(),
            "n_snapshots": len(times),
            "steps": int(system.last_ts_steps),
            "rejected": int(system.last_ts_rejections),
            "snes_failures": int(system.last_ts_snes_failures),
            "final_time": float(system.last_ts_time),
            "state_sha256": fingerprint,
            "power_names": power_names,
            "hamiltonian_shape": list(hamiltonian.shape),
            "snapshot_dtype": str(states.dtype),
        }
        (args.output_dir / f"{args.mode}_physical.json").write_text(
            json.dumps(report, indent=2) + "\n", encoding="utf-8"
        )
        print(json.dumps(report, indent=2), flush=True)


if __name__ == "__main__":
    main()
