# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Very short nonlinear dissipative Saint-Venant 2D output-backend parity.

Only a handful of fixed Crank-Nicolson steps. Reuse the existing
physical SWE model and power/volume diagnostics unchanged. A gentle
start-from-rest spin-up avoids the known inconsistent rotating initial
wall impulse. This test is about IO parity, not about independently
establishing the physical accuracy of the SWE discretization.

Both backends execute identical PETSc/Newton/GetFEM dynamics. Complete
sampled states are materialized only in the TEST to compare them.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path

import numpy as np
from mpi4py import MPI

# Invoked as `python sandbox/verify_async_swe.py`: Python puts
# /workspace/scrimp/sandbox (not its parent) first on sys.path.
# Import the adjacent module, as the acoustic sandbox driver does.
from viscous_swe_2d_mpi import build_system, diagnostics
from scrimp.utils.async_binary_trajectory import LazyTrajectory
from matrix_memory_metrics import snapshot as matrix_memory_snapshot


_DIAGNOSTICS = (
    "t", "kinetic", "potential", "H", "d_grad", "d_div",
    "p_normal", "p_viscous", "supplied", "dissipated",
    "balance", "relative_balance", "relative_balance_exchange",
    "volume", "normal_volume_outflow", "cumulative_outflow",
    "volume_balance", "relative_volume_balance", "lyapunov",
)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--mode", choices=("memory", "async_binary"), required=True
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--mesh", type=float, default=0.35)
    parser.add_argument("--dt", type=float, default=0.001)
    parser.add_argument("--final-time", type=float, default=0.003)
    parser.add_argument("--spinup-duration", type=float, default=0.03)
    args = parser.parse_args()
    if not (args.mesh > 0 and args.dt > 0 and
            args.final_time > 0 and args.spinup_duration > 0):
        parser.error("mesh, dt, T, spin-up duration must be positive")
    expected_steps = round(args.final_time / args.dt)
    if (expected_steps < 2 or
            abs(expected_steps * args.dt - args.final_time) > 1e-12):
        parser.error("T must comprise at least two exact dt increments")

    mpi = MPI.COMM_WORLD
    rank = mpi.Get_rank()
    size = mpi.Get_size()
    os.environ["SCRIMP_TRAJECTORY_BACKEND"] = args.mode
    os.environ["SCRIMP_GETFEM_MATRIX_MODE"] = "distributed"
    os.environ["SCRIMP_MUMPS_RESET_POLICY"] = "on_change"
    os.environ.setdefault("SCRIMP_INIT_FACTOR_SOLVER_TYPE", "mumps")

    swe = build_system(
        gamma=0.0, mesh=args.mesh, final_time=args.final_time,
        save_dt=args.dt, max_step=args.dt, init_step=True,
        ts_type="cn", adapt_type="none", init_dt=1e-9,
        rotation=0.2, initial_profile="rigid",
        wall_velocity_profile="facet_tangent",
        algebraic_initial_guess="zero",
        start_from_rest_duration=args.spinup_duration,
        init_method="bdf", height_order=2,
    )
    swe.solve()
    memory_by_rank = mpi.allgather(matrix_memory_snapshot(swe, rank))
    is_disk = isinstance(swe.solution["z"], LazyTrajectory)
    if is_disk != (args.mode == "async_binary"):
        raise RuntimeError("Incorrect trajectory backend after PETSc solve")
    if swe.linear_mass and swe.linear_stiffness:
        raise RuntimeError("This test must exercise NONLINEAR GetFEM")
    convergence = {
        "reason": int(swe.last_ts_converged_reason),
        "steps": int(swe.last_ts_steps),
        "rejections": int(swe.last_ts_rejections),
        "snes_failures": int(swe.last_ts_snes_failures),
        "t_final": float(swe.last_ts_time),
    }
    if (convergence["reason"] <= 0
            or convergence["steps"] != expected_steps
            or convergence["rejections"] != 0
            or convergence["snes_failures"] != 0
            or abs(convergence["t_final"]-args.final_time) > 1e-12):
        raise RuntimeError(
            "Unsuccessful nonlinear time integration: "
            + json.dumps(convergence)
        )

    times = np.asarray(swe.solution["t"], dtype=np.float64)
    z = swe.solution["z"]
    if (len(times) != len(z) or
            len(times) != expected_steps + 1 or
            abs(times[-1] - args.final_time) > 1e-12):
        raise RuntimeError("Unexpected saved-state time grid")
    # ONLY this tiny comparison driver materializes ALL snapshots.
    states = np.stack([z[k] for k in range(len(z))], axis=0)
    if (states.shape != (len(times), swe.gf_model.nbdof())
            or not np.all(np.isfinite(states))):
        raise RuntimeError("Invalid archived nonlinear state values")
    h_lo, h_size = swe.gf_model.interval_of_variable("h")
    height = states[:, int(h_lo):int(h_lo + h_size)]
    if np.min(height) <= 0.:
        raise RuntimeError("Nonphysical nonpositive water height")
    if not np.any(np.abs(states[-1] - states[0]) > 1e-13):
        raise RuntimeError("Spin-up has no evolving state")

    # Reuse original SWE Hamiltonian, viscous dissipation, boundary
    # power, volume conservation and midpoint CN integration.
    result = diagnostics(swe, gamma=0.0)
    arrays = {key: np.asarray(result[key]) for key in _DIAGNOSTICS}
    for name, values in arrays.items():
        if values.shape != times.shape or not np.all(np.isfinite(values)):
            raise RuntimeError(f"Bad {name} physical diagnostics")

    signature = hashlib.sha256(states.tobytes()).hexdigest()
    signatures = mpi.allgather(signature)
    if len(set(signatures)) != 1:
        raise RuntimeError(
            "Different archived nonlinear states across MPI ranks: "
            + str(signatures)
        )
    if rank == 0:
        args.output_dir.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(
            args.output_dir / f"{args.mode}_swe.npz",
            states=states, **arrays,
        )
        metadata = {
            "mode": args.mode,
            "source_commit": os.environ.get("SCRIMP_BENCH_GIT_SHA", ""),
            "ranks": size,
            "dofs": int(swe.gf_model.nbdof()),
            "snapshots": len(times),
            "convergence": convergence,
            "mesh": args.mesh,
            "dt": args.dt,
            "final_time": args.final_time,
            "spinup_duration": args.spinup_duration,
            "gamma": 0.0,
            "ts_type": "cn",
            "ts_adapt_type": "none",
            "mumps_reset_policy": "on_change",
            "state_sha256": signature,
            "matrix_memory_mode": os.environ.get(
                "SCRIMP_MATRIX_MEMORY_MODE", "legacy"
            ).strip().lower(),
            "matrix_memory_per_rank": memory_by_rank,
            "minimum_water_height": float(np.min(height)),
            "energy_balance_defect_final":
                float(arrays["balance"][-1]),
            "volume_balance_defect_final":
                float(arrays["volume_balance"][-1]),
            "notes": (
                "This tests IO parity with short nonlinear spin-up; "
                "nonzero physical balance defects are NOT automatically "
                "accepted as an accurate published SWE solution."
            ),
        }
        (args.output_dir / f"{args.mode}_swe.json").write_text(
            json.dumps(metadata, indent=2) + "\n",
            encoding="utf-8",
        )
        print(json.dumps(metadata, indent=2), flush=True)


if __name__ == "__main__":
    main()
