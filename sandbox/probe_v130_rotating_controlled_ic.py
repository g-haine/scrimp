# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Native pilot: generate a compatible rotating, dissipative SWE initial state.

A prescribed already-rotating h,p state was proved algebraically
incompatible for the original curved-P2 discretization. Here the real
physical wall is smoothly started from rest (t=0..0.02), with the
original decaying normal-feedback law kept independent. Thus gamma=.25
is allowed without injecting normal feedback into the stationary t0.

The output is an ACTUAL ACCEPTED DAE state checkpoint, never a claim
that the rotating interior has reached an equilibrium. Resume can
continue in bounded short segments, on a new MPI process count, under
the identical absolute physical boundary law. No changes to SCRIMP core.
"""

import argparse
import json
import math
from pathlib import Path
from unittest.mock import patch

import numpy as np
from mpi4py import MPI
from petsc4py import PETSc
import scrimp as S

from viscous_swe_2d_mpi import (
    build_system, diagnostics, ROTATION, H_INITIAL, H_TARGET,
    spinup_control_envelopes,
)
from probe_v130_issue28_one_cn_spinup import (
    native_residual, publish_sandbox_rest,
)
from probe_v130_issue28_structure import indices, layout
from probe_v130_issue28_weak_wall import ensure_real
from scrimp.utils.dae_checkpoint import model_provenance, sha256

RAMP_DURATION = 0.02
RAMP_DECAY_START = 75.0
RAMP_DECAY_END = 200.0
DEFAULT_DT = 0.002
MAX_SEGMENT_STEPS = 100


def checked_grid(t0, tf, dt):
    """Bound memory and refuse an accidental PETSc partial final step."""
    for name, value in (("t0", t0), ("tf", tf), ("dt", dt)):
        if not math.isfinite(value):
            raise ValueError(name + " must be finite")
    if dt <= 0 or t0 < 0 or tf <= t0:
        raise ValueError("Require 0 <= t0 < tf and dt > 0")
    n = int(round((tf - t0) / dt))
    if n < 1 or n > MAX_SEGMENT_STEPS:
        raise ValueError("Limit each segment to 1..100 accepted steps")
    if abs(t0 + n * dt - tf) > 1e-11:
        raise ValueError("Final physical time must lie on the fixed CN grid")
    if tf > RAMP_DECAY_START + 1e-12:
        raise ValueError(
            "Initial-state preparation pilot ends before normal-feedback "
            "activation at t=75. Use a separate controlled evolution "
            "qualification for t>75."
        )
    return n


def model(mesh, gamma, tf, dt):
    return build_system(
        gamma=gamma, mesh=mesh, geometry_order=2, height_order=2,
        final_time=tf, save_dt=dt * (1. - 1e-9), max_step=dt,
        ramp_start=RAMP_DECAY_START, ramp_end=RAMP_DECAY_END,
        init_step=True, ts_type="cn", adapt_type="none",
        init_dt=1e-9, rotation=ROTATION, initial_profile="rigid",
        wall_velocity_profile="facet_tangent",
        algebraic_initial_guess="zero",
        start_from_rest_duration=RAMP_DURATION, init_method="bdf",
    )


def evaluated_envelopes(t):
    """Independent check of the exact strings sent to original GetFEM."""
    wall, feedback = spinup_control_envelopes(
        RAMP_DECAY_START, RAMP_DECAY_END, RAMP_DURATION
    )
    env = {"t": float(t), "min": min, "max": max}
    return (
        float(eval(wall, {"__builtins__": {}}, env)),
        float(eval(feedback, {"__builtins__": {}}, env)),
    )


def run():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--phase", choices=("seed", "resume"), required=True)
    cli.add_argument("--mesh", choices=(0.2, 0.1), type=float, required=True)
    cli.add_argument("--gamma", choices=(0.0, 0.25), type=float, default=0.25)
    cli.add_argument("--final-time", type=float, required=True)
    cli.add_argument("--dt", type=float, default=DEFAULT_DT)
    cli.add_argument(
        "--qualification-mode", choices=("physical", "restart_parity"),
        default="physical",
        help="Physical gate checks model energy defect; restart_parity "
             "records that defect but qualifies only fidelity to uninterrupted TS.",
    )
    cli.add_argument("--checkpoint-in", type=Path)
    cli.add_argument("--checkpoint-out", type=Path, required=True)
    cli.add_argument("--out", type=Path, required=True)
    args = cli.parse_args()
    mpi = PETSc.COMM_WORLD.tompi4py()
    if (args.checkpoint_out.exists()
            or (args.out / "initial_state_report.json").exists()):
        raise FileExistsError("Refusing to overwrite native state evidence")
    if args.phase == "resume" and (
        args.checkpoint_in is None or not args.checkpoint_in.is_file()
    ):
        raise ValueError("Resume requires an existing physical checkpoint")
    if args.phase == "seed" and args.checkpoint_in is not None:
        raise ValueError("Seed must use original analytic compatible rest state")
    checked_grid(0., args.final_time, args.dt) if args.phase == "seed" else None

    system = model(args.mesh, args.gamma, args.final_time, args.dt)
    saved = {}
    if args.phase == "resume":
        header = system.load_checkpoint(args.checkpoint_in)
        t0 = float(system.time_scheme["t_0"])
        if t0 < RAMP_DURATION - 1e-12:
            raise RuntimeError("Cannot resume a pre-ramp preparation")
        nsteps = checked_grid(t0, args.final_time, args.dt)
        extra = header.get("extra", {})
        if (extra.get("experiment") != "controlled_rotating_SWE_preparation"
                or extra.get("gamma") != args.gamma
                or extra.get("nominal_mesh") != args.mesh
                or extra.get("ramp_duration") != RAMP_DURATION
                or extra.get("ramp_decay_start") != RAMP_DECAY_START
                or extra.get("ramp_decay_end") != RAMP_DECAY_END):
            raise RuntimeError("Source checkpoint has incompatible control protocol")
        system.solve()
        start_report = system.last_checkpoint_report
        if (not isinstance(start_report, dict)
                or start_report.get("no_pseudo_BDF_or_projection") is not True
                or start_report.get(
                    "all_GetFEM_and_GlobalVec_coefficients_restored_bitwise"
                ) is not True):
            raise RuntimeError("Production checkpoint cold restoration failed")
    else:
        t0 = 0.
        nsteps = checked_grid(t0, args.final_time, args.dt)
        if args.final_time < RAMP_DURATION - 1e-12:
            raise ValueError("First checkpoint requires completed C1 spinup ramp")

        def exact_rest_initializer(self):
            publish_sandbox_rest(self, mpi, saved)

        # This analytic initial state is ALREADY validated for the open
        # case. The new closed-loop envelope must pass the exact same
        # original native t0 IFunction check (no silent pseudo-BDF).
        with patch.object(S.DPHS, "init_step", exact_rest_initializer):
            system.solve()
        start_report = system.initialization_report
        if (not isinstance(start_report, dict)
                or start_report.get("strategy") !=
                "sandbox_exact_stationary_spinup_only"):
            raise RuntimeError("Exact original t0 startup was not accepted")

    actual = (
        float(system.last_ts_time), int(system.last_ts_steps),
        int(system.last_ts_converged_reason),
        int(system.last_ts_rejections), int(system.last_ts_snes_failures),
    )
    if (abs(actual[0] - args.final_time) > 1e-12
            or actual[1] != nsteps or actual[2] <= 0
            or actual[3] or actual[4]):
        raise RuntimeError("Original PETSc TSCN did not accept full segment: "
                           + repr(actual))
    times = np.asarray(system.solution["t"], dtype=float)
    states = np.asarray([
        ensure_real(z, "original accepted GetFEM FE vector").copy()
        for z in system.solution["z"]
    ])
    dynamic, algebraic, iv, ndof = layout(system)
    if dynamic != ("h", "p") or len(algebraic) == 0:
        raise RuntimeError("Unexpected original SWE DAE variable partition")
    expected_times = np.linspace(t0, args.final_time, nsteps + 1)
    if (states.shape != (nsteps + 1, ndof)
            or times.shape != (nsteps + 1,)
            or not np.allclose(times, expected_times, atol=1e-12, rtol=0.)
            or not np.all(np.isfinite(states))):
        raise RuntimeError("Missing original complete accepted CN trajectory")
    if args.phase == "resume":
        with np.load(args.checkpoint_in, allow_pickle=False) as cp:
            old = np.asarray(cp["state"], dtype=float)
        if not np.array_equal(old, states[0]):
            raise RuntimeError("Cold restore changed full checkpoint initial FE vector")

    di = indices(dynamic, iv)
    ai = indices(algebraic, iv)
    hstart, hcount = iv["h"]
    min_h = float(np.min(states[:, hstart:hstart + hcount]))
    if min_h <= 0:
        raise RuntimeError("Non-positive height FE coefficient in physical trajectory")
    cn_residuals = []
    fa_residuals = []
    derivative_zero = np.zeros(ndof)
    # Check original complete nonlinear GetFEM equations on every
    # accepted endpoint; PETSc success alone is NOT a physics gate.
    for k in range(len(times)):
        tk, zk = float(times[k]), states[k]
        _, fk = native_residual(
            system, mpi, zk, derivative_zero, tk, iv,
            (dynamic, algebraic), with_vector=True
        )
        fa = float(np.linalg.norm(fk[ai]))
        if not math.isfinite(fa) or fa >= 1e-8:
            raise RuntimeError("Original DAE algebraic violation at time "
                               f"{tk:.8g}: {fa:.9g}")
        fa_residuals.append(fa)
        if k == 0:
            continue
        tprev = float(times[k - 1])
        zkprev = states[k - 1]
        _, fres = native_residual(
            system, mpi, zk, 2. * (zk-zkprev)/(tk-tprev), tk,
            iv, (dynamic, algebraic), with_vector=True
        )
        _, fprev = native_residual(
            system, mpi, zkprev, derivative_zero, tprev,
            iv, (dynamic, algebraic), with_vector=True
        )
        cn = float(np.linalg.norm(fres + fprev))
        if not math.isfinite(cn) or cn >= 1e-7:
            raise RuntimeError("Original TSCN residual violation at step "
                               f"{k}: {cn:.9g}")
        cn_residuals.append(cn)

    # Necessary, never sufficient, steady-state metrics.
    stationarity = native_residual(
        system, mpi, states[-1], derivative_zero, float(times[-1]),
        iv, (dynamic, algebraic),
    )
    dzdt = float(np.linalg.norm(
        (states[-1, di] - states[-2, di]) / (times[-1] - times[-2])
    ))
    power = diagnostics(system, gamma=args.gamma)
    names = (
        "H", "volume", "volume_balance", "supplied", "dissipated",
        "p_normal", "p_viscous", "d_grad", "d_div", "balance",
    )
    traces = {name: np.asarray(power[name], dtype=float) for name in names}
    if any(v.shape != times.shape or not np.all(np.isfinite(v))
           for v in traces.values()):
        raise RuntimeError("Invalid original physical power/mass traces")
    dH = traces["H"] - traces["H"][0]
    exchanged = max(
        float(np.max(np.abs(dH))),
        float(np.max(np.abs(traces["supplied"]))),
        float(np.max(np.abs(traces["dissipated"]))),
        1e-14,
    )
    power_defect = float(np.max(np.abs(traces["balance"]))) / exchanged
    mass_defect = float(np.max(np.abs(traces["volume_balance"])))
    # Preserve scientific evidence BEFORE a power gate can reject a
    # numerically accepted PETSc trajectory. A short cold segment can
    # exchange little energy, amplifying its per-segment relative error;
    # only a continuous-vs-split native comparison can diagnose that.
    dH_step = np.diff(traces["H"])
    dW_step = np.diff(traces["supplied"])
    dD_step = np.diff(traces["dissipated"])
    balance_step = dH_step + dD_step - dW_step
    step_exchange = np.maximum.reduce([
        np.abs(dH_step), np.abs(dW_step), np.abs(dD_step),
        np.full_like(dH_step, 1e-14),
    ])
    precheck = {
        "schema": "scrimp-v130-rotating-controlled-precheck-v2",
        "status": "NATIVE_PHYSICAL_DIAGNOSTIC_NOT_A_CHECKPOINT",
        "phase": args.phase, "qualification_mode": args.qualification_mode,
        "source_time": t0, "final_time": float(times[-1]),
        "mesh": args.mesh, "gamma": args.gamma, "ranks": mpi.size,
        "source_checkpoint_sha256": (
            sha256(args.checkpoint_in.read_bytes())
            if args.checkpoint_in is not None else None
        ),
        "source_model_fingerprint": model_provenance(system)["fingerprint_sha256"],
        "PETSc_accepted_steps": actual[1],
        "PETSc_converged_reason": actual[2],
        "PETSc_rejections": actual[3],
        "PETSc_SNES_failures": actual[4],
        "maximum_actual_original_CN_residual": max(cn_residuals),
        "maximum_actual_original_Fa_residual": max(fa_residuals),
        "min_height_FE_coefficient": min_h,
        "segment_energy_exchange_scale": exchanged,
        "maximum_absolute_cumulative_energy_defect": float(
            np.max(np.abs(traces["balance"]))
        ),
        "maximum_absolute_per_step_energy_defect": float(
            np.max(np.abs(balance_step))
        ),
        "maximum_relative_per_step_energy_defect": float(
            np.max(np.abs(balance_step) / step_exchange)
        ),
        "relative_cumulative_defect_over_segment_exchange": power_defect,
        "max_absolute_volume_defect": mass_defect,
        "physical_H_change": float(dH[-1]),
        "physical_boundary_work": float(traces["supplied"][-1]),
        "physical_viscous_dissipation": float(traces["dissipated"][-1]),
        "checkpoint_created": False,
    }
    save_error = None
    if mpi.rank == 0:
        try:
            args.out.mkdir(parents=True, exist_ok=True)
            (args.out / "physical_precheck_report.json").write_text(
                json.dumps(precheck, indent=2) + "\n", encoding="utf-8"
            )
            np.savez_compressed(
                args.out / "physical_precheck_trace.npz",
                times=times, states=states,
                full_state_start=states[0],
                full_state_final=states[-1],
                H=traces["H"], volume=traces["volume"],
                boundary_work=traces["supplied"],
                viscous_dissipation=traces["dissipated"],
                p_normal=traces["p_normal"],
                p_viscous=traces["p_viscous"],
                d_grad=traces["d_grad"], d_div=traces["d_div"],
                power_balance=traces["balance"],
                volume_balance=traces["volume_balance"],
                step_delta_H=dH_step, step_delta_work=dW_step,
                step_delta_dissipation=dD_step,
                step_balance_defect=balance_step,
                step_exchange_scale=step_exchange,
                Fa_residuals=np.asarray(fa_residuals),
                CN_residuals=np.asarray(cn_residuals),
            )
        except Exception as exc:
            save_error = f"{type(exc).__name__}: {exc}"
    save_error = mpi.bcast(save_error, root=0)
    if save_error is not None:
        raise RuntimeError("Physical diagnostic export failed collectively: "
                           + save_error)
    if mpi.rank == 0:
        print(
            "PHYSICAL PRECHECK "
            f"phase={args.phase} source_time={t0:.9g} "
            f"final_time={times[-1]:.9g} "
            f"abs_defect={precheck['maximum_absolute_cumulative_energy_defect']:.9g} "
            f"energy_exchange={exchanged:.9g} "
            f"relative_defect={power_defect:.9g} "
            f"max_step_defect={precheck['maximum_absolute_per_step_energy_defect']:.9g} "
            f"mass_defect={mass_defect:.9g}; "
            "diagnostic trace preserved BEFORE decision", flush=True,
        )

    # Conservative pilot guards (not a claimed general nonlinear-mass
    # CN identity). Historic full ramp observed ~2.8e-5..6.9e-5.
    # The relative Hamiltonian defect is a numerical-model observable,
    # NOT a checkpoint fidelity error. Do not turn an otherwise exact
    # restarted trajectory into a fake checkpoint failure. Preserve the
    # original physical gate as the DEFAULT for scientific model checks;
    # only the explicit restart-parity qualification ignores its magnitude.
    # Mass conservation and native DAE/CN acceptance remain mandatory.
    if (mass_defect > 1e-8 or
            (args.qualification_mode == "physical" and power_defect > 1e-3)):
        raise RuntimeError(
            "Original physical boundary/dissipation/mass balance failed: "
            f"relative power={power_defect:.9g}, mass={mass_defect:.9g}; "
            "NATIVE PETSc trajectory and original per-step power traces "
            "preserved in physical_precheck_{report.json,trace.npz}. "
            "Checkpoint correctly REFUSED; compare with uninterrupted "
            "0->0.024 reference before changing physical acceptance."
        )
    wall0, feedback0 = evaluated_envelopes(0.)
    wallfinal, feedbackfinal = evaluated_envelopes(args.final_time)
    if (abs(wall0) > 1e-14 or abs(feedback0-1.) > 1e-14
            or abs(wallfinal-1.) > 1e-12
            or abs(feedbackfinal-1.) > 1e-14):
        raise RuntimeError("Preparation did not preserve original control clocks")

    # Accepted state is the product. No stationary equilibrium claimed.
    system.save_checkpoint(
        args.checkpoint_out,
        extra={
            "experiment": "controlled_rotating_SWE_preparation",
            "gamma": args.gamma, "nominal_mesh": args.mesh,
            "ramp_duration": RAMP_DURATION,
            "ramp_decay_start": RAMP_DECAY_START,
            "ramp_decay_end": RAMP_DECAY_END,
            "physical_clock": "absolute",
            "previous_checkpoint_sha256": (
                None if args.checkpoint_in is None
                else sha256(args.checkpoint_in.read_bytes())
            ),
            "source_MPI_ranks": int(mpi.size),
            "qualification_mode": args.qualification_mode,
            "max_original_CN_residual": max(cn_residuals),
            "max_original_Fa_residual": max(fa_residuals),
        },
    )
    fingerprint = model_provenance(system)
    report = {
        "schema": "scrimp-v130-rotating-controlled-ic-native-v1",
        "status": "ACCEPTED_COMPATIBLE_ROTATING_DAE_INITIAL_STATE_NOT_STEADY",
        "phase": args.phase, "gamma": args.gamma,
        "qualification_mode": args.qualification_mode,
        "physical_relative_power_gate_applied":
            args.qualification_mode == "physical",
        "physical_energy_defect_retained_not_modified": True,
        "mesh": args.mesh, "MPI": int(mpi.size),
        "physical_time_initial": t0, "physical_time_final": float(times[-1]),
        "dt": args.dt, "accepted_steps": nsteps,
        "PETSc_reason": actual[2], "PETSc_rejected": actual[3],
        "SNES_failures": actual[4],
        "ndof": ndof, "initialization": start_report,
        "original_FE_fingerprint_sha256": fingerprint["fingerprint_sha256"],
        "original_full_state_start_sha256": sha256(states[0].tobytes()),
        "original_full_state_final_sha256": sha256(states[-1].tobytes()),
        "original_Fa_max_l2": max(fa_residuals),
        "original_TSCN_max_l2": max(cn_residuals),
        "minimum_positive_height_FE_coefficient": min_h,
        "original_power_defect_over_exchange": power_defect,
        "original_volume_balance_defect": mass_defect,
        "physical_work": float(traces["supplied"][-1]),
        "physical_dissipation": float(traces["dissipated"][-1]),
        "Hamiltonian_change": float(dH[-1]),
        "final_dynamic_F_with_zero_derivative_l2":
            float(stationarity["dynamic_l2"]),
        "final_dynamic_last_step_dzdt_l2": dzdt,
        "final_wall_envelope": wallfinal,
        "final_normal_feedback_envelope": feedbackfinal,
        "closed_loop_feedback_is_inactive_before_75s": True,
        "stationary_rotating_state_certified": False,
        "full_TS_integrator_history_saved": False,
        "checkpoint_path": str(args.checkpoint_out),
    }
    if mpi.rank == 0:
        args.out.mkdir(parents=True, exist_ok=True)
        (args.out / "initial_state_report.json").write_text(
            json.dumps(report, indent=2) + "\n", encoding="utf-8"
        )
        np.savez_compressed(
            args.out / "initial_state_trace.npz",
            times=times, states=states,
            H=traces["H"], volume=traces["volume"],
            p_normal=traces["p_normal"],
            p_viscous=traces["p_viscous"],
            d_grad=traces["d_grad"], d_div=traces["d_div"],
            boundary_work=traces["supplied"],
            viscous_dissipation=traces["dissipated"],
            power_balance=traces["balance"],
            volume_balance=traces["volume_balance"],
            full_state_start=states[0], full_state_final=states[-1],
        )
        print(
            "ROTATING CONTROLLED IC ACCEPTED "
            f"phase={args.phase} mode={args.qualification_mode} "
            f"gamma={args.gamma:g} "
            f"mesh={args.mesh:g} MPI={mpi.size} "
            f"t={times[-1]:.6g} steps={nsteps} "
            f"Fa={max(fa_residuals):.8g} CN={max(cn_residuals):.8g} "
            f"energy_exchange_error={power_defect:.8g} "
            f"F_d_stationarity={stationarity['dynamic_l2']:.8g} "
            "stationary=NOT_CLAIMED",
            flush=True,
        )
    mpi.Barrier()


if __name__ == "__main__":
    run()
