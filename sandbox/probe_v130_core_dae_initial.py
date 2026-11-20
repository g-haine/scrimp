# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Native qualification of PRODUCTION DPHS.init_step algebraic opt-in.

Unlike the immutable prototype, this script does NOT monkeypatch
init_step. It sets init_step_ts_type=algebraic in SCRIMP's supported
time-scheme options, then compares true GetFEM/PETSc initial and
physical CN results against the frozen independently qualified
diagnostic. Legacy pseudo and tiny negative-time BDF2 are documented
and available as explicit compatibility smoke controls.
"""

import argparse
import json
import os
from pathlib import Path

from mpi4py import MPI
import numpy as np
import scrimp as S
from scrimp.utils.dae_initialization import InconsistentDAEInitialData

from audit_v130_wave_initializers import _build_manufactured_wave
from check_acoustic_energy_mpi import _audit as audit_port_power
from probe_v130_wave_fixed_dynamics import _build_scenario


def _dynamic_intervals(system):
    return {
        state.get_name(): tuple(
            int(x) for x in system.gf_model.interval_of_variable(
                state.get_name()
            )
        )
        for state in system.states.values()
    }


def _clear_benchmark_initializer_for_implicit_default(
    options, *, enabled
):
    """Clear builder-side options so the next call tests SCRIMP's default.

    PETSc.Options.hasName() reports both materialized defaults and
    explicitly selected values: it cannot distinguish their origin.
    Accept the only two known builder outcomes, historical pseudo
    and the newly materialized algebraic default. An unexpected BDF
    or misspelled choice remains an error. Never alter options for
    explicit-mode control runs.
    """
    if not enabled or not options.hasName("init_step_ts_type"):
        return
    inherited = options.getString(
        "init_step_ts_type", ""
    ).strip().lower()
    if inherited not in ("pseudo", "algebraic"):
        raise RuntimeError(
            "Unexpected inherited benchmark initializer mode: "
            f"{inherited!r}"
        )
    options.delValue("init_step_ts_type")
    if options.hasName("init_step_ts_type"):
        raise RuntimeError(
            "Unable to clear builder's PETSc initializer choice"
        )


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--mode", choices=("algebraic", "pseudo", "bdf2"),
        default="algebraic",
    )
    p.add_argument(
        "--scenario", choices=("polynomial_compatible", "published"),
        default="polynomial_compatible",
    )
    p.add_argument("--mesh", type=float, default=.2)
    p.add_argument("--dt", type=float, default=.002)
    p.add_argument("--steps", type=int, default=10)
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--tiny-bdf-dt", type=float, default=1.e-6)
    p.add_argument(
        "--use-default", action="store_true",
        help="Omit init_step_ts_type entirely to exercise the real "
             "DPHS implicit default; only supported with mode=algebraic",
    )
    args = p.parse_args()
    if args.use_default and args.mode != "algebraic":
        p.error("--use-default requires --mode algebraic")
    if args.mode != "algebraic" and args.scenario != "polynomial_compatible":
        p.error("Legacy controls use only polynomial_compatible positive case")
    if (
        args.mesh <= 0 or args.dt <= 0 or args.steps < 2
        or args.tiny_bdf_dt <= 0
    ):
        p.error("Expected strictly positive mesh, dt and tiny BDF2 step")
    mpi = MPI.COMM_WORLD
    if args.scenario == "polynomial_compatible":
        system = _build_manufactured_wave(
            args.mesh, "0.08*sin(20*t)"
        )
    else:
        system = _build_scenario(args.mesh, "published")

    # These diagnostic builders construct a wave by calling an
    # existing example's set_time_scheme(), which already populates
    # PETSc.Options with either an explicitly selected historical
    # pseudo initializer or SCRIMP's materialized algebraic default.
    # Neither is a keyword in the *new* set_time_scheme() call below.
    # Clear the inherited builder configuration only for the true
    # no-keyword diagnostic. Explicit pseudo/BDF controls are unchanged.
    _clear_benchmark_initializer_for_implicit_default(
        system.time_scheme, enabled=args.use_default
    )

    original = np.asarray(
        system.gf_model.from_variables(), dtype=np.float64
    ).copy()
    dynamic = _dynamic_intervals(system)
    option = "bdf" if args.mode == "bdf2" else args.mode
    kwargs = dict(
        ts_type="cn", ts_adapt_type="none",
        ts_exact_final_time="matchstep",
        dt=args.dt, dt_save=args.dt,
        t_f=args.dt * args.steps,
        ksp_type="preonly", pc_type="lu",
        pc_factor_mat_solver_type="mumps",
        init_step=True,
    )
    if not args.use_default:
        kwargs["init_step_ts_type"] = option
    if args.mode == "algebraic":
        kwargs.update(
            init_algebraic_atol=1.e-10,
            init_algebraic_rtol=1.e-8,
            init_algebraic_accuracy=1.e-11,
            init_algebraic_max_passes=3,
        )
    elif args.mode == "bdf2":
        kwargs.update(
            init_step_bdf_order=2,
            init_step_dt=args.tiny_bdf_dt,
            init_step_from_negative_time=True,
        )
    # Test which options the *caller supplied*; PETSc.Options.hasName()
    # is true for both caller-provided and SCRIMP-materialized defaults.
    explicit_present = "init_step_ts_type" in kwargs
    if args.use_default and explicit_present:
        raise RuntimeError(
            "Implicit-default regression accidentally supplied "
            "an init_step_ts_type keyword argument"
        )
    system.set_time_scheme(**kwargs)
    effective_strategy = system.time_scheme.getString(
        "init_step_ts_type", ""
    ).strip().lower()
    if effective_strategy != option:
        raise RuntimeError(
            f"Production PETSc.Options selected {effective_strategy!r}"
            f", expected {option!r}"
        )
    if not system.time_scheme.hasName("init_step_ts_type"):
        raise RuntimeError(
            "Production failed to register its effective initializer "
            "in PETSc.Options"
        )

    report = {
        "schema": "scrimp-v130-production-dae-initial-native-v1",
        "source_revision": os.environ.get(
            "SCRIMP_BENCH_GIT_SHA", ""
        ),
        "frozen_qualified_image_source": os.environ.get(
            "SCRIMP_IMAGE_REVISION", ""
        ),
        "mpi_ranks": mpi.Get_size(),
        "scenario": args.scenario,
        "mode": args.mode,
        "implicit_default_requested": bool(args.use_default),
        # Caller-supplied kwargs, not PETSc database presence.
        "explicit_init_step_ts_type_present": bool(explicit_present),
        "effective_PETSc_initial_strategy": effective_strategy,
        "time_initial": 0.,
        "physical_steps_requested": args.steps,
        "main_dt": args.dt,
        "tiny_initial_BDF2_dt": (
            args.tiny_bdf_dt if args.mode == "bdf2" else None
        ),
        "initialization_options": {
            k: str(v) for k, v in kwargs.items()
            if k.startswith("init_")
        },
    }
    try:
        system.solve()
    except InconsistentDAEInitialData as exc:
        if (
            args.mode != "algebraic" or
            args.scenario != "published"
        ):
            raise
        info = getattr(system, "initialization_report", None)
        if not info or info["status"] != "UNRESOLVED":
            raise RuntimeError(
                "Production rejected Gaussian without native residual report"
            ) from exc
        report.update({
            "status": "EXPECTED_REJECTION",
            "exception": str(exc),
            "init_report": info,
            "physical_cn_executed": False,
            "initial_data_preserved": bool(
                np.array_equal(
                    np.asarray(system.gf_model.from_variables()),
                    original,
                )
            ),
        })
        if not report["initial_data_preserved"]:
            raise RuntimeError(
                "Rejected Gaussian left modified GetFEM initial state"
            )
    else:
        if args.scenario == "published":
            raise RuntimeError(
                "Gaussian unexpectedly passed; investigate physical "
                "feasibility before updating the baseline"
            )
        times = np.asarray(system.solution["t"], dtype=float)
        states = np.asarray(system.solution["z"], dtype=np.float64)
        expected = np.arange(args.steps + 1) * args.dt
        if (
            times.shape != expected.shape
            or not np.allclose(
                times, expected, atol=1e-10, rtol=0
            )
            or states.shape != (len(expected), len(original))
        ):
            raise RuntimeError("Incomplete production CN snapshots")
        if (
            system.last_ts_converged_reason <= 0
            or system.last_ts_steps != args.steps
            or system.last_ts_rejections
            or system.last_ts_snes_failures
        ):
            raise RuntimeError(
                "Production CN TS did not converge to the expected grid"
            )
        drift = {
            name: float(np.linalg.norm(
                states[0, start:start+count]
                - original[start:start+count]
            ))
            for name, (start, count) in dynamic.items()
        }
        if args.mode == "algebraic" and any(drift.values()):
            raise RuntimeError(
                "Production algebraic initializer changed prescribed "
                f"dynamic coefficients: {drift}"
            )
        power = audit_port_power(system, mpi)
        tol = float(
            1e-9 + 1e-5 * power["balance_reference_scale"]
        )
        if (
            power["absolute_port_work"] <= 1e-9
            or power["maximum_absolute_balance_defect"] > tol
        ):
            raise RuntimeError(
                "Production initial state failed non-vacuous CN "
                "Hamiltonian+signed port work gate"
            )
        info = getattr(system, "initialization_report", None)
        if args.mode == "algebraic":
            if not info or info["status"] != "PASS":
                raise RuntimeError(
                    "Production init report missing or not qualified"
                )
            if info["native_final_residual_l2"] > info[
                "feasibility_tolerance"
            ]:
                raise RuntimeError("True initial IFunction residual fails")
            if not info["accuracy_target_met"]:
                raise RuntimeError(
                    "Production projected t0 residual misses strict "
                    "separate accuracy goal"
                )
        report.update({
            "status": "PASS",
            "init_report": info,
            "physical_cn_executed": True,
            "initial_dynamic_drift_l2": drift,
            "full_differential_initial_exact": (
                args.mode == "algebraic"
            ),
            "physical_times": times.tolist(),
            "physical_power": power,
            "power_balance_tolerance": tol,
        })
        if mpi.Get_rank() == 0:
            args.output_dir.mkdir(parents=True, exist_ok=True)
            stem = f"core_{args.mode}_{args.scenario}_{mpi.Get_size()}r"
            np.savez_compressed(
                args.output_dir / f"{stem}.npz",
                times=times, states=states, initial=original,
                projected=states[0].copy(),
            )
    if mpi.Get_rank() == 0:
        args.output_dir.mkdir(parents=True, exist_ok=True)
        stem = f"core_{args.mode}_{args.scenario}_{mpi.Get_size()}r"
        (args.output_dir / f"{stem}.json").write_text(
            json.dumps(report, indent=2) + "\n", encoding="utf-8"
        )
        if report["status"] == "PASS":
            print(
                f"CORE DAE {args.mode} MPI={mpi.Get_size()} "
                f"{args.scenario}: PASS; "
                f"drift={report['initial_dynamic_drift_l2']}, "
                f"max power defect={report['physical_power']['maximum_absolute_balance_defect']:.9g}, "
                f"work={report['physical_power']['port_work_integral'][-1]:.9g}",
                flush=True,
            )
        else:
            print(
                f"CORE DAE {args.mode} MPI={mpi.Get_size()} "
                f"{args.scenario}: EXPECTED_REJECTION; "
                f"F_a={report['init_report']['native_final_residual_l2']:.9g}, "
                f"tolerance={report['init_report']['feasibility_tolerance']:.9g}",
                flush=True,
            )


if __name__ == "__main__":
    main()
