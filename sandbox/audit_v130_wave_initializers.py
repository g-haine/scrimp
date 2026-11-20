# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Diagnostic-only comparison of physical CN starts from two DAE initializers.

The PDE, FE model, nonzero manufactured p0, boundary drive, fixed-CN
scheme and rc3 image are IDENTICAL across modes and MPI rank counts.

Modes:
 * pseudo: unchanged production accepted pseudo-TS is the main CN start.
 * algebraic: run the SAME production pseudo-TS first, but replace
   its accepted GetFEM initial state with a diagnostic algebraic-only
   projection of the original q/p, before the main CN TS is constructed.

This is a controlled physics comparison, NOT an implementation of a new
production initializer. It may only qualify the specific linear wave
model. Report the full states and actual midpoint Hamiltonian/port power
with the previously validated native GetFEM assembly routine.
"""

import argparse
import json
import os
from pathlib import Path
import types
from unittest.mock import patch

from mpi4py import MPI
import numpy as np
import scrimp as S

from check_acoustic_energy_mpi import _audit as audit_port_power
from check_acoustic_energy_mpi import _build_with_left_drive
from audit_v130_wave_dae_initial import _intervals
from probe_v130_wave_fixed_dynamics import (
    EXACT_FE_COMPATIBLE_P0,
    _algebraic_indices,
    _owned_residual_and_block,
)
from solve_v130_dae_algebraic_block import solve_algebraic_block


def _build_manufactured_wave(mesh, drive):
    original = S.DPHS.set_initial_value
    seen = []

    def override_initial(system, name, expression):
        if name == "p":
            seen.append(expression)
            expression = EXACT_FE_COMPATIBLE_P0
        return original(system, name, expression)

    with patch.object(S.DPHS, "set_initial_value", override_initial):
        model = _build_with_left_drive(mesh, drive)
    if len(seen) != 1:
        raise RuntimeError("Expected exactly one published p0 override")
    return model


def _positive_projection(system, initial, dynamic, algebraic, mpi):
    """Independent real PETSc algebraic block; full q/p coefficients frozen."""
    idx = _algebraic_indices(algebraic)
    residual, Kaa = _owned_residual_and_block(
        system, initial, idx, mpi, extract_matrix=True
    )
    if mpi.Get_rank() == 0:
        try:
            correction, details = solve_algebraic_block(Kaa, residual)
            if np.iscomplexobj(correction):
                if np.max(np.abs(np.imag(correction))) > 1e-10:
                    raise ValueError("Complex correction for real wave")
                correction = np.real(correction)
            out = {
                "error": None,
                "correction": np.asarray(
                    correction, dtype=np.float64
                ),
                "details": details,
            }
        except Exception as exc:
            out = {"error": f"{type(exc).__name__}: {exc}"}
    else:
        out = None
    out = mpi.bcast(out, root=0)
    if out["error"]:
        raise RuntimeError(out["error"])

    candidate = np.array(initial, dtype=np.float64, copy=True)
    candidate[idx] += out["correction"]
    real_residual, _ = _owned_residual_and_block(
        system, candidate, idx, mpi
    )
    if mpi.Get_rank() == 0:
        residual_norm = float(np.linalg.norm(real_residual))
        threshold = float(out["details"]["requested_tolerance"])
        validity_error = None
        if residual_norm > threshold:
            validity_error = (
                "Algebraic-only projected model failed its native "
                f"IFunction residual gate: {residual_norm} > {threshold}"
            )
        elif not out["details"]["linearized_compatible_within_tolerance"]:
            validity_error = (
                "LSMR linearized algebraic-only correction did not pass"
            )
        decision = (residual_norm, threshold, validity_error)
    else:
        decision = None
    residual_norm, threshold, validity_error = mpi.bcast(
        decision, root=0
    )
    if validity_error:
        raise RuntimeError(validity_error)

    changes = {}
    for name, (start, count) in dynamic.items():
        diff = candidate[start:start + count] - initial[
            start:start + count
        ]
        changes[name] = float(np.linalg.norm(diff))
        if not np.array_equal(
            candidate[start:start + count],
            initial[start:start + count],
        ):
            raise RuntimeError(
                f"Projected initializer changed prescribed {name} coefficients"
            )

    return candidate, {
        "algebraic_residual_after_project_l2": residual_norm,
        "algebraic_residual_before_project_l2": float(
            out["details"]["raw_constraint_residual_l2"]
        ),
        "native_algebraic_tolerance": threshold,
        "projection_correction_l2": float(
            np.linalg.norm(out["correction"])
        ),
        "dynamic_fields_exactly_preserved": changes,
    }


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--mode", choices=("pseudo", "algebraic"), required=True)
    p.add_argument("--mesh", type=float, default=0.2)
    p.add_argument("--dt", type=float, default=0.002)
    p.add_argument("--steps", type=int, default=10)
    p.add_argument(
        "--drive-left", default="0.08*sin(20*t)",
        help="Nonzero left boundary drive with U_L(0)=0",
    )
    p.add_argument("--output-dir", type=Path, required=True)
    args = p.parse_args()
    if args.mesh <= 0 or args.dt <= 0 or args.steps < 2:
        p.error("Positive mesh/dt and at least two CN steps required")
    if args.drive_left != "0.08*sin(20*t)":
        p.error(
            "This acceptance diagnostic is restricted to the specified "
            "nonzero, initially compatible boundary drive"
        )
    mpi = MPI.COMM_WORLD
    os.environ.setdefault("SCRIMP_INIT_FACTOR_SOLVER_TYPE", "mumps")
    system = _build_manufactured_wave(args.mesh, args.drive_left)
    t_final = args.dt * args.steps
    system.set_time_scheme(
        ts_type="cn", ts_adapt_type="none",
        ts_exact_final_time="matchstep",
        dt=args.dt, dt_save=args.dt, t_f=t_final,
        ksp_type="preonly", pc_type="lu",
        pc_factor_mat_solver_type="mumps",
    )
    original_init = system.init_step
    capture = {}

    def instrumented_initial_step(self):
        _, dynamic, algebraic = _intervals(self)
        initial = np.array(
            self.gf_model.from_variables(), dtype=np.float64, copy=True
        )
        original_init()
        accepted_pseudo = np.array(
            self.gf_model.from_variables(), dtype=np.float64, copy=True
        )
        if args.mode == "algebraic":
            candidate, projection = _positive_projection(
                self, initial, dynamic, algebraic, mpi
            )
        else:
            candidate = accepted_pseudo.copy()
            projection = None

        self.gf_model.set_time(float(self.time_scheme["t_0"]))
        self.gf_model.to_variables(
            np.array(candidate, dtype=np.float64, copy=True)
        )
        drift = {
            name: float(np.linalg.norm(
                accepted_pseudo[start:start + count]
                - initial[start:start + count]
            ))
            for name, (start, count) in dynamic.items()
        }
        capture.update({
            "initial": initial,
            "candidate": candidate,
            "pseudo": accepted_pseudo,
            "projection": projection,
            "pseudo_dynamic_drift_l2": drift,
            "dynamic_intervals": dynamic,
        })

    system.init_step = types.MethodType(
        instrumented_initial_step, system
    )
    system.solve()
    if not capture:
        raise RuntimeError("Initial consistency diagnostic was not invoked")
    if (
        system.last_ts_converged_reason <= 0
        or system.last_ts_steps != args.steps
        or system.last_ts_rejections != 0
        or system.last_ts_snes_failures != 0
    ):
        raise RuntimeError("Main physical CN TS did not meet convergence gate")

    times = np.asarray(system.solution["t"], dtype=float)
    trajectory = np.asarray(system.solution["z"])
    if np.iscomplexobj(trajectory):
        if np.max(np.abs(trajectory.imag)) > 1e-10:
            raise RuntimeError("Unexpected complex wave trajectory")
        trajectory = trajectory.real
    trajectory = np.asarray(trajectory, dtype=np.float64)
    expected = np.arange(args.steps + 1) * args.dt
    if (
        times.shape != expected.shape
        or not np.allclose(times, expected, atol=1e-10, rtol=0)
        or trajectory.shape != (len(expected), len(capture["candidate"]))
    ):
        raise RuntimeError(
            "Missing/nonuniform complete initial+CN trajectory: "
            f"t={times}, z.shape={trajectory.shape}"
        )
    main_start = float(np.linalg.norm(
        trajectory[0] - capture["candidate"]
    ))
    if main_start > 5e-10:
        raise RuntimeError(
            "Main physical TS did not start from intended accepted "
            f"initializer: L2={main_start}"
        )
    dynamic_start_drifts = {}
    for name, (start, count) in capture["dynamic_intervals"].items():
        dynamic_start_drifts[name] = float(np.linalg.norm(
            trajectory[0, start:start + count]
            - capture["initial"][start:start + count]
        ))
    if args.mode == "algebraic":
        if any(x > 1e-12 for x in dynamic_start_drifts.values()):
            raise RuntimeError(
                "Projected main CN did not start from the exact "
                f"prescribed q0/p0: {dynamic_start_drifts}"
            )

    # Same tested Hamiltonian and per-boundary midpoint-power
    # implementation as the existing acoustic physical gate.
    power = audit_port_power(system, mpi)
    if len(power["times"]) != args.steps + 1:
        raise RuntimeError("Port-power audit lost CN snapshots")
    nonvacuous = power["absolute_port_work"] > 1e-9
    if not nonvacuous and mpi.Get_rank() == 0:
        print(
            "SCIENTIFIC WARNING: non-vacuous port-power threshold NOT "
            "MET; the final strict comparator will reject this case.",
            flush=True,
        )
    energy_tolerance = float(
        1e-9 + 1e-5 * power["balance_reference_scale"]
    )
    energy_pass = (
        power["maximum_absolute_balance_defect"]
        <= energy_tolerance
    )
    report = {
        "schema": "scrimp-v130-wave-initializer-physics-v1",
        "image_source_revision": os.environ.get(
            "SCRIMP_IMAGE_REVISION", ""
        ),
        "initializer": args.mode,
        "mpi_ranks": mpi.Get_size(),
        "mesh_h": args.mesh,
        "dt": args.dt,
        "steps": args.steps,
        "final_time": t_final,
        "drive_left": args.drive_left,
        "initial_p0_expression": EXACT_FE_COMPATIBLE_P0,
        "dofs": len(capture["candidate"]),
        "initial_state_mode": (
            "unchanged production accepted pseudo-TS"
            if args.mode == "pseudo"
            else "diagnostic projection of prescribed q0,p0; "
                 "production pseudo-TS still ran but state replaced "
                 "before main CN"
        ),
        "pseudo_dynamic_drift_l2": capture[
            "pseudo_dynamic_drift_l2"
        ],
        "main_initial_dynamic_drift_l2": dynamic_start_drifts,
        "main_start_vs_initializer_l2": main_start,
        "algebraic_projection": capture["projection"],
        "main_ts_converged_reason": int(
            system.last_ts_converged_reason
        ),
        "main_ts_rejected_steps": int(system.last_ts_rejections),
        "main_ts_snes_failures": int(system.last_ts_snes_failures),
        "power": power,
        "absolute_balance_tolerance": energy_tolerance,
        "nonvacuous_port_power": nonvacuous,
        "physical_power_balance_within_tolerance": energy_pass,
        "caution": (
            "This is a diagnostic-only initializer override on a "
            "linear compatible wave. Not a production DAE strategy, "
            "nonlinear proof or demonstration of all DAE index classes."
        ),
    }
    if mpi.Get_rank() == 0:
        args.output_dir.mkdir(parents=True, exist_ok=True)
        stem = f"wave_{args.mode}_{mpi.Get_size()}r"
        (args.output_dir / f"{stem}.json").write_text(
            json.dumps(report, indent=2) + "\n",
            encoding="utf-8",
        )
        np.savez_compressed(
            args.output_dir / f"{stem}.npz",
            times=times,
            trajectory=trajectory,
            prescribed_initial=capture["initial"],
            pseudo_accepted=capture["pseudo"],
            main_initializer=capture["candidate"],
        )
        print(
            f"SCIENCE {args.mode} MPI={mpi.Get_size()}: "
            f"q0 drift={dynamic_start_drifts['q']:.9g}, "
            f"p0 drift={dynamic_start_drifts['p']:.9g}, "
            f"max |H-H0+W|={power['maximum_absolute_balance_defect']:.9g}, "
            f"tolerance={energy_tolerance:.9g}, "
            f"integrated work={power['port_work_integral'][-1]:.9g}, "
            f"absolute work={power['absolute_port_work']:.9g}, "
            f"power balance={'PASS' if energy_pass else 'FAIL'}",
            flush=True,
        )
    # Keep all four modes/ranks available even if one physical power
    # check fails: the independent final comparator applies the strict
    # scientific gate after collecting complete native evidence.
    # This does NOT mark a failed power balance as accepted.
    if not energy_pass and mpi.Get_rank() == 0:
        print(
            "SCIENTIFIC WARNING: physical CN power balance FAILED; "
            "comparison will reject this run after collecting all "
            "MPI and initializer evidence.",
            flush=True,
        )


if __name__ == "__main__":
    main()
