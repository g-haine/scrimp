# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""#28: native short physical TS after production nonlinear-affine initialization.

Load the published original examples/shallow_water.py, suppress ONLY its
original long TS and figures during construction; retain its original
rectangle, physical FEM, weak forms, initial h0,p0, four zero normal
controls, and explicit gyro. Then opt into the now native-qualified
production nonlinear_affine initialization and run a configurable
number of short genuine BDF2 steps at a fixed physical final time.
The default remains the original 3-step smoke (dt=1e-4, tf=3e-4).

This is a first actual physical-time diagnostic, not an exact discrete
Hamiltonian-conservation theorem and not qualification of the higher-index
viscous rotating SWE traction Y0. No modifications to the PDE.
"""

import argparse
import importlib.util
import json
from pathlib import Path
from unittest.mock import patch

import numpy as np
from mpi4py import MPI
from petsc4py import PETSc

import scrimp as S


TF = 3e-4
DT = 1e-4
DOMAIN_REGION = 1
BORDER_REGIONS = (10, 11, 12, 13)


def build_original_without_published_long_run():
    """Use the actual shipped model, skip only 5000-step solve and plots."""
    original = Path(__file__).resolve().parents[1] / "examples/shallow_water.py"
    spec = importlib.util.spec_from_file_location(
        "scrimp_v130_original_shallow_water_time_gate", original
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    def skip_solve(self):
        return None
    def skip_plot(self, *args, **kwargs):
        return None
    def skip_export(self, *args, **kwargs):
        return None
    with (
        patch.object(S.DPHS, "solve", skip_solve),
        patch.object(S.DPHS, "plot_Hamiltonian", skip_plot),
        patch.object(S.DPHS, "export_to_pv", skip_export),
    ):
        model = module.shallow_water_eq()
    return model


def field_indices(model, name):
    start, width = map(int, model.gf_model.interval_of_variable(name))
    return np.arange(start, start+width, dtype=np.int64)


def quantities(model, expression, region=DOMAIN_REGION):
    out = np.asarray(model.get_quantity(expression, region=region), dtype=float).reshape(-1)
    if len(out) != len(model.solution["t"]):
        raise RuntimeError(
            f"Quantity {expression!r} has {len(out)} values for "
            f"{len(model.solution['t'])} times"
        )
    if not np.all(np.isfinite(out)):
        raise FloatingPointError(f"Nonfinite physical quantity {expression!r}")
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--dt", type=float, default=DT)
    parser.add_argument("--final-time", type=float, default=TF)
    args = parser.parse_args()
    dt = float(args.dt)
    tf = float(args.final_time)
    if dt <= 0. or tf <= 0.:
        parser.error("dt and final-time must both be positive")
    expected_steps = int(round(tf / dt))
    if (expected_steps < 3
            or not np.isclose(expected_steps*dt, tf, atol=1e-13, rtol=0.)):
        parser.error(
            "Need at least 3 exact equal-size time steps; "
            "final-time must be an integer multiple of dt"
        )
    mpi = PETSc.COMM_WORLD.tompi4py()

    model = build_original_without_published_long_run()
    original = np.asarray(
        model.gf_model.from_variables(), dtype=float
    ).copy()
    initial_h = original[field_indices(model, "h")].copy()
    initial_p = original[field_indices(model, "p")].copy()
    if not np.array_equal(initial_p, np.zeros_like(initial_p)):
        raise RuntimeError("Original prescribed p0 is no longer zero")

    model.set_time_scheme(
        # The ORIGINAL published example has a deliberate legacy pseudo
        # initializer. Opt in ONLY to the now independently tested
        # production nonlinear-state / affine-algebraic initializer.
        init_step=True,
        init_step_ts_type="nonlinear_affine",
        init_algebraic_accuracy=1e-11,
        # Preserve original BDF2 family, physical FEM and zero controls.
        # Only shorten the temporal horizon to a small, recorded pilot.
        ts_type="bdf", ts_bdf_order=2,
        t_0=0., t_f=tf, dt=dt, dt_save=dt,
        ts_adapt_type="none",
        ts_adapt_dt_min=1e-12,
        ts_adapt_dt_max=dt,
        ksp_type="preonly", pc_type="lu",
        pc_factor_mat_solver_type="mumps",
    )
    model.solve()

    report = getattr(model, "initialization_report", None)
    if not isinstance(report, dict) or report.get("status") != "PASS":
        raise RuntimeError("Production nonlinear-affine initializer did not pass")
    if report.get("strategy") != "nonlinear_affine":
        raise RuntimeError("Unexpected initial strategy")
    if not report.get("differential_initial_exact"):
        raise RuntimeError("Production initializer did not preserve h0,p0")
    if report["native_final_algebraic_residual_l2"] > report["accuracy_target"]:
        raise RuntimeError("Actual production initial algebraic residual is invalid")

    reason = int(model.last_ts_converged_reason)
    steps = int(model.last_ts_steps)
    rejected = int(model.last_ts_rejections)
    snes_failures = int(model.last_ts_snes_failures)
    actual_tf = float(model.last_ts_time)
    times = np.asarray(model.solution["t"], dtype=float)
    if (reason <= 0 or steps != expected_steps or rejected or snes_failures
            or abs(actual_tf-tf) > 1e-12):
        raise RuntimeError(
            "Physical BDF2 TS did not complete required short-run gate: "
            f"reason={reason}, steps={steps}, rejected={rejected}, "
            f"SNESfail={snes_failures}, t={actual_tf:.12g}"
        )
    if (len(times) < 3 or not np.all(np.isfinite(times))
            or abs(times[0]) > 1e-12
            or abs(times[-1]-tf) > 1e-12
            or not np.all(np.diff(times) > 0)):
        raise RuntimeError(f"Bad initial/final physical trajectory timestamps: {times}")

    zfirst = np.asarray(model.solution["z"][0], dtype=float)
    same_h = bool(np.array_equal(zfirst[field_indices(model,"h")], initial_h))
    same_p = bool(np.array_equal(zfirst[field_indices(model,"p")], initial_p))
    if not mpi.allreduce(same_h and same_p, op=MPI.LAND):
        raise RuntimeError("First physically published state differs from original h0,p0")

    mass = quantities(model, "h")
    kinetic = quantities(model, "0.5*h*(p.p)/rho")
    potential = quantities(model, "0.5*rho*g*h*h")
    H = kinetic + potential
    powers = {}
    for i, region in enumerate(BORDER_REGIONS):
        powers[f"Gamma_{i}"] = quantities(
            model, f"U_{i}*Y_{i}", region=region
        )
    power_sum = sum(powers.values())
    # Independently exercise the actual SCRIMP Hamiltonian term registry,
    # which is configured by the original published example.
    model.compute_Hamiltonian()
    H_registered = np.asarray(model.get_Hamiltonian(),dtype=float).reshape(-1)
    if H_registered.size != H.size or not np.allclose(
        H_registered, H, rtol=1e-11, atol=1e-7
    ):
        raise RuntimeError("Registered Hamiltonian differs from original weak-form energy")

    m0 = float(mass[0])
    H0 = float(H[0])
    if m0 <= 0 or H0 <= 0:
        raise RuntimeError("Unexpected nonphysical original initial integrals")
    relative_mass_drift = float(np.max(np.abs(mass-m0))/max(abs(m0),1e-15))
    relative_H_drift = float(np.max(np.abs(H-H0))/max(abs(H0),1e-15))
    max_abs_power = float(np.max(np.abs(power_sum)))
    cumulative_power = np.zeros(len(times), dtype=float)
    cumulative_power[1:] = np.cumsum(
        0.5*(power_sum[1:]+power_sum[:-1])*np.diff(times)
    )
    # No viscosity in this original inviscid model. Total input power
    # is zero under homogeneous boundary controls. BDF2 is not exactly
    # energy preserving, so record its *actual* discrete energy defect,
    # including its normalization by energy CHANGES rather than only H0.
    energy_balance_defect = H - H0 - cumulative_power
    exchanged_energy_scale = max(
        float(np.max(np.abs(H-H0))),
        float(np.max(np.abs(cumulative_power))),
        1e-14,
    )
    relative_balance_to_exchange = float(
        np.max(np.abs(energy_balance_defect))/exchanged_energy_scale
    )
    relative_balance_to_H0 = float(
        np.max(np.abs(energy_balance_defect))/max(abs(H0), 1e-14)
    )
    physical_gate_failures = []
    # BDF2 is NOT exactly energy preserving. These are initial SHORT
    # physical sanity gates, NOT an exact power-balance certification.
    if relative_mass_drift > 1e-5:
        physical_gate_failures.append(
            f"relative mass drift {relative_mass_drift:.9g} > 1e-5"
        )
    if relative_H_drift > 5e-3:
        physical_gate_failures.append(
            f"relative H drift {relative_H_drift:.9g} > 5e-3"
        )
    if max_abs_power > 1e-7:
        physical_gate_failures.append(
            f"max abs homogeneous boundary power {max_abs_power:.9g} > 1e-7"
        )

    result = {
        "schema": "scrimp-v130-issue28-original-shallow-short-physical-v1",
        "status": (
            "FAILED_SHORT_PHYSICAL_SANITY" if physical_gate_failures
            else "PASS_SHORT_PHYSICAL_TIME_GATE"
        ),
        "model": "original shipped examples/shallow_water.py",
        "mpi_ranks": int(mpi.size),
        "mesh_h": 0.1,
        "t_final_requested": tf,
        "time_step_requested": dt,
        "expected_accepted_steps": expected_steps,
        "t_final_actual": actual_tf,
        "saved_times": times.tolist(),
        "ts_converged_reason": reason,
        "accepted_steps": steps,
        "rejected_steps": rejected,
        "nonlinear_solver_failures": snes_failures,
        "production_initializer_report": report,
        "original_prescribed_h0_p0_bitwise": True,
        "original_boundary_controls_unchanged_zero": True,
        "original_explicit_gyro_unchanged": bool(
            model.bricks["G"].get_explicit()
        ),
        "mass": mass.tolist(),
        "kinetic_energy": kinetic.tolist(),
        "potential_energy": potential.tolist(),
        "hamiltonian": H.tolist(),
        "boundary_powers": {k: v.tolist() for k,v in powers.items()},
        "total_boundary_power": power_sum.tolist(),
        "relative_mass_drift": relative_mass_drift,
        "relative_hamiltonian_drift": relative_H_drift,
        "maximum_abs_boundary_power": max_abs_power,
        "cumulative_boundary_power": cumulative_power.tolist(),
        "energy_balance_defect": energy_balance_defect.tolist(),
        "signed_final_energy_balance_defect": float(energy_balance_defect[-1]),
        "relative_abs_final_energy_balance_defect_to_H0": float(
            abs(energy_balance_defect[-1]) / max(abs(H0), 1e-14)
        ),
        "saved_sample_count": int(len(times)),
        "max_relative_energy_balance_to_H0": relative_balance_to_H0,
        "max_relative_energy_balance_to_exchanged_energy": relative_balance_to_exchange,
        "physical_sanity_gate_failures": physical_gate_failures,
        "Hamiltonian_GetFEM_and_registered_agree": True,
        "BDF2_exact_power_conservation_claimed": False,
        "physical_spatial_convergence_qualified": False,
        "original_rotating_viscous_Y0_qualified": False,
    }
    if mpi.rank == 0:
        args.out.mkdir(parents=True, exist_ok=True)
        (args.out/"shallow_short_physical.json").write_text(
            json.dumps(result,indent=2)+"\n", encoding="utf-8"
        )
        print(
            f"ISSUE28 SHALLOW SHORT PHYSICAL MPI={mpi.size} "
            f"{result['status']} steps={steps} t={actual_tf:.9g} "
            f"mass_drift={relative_mass_drift:.9g} "
            f"H_drift={relative_H_drift:.9g} "
            f"relative_energy_defect_exchanged={relative_balance_to_exchange:.9g} "
            f"max_boundary_power={max_abs_power:.9g}",
            flush=True,
        )
    mpi.Barrier()
    if physical_gate_failures:
        raise AssertionError(
            "Native short-time physical sanity gate failed (report retained): "
            + "; ".join(physical_gate_failures)
        )


if __name__ == "__main__":
    main()
