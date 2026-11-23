# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""#28/#6 physical C1 spinup: FULL RAMP PETSc TSCN refinement pilot.

Independent sandbox study, never a production high-index initializer:
- original curved P2 viscous SWE and h-weighted physical wall powers,
- explicit NEW consistent stationary q0 at t0, unchanged physical bricks,
- t_final=0.02 fixed: full C1 wall acceleration, steps=10/20/40,
- truly accepted PETSc TSCN steps and real GetFEM IFunction checked
  at *each* endpoint, not an incorrect implicit midpoint residual,
- endpoint algebraic constraints + GetFEM Hamiltonian, wall work,
  grad/div viscous dissipation and mass law at each saved time,
- store actual physical FE coefficients and phase-resolved power traces.

PETSc TSCN endpoint assumes constant linear mass, whereas SWE M(h)
depends on h. Even good TSCN/energy residuals do not automatically
qualify time accuracy or a general nonlinear-mass integrator.
"""

import argparse
import json
from pathlib import Path
from unittest.mock import patch

import numpy as np
from mpi4py import MPI
from petsc4py import PETSc
import scrimp as S

from viscous_swe_2d_mpi import (
    build_system, diagnostics, ROTATION, H_INITIAL,
)
from probe_v130_issue28_structure import indices, layout
from probe_v130_issue28_one_cn_spinup import (
    native_residual, publish_sandbox_rest,
)
from probe_v130_issue28_weak_wall import ensure_real

FINAL_TIME = 0.02
COUNTS = (10, 20, 40)
RAMP_DURATION = 0.02


def numeric_trace(seq, n, name):
    arr = np.asarray(seq, dtype=float)
    if arr.shape != (n,) or not np.all(np.isfinite(arr)):
        raise RuntimeError(f"Native physical {name}: wrong shape/nonfinite values")
    return arr


def run():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mesh", type=float, required=True, choices=(.2, .1))
    parser.add_argument("--steps", type=int, required=True, choices=COUNTS)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    mpi = PETSc.COMM_WORLD.tompi4py()
    dt = FINAL_TIME / args.steps

    # All runs model the SAME physical experiment and stop at SAME t_f.
    # Only the temporal grid changes. P2 geometry, viscous and
    # gyroscopic weak bricks remain the original ones.
    system = build_system(
        gamma=0., mesh=args.mesh, geometry_order=2, height_order=2,
        final_time=FINAL_TIME, save_dt=dt*(1.-1e-9), max_step=dt,
        ramp_start=.03, ramp_end=.06,
        init_step=True, ts_type="cn", adapt_type="none",
        init_dt=1e-9, rotation=ROTATION, initial_profile="rigid",
        wall_velocity_profile="facet_tangent",
        algebraic_initial_guess="zero",
        start_from_rest_duration=RAMP_DURATION, init_method="bdf",
    )
    published = {}

    def sandbox_initializer(self):
        publish_sandbox_rest(self, mpi, published)

    with patch.object(S.DPHS, "init_step", sandbox_initializer):
        system.solve()

    initialization = getattr(system, "initialization_report", {})
    if (not isinstance(initialization, dict)
            or initialization.get("getfem_globalvec_published_bitwise") is not True):
        raise RuntimeError("Native TS did not publish exact analytic sandbox t0")
    reason = int(getattr(system, "last_ts_converged_reason", 0))
    steps = int(getattr(system, "last_ts_steps", -1))
    rejected = int(getattr(system, "last_ts_rejections", -1))
    failures = int(getattr(system, "last_ts_snes_failures", -1))
    reached = float(getattr(system, "last_ts_time", float("nan")))
    if (reason <= 0 or steps != args.steps or rejected != 0
            or failures != 0 or not np.isfinite(reached)
            or abs(reached-FINAL_TIME) > 1e-12):
        raise RuntimeError(
            f"Physical TSCN fails accepted-step gate: reason={reason}, "
            f"steps={steps}/{args.steps}, rejected={rejected}, "
            f"SNES={failures}, reached={reached:.17g}"
        )

    times = np.asarray(system.solution["t"], dtype=float)
    states = np.asarray([
        ensure_real(z, "archived native physical state").copy()
        for z in system.solution["z"]
    ])
    _, _, iv, n = layout(system)
    dnames, anames, _, _ = layout(system)
    if (times.shape != (args.steps+1,)
            or states.shape != (args.steps+1, n)
            or not np.all(np.isfinite(states))
            or not np.allclose(
                times, np.linspace(0., FINAL_TIME, args.steps+1),
                rtol=0., atol=5e-13
            )
            or not np.array_equal(states[0], published["candidate"])):
        raise RuntimeError("Native physical TSCN trajectory/archive incomplete")

    d = indices(dnames, iv)
    a = indices(anames, iv)
    hfirst, nh = iv["h"]
    stages = []
    endpoints = []
    out_of_scheme_midpoints = []
    zero = np.zeros(n, dtype=float)
    for k in range(args.steps):
        t0, t1 = float(times[k]), float(times[k+1])
        h = t1-t0
        if not abs(h-dt) <= 1e-12:
            raise RuntimeError("Physical timestep changed without adaptivity")
        z0, z1 = states[k], states[k+1]
        z_dot = (z1-z0)/h
        _, prev = native_residual(
            system, mpi, z0, zero, t0, iv, (dnames, anames),
            with_vector=True
        )
        _, next_residual = native_residual(
            system, mpi, z1, 2*z_dot, t1, iv, (dnames, anames),
            with_vector=True
        )
        # Exact PETSc endpoint theta=0.5 discrete equation. NOTE PLUS.
        discrepancy = prev+next_residual
        stage = {
            "full_l2": float(np.linalg.norm(discrepancy)),
            "dynamic_l2": float(np.linalg.norm(discrepancy[d])),
            "algebraic_l2": float(np.linalg.norm(discrepancy[a])),
        }
        endpoint = native_residual(
            system, mpi, z1, zero, t1, iv, (dnames, anames)
        )
        mid = native_residual(
            system, mpi, .5*(z0+z1), z_dot,
            .5*(t0+t1), iv, (dnames, anames)
        )
        if (not np.isfinite(stage["full_l2"])
                or stage["full_l2"] >= 1e-7
                or endpoint["algebraic_l2"] >= 1e-7):
            raise RuntimeError(
                f"Native PETSc TSCN step {k+1}/{args.steps} violates "
                "the ORIGINAL GetFEM endpoint equations: "
                f"TSCN={stage['full_l2']:.9g}, "
                f"Fa(t1)={endpoint['algebraic_l2']:.9g}"
            )
        if (float(np.min(z1[hfirst:hfirst+nh])) <= 0.
                or not np.all(np.isfinite(z1[hfirst:hfirst+nh]))):
            raise RuntimeError("Physical Saint-Venant height coefficients nonpositive")
        stages.append(stage)
        endpoints.append(endpoint)
        out_of_scheme_midpoints.append(mid["full_l2"])

    # The wall is required to have reached its full C1 target omega=.2
    # at t=.02. The original decreasing envelope remains 1 until .03;
    # the original native port-U control must be at its plateau.
    full_control = native_residual(
        system, mpi, states[0], zero, FINAL_TIME, iv, (dnames, anames)
    )["U_l2"]
    slightly_later_control = native_residual(
        system, mpi, states[0], zero, FINAL_TIME+1e-5,
        iv, (dnames, anames)
    )["U_l2"]
    if (not np.isfinite(full_control) or full_control <= 1e-12
            or not np.isfinite(slightly_later_control)
            or abs(slightly_later_control/full_control-1.) > 1e-10):
        raise RuntimeError(
            "Original physical facet-tangent control has not reached "
            f"the full C1 wall-rotation plateau at t=.02: "
            f"U={full_control:.9g}, U(t+.00001)={slightly_later_control:.9g}"
        )

    # Original h-weighted ports, true SWE dissipative terms, and
    # GetFEM Hamiltonian evaluated on ALL ranks, not rank-0 only.
    powers = diagnostics(system, gamma=0.)
    ntime = args.steps+1
    fields = (
        "H", "kinetic", "potential", "p_normal", "p_viscous",
        "d_grad", "d_div", "volume", "volume_balance",
        "supplied", "dissipated", "balance",
        "relative_balance_exchange",
    )
    traces = {key:numeric_trace(powers[key], ntime, key) for key in fields}
    H = traces["H"]
    W = traces["supplied"]
    D = traces["dissipated"]
    balance = H-H[0]+D-W
    if np.max(np.abs(balance-traces["balance"])) > 1e-10:
        raise RuntimeError("Native SWE physical power trace disagrees with own components")
    if np.min(traces["d_grad"][1:]+traces["d_div"][1:]) < -1e-12:
        raise RuntimeError("Viscous SWE physical dissipated power became negative")

    change = float(H[-1]-H[0])
    work = float(W[-1])
    diss = float(D[-1])
    defect = float(balance[-1])
    volume_defect = float(traces["volume_balance"][-1])
    scale = max(abs(change), abs(work), abs(diss), 1e-14)
    exchanged_ratio = defect/scale
    final_d = states[-1, d]-states[0, d]
    if float(np.linalg.norm(final_d)) <= 1e-12:
        raise RuntimeError("Physical state has not evolved despite nonzero wall ramp")

    max_tscn = max(x["full_l2"] for x in stages)
    max_fa = max(x["algebraic_l2"] for x in endpoints)
    result = {
        "schema": "scrimp-v130-issue28-full-spinup-v1",
        "status": "NATIVE_TSCN_FULL_C1_RAMP_DIAGNOSTIC",
        "nominal_mesh": args.mesh, "mpi_ranks": mpi.size,
        "target_rotation": ROTATION,
        "initial_height": H_INITIAL,
        "initial_momentum_identically_zero": True,
        "geometry_order": 2, "height_order": 2,
        "original_physical_facet_tangent_wall": True,
        "original_viscous_and_gyro_weak_bricks": True,
        "initial_analytic_publication": initialization,
        "original_incompatible_rotating_initial_q0_untouched": True,
        "end_of_C1_ramp_reached": True,
        "native_full_target_wall_control_U_l2": float(full_control),
        "native_post_ramp_plateau_wall_control_U_l2": float(slightly_later_control),
        "original_ramp_start_of_decreasing_envelope": .03,
        "physical_scheme": "petsc_tscn_endpoint_theta_half",
        "nonlinear_mass_PETSc_endpoint_assumption_violated": True,
        "nonlinear_mass_time_accuracy_qualified": False,
        "final_time_fixed": FINAL_TIME,
        "fixed_step": dt,
        "requested_steps": args.steps,
        "actual_total_dofs": n,
        "differential_dofs": int(len(d)),
        "dynamic_variable_intervals": {
            name: [int(iv[name][0]), int(iv[name][1])]
            for name in dnames
        },
        "actual_TS_time": reached,
        "ts_reason": reason,
        "ts_accepted_steps": steps,
        "ts_rejected_steps": rejected,
        "ts_snes_failures": failures,
        "all_full_native_IFunction_TSCN_step_residuals": stages,
        "all_endpoint_original_algebraic_residuals": [
            x["algebraic_l2"] for x in endpoints
        ],
        "all_implicit_midpoint_OUT_OF_SCHEME_residuals":
            out_of_scheme_midpoints,
        "maximum_native_TSCN_equation_residual_l2": max_tscn,
        "maximum_endpoint_native_algebraic_residual_l2": max_fa,
        "maximum_implicit_midpoint_OUT_OF_SCHEME_residual_l2":
            float(max(out_of_scheme_midpoints)),
        "final_dynamic_state_change_l2": float(np.linalg.norm(final_d)),
        "minimum_height_coefficient_over_trajectory": float(
            np.min(states[:, hfirst:hfirst+nh])
        ),
        "initial_hamiltonian": float(H[0]),
        "final_hamiltonian": float(H[-1]),
        "final_hamiltonian_change": change,
        "integrated_physical_wall_work": work,
        "integrated_viscous_dissipation": diss,
        "integrated_viscous_grad": float(np.sum(
            np.diff(times)*traces["d_grad"][1:]
        )),
        "integrated_viscous_div": float(np.sum(
            np.diff(times)*traces["d_div"][1:]
        )),
        "integrated_normal_port_work": float(np.sum(
            np.diff(times)*traces["p_normal"][1:]
        )),
        "integrated_tangent_port_work": float(np.sum(
            np.diff(times)*traces["p_viscous"][1:]
        )),
        "final_physical_power_balance_defect": defect,
        "maximum_physical_power_balance_defect_l2":
            float(np.max(np.abs(balance))),
        "final_defect_over_actual_energy_exchange": exchanged_ratio,
        "final_volume_balance_defect": volume_defect,
        "maximum_volume_balance_defect": float(
            np.max(np.abs(traces["volume_balance"]))
        ),
        "time_order_qualified": False,
        "physical_energy_balance_convergence_qualified": False,
        "multistep_rotating_trajectory_qualified": False,
        "general_high_index_initializer_qualified": False,
        "physical_power_audit_used_h_weighted_boundary_terms": True,
        "original_MatNest_first_derivative_gate_previously_qualified": True,
        "interpretation": (
            "True original GetFEM endpoint TSCN equations and algebraic "
            "constraints checked at EACH accepted step on one "
            "fixed-final-time refinement run. Original physical weighted "
            "Hamiltonian/port-work/viscous dissipation and endpoint "
            "state saved for a separate strict MPI and time-refinement "
            "analysis. PETSc endpoint nonlinear-mass accuracy and pH "
            "physical energy time convergence remain unqualified."
        ),
    }
    if mpi.rank == 0:
        args.out.mkdir(parents=True, exist_ok=True)
        (args.out/"full_spinup.json").write_text(
            json.dumps(result, indent=2)+"\n", encoding="utf-8"
        )
        np.savez_compressed(
            args.out/"full_spinup.npz",
            times=times, states=states, **traces,
        )
        print(
            f"ISSUE28 FULL SPINUP h={args.mesh:g} MPI={mpi.size} "
            f"n={args.steps} dt={dt:.9g} "
            f"TSCN_max={max_tscn:.9g} Fa_max={max_fa:.9g} "
            f"dH={change:.9g} W={work:.9g} D={diss:.9g} "
            f"balance={defect:.9g} balance/exchange={exchanged_ratio:.9g} "
            f"mass_defect={volume_defect:.9g}",
            flush=True
        )
    mpi.Barrier()


if __name__ == "__main__":
    run()
