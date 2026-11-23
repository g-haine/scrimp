# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""#28/#6: sandbox-only analytic startup and one real PETSc TSCN step.

New compatible rest-to-spinup experiment, original viscous SWE bricks.
Publish already native-qualified rest coefficients in the isolated
sandbox init hook, preserving differential q0 and production source.

PETSc ts_type=cn is ENDPOINT theta=0.5 TRAPEZOIDAL, not implicit
midpoint. Exact PETSc equation at one step:
F(t1,z1,2*(z1-z0)/dt) + F(t0,z0,0) = 0.
Check this equation with the original unmodified IFunction; retain
the midpoint residual as OUT-OF-SCHEME diagnostic only.

PETSc explicitly cautions that TSCN endpoint assumes linear
time-independent mass, unlike nonlinear SWE M(h). A small
TSCN equation residual is a solver qualification ONLY, not
a proof of physical order, Hamiltonian/power balance or
general high-index DAE initialisation.
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
    build_system, diagnostics, H_INITIAL, RHO, GRAVITY,
    ROTATION,
)
from probe_v130_issue28_structure import indices, layout
from probe_v130_issue28_weak_wall import ensure_real

STEP = 0.002
RAMP_DURATION = 0.02


def gathered(vec, mpi):
    return np.concatenate(mpi.allgather(
        ensure_real(vec.getArray(readonly=True), "native PETSc vector").copy()
    ))


def native_residual(system, mpi, values, derivative, time, iv, names, *, with_vector=False):
    """Evaluate original IFunction; optionally return full native PETSc coefficients."""
    n = int(system.gf_model.nbdof())
    if values.shape != (n,) or derivative.shape != (n,):
        raise RuntimeError("Native physical residual wrong FE dimensions")
    created = []
    try:
        z = system.tangent_stiffness.createVecRight()
        created.append(z)
        zd = z.duplicate()
        created.append(zd)
        f = z.duplicate()
        created.append(f)
        lo, hi = z.getOwnershipRange()
        if hi > lo:
            ids = np.arange(lo, hi, dtype=PETSc.IntType)
            z.setValues(ids, values[lo:hi])
            zd.setValues(ids, derivative[lo:hi])
        z.assemble()
        zd.assemble()
        system.IFunction(None, float(time), z, zd, f)
        full = gathered(f, mpi)
        d = indices(names[0], iv)
        a = indices(names[1], iv)
        stats = {
            "full_l2": float(np.linalg.norm(full)),
            "dynamic_l2": float(np.linalg.norm(full[d])),
            "algebraic_l2": float(np.linalg.norm(full[a])),
            "U_l2": float(np.linalg.norm(full[
                slice(iv["U"][0], sum(iv["U"]))
            ])),
        }
        return (stats, full) if with_vector else stats
    finally:
        for obj in reversed(created):
            obj.destroy()


def publish_sandbox_rest(system, mpi, stored):
    dynamic, algebraic, iv, n = layout(system)
    if dynamic != ("h", "p"):
        raise RuntimeError("Unexpected original SWE dynamic partition")
    original = ensure_real(
        system.gf_model.from_variables(), "original prescribed initial state"
    ).copy()
    orig_global = ensure_real(
        system.GlobalVec.getArray(readonly=True), "GetFEM replicated GlobalVec"
    ).copy()
    if original.shape != (n,) or orig_global.shape != (n,):
        raise RuntimeError("Unexpected original replicated P2 state size")
    d = indices(dynamic, iv)
    pfirst, psize = iv["p"]
    hfirst, hsize = iv["h"]
    if not np.array_equal(original[pfirst:pfirst+psize],
                          np.zeros(psize)):
        raise RuntimeError("Physical spinup p0 must be identically zero")
    if not np.array_equal(original[hfirst:hfirst+hsize],
                          np.full(hsize, H_INITIAL)):
        raise RuntimeError("Physical spinup h0 must equal constant 5")
    candidate = original.copy()
    for name, value in (
        ("e_h", RHO*GRAVITY*H_INITIAL),
        ("Y_n", -RHO*GRAVITY*H_INITIAL),
    ):
        first, width = iv[name]
        candidate[first:first+width] = value
    if not np.array_equal(original[d], candidate[d]):
        raise RuntimeError("Sandbox analytic startup modified dynamic q0")
    if np.any(candidate[indices(("Y", "e_p", "U", "U_n"), iv)]):
        raise RuntimeError("Unexpected initially nonzero spinup port state")

    initial = native_residual(
        system, mpi, candidate, np.zeros(n), 0., iv, (dynamic, algebraic)
    )
    if initial["full_l2"] > 1e-10 or not np.isfinite(initial["full_l2"]):
        raise RuntimeError("Candidate analytic t0 is not native DAE-compatible")

    system.gf_model.set_time(0.)
    system.gf_model.to_variables(candidate.copy())
    system.GlobalVec.getArray()[:] = candidate
    gf_published = ensure_real(
        system.gf_model.from_variables(), "published native GetFEM"
    ).copy()
    global_published = ensure_real(
        system.GlobalVec.getArray(readonly=True), "published GlobalVec"
    ).copy()
    if not mpi.allreduce(
        bool(np.array_equal(gf_published, candidate)
             and np.array_equal(global_published, candidate)
             and np.array_equal(gf_published[d], original[d])),
        op=MPI.LAND,
    ):
        raise RuntimeError("Sandbox startup publication is not bitwise exact")
    system.initialization_report = {
        "strategy": "sandbox_exact_stationary_spinup_only",
        "status": "SANDBOX_ONLY_PUBLISHED_NOT_GENERAL_INITIALIZER",
        "initial_native_IFunction": initial,
        "differential_q0_preserved_bitwise": True,
        "getfem_globalvec_published_bitwise": True,
        "no_pseudo_or_BDF_initialization": True,
        "analytic_gprime_t0_exact_zero": True,
        "physical_time_step_not_yet_accepted": True,
    }
    stored["original"] = original
    stored["candidate"] = candidate
    stored["initial_residual"] = initial
    if mpi.rank == 0:
        print(
            "ISSUE28 ONE CN STARTUP "
            f"MPI={mpi.size} native_Fa={initial['algebraic_l2']:.9g} "
            f"Fd={initial['dynamic_l2']:.9g} "
            "q0_preserved=PASS published=PASS",
            flush=True
        )


def _series(values, size, field):
    v = np.asarray(values, dtype=float)
    if v.shape != (size,) or not np.all(np.isfinite(v)):
        raise RuntimeError(f"Invalid physical {field} trace")
    return v


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--mesh", type=float, choices=(.2, .1), required=True)
    p.add_argument("--out", type=Path, required=True)
    args = p.parse_args()
    mpi = PETSc.COMM_WORLD.tompi4py()
    # Use the same original physical weak model and geometry. Only
    # the *declared* physical initial experiment is different.
    system = build_system(
        gamma=0., mesh=args.mesh, geometry_order=2, height_order=2,
        final_time=STEP, save_dt=STEP, max_step=STEP,
        ramp_start=.03, ramp_end=.06,
        init_step=True, ts_type="cn", adapt_type="none",
        init_dt=1e-9, rotation=ROTATION, initial_profile="rigid",
        wall_velocity_profile="facet_tangent",
        algebraic_initial_guess="zero",
        start_from_rest_duration=RAMP_DURATION,
        init_method="bdf",
    )
    published = {}

    def sandbox_initializer(self):
        publish_sandbox_rest(self, mpi, published)

    # Patching applies ONLY within this isolated diagnostic Python
    # process. The production source is mounted read-only by the runner.
    with patch.object(S.DPHS, "init_step", sandbox_initializer):
        system.solve()

    report = getattr(system, "initialization_report", None)
    if not isinstance(report, dict) or not report.get(
        "getfem_globalvec_published_bitwise", False
    ):
        raise RuntimeError("CN run lost its exact sandbox startup report")
    reason = int(getattr(system, "last_ts_converged_reason", 0))
    accepted_steps = int(getattr(system, "last_ts_steps", -1))
    rejected_steps = int(getattr(system, "last_ts_rejections", -1))
    snes_failures = int(getattr(system, "last_ts_snes_failures", -1))
    actual_end = float(getattr(system, "last_ts_time", float("nan")))
    if (reason <= 0 or accepted_steps != 1
            or rejected_steps != 0 or snes_failures != 0
            or not np.isfinite(actual_end) or abs(actual_end-STEP) > 1e-12):
        raise RuntimeError(
            f"Real CN did not accept exactly one physical step: "
            f"reason={reason}, steps={accepted_steps}, "
            f"rejections={rejected_steps}, SNESfail={snes_failures}, "
            f"actual_t={actual_end:.12g}"
        )
    times = np.asarray(system.solution["t"], dtype=float)
    samples = [ensure_real(z, "stored physical state").copy()
               for z in system.solution["z"]]
    if (times.shape != (2,) or len(samples) != 2
            or abs(times[0]) > 1e-12
            or abs(times[1] - STEP) > 1e-12):
        raise RuntimeError(
            f"One-step CN did not archive two physical endpoints: {times}"
        )
    initial, final = samples
    candidate = published["candidate"]
    if not mpi.allreduce(bool(np.array_equal(initial, candidate)), op=MPI.LAND):
        raise RuntimeError("First stored physical state differs from published t0")

    dnames, anames, iv, n = layout(system)
    d = indices(dnames, iv)
    qdot = (final-initial)/STEP
    mid = (final+initial)/2
    # PETSc TSCN is the ENDPOINT trapezoidal theta method (theta=.5).
    # From PETSc's theta.c: shift=1/(theta*dt)=2/dt, and
    # SNES receives affine RHS ((theta-1)/theta)*F(t0,z0,0)=-F0.
    # The TRUE equation is F(t1,z1,2*(z1-z0)/dt) + F(t0,z0,0)=0.
    # This is NOT the implicit midpoint scheme; PETSc additionally
    # warns that its endpoint method assumes constant linear mass,
    # while our actual Saint-Venant M(h) is state dependent.
    zero = np.zeros(n, dtype=float)
    initial_zero, v_initial_zero = native_residual(
        system, mpi, initial, zero, 0., iv, (dnames, anames),
        with_vector=True,
    )
    final_twice_dot, v_final_twice_dot = native_residual(
        system, mpi, final, 2.*qdot, STEP, iv, (dnames, anames),
        with_vector=True,
    )
    cn_vec = v_final_twice_dot + v_initial_zero
    a = indices(anames, iv)
    cn_equation = {
        "full_l2": float(np.linalg.norm(cn_vec)),
        "dynamic_l2": float(np.linalg.norm(cn_vec[d])),
        "algebraic_l2": float(np.linalg.norm(cn_vec[a])),
    }
    # OUT-OF-SCHEME midpoint is information, never a TSCN gate.
    midpoint_diagnostic = native_residual(
        system, mpi, mid, qdot, STEP/2, iv, (dnames, anames)
    )
    endpoint_zero_dot = native_residual(
        system, mpi, final, zero, STEP, iv, (dnames, anames)
    )
    if (not np.isfinite(cn_equation["full_l2"])
            or cn_equation["full_l2"] > 1e-7):
        raise RuntimeError(
            "Independent original IFunction does not satisfy actual PETSc "
            "TSCN endpoint equation: "
            f"norm={cn_equation['full_l2']:.9g} "
            f"Fd={cn_equation['dynamic_l2']:.9g} "
            f"Fa={cn_equation['algebraic_l2']:.9g}"
        )
    if (not np.isfinite(endpoint_zero_dot["algebraic_l2"])
            or endpoint_zero_dot["algebraic_l2"] > 1e-7):
        raise RuntimeError(
            "Endpoint algebraic constraint fails after PETSc TSCN: "
            f"Fa={endpoint_zero_dot['algebraic_l2']:.9g}"
        )

    # The independent original SWE diagnostic computes the weighted
    # physical wall powers, NOT generic (unweighted) port power.
    powers = diagnostics(system, gamma=0.)
    for key in (
        "H", "kinetic", "potential", "p_normal", "p_viscous",
        "d_grad", "d_div", "volume", "volume_balance",
        "supplied", "dissipated", "balance",
        "relative_balance_exchange",
    ):
        powers[key] = _series(powers[key], 2, key)
    h = powers["H"]
    change = float(h[1]-h[0])
    work = float(powers["supplied"][1])
    diss = float(powers["dissipated"][1])
    defect = float(change + diss-work)
    scale = max(abs(change), abs(work), abs(diss), 1e-14)
    if not np.all(np.isfinite(
        [change, work, diss, defect, scale]
    )) or diss < -1e-12:
        raise RuntimeError("Nonfinite or negative physical viscous dissipation")
    hstart, nh = iv["h"]
    height = final[hstart:hstart+nh]
    if not np.isfinite(height).all() or float(height.min()) <= 0:
        raise RuntimeError("Nonpositive/invalid physical SWE height coefficients")

    # DO NOT claim physical power-balance convergence from one dt.
    # The defect is measured relative to actual exchanged energy;
    # h-weighted wall work and separate grad/div dissipation are saved.
    result = {
        "schema": "scrimp-v130-issue28-one-cn-spinup-v1",
        "status": "NATIVE_SINGLE_PHYSICAL_CN_STEP_DIAGNOSTIC_ONLY",
        "source_mesh": args.mesh, "mpi_ranks": mpi.size,
        "new_physical_rest_spinup_experiment": True,
        "original_incompatible_rotation_initial_state_unchanged": True,
        "geometry_order": 2, "height_order": 2,
        "physical_facet_tangent_wall": True, "rotation_after_spinup": ROTATION,
        "spinup_duration": RAMP_DURATION,
        "physical_scheme": "cn", "physical_time_step": STEP,
        "initial": report,
        "ts_converged_reason": reason,
        "ts_accepted_steps": accepted_steps,
        "ts_rejected_steps": rejected_steps,
        "ts_snes_failures": snes_failures,
        "ts_actual_final_time": actual_end,
        "petsc_TSCN_form": "theta_0p5_endpoint_trapezoidal",
        "petsc_TSCN_shift": 2.0/STEP,
        "petsc_TSCN_affine_old_F_coefficient": -1.0,
        "nonlinear_mass_PETSc_endpoint_assumption_violated": True,
        "nonlinear_mass_time_accuracy_qualified": False,
        "original_IFunction_t0_zero_dot_l2": initial_zero,
        "original_IFunction_t1_twice_difference_dot_l2": final_twice_dot,
        "original_IFunction_TSCN_endpoint_equation_l2": cn_equation,
        "original_IFunction_midpoint_OUT_OF_SCHEME_l2":
            midpoint_diagnostic,
        "original_IFunction_endpoint_zero_dot_l2": endpoint_zero_dot,
        "first_physical_state_changed": bool(
            not np.array_equal(final[d], initial[d])
        ),
        "dynamic_step_change_l2": float(np.linalg.norm(final[d]-initial[d])),
        "first_step_minimum_height_coefficient": float(height.min()),
        "hamiltonian_initial": float(h[0]),
        "hamiltonian_final": float(h[1]),
        "hamiltonian_change": change,
        "physical_boundary_normal_work": float(
            STEP * powers["p_normal"][1]
        ),
        "physical_boundary_viscous_work": float(
            STEP * powers["p_viscous"][1]
        ),
        "integrated_boundary_work": work,
        "integrated_viscous_dissipation": diss,
        "integrated_dissipation_grad": float(
            STEP * powers["d_grad"][1]
        ),
        "integrated_dissipation_div": float(
            STEP * powers["d_div"][1]
        ),
        "physical_power_balance_defect": defect,
        "physical_power_balance_defect_over_actual_exchange":
            float(defect/scale),
        "physical_volume_balance_defect":
            float(powers["volume_balance"][1]),
        "two_endpoint_states_archived": True,
        "native_original_TSCN_endpoint_equation_verified": True,
        "implicit_midpoint_scheme_used": False,
        "energy_balance_convergence_qualified": False,
        "physical_multistep_trajectory_qualified": False,
        "production_high_index_initializer_qualified": False,
        "interpretation": (
            "The exact t0 analytic physical state was published by "
            "an isolated sandbox patch, and one actual PETSc TSCN "
            "endpoint trapezoidal step was accepted. Original IFunction "
            "was tested against F(t1,z1,2*delta/dt)+F(t0,z0,0)=0. "
            "The implicit-midpoint evaluation is OUT-OF-SCHEME only. "
            "PETSc's TSCN endpoint assumes time-independent linear mass, "
            "not satisfied by the nonlinear SWE M(h), so no order, "
            "physical power or high-index production claim is made."
        ),
    }
    if mpi.rank == 0:
        args.out.mkdir(parents=True, exist_ok=True)
        (args.out/"one_cn_spinup.json").write_text(
            json.dumps(result, indent=2)+"\n", encoding="utf-8"
        )
        np.savez_compressed(
            args.out/"one_cn_spinup_endpoints.npz",
            times=times, initial=initial, final=final,
            H=powers["H"], p_normal=powers["p_normal"],
            p_viscous=powers["p_viscous"],
            d_grad=powers["d_grad"], d_div=powers["d_div"],
        )
        print(
            f"ISSUE28 ONE CN PHYSICAL mesh={args.mesh:g} MPI={mpi.size} "
            f"TS_reason={reason} accepted=1 "
            f"TSCN_F={cn_equation['full_l2']:.9g} "
            f"mid_OUT_OF_SCHEME={midpoint_diagnostic['full_l2']:.9g} "
            f"end_Fa={endpoint_zero_dot['algebraic_l2']:.9g} "
            f"dH={change:.9g} W={work:.9g} D={diss:.9g} "
            f"balance={defect:.9g} "
            f"balance/exchange={defect/scale:.9g}",
            flush=True,
        )
    mpi.Barrier()


if __name__ == "__main__":
    main()
