# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Diagnostic-only physical-initial-data and algebraic-DAE wave audit.

Uses the unchanged published examples/wave.py equations through
bench_acoustic_2d._build_wave. Intercepts DPHS.init_step only to read
pre/post GetFEM state and independently evaluate IFunction's algebraic
rows. It does NOT modify the accepted initial vector or solve policy.

A pseudo-time step can reduce algebraic constraint residuals while
drifting away from explicitly prescribed differential initial data.
Report the two properties separately; do NOT label a historical
v1.2.0 trajectory match as proof of physical initialization.

Run via scripts/run-v130-dae-initial-audit.sh; outputs are JSON and
full time-series NPZ for an independent 1-/2-rank comparison.
"""

import argparse
import json
import os
from pathlib import Path
import types

from mpi4py import MPI
import numpy as np
import scrimp as S

from bench_acoustic_2d import _build_wave


def _intervals(system):
    """Obtain all named unknown blocks and their mathematical roles."""
    names = {}
    dynamic = {}
    algebraic = {}
    for state in system.states.values():
        name = state.get_name()
        start, count = map(int, system.gf_model.interval_of_variable(name))
        dynamic[name] = (start, count)
    for port in system.ports.values():
        fields = []
        if port.get_algebraic():
            fields.append(port.get_flow())
            if not port.get_substituted():
                fields.append(port.get_effort())
        elif not port.get_substituted():
            fields.append(port.get_effort())
        for name in fields:
            if name:
                start, count = map(
                    int, system.gf_model.interval_of_variable(name)
                )
                algebraic[name] = (start, count)
    for name, interval in list(dynamic.items()) + list(algebraic.items()):
        if name in names and names[name] != interval:
            raise RuntimeError(f"Inconsistent variable interval: {name}")
        names[name] = interval
    total = int(system.gf_model.nbdof())
    coverage = np.zeros(total, dtype=np.int8)
    for start, count in set(names.values()):
        coverage[start:start + count] += 1
    if np.any(coverage != 1):
        raise RuntimeError(
            "Unclassified/overlapping DAE unknowns; do not infer rows"
        )
    if set(dynamic) & set(algebraic):
        raise RuntimeError("Dynamical and algebraic unknowns overlap")
    return names, dynamic, algebraic


def _sync_global_across_ranks(state, mpi, description):
    """The small wave audit must see identical GetFEM global arrays."""
    gathered = mpi.allgather(np.asarray(state, dtype=np.float64))
    ref = gathered[0]
    differences = [
        np.linalg.norm(other - ref) / max(1.0, np.linalg.norm(ref))
        for other in gathered[1:]
    ]
    deviation = float(max(differences, default=0.0))
    if deviation > 5e-8:
        raise RuntimeError(
            f"{description}: rank disagreement {deviation:.6g}"
        )
    return deviation


def _algebraic_residual(system, state, algebraic, mpi):
    """Evaluate F(t0,z,zd=0) ONLY on owned algebraic residual rows.

    Algebraic mass rows have zero derivatives in the published model;
    no claim is made about a norm of the full differential residual.
    This invokes the existing production IFunction, not an invented
    reconstruction of constitutive and boundary equations.
    """
    t0 = float(system.time_scheme["t_0"])
    vec = system.tangent_stiffness.createVecRight()
    lo, hi = map(int, vec.getOwnershipRange())
    if hi > lo:
        vec.setValues(
            np.arange(lo, hi, dtype=np.int32),
            np.asarray(state[lo:hi], dtype=np.float64),
        )
    vec.assemble()
    zd = vec.duplicate()
    zd.set(0.0)
    residual = vec.duplicate()
    residual.set(0.0)
    # The linear example's GetFEM model is explicitly synchronized here,
    # so time-dependent controls and source RHS see the requested state.
    system.gf_model.set_time(t0)
    system.gf_model.to_variables(
        np.array(state, dtype=np.float64, copy=True)
    )
    try:
        system.IFunction(None, t0, vec, zd, residual)
        entries = residual.getArray(readonly=True)
        alg_indices = np.concatenate([
            np.arange(start, start + size, dtype=np.int64)
            for start, size in algebraic.values()
        ])
        owned = alg_indices[(alg_indices >= lo) & (alg_indices < hi)]
        local = entries[owned - lo]
        squared = float(np.vdot(local, local).real)
        full_squared = float(np.vdot(entries, entries).real)
        norm = float(np.sqrt(mpi.allreduce(squared, op=MPI.SUM)))
        full_norm = float(
            np.sqrt(mpi.allreduce(full_squared, op=MPI.SUM))
        )
        # Structural check: a residual selected as algebraic must have
        # zero rows in the actual distributed mass matrix.
        ptr, _cols, data = system.mass.getValuesCSR()
        local_mass_max = 0.0
        for index in owned:
            row = int(index - lo)
            if ptr[row + 1] > ptr[row]:
                local_mass_max = max(
                    local_mass_max,
                    float(np.max(np.abs(data[ptr[row]:ptr[row + 1]]))),
                )
        max_mass = float(mpi.allreduce(local_mass_max, op=MPI.MAX))
        if max_mass > 1e-12:
            raise RuntimeError(
                "A purported algebraic row carries mass: "
                f"max={max_mass:.6g}"
            )
        return {
            "algebraic_residual_l2": norm,
            "full_residual_at_zero_derivative_l2": full_norm,
            "max_algebraic_mass_entry": max_mass,
            "algebraic_row_count": int(len(alg_indices)),
        }
    finally:
        vec.destroy()
        zd.destroy()
        residual.destroy()


def _block_norms(z, intervals):
    return {
        name: float(np.linalg.norm(z[start:start + size]))
        for name, (start, size) in intervals.items()
    }


def _audit(args, mpi):
    system = _build_wave(args.mesh)
    if not system.initial_value_set.get("q"):
        raise RuntimeError("Published wave q(0) was not prescribed")
    if not system.initial_value_set.get("p"):
        raise RuntimeError("Published wave p(0) was not prescribed")
    system.set_time_scheme(
        ts_type="cn", ts_adapt_type="none",
        ts_exact_final_time="matchstep",
        dt=args.dt, dt_save=args.dt, t_f=args.final_time,
        ksp_type="preonly", pc_type="lu",
        pc_factor_mat_solver_type="mumps",
    )
    captures = {}
    original_init = system.init_step

    def monitored_initial_step(self):
        names, dynamic, algebraic = _intervals(self)
        before = np.asarray(
            self.gf_model.from_variables(), dtype=np.float64
        ).copy()
        before_rank_delta = _sync_global_across_ranks(
            before, mpi, "before pseudo-TS"
        )
        original_init()
        after = np.asarray(
            self.gf_model.from_variables(), dtype=np.float64
        ).copy()
        after_rank_delta = _sync_global_across_ranks(
            after, mpi, "after pseudo-TS"
        )
        before_res = _algebraic_residual(
            self, before, algebraic, mpi
        )
        after_res = _algebraic_residual(
            self, after, algebraic, mpi
        )
        # The IFunction audit is read-only with respect to accepted
        # time integration: restore the genuine pseudo-TS accepted
        # GetFEM state before the main physical TS is constructed.
        self.gf_model.set_time(float(self.time_scheme["t_0"]))
        self.gf_model.to_variables(
            np.array(after, dtype=np.float64, copy=True)
        )
        captures.update({
            "before": before, "after": after,
            "names": names, "dynamic": dynamic,
            "algebraic": algebraic,
            "before_residual": before_res,
            "after_residual": after_res,
            "before_rank_delta": before_rank_delta,
            "after_rank_delta": after_rank_delta,
        })

    system.init_step = types.MethodType(monitored_initial_step, system)
    system.solve()

    if not captures:
        raise RuntimeError("Expected pseudo-TS initialization was skipped")
    if (int(system.last_ts_converged_reason) <= 0
            or int(system.last_ts_steps) != args.steps
            or int(system.last_ts_rejections) != 0
            or int(system.last_ts_snes_failures) != 0):
        raise RuntimeError("Main fixed-CN solver did not pass its step gate")
    times = np.asarray(system.solution["t"], dtype=float)
    trajectory = np.asarray(system.solution["z"])
    if np.iscomplexobj(trajectory):
        imaginary = float(np.max(np.abs(trajectory.imag)))
        if imaginary > 1e-10:
            raise RuntimeError(f"Unphysical imaginary state {imaginary}")
        trajectory = trajectory.real
    trajectory = np.asarray(trajectory, dtype=np.float64)
    expected_times = np.arange(args.steps + 1, dtype=float) * args.dt
    if (trajectory.shape != (args.steps + 1, len(captures["after"]))
            or times.shape != expected_times.shape
            or not np.allclose(times, expected_times, atol=1e-10, rtol=0)):
        raise RuntimeError(
            f"Expected {args.steps + 1} complete fixed-CN snapshots; "
            f"got t={times}, shape={trajectory.shape}"
        )
    after = captures["after"]
    start_difference = float(
        np.linalg.norm(trajectory[0] - after)
        / max(1.0, np.linalg.norm(after))
    )
    if start_difference > 5e-8:
        raise RuntimeError(
            f"Main physical TS did not start at accepted pseudo state: "
            f"{start_difference:.6g}"
        )
    before = captures["before"]
    changes = {}
    physical_fidelity = True
    for name, (start, size) in captures["dynamic"].items():
        initial = before[start:start + size]
        accepted = after[start:start + size]
        abs_change = float(np.linalg.norm(accepted - initial))
        baseline = float(np.linalg.norm(initial))
        scaled_change = abs_change / max(1.0, baseline)
        faithful = bool(
            abs_change <= args.fidelity_atol
            + args.fidelity_rtol * baseline
        )
        physical_fidelity = physical_fidelity and faithful
        changes[name] = {
            "prescribed_l2": baseline,
            "accepted_l2": float(np.linalg.norm(accepted)),
            "absolute_change_l2": abs_change,
            "scaled_change_l2": scaled_change,
            "preserved_within_tolerance": faithful,
        }
    if changes["q"]["prescribed_l2"] > 1e-10:
        raise RuntimeError("Expected exact q(0)=0 of published wave")
    report = {
        "schema": "scrimp-v130-wave-dae-initial-v1",
        "source_image_revision": os.environ.get("SCRIMP_IMAGE_REVISION", ""),
        "mpi_ranks": mpi.Get_size(),
        "mesh_h": args.mesh,
        "dt": args.dt,
        "final_time": args.final_time,
        "cn_steps": args.steps,
        "dofs": len(after),
        "initial_consistency_type": "pseudo",
        "initial_consistency_dt": args.dt ** 2,
        "variable_intervals": {
            k: [v[0], v[0] + v[1]] for k, v in captures["names"].items()
        },
        "dynamic_names": list(captures["dynamic"]),
        "algebraic_names": list(captures["algebraic"]),
        "dynamic_initial_changes": changes,
        "physical_initial_value_fidelity": (
            "PASS" if physical_fidelity else "FAIL"
        ),
        "algebraic_residual_before": captures["before_residual"],
        "algebraic_residual_after": captures["after_residual"],
        "getfem_rank_delta_before": captures["before_rank_delta"],
        "getfem_rank_delta_after": captures["after_rank_delta"],
        "accepted_pseudo_vs_main_start_scaled_l2": start_difference,
        "saved_times": times.tolist(),
        "saved_variable_norms": {
            name: [
                float(np.linalg.norm(trajectory[:, start:start + size][i]))
                for i in range(len(times))
            ]
            for name, (start, size) in captures["names"].items()
        },
        "fidelity_atol": args.fidelity_atol,
        "fidelity_rtol": args.fidelity_rtol,
        "caution": (
            "A small algebraic residual is distinct from physical "
            "initial-value fidelity. No production solver is changed. "
            "The full port-Hamiltonian power balance is a separate gate."
        ),
    }
    if mpi.Get_rank() == 0:
        args.output_dir.mkdir(parents=True, exist_ok=True)
        stem = f"wave_initial_{mpi.Get_size()}r"
        (args.output_dir / f"{stem}.json").write_text(
            json.dumps(report, indent=2) + "\n", encoding="utf-8"
        )
        np.savez_compressed(
            args.output_dir / f"{stem}.npz",
            before=before, after=after,
            times=times, trajectory=trajectory,
        )
        print(
            f"DIAGNOSTIC DAE {mpi.Get_size()} rank(s): "
            f"q(0) drift L2={changes['q']['absolute_change_l2']:.9g}; "
            f"p(0) drift L2={changes['p']['absolute_change_l2']:.9g}; "
            f"algebraic residual before/after "
            f"{captures['before_residual']['algebraic_residual_l2']:.9g}/"
            f"{captures['after_residual']['algebraic_residual_l2']:.9g}; "
            f"dynamic fidelity={report['physical_initial_value_fidelity']}",
            flush=True,
        )
        print(
            "Initial-state physics gate is a DIAGNOSTIC; no automatic "
            "production initial-condition correction is applied.",
            flush=True,
        )
    if args.fail_on_dynamic_drift and not physical_fidelity:
        raise RuntimeError(
            "Pseudo-TS changed explicitly prescribed dynamic initial data"
        )


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--mesh", type=float, default=0.2)
    p.add_argument("--dt", type=float, default=0.002)
    p.add_argument("--steps", type=int, default=5)
    p.add_argument("--final-time", type=float, default=0.01)
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--fidelity-atol", type=float, default=1e-8)
    p.add_argument("--fidelity-rtol", type=float, default=1e-8)
    p.add_argument(
        "--fail-on-dynamic-drift", action="store_true",
        help="Otherwise report expected physical fidelity failure "
             "without hiding numerical DAE evidence",
    )
    args = p.parse_args()
    if not (
        0 < args.mesh <= 1
        and args.dt > 0
        and args.steps >= 2
        and abs(args.steps * args.dt - args.final_time) <= 1e-12
        and args.fidelity_atol >= 0 and args.fidelity_rtol >= 0
    ):
        p.error("Need valid h, dt, steps, final-time and tolerances")
    os.environ.setdefault("SCRIMP_INIT_FACTOR_SOLVER_TYPE", "mumps")
    _audit(args, MPI.COMM_WORLD)


if __name__ == "__main__":
    main()
