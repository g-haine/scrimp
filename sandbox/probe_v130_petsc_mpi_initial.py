# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Prototype MPI-distributed PETSc algebraic-only DAE initialization.

Diagnostic only; it patches a single freshly built DPHS instance's
init_step. The production rc3 source, default initializer and Docker
image are unchanged, and NO pseudo-time TS runs.

For a *linear* state/time-independent stiffness and mass wave,
extract algebraic rows/columns using distributed PETSc IS, then use
PETSc KSP LSQR on MatCreateSubMatrix(K, I_alg, I_alg) to correct ONLY
algebraic unknowns. It is essential that the full algebraic equation
is re-evaluated with the original IFunction after solving, and that
the original differential state remains bitwise identical.

A nonzero residual is explicitly reported as UNRESOLVED; do not
equate KSP convergence (possibly least squares) with DAE consistency.
The published Gaussian is an expected negative feasibility control;
the manufactured polynomial p0=x*(3-x) is a positive control.
"""

import argparse
import json
import os
from pathlib import Path
import types

from mpi4py import MPI
import numpy as np
from petsc4py import PETSc

from audit_v130_wave_dae_initial import _intervals
from audit_v130_wave_initializers import _build_manufactured_wave
from probe_v130_wave_fixed_dynamics import (
    _algebraic_indices, _build_scenario,
)
from check_acoustic_energy_mpi import _audit as audit_port_power


def _layout(system, mpi):
    """Return complete, disjoint dynamic and algebraic index sets."""
    _, dynamic, algebraic = _intervals(system)
    n = int(system.gf_model.nbdof())
    aidx = _algebraic_indices(algebraic)
    didx = np.concatenate([
        np.arange(start, start+count, dtype=np.int64)
        for start, count in dynamic.values()
    ])
    all_ids = np.concatenate((aidx, didx))
    if not np.array_equal(
        np.sort(all_ids), np.arange(n, dtype=np.int64)
    ):
        raise ValueError(
            "Unsupported layout: differential and algebraic intervals "
            "must partition all GetFEM DOFs exactly"
        )
    if not (aidx.size and didx.size):
        raise ValueError("This diagnostic requires both DAE variable types")
    if not (system.linear_mass and system.linear_stiffness):
        raise ValueError(
            "Linear DAE gate: nonlinear or state-dependent matrices "
            "are unsupported, not silently approximated"
        )
    if system.lumped_variables:
        raise ValueError("Lumped DAE states are out of diagnostic scope")
    if any(
        port.get_substituted() for port in system.ports.values()
    ):
        raise ValueError("Substituted variables unsupported")
    rows_lo, rows_hi = map(
        int, system.stiffness.getOwnershipRange()
    )
    if system.mass.getOwnershipRange() != (rows_lo, rows_hi):
        raise ValueError("Mass and stiffness PETSc row layout mismatch")

    # All ranks hold the same global GetFEM ordering metadata;
    # only the owned algebraic coordinates enter the PETSc IS.
    local = aidx[
        (aidx >= rows_lo) & (aidx < rows_hi)
    ]
    ptr, cols, data = system.mass.getValuesCSR()
    if len(ptr) != rows_hi-rows_lo+1:
        raise ValueError("Unexpected local mass CSR shape")
    bad_mass = max(
        (
            float(np.max(np.abs(data[ptr[i-rows_lo]:ptr[i-rows_lo+1]])))
            for i in local
            if ptr[i-rows_lo+1] > ptr[i-rows_lo]
        ),
        default=0.0,
    )
    global_max_mass = mpi.allreduce(bad_mass, op=MPI.MAX)
    if global_max_mass > 1e-12:
        raise ValueError(
            "Algebraic rows contain time-derivative mass: "
            f"max coefficient {global_max_mass:.9g}"
        )
    return n, dynamic, algebraic, aidx, local, rows_lo, rows_hi, global_max_mass


def _petsc_initial_vec(system, initial, lo, hi):
    z = system.stiffness.createVecRight()
    if hi>lo:
        z.setValues(
            range(lo,hi), initial[lo:hi],
            addv=PETSc.InsertMode.INSERT_VALUES,
        )
    z.assemble()
    return z


def _native_algebraic_residual(system, z, zero_dot, out, isalg):
    """Re-evaluate original production IFunction (including RHS)."""
    system.IFunction(
        None, float(system.time_scheme["t_0"]),
        z, zero_dot, out,
    )
    selected = out.getSubVector(isalg)
    try:
        return float(selected.norm(PETSc.NormType.NORM_2))
    finally:
        out.restoreSubVector(isalg, selected)


def _project_distributed(system, mpi, *, atol=1e-10, rtol=1e-8):
    n, dynamic, algebraic, aidx, local, lo, hi, algebraic_mass_max = _layout(
        system, mpi
    )
    initial = np.array(
        system.gf_model.from_variables(), dtype=np.float64, copy=True
    )
    if initial.shape != (n,):
        raise ValueError("GetFEM initial coefficient layout changed")
    t0 = float(system.time_scheme["t_0"])
    z = _petsc_initial_vec(system, initial, lo, hi)
    z_dot = z.duplicate()
    z_dot.set(0)
    residual = z.duplicate()
    residual.set(0)

    # Each rank contributes only algebraic indices from its PETSc
    # owned rows. MatCreateSubMatrix is collective, never gathered
    # to rank 0; the resulting K_aa remains distributed.
    isalg = PETSc.IS().createGeneral(
        np.asarray(local, dtype=PETSc.IntType),
        comm=system.stiffness.getComm(),
    )
    kaa = None
    ksp = None
    rhs = None
    correction = None
    try:
        system.gf_model.set_time(t0)
        system.gf_model.to_variables(initial.copy())
        norm_initial = _native_algebraic_residual(
            system, z, z_dot, residual, isalg
        )
        kaa = system.stiffness.createSubMatrix(isalg, isalg)
        rhs = kaa.createVecRight()
        correction = rhs.duplicate()
        correction.set(0)
        # The local PETSc subvector and submatrix rows are in the
        # index set's local-order position; reject any mismatch.
        sublo, subhi = rhs.getOwnershipRange()
        if subhi-sublo != len(local):
            raise ValueError(
                "PETSc extracted K_aa ownership differs from local "
                "algebraic index-set selection"
            )
        ksp = PETSc.KSP().create(
            comm=system.stiffness.getComm()
        )
        ksp.setOptionsPrefix("scrimp_dae_trial_")
        ksp.setOperators(kaa)
        ksp.setType(PETSc.KSP.Type.LSQR)
        ksp.getPC().setType(PETSc.PC.Type.NONE)
        # Do not relax the independently specified native
        # algebraic IFunction acceptance tolerance.
        ksp.setTolerances(
            rtol=1e-13, atol=1e-13,
            max_it=max(1000, 5*int(aidx.size)),
        )
        ksp.setFromOptions()

        tolerance = float(atol + rtol*norm_initial)
        # SCIENTIFIC distinction: feasibility of the algebraic state
        # is governed ONLY by the unchanged 'tolerance'.  The full-CN
        # MPI trajectory comparison is more sensitive than this
        # residual tolerance.  In particular, MPI=2 previously stopped
        # after one pass at 1.7163317e-9 while MPI=1 refined to
        # 9.8141843e-14.  Continue improving a *compatible* positive
        # state toward this smaller, diagnostic-only fidelity target
        # even after the physical feasibility gate has passed.
        # Never weaken either the feasibility or full-CN MPI criteria.
        diagnostic_accuracy_target = (
            min(tolerance, 1.e-11)
            if getattr(system, "_v130_fe_polynomial_diagnostic", False)
            else tolerance
        )
        # PETSc LSQR may stop when its normal-equation/iterative
        # criterion has converged, even when the original weak
        # constraint rows are just above our stricter physical gate.
        # Solve at most two additional distributed defect-correction
        # problems Kaa * d = -F_a(t0, z, 0), re-evaluating the real
        # native IFunction after EACH pass. This is NOT a nonlinear
        # solver or a relaxation of the compatibility requirement.
        max_refinement_passes = 3
        history = []
        best_norm = norm_initial
        best_candidate = initial.copy()
        previous_norm = norm_initial
        reason = 0
        iters = 0
        ksp_final_norm = float("nan")
        # Each rank inserts ONLY its original owned algebraic DOFs.
        # GetFEM/NumPy indices are int64 while the qualified rc3
        # PETSc build uses 32-bit PetscInt.
        if len(local):
            if (
                np.min(local) < np.iinfo(PETSc.IntType).min
                or np.max(local) > np.iinfo(PETSc.IntType).max
            ):
                raise OverflowError("DAE index outside PETSc.IntType range")
        owned_petsc_indices = np.asarray(
            local, dtype=PETSc.IntType
        )

        for pass_number in range(1, max_refinement_passes+1):
            # PETSc submatrix row ownership must match the original
            # local algebraic ordering selected by the distributed IS.
            owned_f = np.asarray(residual.getArray(readonly=True))
            rhs.getArray()[:] = -owned_f[local-lo]
            rhs.assemble()
            correction.set(0)
            ksp.solve(rhs, correction)
            reason = int(ksp.getConvergedReason())
            iters = int(ksp.getIterationNumber())
            ksp_final_norm = float(ksp.getResidualNorm())
            if not np.all(np.isfinite(correction.getArray(readonly=True))):
                raise RuntimeError("Non-finite PETSc algebraic correction")
            if len(local):
                z.setValues(
                    owned_petsc_indices,
                    correction.getArray(readonly=True).copy(),
                    addv=PETSc.InsertMode.ADD_VALUES,
                )
            z.assemble()
            system.Scatter[0].scatter(z, system.GlobalVec)
            candidate = np.array(
                system.GlobalVec.getArray(readonly=True).real,
                dtype=np.float64, copy=True,
            )
            for name, (start, count) in dynamic.items():
                if not np.array_equal(
                    initial[start:start+count],
                    candidate[start:start+count],
                ):
                    raise RuntimeError(
                        f"Distributed projected DAE initializer changed {name}"
                    )
            system.gf_model.set_time(t0)
            system.gf_model.to_variables(candidate.copy())
            norm_after = _native_algebraic_residual(
                system, z, z_dot, residual, isalg
            )
            history.append({
                "pass_number": pass_number,
                "native_algebraic_residual_l2": norm_after,
                "petsc_ksp_convergence_reason": reason,
                "petsc_ksp_iterations": iters,
                "petsc_ksp_reported_residual_norm": ksp_final_norm,
                "petsc_correction_l2": float(
                    correction.norm(PETSc.NormType.NORM_2)
                ),
            })
            if not np.isfinite(norm_after):
                raise RuntimeError(
                    "Non-finite actual DAE algebraic IFunction residual"
                )
            if norm_after < best_norm:
                best_norm = norm_after
                best_candidate = candidate.copy()
            # The physical feasibility threshold is unchanged.
            # Do not truncate a positive-model refinement merely
            # because it passed the looser DAE feasibility gate:
            # the MPI full-state comparison remains independently
            # strict and needs a reproducibly accurate initial state.
            if norm_after <= diagnostic_accuracy_target:
                break
            if reason <= 0:
                # An iterative KSP failure is not proof that the
                # prescribed initial data are incompatible.
                break
            if norm_after >= previous_norm * (1. - 1.e-6):
                # Abort stagnating/oscillating refinement rather than
                # manufacturing a PASS by unbounded iterations.
                break
            previous_norm = norm_after

        # Only publish the best actually evaluated native residual.
        # In particular, a final unsuccessful refinement must never
        # replace a more accurate preceding projected state.
        if not np.array_equal(candidate, best_candidate):
            z.getArray()[:] = best_candidate[lo:hi]
            z.assemble()
            system.Scatter[0].scatter(z, system.GlobalVec)
            candidate = np.array(
                system.GlobalVec.getArray(readonly=True).real,
                dtype=np.float64, copy=True,
            )
            system.gf_model.set_time(t0)
            system.gf_model.to_variables(candidate.copy())
        norm_final = _native_algebraic_residual(
            system, z, z_dot, residual, isalg
        )
        if abs(norm_final - best_norm) > max(1e-13, 1e-6*best_norm):
            raise RuntimeError(
                "Restored PETSc state does not reproduce the best "
                "original IFunction residual"
            )
        # Check consistent across ranks, not just the printed root.
        verdict = {
            "global_dofs": n,
            "algebraic_dofs": int(len(aidx)),
            "algebraic_mass_row_max_abs": float(algebraic_mass_max),
            "native_algebraic_before_l2": norm_initial,
            "native_algebraic_after_l2": norm_final,
            "native_residual_tolerance": tolerance,
            "diagnostic_accuracy_target": diagnostic_accuracy_target,
            "diagnostic_accuracy_target_met": bool(
                norm_final <= diagnostic_accuracy_target
            ),
            "ksp_type": str(ksp.getType()),
            "ksp_pc_type": str(ksp.getPC().getType()),
            "ksp_iterations": iters,
            "ksp_convergence_reason": reason,
            "ksp_reported_residual_norm": ksp_final_norm,
            "lsqr_refinement_passes": history,
            "lsqr_refinement_pass_count": len(history),
            "lsqr_max_refinement_passes": max_refinement_passes,
            "lsqr_success_requires_native_ifunction": True,
            "differential_state_exactly_preserved": True,
            "algebraic_state_correction_l2": float(
                np.linalg.norm(candidate[aidx] - initial[aidx])
            ),
            "algebraic_feasibility": (
                "PASS" if norm_final <= tolerance else "UNRESOLVED"
            ),
        }
        if reason <= 0 and norm_final > tolerance:
            verdict["interpretation"] = (
                "PETSc LSQR did not meet its iterative stopping rule; "
                "this is NOT a proof of mathematical incompatibility"
            )
        elif norm_final > tolerance:
            verdict["interpretation"] = (
                "PETSc LSQR stopped but true DAE algebraic residual "
                "exceeds tolerance; no valid projected state exists "
                "among the computed corrections"
            )
        else:
            verdict["interpretation"] = (
                "Actual IFunction algebraic rows meet tolerance "
                "with prescribed q/p exactly frozen on this mesh"
            )
        return initial, candidate, dynamic, verdict
    finally:
        if ksp is not None:
            ksp.destroy()
        if correction is not None:
            correction.destroy()
        if rhs is not None:
            rhs.destroy()
        if kaa is not None:
            kaa.destroy()
        isalg.destroy()
        residual.destroy()
        z_dot.destroy()
        z.destroy()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--scenario", required=True,
        choices=("polynomial_compatible", "published"),
    )
    p.add_argument("--mesh", type=float, default=0.2)
    p.add_argument("--dt", type=float, default=0.002)
    p.add_argument("--steps", type=int, default=10)
    p.add_argument("--output-dir", type=Path, required=True)
    args = p.parse_args()
    if args.mesh <= 0 or args.dt <= 0 or args.steps < 2:
        p.error("Positive mesh/dt and at least two steps are required")
    mpi = MPI.COMM_WORLD
    os.environ.setdefault("SCRIMP_INIT_FACTOR_SOLVER_TYPE", "mumps")
    if args.scenario == "polynomial_compatible":
        system = _build_manufactured_wave(
            args.mesh, "0.08*sin(20*t)"
        )
    else:
        system = _build_scenario(args.mesh, "published")
    # This marker is confined to this freshly built experimental
    # DPHS instance; never change the published example or production
    # solver.  The tighter iteration goal is a diagnostic fidelity
    # goal, NOT an acceptance-tolerance relaxation.
    system._v130_fe_polynomial_diagnostic = (
        args.scenario == "polynomial_compatible"
    )
    system.set_time_scheme(
        ts_type="cn", ts_adapt_type="none",
        ts_exact_final_time="matchstep",
        dt=args.dt, dt_save=args.dt,
        t_f=args.dt*args.steps,
        ksp_type="preonly", pc_type="lu",
        pc_factor_mat_solver_type="mumps",
    )
    outcome = {}

    def distributed_init(self):
        initial, candidate, dynamic, verdict = _project_distributed(
            self, mpi
        )
        outcome.update(
            initial=initial, candidate=candidate,
            dynamic=dynamic, verdict=verdict,
        )
        if verdict["algebraic_feasibility"] != "PASS":
            # Reject atomically: restore the untouched prescribed
            # GetFEM variables before aborting the physical TS.
            # In particular, never leak an inconsistent candidate
            # into a later solve after caller-side exception handling.
            self.gf_model.set_time(float(self.time_scheme["t_0"]))
            self.gf_model.to_variables(initial.copy())
            # Do NOT give physical TS a silently inconsistent start;
            # use a distinct exception/status to distinguish this
            # expected Gaussian feasibility failure from a KSP crash.
            raise RuntimeError(
                "DAE_ALGEBRAIC_INITIAL_DATA_UNRESOLVED "
                f"residual={verdict['native_algebraic_after_l2']:.9g} "
                f"threshold={verdict['native_residual_tolerance']:.9g}"
            )

    system.init_step = types.MethodType(distributed_init, system)
    report = {
        "schema": "scrimp-v130-distributed-petsc-dae-probe-v1",
        "image_source_revision": os.environ.get(
            "SCRIMP_IMAGE_REVISION", ""
        ),
        "mpi_ranks": mpi.Get_size(),
        "scenario": args.scenario,
        "mesh_h": args.mesh,
        "dt": args.dt,
        "steps": args.steps,
        "pseudo_ts_executed": False,
        "getfem_stiffness_gathered_to_rank_0": False,
    }
    try:
        system.solve()
    except RuntimeError as exc:
        if (
            "DAE_ALGEBRAIC_INITIAL_DATA_UNRESOLVED" not in str(exc)
            or not outcome
        ):
            raise
        # Always preserve the complete PETSc iteration history and
        # actual IFunction results on both a positive-control FAILURE
        # and an expected negative-control REJECTION.  This avoids
        # losing the decisive diagnostics when the physical TS is
        # correctly prevented from starting.
        report.update(outcome["verdict"])
        report["expected_incompatible_gaussian_rejected"] = (
            args.scenario == "published"
        )
        report["physical_cn_executed"] = False
        report["scientific_status"] = (
            "EXPECTED_REJECTION"
            if args.scenario == "published" else "UNRESOLVED"
        )
    else:
        if args.scenario != "polynomial_compatible":
            raise RuntimeError(
                "Previously obstructed Gaussian surprisingly passed "
                "distributed PETSc projection: investigate, do not "
                "treat as predefined expected success"
            )
        trajectory = np.asarray(
            system.solution["z"], dtype=np.float64
        )
        times = np.asarray(system.solution["t"], dtype=float)
        expected_times = np.arange(args.steps+1)*args.dt
        if (
            times.shape != expected_times.shape
            or not np.allclose(
                times, expected_times, atol=1e-10, rtol=0
            )
            or trajectory.shape != (
                len(expected_times), len(outcome["candidate"])
            )
        ):
            raise RuntimeError(
                "Distributed PETSc DAE trial: missing CN snapshots"
            )
        start_gap = float(np.linalg.norm(
            trajectory[0]-outcome["candidate"]
        ))
        if start_gap > 5e-10:
            raise RuntimeError(
                "Main physical CN TS did not start from projected "
                f"PETSc DAE state: {start_gap}"
            )
        if (
            system.last_ts_converged_reason <= 0
            or system.last_ts_steps != args.steps
            or system.last_ts_rejections
            or system.last_ts_snes_failures
        ):
            raise RuntimeError("Projected native CN TS did not converge")
        power = audit_port_power(system, mpi)
        tolerance = float(
            1e-9 + 1e-5*power["balance_reference_scale"]
        )
        if power["absolute_port_work"] <= 1e-9:
            raise RuntimeError(
                "Distributed PETSc projection CN has vacuous port work"
            )
        report.update(outcome["verdict"])
        report.update({
            "expected_incompatible_gaussian_rejected": False,
            "physical_cn_executed": True,
            "scientific_status": (
                "PASS" if (
                    power["maximum_absolute_balance_defect"] <= tolerance
                ) else "POWER_BALANCE_FAIL"
            ),
            "main_initial_vs_projected_l2": start_gap,
            "power_balance_tolerance": tolerance,
            "power": power,
            "physical_times": times.tolist(),
            "main_ts_converged_reason": int(
                system.last_ts_converged_reason
            ),
        })
        if mpi.Get_rank() == 0:
            args.output_dir.mkdir(parents=True, exist_ok=True)
            np.savez_compressed(
                args.output_dir / f"distributed_{args.scenario}_{mpi.Get_size()}r.npz",
                times=times, states=trajectory,
                initial=outcome["initial"],
                projected=outcome["candidate"],
            )
    if mpi.Get_rank() == 0:
        args.output_dir.mkdir(parents=True, exist_ok=True)
        dest = args.output_dir / (
            f"distributed_{args.scenario}_{mpi.Get_size()}r.json"
        )
        dest.write_text(
            json.dumps(report, indent=2)+"\n",
            encoding="utf-8",
        )
        print(
            f"DISTRIBUTED PETSC DAE MPI={mpi.Get_size()} "
            f"scenario={args.scenario}: "
            f"result={report['scientific_status']}, "
            f"F_a={report['native_algebraic_before_l2']:.8g} -> "
            f"{report['native_algebraic_after_l2']:.8g}, "
            f"tol={report['native_residual_tolerance']:.8g}, "
            f"LSQR passes={report['lsqr_refinement_pass_count']}, "
            f"accuracy_goal={report['diagnostic_accuracy_target']:.8g}, "
            f"accuracy_goal_met={report['diagnostic_accuracy_target_met']}, "
            f"LSQR iterations={report['ksp_iterations']}, "
            f"KSP reason={report['ksp_convergence_reason']}, "
            f"pseudo=NEVER, global Kaa gathering=NEVER",
            flush=True,
        )
    if report["scientific_status"] not in (
        "PASS", "EXPECTED_REJECTION"
    ):
        raise RuntimeError(
            "PETSc distributed DAE positive-control initial "
            "feasibility UNRESOLVED: diagnostic JSON archived "
            "before main CN; acceptance threshold unchanged"
        )


if __name__ == "__main__":
    main()
