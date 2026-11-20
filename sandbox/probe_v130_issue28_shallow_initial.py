# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Issue #28: original rectangular nonlinear shallow-water DAE initial audit.

Loads the shipped examples/shallow_water.py unchanged, intercepting only
its long solve/output as in the qualified issue-27 mass-Jacobian probe.
Preserves all published weak bricks, control laws, fields and initial
h0,p0 (the diagnostic mesh parameter defaults to original 0.1).

At t0 it freezes the existing explicit=True gyroscopic weak brick at
the accepted INITIAL state. The original model has nonlinear mass
and constitutive dependence on the differential states, but its
algebraic rows are AFFINE in algebraic unknowns at fixed h,p.
This test attempts exactly ONE algebraic Newton/direct Kaa solve;
checks native IFunction, analytic GetFEM tangent versus independent
FD in algebraic directions, zero algebraic mass rows, preserved h0,p0,
exact GetFEM rollback, and MPI 1/2 parity.

This is a diagnostic, NOT a production nonlinear initializer,
not a general nonlinear-algebraic Newton, and not a physical TS run.
"""

import argparse
import json
from pathlib import Path

import numpy as np
from mpi4py import MPI
from petsc4py import PETSc

from check_v130_shallow_water_mass_jacobian import _build_original
TOL_RESIDUAL = 1e-9
TOL_DERIVATIVE = 1e-7


def indices(names, iv):
    """Global GetFEM variable DOF indices, independent of disk probes."""
    return np.concatenate([
        np.arange(iv[name][0], iv[name][0] + iv[name][1], dtype=np.int64)
        for name in names
    ])


def indexset(global_rows, lo, hi, comm):
    owned = global_rows[
        (global_rows >= lo) & (global_rows < hi)
    ]
    return PETSc.IS().createGeneral(
        np.asarray(owned, dtype=PETSc.IntType), comm=comm
    )


def rowmax(mat, selected_rows, mpi):
    pointer, _, values = mat.getValuesCSR()
    lo, hi = mat.getOwnershipRange()
    mx = 0.
    for row in selected_rows[
        (selected_rows >= lo) & (selected_rows < hi)
    ]:
        segment = values[pointer[row-lo]:pointer[row-lo+1]]
        if len(segment):
            mx = max(mx, float(np.max(np.abs(segment))))
    return float(mpi.allreduce(mx, op=MPI.MAX))


def norm_sub(v, iset):
    subv = v.getSubVector(iset)
    try:
        return float(subv.norm())
    finally:
        v.restoreSubVector(iset, subv)


def layout(system):
    dynamic = tuple(s.get_name() for s in system.states.values())
    if dynamic != ("h", "p"):
        raise ValueError("Unexpected original shallow-water states")
    names = list(dynamic)
    for port in system.ports.values():
        if port.get_substituted():
            raise ValueError("Substituted port unexpected")
        if port.get_algebraic():
            names.extend(x for x in (port.get_flow(), port.get_effort()) if x)
        else:
            names.append(port.get_effort())
    if len(set(names)) != len(names):
        raise ValueError("Duplicate original GetFEM variable")
    n = int(system.gf_model.nbdof())
    iv = {
        key: tuple(map(int, system.gf_model.interval_of_variable(key)))
        for key in names
    }
    covered = np.zeros(n, dtype=np.int8)
    for start, width in iv.values():
        if width <= 0 or start < 0 or start + width > n:
            raise ValueError("Variable interval out of bounds")
        covered[start:start+width] += 1
    if not np.all(covered == 1):
        raise ValueError("GetFEM variable intervals do not partition the state")
    return dynamic, tuple(x for x in names if x not in dynamic), iv, n


def blocks(F, names, iv, lo, hi):
    result = {}
    for name in names:
        iset = indexset(indices((name,), iv), lo, hi, F.getComm())
        try:
            result[name] = norm_sub(F, iset)
        finally:
            iset.destroy()
    return result


def run_audit(system, mpi):
    dnames, anames, iv, n = layout(system)
    if system.linear_mass or system.linear_stiffness:
        raise ValueError("The original shallow-water DAE must have nonlinear mass and stiffness")
    if any(name.startswith("Y") for name in dnames):
        raise ValueError("Boundary multipliers incorrectly classified as differential")
    if not system.bricks["G"].get_explicit():
        raise RuntimeError("Published shallow-water gyro G is not explicit")
    if system.bricks["G"].get_form() != "(Gyro(p) * e_p) . Test_p":
        raise RuntimeError("Published gyroscopic form has changed")

    d = np.sort(indices(dnames, iv))
    a = np.sort(indices(anames, iv))
    original = np.asarray(system.gf_model.from_variables(), dtype=float).copy()
    if original.shape != (n,) or not np.isfinite(original).all():
        raise ValueError("Invalid original initial FE vector")
    lo, hi = system.tangent_stiffness.getOwnershipRange()
    objects = []
    def keep(o):
        objects.append(o)
        return o

    try:
        z = keep(system.J.createVecRight())
        owned = np.arange(lo, hi, dtype=PETSc.IntType)
        if len(owned):
            z.setValues(owned, original[lo:hi])
        z.assemble()
        zd = keep(z.duplicate())
        zd.set(0.)
        F0 = keep(z.duplicate())
        Fa = keep(z.duplicate())
        ia = keep(indexset(a, lo, hi, z.getComm()))
        idyn = keep(indexset(d, lo, hi, z.getComm()))

        # The original example's gyroscopic term is fully explicit.
        # It must remain frozen at the same accepted INITIAL state
        # across the algebraic trial and both finite differences.
        system._freeze_explicit_rhs(0., z)
        if not system._explicit_rhs_ready:
            raise RuntimeError("Original explicit G was not frozen")
        gyro_norm = float(system._explicit_rhs.norm())
        frozen = keep(system._explicit_rhs.copy())

        system.IFunction(None, 0., z, zd, F0)
        initial = norm_sub(F0, ia)
        initial_blocks = blocks(F0, anames, iv, lo, hi)
        mass_a_max = rowmax(system.tangent_mass, a, mpi)
        if mass_a_max > 1e-12:
            raise RuntimeError("Nonzero algebraic mass rows")

        system._sync_lumped_trial_state(0., z)
        system.assemble_nl_stiffness()
        A = keep(system.tangent_stiffness.createSubMatrix(ia, ia))
        if A.getSize() != (len(a), len(a)):
            raise RuntimeError("Non-square shallow-water algebraic block")
        rhs = keep(A.createVecLeft())
        original_fa = F0.getSubVector(ia)
        try:
            rhs.getArray()[:] = -np.asarray(
                original_fa.getArray(readonly=True)
            )
        finally:
            F0.restoreSubVector(ia, original_fa)
        rhs.assemble()
        delta = keep(A.createVecRight())
        delta.set(0.)

        ksp = keep(PETSc.KSP().create(comm=A.getComm()))
        ksp.setOperators(A)
        ksp.setType("preonly")
        ksp.getPC().setType("lu")
        ksp.getPC().setFactorSolverType("mumps")
        ksp.setUp()
        ksp.solve(rhs, delta)
        reason = int(ksp.getConvergedReason())
        if reason <= 0 or not np.isfinite(delta.norm()):
            raise RuntimeError(
                f"Native algebraic direct KSP failed, reason={reason}"
            )

        candidate = keep(z.copy())
        owned_algebraic = a[(a >= lo) & (a < hi)]
        values = np.asarray(delta.getArray(readonly=True)).copy()
        if len(owned_algebraic) != len(values):
            raise RuntimeError("Nonmatching PETSc algebraic ownership")
        if len(values):
            candidate.setValues(
                np.asarray(owned_algebraic, dtype=PETSc.IntType),
                values,
                addv=PETSc.InsertMode.ADD_VALUES,
            )
        candidate.assemble()
        system.IFunction(None, 0., candidate, zd, Fa)
        final = norm_sub(Fa, ia)
        final_blocks = blocks(Fa, anames, iv, lo, hi)

        # The residual is affine in algebraic variables at fixed h,p.
        predicted = keep(A.createVecLeft())
        A.mult(delta, predicted)
        predicted.axpy(-1., rhs)
        full_fa = Fa.getSubVector(ia)
        try:
            predicted.axpy(-1., full_fa)
            affine_error = float(predicted.norm())
        finally:
            Fa.restoreSubVector(ia, full_fa)

        # Independent finite difference, algebraic DOFs only.
        direction = keep(A.createVecRight())
        local_a = np.arange(*direction.getOwnershipRange())
        direction.getArray()[:] = np.sin(1. + 0.017 * local_a)
        nrm = float(direction.norm())
        if nrm == 0:
            raise RuntimeError("Empty algebraic test direction")
        direction.scale(1. / nrm)
        global_dir = keep(z.duplicate())
        global_dir.set(0.)
        owned_a = a[(a >= lo) & (a < hi)]
        if len(owned_a):
            global_dir.setValues(
                np.asarray(owned_a, dtype=PETSc.IntType),
                np.asarray(direction.getArray(readonly=True)),
            )
        global_dir.assemble()
        # Numerical precision: e_h is O(rho*g*h) = O(1e4) in the
        # original SI-scale model. A normalized algebraic direction with
        # a 1e-6 finite-difference step creates perturbations of order
        # 1e-8 per DOF, so subtracting two assembled residuals can
        # amplify floating-point cancellation. For an AFFINE Fa, finite
        # differences should be stable at several larger step sizes.
        # Do NOT relax the derivative threshold on the basis of a
        # single bad small-epsilon experiment.
        steps = (1e-1, 1e-2, 1e-3, 1e-4, 1e-6)
        zp = keep(candidate.copy())
        zm = keep(candidate.copy())
        Fp = keep(z.duplicate())
        Fm = keep(z.duplicate())
        exact_dir = keep(A.createVecLeft())
        A.mult(direction, exact_dir)
        directional_reference_norm = float(exact_dir.norm())
        fd_scan = []
        for eps in steps:
            candidate.copy(zp)
            candidate.copy(zm)
            zp.axpy(eps, global_dir)
            zm.axpy(-eps, global_dir)
            system.IFunction(None, 0., zp, zd, Fp)
            system.IFunction(None, 0., zm, zd, Fm)
            Fp.axpy(-1., Fm)
            Fp.scale(1./(2.*eps))
            fd = Fp.getSubVector(ia)
            try:
                error = keep(exact_dir.copy())
                error.axpy(-1., fd)
                absolute_error = float(error.norm())
                fd_norm = float(fd.norm())
                error.destroy()
                objects.pop()
            finally:
                Fp.restoreSubVector(ia, fd)
            relative_error = absolute_error / max(
                1e-12, directional_reference_norm, fd_norm
            )
            fd_scan.append({
                "epsilon": float(eps),
                "fd_norm": fd_norm,
                "jacobian_action_norm": directional_reference_norm,
                "absolute_error": absolute_error,
                "relative_error": relative_error,
            })
        # A native FD check is qualified only when at least two
        # independently chosen amplitudes achieve the SAME strict
        # relative tolerance, with the native analytic affine identity
        # checked separately. The smallest step is diagnostic-only.
        stable_large_steps = (
            fd_scan[0]["relative_error"] <= TOL_DERIVATIVE
            and fd_scan[1]["relative_error"] <= TOL_DERIVATIVE
        )
        best_fd = min(fd_scan, key=lambda row: row["relative_error"])
        relative_fd = best_fd["relative_error"]
        print(
            f"ISSUE28 SHALLOW FD SCAN MPI={mpi.size}: "
            f"initial |Fa|={initial:.9g}, candidate |Fa|={final:.9g}, "
            f"affine={affine_error:.9g}, "
            f"steps={[(q['epsilon'], q['relative_error']) for q in fd_scan]}",
            flush=True,
        )

        # Trial h,p remain bitwise original; return ALL GetFEM
        # coefficients to their exact original values on all ranks.
        system._sync_lumped_trial_state(0., candidate)
        model_candidate = np.asarray(system.gf_model.from_variables(), dtype=float)
        h_p_preserved = bool(np.array_equal(model_candidate[d], original[d]))
        gyro_drift = keep(frozen.copy())
        gyro_drift.axpy(-1., system._explicit_rhs)
        gyro_drift_norm = float(gyro_drift.norm())
        system.gf_model.set_time(0.)
        system.gf_model.to_variables(original.copy())
        state_restored = bool(np.array_equal(
            np.asarray(system.gf_model.from_variables(), dtype=float),
            original,
        ))
        h_p_preserved = bool(mpi.allreduce(h_p_preserved, op=MPI.LAND))
        state_restored = bool(mpi.allreduce(state_restored, op=MPI.LAND))
        if not h_p_preserved or not state_restored or gyro_drift_norm != 0.:
            raise RuntimeError("State or frozen explicit gyro changed in diagnostic")
        # Report a FAILED numerical qualification rather than throwing
        # before writing the evidence JSON; native runner still fails.
        affine_ok = affine_error <= 1e-8 * max(1., initial)
        fd_ok = stable_large_steps
        residual_ok = final <= TOL_RESIDUAL
        status = (
            "PASS_DIAGNOSTIC" if affine_ok and fd_ok and residual_ok
            else "FAILED_AFFINITY" if not affine_ok
            else "FAILED_JACOBIAN_FD" if not fd_ok
            else "FAILED_ALGEBRAIC_RESIDUAL"
        )
        report = {
            "schema": "scrimp-v130-issue28-shallow-fixed-dynamics-v1",
            "status": status,
            "model": "unchanged examples/shallow_water.py weak forms, original h0,p0",
            "mpi_ranks": int(mpi.size),
            "dofs": int(n),
            "dynamic_dofs": int(len(d)),
            "algebraic_dofs": int(len(a)),
            "algebraic_names": list(anames),
            "nonlinear_mass": True,
            "nonlinear_stiffness": True,
            "explicit_gyroscopic_G_frozen_at_original_initial_state": True,
            "frozen_gyro_rhs_norm": gyro_norm,
            "frozen_gyro_rhs_drift": gyro_drift_norm,
            "max_abs_algebraic_mass_rows": mass_a_max,
            "initial_Fa_l2": initial,
            "initial_Fa_block_l2": initial_blocks,
            "after_one_direct_algebraic_Newton_Fa_l2": final,
            "after_Fa_block_l2": final_blocks,
            "one_step_affine_model_error_l2": affine_error,
            "relative_algebraic_Jacobian_FD_error": relative_fd,
            "finite_difference_scan": fd_scan,
            "fd_large_steps_both_pass": bool(stable_large_steps),
            "algebraic_residual_threshold": TOL_RESIDUAL,
            "algebraic_fd_relative_threshold": TOL_DERIVATIVE,
            "algebraic_direct_KSP_reason": reason,
            "h0_p0_bitwise_preserved": h_p_preserved,
            "getfem_initial_state_rollback_bitwise": state_restored,
            "main_TS_executed": False,
            "initializer_TS_executed": False,
            "genuine_nonlinear_algebraic_Newton_qualified": False,
            "production_nonlinear_initializer_qualified": False,
            "interpretation": (
                "Native original nonlinear shallow-water DAE with "
                "state-dependent mass has affine algebraic constraints "
                "at fixed dynamic h,p; one algebraic direct solve is "
                "diagnostic only and not full time integration."
            ),
        }
        if mpi.rank == 0:
            print(
                f"ISSUE28 SHALLOW FIXED DYNAMICS {status} MPI={mpi.size}: "
                f"|Fa(initial)|={initial:.9g}, |Fa(after)|={final:.9g}, "
                f"FDrel={relative_fd:.9g}, affine={affine_error:.9g}, "
                f"gyro_frozen={gyro_norm:.9g}",
                flush=True,
            )
        return report
    finally:
        system.gf_model.set_time(0.)
        system.gf_model.to_variables(original.copy())
        for obj in reversed(objects):
            obj.destroy()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--mesh", type=float, default=.1)
    p.add_argument("--out", type=Path, required=True)
    args = p.parse_args()
    if args.mesh <= 0:
        p.error("mesh must be positive")
    mpi = PETSc.COMM_WORLD.tompi4py()
    system = _build_original(args.mesh)
    result = run_audit(system, mpi)
    if mpi.rank == 0:
        args.out.mkdir(parents=True, exist_ok=True)
        (args.out / "shallow_fixed_dynamics.json").write_text(
            json.dumps(result, indent=2) + "\n", encoding="utf-8"
        )
    mpi.Barrier()
    if result["status"] != "PASS_DIAGNOSTIC":
        raise RuntimeError(
            "Native shallow-water algebraic diagnostic rejected: "
            f"{result['status']}; full JSON evidence saved. "
            f"Final |Fa|={result['after_one_direct_algebraic_Newton_Fa_l2']:.9g}; "
            f"FD scan={result['finite_difference_scan']}"
        )


if __name__ == "__main__":
    main()
