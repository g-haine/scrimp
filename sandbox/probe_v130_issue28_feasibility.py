# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Issue #28.1: fixed-differential-data ALGEBRAIC FEASIBILITY audit.

Unlike an initializer, this is read-only evidence: intercept the normal
SCRIMP setup before init_step, hold original h0,p0 exactly fixed, remove
only the known null-traction Y *columns* from the algebraic Jacobian,
and solve the resulting (possibly rectangular) linear least-squares
problem on true PETSc-assembled algebraic rows using distributed LSQR.

The actual original IFunction is then recomputed at the diagnostic
least-squares candidate. Confirm algebraic linearity with A*delta,
assess A^T residual, report individual constraint blocks, and roll
back ALL GetFEM variables. No accepted state and no TS time step.
Nonzero residual with a certified stationary least-squares solution
shows a numerical compatibility obstruction in this discretization.
This does not determine DAE differential index or implement Newton.
"""

import argparse
import json
from pathlib import Path
from unittest.mock import patch

import numpy as np
from mpi4py import MPI
from petsc4py import PETSc
import scrimp as S
from viscous_swe_2d_mpi import build_system
from probe_v130_issue28_structure import (
    AuditComplete, indices, indexset, layout, norm_sub, rowmax,
)

ATOL_FEAS = 1e-10
RTOL_FEAS = 1e-8
LINEARITY_REL = 1e-8
NORMAL_REL = 1e-7


def _block_norms(v, names, iv, lo, hi):
    data = {}
    for name in names:
        ix = indexset(indices((name,), iv), lo, hi, v.getComm())
        try:
            data[name] = norm_sub(v, ix)
        finally:
            ix.destroy()
    return data


def _audit(system, mpi):
    dynamic, algebraic, iv, n = layout(system)
    if "Y" not in algebraic or tuple(dynamic) != ("h", "p"):
        raise ValueError("Unsupported SWE variable layout")
    rows = np.sort(indices(algebraic, iv))
    unknown_names = tuple(x for x in algebraic if x != "Y")
    unknowns = np.sort(indices(unknown_names, iv))
    Y = indices(("Y",), iv)
    fixed = indices(dynamic, iv)
    start = np.asarray(system.gf_model.from_variables(), dtype=float).copy()
    if start.shape != (n,) or not np.all(np.isfinite(start)):
        raise ValueError("Invalid original initial GetFEM vector")
    lo, hi = system.tangent_stiffness.getOwnershipRange()
    if system.tangent_mass.getOwnershipRange() != (lo, hi):
        raise ValueError("Distributed mass/tangent ownership differs")
    created = []

    def keep(obj):
        created.append(obj)
        return obj

    try:
        z0 = keep(system.tangent_stiffness.createVecRight())
        if hi > lo:
            z0.setValues(
                np.arange(lo, hi, dtype=PETSc.IntType),
                start[lo:hi],
                addv=PETSc.InsertMode.INSERT_VALUES,
            )
        z0.assemble()
        zdot = keep(z0.duplicate())
        zdot.set(0)
        F = keep(z0.duplicate())
        ia = keep(indexset(rows, lo, hi, z0.getComm()))
        iu = keep(indexset(unknowns, lo, hi, z0.getComm()))
        baseline_blocks = None

        # The original (zero algebraic initial guess) residual.
        system.IFunction(None, 0., z0, zdot, F)
        initial_norm = norm_sub(F, ia)
        baseline_blocks = _block_norms(F, algebraic, iv, lo, hi)
        mass_rows = rowmax(system.tangent_mass, rows, mpi)
        if mass_rows > 1e-12:
            raise ValueError("Algebraic tangent mass rows are nonzero")

        # z0 is held fixed in h,p. Here the open-loop SWE algebraic
        # equations are linear in the remaining algebraic unknowns,
        # although the overall PDE is nonlinear in h,p.
        system._sync_lumped_trial_state(0., z0)
        system.assemble_nl_stiffness()
        A = keep(system.tangent_stiffness.createSubMatrix(ia, iu))
        m, nfree = A.getSize()
        if (m != len(rows) or nfree != len(unknowns)
                or len(rows) != len(unknowns) + len(Y)):
            raise RuntimeError("Unexpected reduced rectangular Jacobian shape")

        rhs = keep(A.createVecLeft())
        fa = F.getSubVector(ia)
        try:
            rhs.getArray()[:] = -np.asarray(fa.getArray(readonly=True))
        finally:
            F.restoreSubVector(ia, fa)
        rhs.assemble()
        delta = keep(A.createVecRight())
        delta.set(0)

        ksp = keep(PETSc.KSP().create(comm=A.getComm()))
        ksp.setOptionsPrefix("scrimp_issue28_feas_")
        ksp.setOperators(A)
        ksp.setType("lsqr")
        ksp.getPC().setType("none")
        ksp.setTolerances(rtol=1e-12, atol=1e-13, max_it=10000)
        ksp.setFromOptions()
        ksp.solve(rhs, delta)
        reason = int(ksp.getConvergedReason())
        iterations = int(ksp.getIterationNumber())
        correction_norm = float(delta.norm())
        correction_local = np.asarray(delta.getArray(readonly=True)).copy()
        if not mpi.allreduce(bool(np.all(np.isfinite(correction_local))), op=MPI.LAND):
            raise FloatingPointError("LSQR correction is not finite")
        owned_unknowns = unknowns[(unknowns >= lo) & (unknowns < hi)]
        if len(owned_unknowns) != len(correction_local):
            raise RuntimeError("Reduced PETSc correction ownership/index mismatch")

        # Verify that original F_a is affine in these algebraic
        # variables by comparing F_a(z0+delta) to F_a(z0)+A*delta.
        predicted = keep(A.createVecLeft())
        A.mult(delta, predicted)
        predicted.axpy(-1, rhs)  # F_a(z0) == -rhs
        candidate = keep(z0.copy())
        if len(owned_unknowns):
            candidate.setValues(
                np.asarray(owned_unknowns, dtype=PETSc.IntType),
                correction_local,
                addv=PETSc.InsertMode.ADD_VALUES,
            )
        candidate.assemble()
        # No candidate is accepted into the physical state.
        system.IFunction(None, 0., candidate, zdot, F)
        final_norm = norm_sub(F, ia)
        final_blocks = _block_norms(F, algebraic, iv, lo, hi)
        actual = F.getSubVector(ia)
        try:
            discrepancy = keep(predicted.copy())
            discrepancy.axpy(-1, actual)
            affine_error = float(discrepancy.norm())
            gradient = keep(A.createVecRight())
            A.multTranspose(actual, gradient)
            normal_gradient_norm = float(gradient.norm())
        finally:
            F.restoreSubVector(ia, actual)

        # Observe the FULL state after the diagnostic LSQR update,
        # then roll back the original GetFEM values immediately.
        system._sync_lumped_trial_state(0., candidate)
        candidate_full = np.asarray(system.gf_model.from_variables(), dtype=float)
        exact_fixed_during_trial = bool(np.array_equal(
            candidate_full[fixed], start[fixed]
        ))
        exact_Y_unchanged = bool(np.array_equal(candidate_full[Y], start[Y]))
        all_exact = bool(mpi.allreduce(
            exact_fixed_during_trial and exact_Y_unchanged, op=MPI.LAND
        ))
        system._sync_lumped_trial_state(0., z0)
        restored = np.asarray(system.gf_model.from_variables(), dtype=float)
        rollback_exact = bool(mpi.allreduce(
            np.array_equal(restored, start), op=MPI.LAND
        ))
        if not all_exact or not rollback_exact:
            raise RuntimeError("Diagnostic changed fixed h/p or Y, or rollback failed")

        feasibility_tol = ATOL_FEAS + RTOL_FEAS * initial_norm
        affine_tol = LINEARITY_REL * max(1., initial_norm)
        gradient_tol = NORMAL_REL * max(1., initial_norm)
        certifiable = (
            reason > 0
            and affine_error <= affine_tol
            and normal_gradient_norm <= gradient_tol
        )
        if not all(np.isfinite(x) for x in (
            initial_norm, correction_norm, final_norm, affine_error,
            normal_gradient_norm,
        )):
            raise FloatingPointError("Nonfinite least-squares feasibility evidence")
        if affine_error > affine_tol:
            raise RuntimeError(
                f"Algebraic dependence is not affine at fixed h,p: "
                f"linearization error {affine_error:.9g} > {affine_tol:.9g}"
            )
        status = (
            "FEASIBLE_AT_TOLERANCE" if certifiable and final_norm <= feasibility_tol
            else "INCOMPATIBLE_AT_TOLERANCE" if certifiable
            else "INCONCLUSIVE_LSQR"
        )
        # Independent confirmation of the actual GetFEM geometric
        # mapping, NOT just the number of nodes in a Gmsh .msh file.
        # A curved mesh must retain a nonlinear GeoTrans after import.
        gt_raw = system.domain.get_mesh()[0].geotrans()
        geotrans = gt_raw[0] if isinstance(gt_raw, tuple) else gt_raw
        if not isinstance(geotrans, (list, tuple, np.ndarray)):
            geotrans = [geotrans]
        gt_names = [gt.char() for gt in geotrans]
        gt_has_curve = any(not bool(gt.is_linear()) for gt in geotrans)
        report = {
            "schema": "scrimp-v130-issue28-fixed-state-feasibility-v1",
            "status": status,
            "mpi_ranks": int(mpi.size),
            "rotation": float(getattr(system, "_issue28_rotation", 0.)),
            "original_tangent_profile": True,
            "original_polygon_facet_boundary": int(system._issue28_geometry_order) == 1,
            "geometry_order_requested": int(system._issue28_geometry_order),
            "getfem_imported_geometry_nodes": int(
                system.domain.get_mesh()[0].nbpts()
            ),
            "getfem_imported_element_count": int(
                system.domain.get_mesh()[0].nbcvs()
            ),
            "getfem_imported_geotransform_descriptors": gt_names,
            "getfem_import_has_nonlinear_geometric_transform":
                bool(gt_has_curve),
            "compatibility_advice": [
                "The field FEM degree and the mesh GEOMETRIC order are independent.",
                "For rotating data on a circular physical boundary, first test "
                "curved geometry order 2 rather than changing the normal wall control.",
                "An isoparametric P2 boundary still only approximates the exact circle.",
                "If the native weak boundary residual remains inconsistent, "
                "choose an initial momentum field that satisfies the actual "
                "discrete no-penetration and tangential constraints, or drive "
                "the wall smoothly from a compatible rest state.",
                "Never silently change prescribed p0 or boundary controls, "
                "and never use a nonphysical rigid_trace wall as a PASS."
            ],
            "dofs": n,
            "algebraic_rows": len(rows),
            "adjustable_algebraic_dofs": len(unknowns),
            "excluded_null_traction_Y_dofs": len(Y),
            "reduced_matrix_shape": [m, nfree],
            "initial_Fa_l2": initial_norm,
            "least_squares_final_Fa_l2": final_norm,
            "initial_residual_by_block": baseline_blocks,
            "projected_residual_by_block": final_blocks,
            "numerical_feasibility_tol": feasibility_tol,
            "linear_dependence_verified": bool(affine_error <= affine_tol),
            "affine_model_error_l2": affine_error,
            "affine_model_error_limit": affine_tol,
            "normal_equation_AT_Fa_l2": normal_gradient_norm,
            "normal_equation_limit": gradient_tol,
            "distributed_lsqr_converged_reason": reason,
            "distributed_lsqr_iterations": iterations,
            "distributed_lsqr_reported_residual": float(ksp.getResidualNorm()),
            "algebraic_correction_l2": correction_norm,
            "differential_h0_p0_bitwise_preserved": all_exact,
            "traction_Y_unchanged": all_exact,
            "entire_initial_state_rollback_bitwise": rollback_exact,
            "no_accepted_initial_state": True,
            "initialization_TS_executed": False,
            "main_TS_executed": False,
            "nonlinear_initializer_qualified": False,
            "index_classification_proven": False,
            "interpretation": (
                "Independent fixed-dynamic-data native LSQR feasibility "
                "diagnostic; classification only when both KSP convergence "
                "and A^T residual certify the affine least-squares minimum."
            ),
        }
        if mpi.rank == 0:
            print("ISSUE28 FEASIBILITY:", status, flush=True)
            print(json.dumps({
                "mpi_ranks": mpi.size,
                "initial_Fa_l2": initial_norm,
                "least_squares_final_Fa_l2": final_norm,
                "residual_by_block": final_blocks,
                "lsqr_reason": reason,
                "lsqr_iterations": iterations,
                "normal_equation_AT_Fa_l2": normal_gradient_norm,
                "affine_error": affine_error,
                "fixed_h_p_preserved": all_exact,
                "rollback_exact": rollback_exact,
            }, indent=2), flush=True)
        return report
    finally:
        # Restore even after an unexpected failure or a rejected LSQR
        # trial: never leave the global GetFEM model at a candidate.
        system.gf_model.set_time(0.)
        system.gf_model.to_variables(start.copy())
        for obj in reversed(created):
            obj.destroy()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--mesh", type=float, default=.2)
    p.add_argument("--rotation", type=float, choices=None, default=.2)
    p.add_argument("--geometry-order", type=int, choices=(1, 2), default=1)
    p.add_argument("--out", type=Path, required=True)
    args = p.parse_args()
    if not np.isfinite(args.rotation) or args.rotation not in (0., .2):
        p.error("Only original 0.2 rotation and zero-rotation control are accepted")
    if args.mesh != .2:
        p.error("Qualification audit pins exactly mesh=0.2")
    mpi = PETSc.COMM_WORLD.tompi4py()
    system = build_system(
        gamma=0., mesh=args.mesh,
        final_time=.02, save_dt=.002, max_step=.002,
        ramp_start=.003, ramp_end=.009,
        init_step=True, ts_type="bdf", adapt_type="none",
        init_dt=1e-9, rotation=args.rotation,
        initial_profile="tangent",
        wall_velocity_profile="facet_tangent",
        algebraic_initial_guess="zero", init_method="bdf",
        height_order=2, geometry_order=args.geometry_order,
    )
    system._issue28_rotation = args.rotation
    system._issue28_geometry_order = args.geometry_order
    result = {}

    def intercept(self):
        result["audit"] = _audit(self, mpi)
        raise AuditComplete

    with patch.object(S.DPHS, "init_step", intercept):
        try:
            system.solve()
        except AuditComplete:
            pass
        else:
            raise RuntimeError("Unexpected physical time integration")
    if "audit" not in result:
        raise RuntimeError("No feasibility result")
    if mpi.rank == 0:
        args.out.mkdir(parents=True, exist_ok=True)
        (args.out / "feasibility.json").write_text(
            json.dumps(result["audit"], indent=2) + "\n",
            encoding="utf-8",
        )
    mpi.Barrier()


if __name__ == "__main__":
    main()
