# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""#28: independently eliminate the exact interior and wall-control blocks.

Read-only, original P2 rotating shallow-water problem at prescribed q0.
At constant h0, the e_p and p FEM bases coincide: impose e_p=p0/rho
by *copying the actual original FE coefficients*, not by LSQR.
Then solve only the original F_U=0 wall-control block with its native
Jacobian -M_U. Verify native F_ep and F_U vanish. Measure the remaining
native F_Y and compare coefficient by coefficient with +/-h0 times the
independently GetFEM-assembled weak wall moment b_i.

A nonzero F_Y then demonstrates the joint incompatibility of these
specific original algebraic blocks at the prescribed q0, under
nondegenerate same-basis FE mass maps. It does not prove a general
high-index initializer or an index classification. No physical TS,
initialization, change to q0, artificial wall control, or accepted Y0.
"""

import argparse
import json
from pathlib import Path
from unittest.mock import patch

import numpy as np
from mpi4py import MPI
from petsc4py import PETSc
import scrimp as S

from viscous_swe_2d_mpi import build_system, H_INITIAL, RHO, ROTATION
from probe_v130_issue28_structure import (
    AuditComplete, indices, indexset, layout, norm_sub,
)
from probe_v130_issue28_weak_wall import (
    ensure_real, assemble_wall_vector,
)


def global_array(vec, mpi):
    """Gather a small diagnostic vector without truncating complex values."""
    return np.concatenate(mpi.allgather(
        ensure_real(vec.getArray(readonly=True), "PETSc diagnostic vector").copy()
    ))


def audit(s, mpi, mesh, omega):
    differential, algebraic, iv, n = layout(s)
    for name in ("h", "p", "e_p", "U", "Y", "U_n"):
        if name not in iv:
            raise RuntimeError(f"Missing original model variable {name}")
    p_first, p_size = iv["p"]
    ep_first, ep_size = iv["e_p"]
    u_first, u_size = iv["U"]
    y_first, y_size = iv["Y"]
    un_first, un_size = iv["U_n"]
    if p_size != ep_size or u_size != y_size:
        raise RuntimeError("Original p/e_p or U/Y FE sizes differ")
    original = ensure_real(s.gf_model.from_variables(), "original initial state")
    if original.size != n:
        raise RuntimeError("Incorrect native FE state size")
    if np.linalg.norm(original[un_first:un_first+un_size]) != 0.:
        raise RuntimeError("Original zero normal wall control was changed")
    lo, hi = s.tangent_stiffness.getOwnershipRange()
    created = []

    def keep(obj):
        created.append(obj)
        return obj

    snapshot_global = s.GlobalVec.copy()
    snapshot_ready = getattr(s, "_explicit_rhs_ready", False)
    snapshot_time = getattr(s, "_explicit_rhs_time", None)
    rhs = getattr(s, "_explicit_rhs", None)
    rhs_copy = rhs.copy() if rhs is not None else None

    try:
        # Construct the *exact* constitutive FE candidate. Because p and
        # e_p share the original P2 mesh FEM, this is a coefficient identity.
        constitutive = original.copy()
        constitutive[ep_first:ep_first+ep_size] = (
            original[p_first:p_first+p_size] / RHO
        )
        z = keep(s.tangent_stiffness.createVecRight())
        if hi > lo:
            z.setValues(np.arange(lo, hi, dtype=PETSc.IntType),
                        constitutive[lo:hi])
        z.assemble()
        zdot = keep(z.duplicate())
        zdot.set(0.)
        f = keep(z.duplicate())
        iep = keep(indexset(indices(("e_p",), iv), lo, hi, z.getComm()))
        iu = keep(indexset(indices(("U",), iv), lo, hi, z.getComm()))
        iy = keep(indexset(indices(("Y",), iv), lo, hi, z.getComm()))

        s.IFunction(None, 0., z, zdot, f)
        ep_before = norm_sub(f, iep)
        u_before = norm_sub(f, iu)
        if ep_before > 1e-10:
            raise RuntimeError(
                f"Exact p/e_p coefficient copy does not satisfy native F_ep: "
                f"{ep_before:.9g}"
            )

        # For the untouched wall-control brick,
        # F_U(U+delta) = F_U(U) + J_UU delta, J_UU = -M_U.
        s._sync_lumped_trial_state(0., z)
        s.assemble_nl_stiffness()
        M = keep(s.tangent_stiffness.createSubMatrix(iu, iu))
        if M.getSize() != (u_size, u_size):
            raise RuntimeError("Unexpected native wall U mass size")
        M.scale(-1.)

        b = keep(M.createVecLeft())
        sub = f.getSubVector(iu)
        try:
            b.getArray()[:] = ensure_real(
                sub.getArray(readonly=True), "original native F_U"
            )
        finally:
            f.restoreSubVector(iu, sub)
        b.assemble()

        wall = keep(M.createVecRight())
        wall.set(0.)
        ksp = keep(PETSc.KSP().create(comm=M.getComm()))
        ksp.setOptionsPrefix("scrimp_issue28_eliminate_wall_")
        ksp.setOperators(M)
        ksp.setType(PETSc.KSP.Type.PREONLY)
        ksp.getPC().setType(PETSc.PC.Type.LU)
        ksp.getPC().setFactorSolverType("mumps")
        if b.norm() > 0.:
            ksp.solve(b, wall)
            if ksp.getConvergedReason() <= 0:
                raise RuntimeError("Native U wall mass elimination failed")
        wall_coeff = global_array(wall, mpi)
        if wall_coeff.size != u_size:
            raise RuntimeError("Wrong distributed wall coefficient layout")
        trial = constitutive.copy()
        trial[u_first:u_first+u_size] = (
            constitutive[u_first:u_first+u_size] + wall_coeff
        )
        z_trial = keep(z.copy())
        if hi > lo:
            z_trial.setValues(np.arange(lo, hi, dtype=PETSc.IntType),
                              trial[lo:hi])
        z_trial.assemble()
        s.IFunction(None, 0., z_trial, zdot, f)
        full_f = global_array(f, mpi)
        ep_final = norm_sub(f, iep)
        u_final = norm_sub(f, iu)
        y_final = norm_sub(f, iy)
        Fy = full_f[y_first:y_first+y_size]
        if max(ep_final, u_final) > 1e-10:
            raise RuntimeError(
                f"Native exact elimination failed: F_ep={ep_final:.9g}, "
                f"F_U={u_final:.9g}"
            )
        s._sync_lumped_trial_state(0., z_trial)
        expr = f"((p/rho)-({omega:.17g})*Tangent).Test_U"
        moment, metadata = assemble_wall_vector(
            s, mpi, iv, n, "direct-wall-mismatch", expr
        )
        prediction = H_INITIAL * moment
        error_plus = float(np.linalg.norm(Fy-prediction))
        error_minus = float(np.linalg.norm(Fy+prediction))
        sign = 1 if error_plus <= error_minus else -1
        error = min(error_plus, error_minus)
        residual_scale = float(np.linalg.norm(prediction))
        # Different GetFEM assembly paths accumulate O(machine epsilon)
        # absolute coefficient rounding. The former
        # 1e-10*max(1e-6, ||5b||) accidentally rejected 1.0994e-16
        # at h=.1 although F_Y was ~6.38e-7. Compare with both a
        # strict 5e-15 absolute rounding budget and 2e-9 *physical*
        # relative coefficient tolerance; do not scale by a made-up
        # floor on ||5b||. Gross moment/sign errors remain rejected.
        identity_tolerance = 5e-15 + 2e-9 * residual_scale
        identity_relative_error = error / max(residual_scale, 1e-30)
        if not np.isfinite(error) or error > identity_tolerance:
            raise RuntimeError(
                "Exact native F_Y differs from +/-h0*independent weak "
                f"U wall moments: error={error:.9g}, "
                f"relative={identity_relative_error:.9g}, "
                f"limit={identity_tolerance:.9g}, scale={residual_scale:.9g}"
            )
        if omega == 0.:
            if y_final > 1e-10 or residual_scale > 1e-10:
                raise RuntimeError("Exact rest positive control failed")
        elif residual_scale < 1e-10 or y_final < 1e-10:
            raise RuntimeError(
                "Original rotating weak incompatibility unexpectedly vanished"
            )
        if not np.array_equal(trial[indices(differential, iv)],
                              original[indices(differential, iv)]):
            raise RuntimeError("Original physical differential q0 changed")
        if not np.array_equal(trial[y_first:y_first+y_size],
                              original[y_first:y_first+y_size]):
            raise RuntimeError("Original traction Y changed")

        report = {
            "schema": "scrimp-v130-issue28-exact-weak-elimination-v1",
            "status": "EXACT_BLOCK_ELIMINATION_DIAGNOSTIC_ONLY",
            "nominal_mesh": mesh, "mpi_ranks": mpi.size,
            "rotation": omega, "dofs": n,
            "original_geometry_P2": True,
            "original_p_e_p_same_FE_dimension": p_size == ep_size,
            "original_U_Y_same_FE_dimension": u_size == y_size,
            "native_Fep_after_exact_copy_l2": ep_before,
            "native_Fep_after_wall_elimination_l2": ep_final,
            "native_FU_before_wall_elimination_l2": u_before,
            "native_FU_after_wall_elimination_l2": u_final,
            "native_FY_after_both_exact_eliminations_l2": y_final,
            "native_FY_coefficients": Fy.tolist(),
            "independent_native_weak_U_moment_coefficients": moment.tolist(),
            "native_FY_equals_signed_h0_wall_moment_sign": sign,
            "native_FY_weak_moment_identity_error_l2": error,
            "native_FY_weak_moment_identity_relative_error":
                identity_relative_error,
            "native_FY_weak_moment_identity_tolerance_l2":
                identity_tolerance,
            "expected_signed_h0_wall_moment_l2": residual_scale,
            "native_U_mass_KSP_reason": (
                int(ksp.getConvergedReason()) if b.norm() > 0. else 0
            ),
            "original_fixed_h0_p0_preserved_bitwise": True,
            "original_traction_Y_preserved_bitwise": True,
            "entire_state_rollback_bitwise": True,
            "no_physical_TS": True,
            "no_pseudo_or_BDF_initializer": True,
            "consistent_rotating_initial_Y0_accepted": False,
            "differential_index_proven": False,
            "generic_initialization_qualified": False,
            "interpretation": (
                "At fixed original h0,p0, original P2 e_p=p/rho "
                "satisfies the interior constitutive equation. Solving "
                "the unmodified U control block exactly leaves a nonzero "
                "native Y algebraic residual equal to the independently "
                "integrated original weak wall mismatch times h0 (up to "
                "the native flow/effort sign). This isolates these "
                "algebraic blocks; no DAE state has been accepted."
            ),
        }
        if mpi.rank == 0:
            print(
                f"ISSUE28 EXACT WEAK ELIMINATION h={mesh:g} omega={omega:g} "
                f"MPI={mpi.size} F_ep={ep_final:.9g} F_U={u_final:.9g} "
                f"F_Y={y_final:.9g} expected_h_b={residual_scale:.9g} "
                f"identity_error={error:.9g}",
                flush=True,
            )
        return report
    finally:
        s.gf_model.set_time(0.)
        s.gf_model.to_variables(original.copy())
        snapshot_global.copy(s.GlobalVec)
        s._explicit_rhs_ready = snapshot_ready
        s._explicit_rhs_time = snapshot_time
        if rhs_copy is not None and rhs is not None:
            rhs_copy.copy(rhs)
        restored = ensure_real(
            s.gf_model.from_variables(), "GetFEM exact elimination rollback"
        )
        if not mpi.allreduce(bool(np.array_equal(restored, original)),
                             op=MPI.LAND):
            raise RuntimeError("Exact-block diagnostic failed bitwise rollback")
        for obj in reversed(created):
            obj.destroy()
        snapshot_global.destroy()
        if rhs_copy is not None:
            rhs_copy.destroy()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mesh", type=float, choices=(.2, .1), required=True)
    parser.add_argument("--case", choices=("rest", "rotating"), required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    mpi = PETSc.COMM_WORLD.tompi4py()
    omega = 0. if args.case == "rest" else ROTATION
    system = build_system(
        gamma=0., mesh=args.mesh, geometry_order=2, height_order=2,
        final_time=.02, save_dt=.002, max_step=.002,
        ramp_start=.003, ramp_end=.009, init_step=True,
        ts_type="bdf", adapt_type="none", init_dt=1e-9,
        rotation=omega, initial_profile="tangent",
        wall_velocity_profile="facet_tangent",
        algebraic_initial_guess="zero", init_method="bdf",
    )
    result = {}

    def intercept(self):
        result["report"] = audit(self, mpi, args.mesh, omega)
        raise AuditComplete

    with patch.object(S.DPHS, "init_step", intercept):
        try:
            system.solve()
        except AuditComplete:
            pass
        else:
            raise RuntimeError("Unexpected physical time integration")
    if "report" not in result:
        raise RuntimeError("No original algebraic elimination audit ran")
    if mpi.rank == 0:
        args.out.mkdir(parents=True, exist_ok=True)
        (args.out/"exact_weak_elimination.json").write_text(
            json.dumps(result["report"], indent=2)+"\n",
            encoding="utf-8",
        )
    mpi.Barrier()


if __name__ == "__main__":
    main()
