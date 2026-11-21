# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""#28/#6 compatible physical rest-to-spinup first differentiated pilot.

Use original nonlinear dissipative Saint-Venant curved-P2 physical bricks.
Explicitly new initial data: h0=5, p0=0, physical facet-tangent wall,
omega=.2 applied via an opt-in C1 rising ramp of duration .02.
e_h=.5, Y_n=-.5, all other initial algebraic fields zero.

At t=0 g'(0)=0 EXACTLY, so partial_t F_a=0. Independently verify
time-dependent native GetFEM U residual at delta and 2delta follows
exact smoothstep g(2delta)/g(delta). Never use a spurious fixed central
difference at the clipped start. Then solve original native MatNest
first-differentiated system with true physical -F_d and exact -F_a,t.
Verify qdot0, diagnostic Y0 and original IFunction, MPI1/2, and rollback.

Diagnostic ONLY, no production initializer, no physical TS.
"""

import argparse"""

import argparse
import json
from pathlib import Path
from unittest.mock import patch

import numpy as np
from mpi4py import MPI
from petsc4py import PETSc
import scrimp as S

from viscous_swe_2d_mpi import (
    build_system, H_INITIAL, RHO, GRAVITY, ROTATION,
)
from probe_v130_issue28_structure import (
    AuditComplete, indices, indexset, layout, matmax, norm_sub, rowmax,
)
from probe_v130_issue28_weak_wall import ensure_real


def full_vector(vector, mpi):
    """Collect only the small diagnostic (3300 DOF) PETSc distributed vector."""
    return np.concatenate(mpi.allgather(
        ensure_real(vector.getArray(readonly=True), "PETSc native vector").copy()
    ))


def set_full(vector, values):
    lo, hi = vector.getOwnershipRange()
    if lo < hi:
        vector.setValues(np.arange(lo, hi, dtype=PETSc.IntType),
                         np.asarray(values[lo:hi], dtype=float))
    vector.assemble()


def report_block_norms(vector, names, intervals, lo, hi):
    result = {}
    for name in names:
        iset = indexset(indices((name,), intervals), lo, hi, vector.getComm())
        try:
            result[name] = norm_sub(vector, iset)
        finally:
            iset.destroy()
    return result


def audit(s, mpi, mesh):
    label = "spinup"
    dnames, anames, iv, n = layout(s)
    if tuple(dnames) != ("h", "p") or "Y" not in anames:
        raise RuntimeError("Original rotating SWE variable structure changed")
    d = np.sort(indices(dnames, iv))
    a = np.sort(indices(anames, iv))
    y = np.sort(indices(("Y",), iv))
    regular_names = tuple(name for name in anames if name != "Y")
    regular = np.sort(indices(regular_names, iv))
    expected = {
        .2: (1452, 1848, 1716, 132, 3300),
        .1: (4758, 5514, 5262, 252, 10272),
    }
    if (len(d), len(a), len(regular), len(y), n) != expected[mesh]:
        raise RuntimeError("Unexpected original curved-P2 mesh dimensions")
    initial = ensure_real(s.gf_model.from_variables(), "original physical q0").copy()
    if initial.shape != (n,) or not np.isfinite(initial).all():
        raise RuntimeError("Invalid original prescribed GetFEM state")
    start_d = initial[d].copy()
    start_Y = initial[y].copy()
    if not np.array_equal(start_Y, np.zeros_like(start_Y)):
        raise RuntimeError("Original Y0 initial guess was changed")
    if label == "spinup":
        start = initial.copy()
        for name, value in (
            ("e_h", float(RHO*GRAVITY*H_INITIAL)),
            ("Y_n", -float(RHO*GRAVITY*H_INITIAL)),
        ):
            first, count = iv[name]
            start[first:first+count] = value
    elif label == "rotating":
        start = initial.copy()
    else:
        raise ValueError("Unknown physical case")

    snapshot_ready = getattr(s, "_explicit_rhs_ready", False)
    snapshot_time = getattr(s, "_explicit_rhs_time", None)
    snap_rhs = getattr(s, "_explicit_rhs", None)
    snap_rhs_copy = snap_rhs.copy() if snap_rhs is not None else None
    snap_global = s.GlobalVec.copy()
    owned_lo, owned_hi = s.tangent_stiffness.getOwnershipRange()
    objects = []

    def keep(obj):
        objects.append(obj)
        return obj

    try:
        z = keep(s.tangent_stiffness.createVecRight())
        set_full(z, start)
        zero_dot = keep(z.duplicate())
        zero_dot.set(0.)
        f = keep(z.duplicate())
        idyn = keep(indexset(d, owned_lo, owned_hi, z.getComm()))
        ialg = keep(indexset(a, owned_lo, owned_hi, z.getComm()))
        ireg = keep(indexset(regular, owned_lo, owned_hi, z.getComm()))
        iy = keep(indexset(y, owned_lo, owned_hi, z.getComm()))
        s.IFunction(None, 0., z, zero_dot, f)
        initial_fa = norm_sub(f, ialg)
        projection = {
            "type": "known-exact-analytical-spinup-rest",
            "original_zero_algebraic_Fa_l2": initial_fa,
        }

        # No LSQR or algebraic fitting for the exact positive rest state.
        candidate = start.copy()

        s.IFunction(None, 0., z, zero_dot, f)
        fa_candidate = norm_sub(f, ialg)
        fd_candidate = norm_sub(f, idyn)
        projection["native_final_Fa_l2"] = fa_candidate
        projection["native_final_algebraic_blocks_l2"] = report_block_norms(
            f, anames, iv, owned_lo, owned_hi
        )
        projection["exact_discrete_spinup_rest_t0_candidate"] = True
        projection["raw_fixed_mesh_Fa_below_1e-9"] = bool(fa_candidate <= 1e-9)

        if label == "spinup" and (fa_candidate > 1e-11
                                or fd_candidate > 1e-11):
            raise RuntimeError("Known exact stationary reference failed")

        # The actual dynamic source at the projected effort/control
        # fields, not the misleading zero-effort original guess.
        original_f = full_vector(f, mpi)
        dynamic_rhs = -original_f[d]

        # The new compatible rest-to-spinup wall has g'(0)=0 exactly,
        # whereas the ORIGINAL decreasing envelope remains 1 for
        # 0 <= t < .03. Thus true partial_t F_a(0) = 0 analytically.
        # Verify the smoothstep with real GetFEM off-t0 evaluations
        # instead of taking a spurious fixed central FD at t=0.
        time_eps = 1e-4
        tau = .02
        q1 = time_eps/tau
        q2 = 2*time_eps/tau
        g1 = 3*q1*q1-2*q1*q1*q1
        g2 = 3*q2*q2-2*q2*q2*q2
        analytic_ratio = g2/g1
        one = keep(z.duplicate())
        two = keep(z.duplicate())
        s.IFunction(None, time_eps, z, zero_dot, one)
        s.IFunction(None, 2*time_eps, z, zero_dot, two)
        time_diff_one = full_vector(one, mpi)-original_f
        time_diff_two = full_vector(two, mpi)-original_f
        u = indices(("U",), iv)
        not_u = np.setdiff1d(np.arange(n), u)
        other_time_change_l2 = max(
            np.linalg.norm(time_diff_one[not_u]),
            np.linalg.norm(time_diff_two[not_u]),
        )
        u_time_change_one_l2 = np.linalg.norm(time_diff_one[u])
        u_time_change_two_l2 = np.linalg.norm(time_diff_two[u])
        observed_ratio = (
            u_time_change_two_l2/u_time_change_one_l2
            if u_time_change_one_l2 else float("nan")
        )
        if (u_time_change_one_l2 <= 1e-15
                or other_time_change_l2 > 1e-11
                or not np.isfinite(observed_ratio)
                or abs(observed_ratio/analytic_ratio-1.) > 2e-6):
            raise RuntimeError(
                f"Original GetFEM C1 wall time law invalid: "
                f"observed={observed_ratio:.9g}, "
                f"analytic={analytic_ratio:.9g}, "
                f"other={other_time_change_l2:.9g}"
            )
        # Analytically EXACT, after independently verifying the
        # native original smoothstep and its quadratic onset.
        algebraic_time_rhs = np.zeros(len(a), dtype=float)

        # Assemble the *real* mass and nonlinear residual Jacobian at
        # the chosen regular-algebraic state, not zero algebraic guesses.
        s._sync_lumped_trial_state(0., z)
        s.assemble_nl_stiffness()
        Mdd = keep(s.tangent_mass.createSubMatrix(idyn, idyn))
        JdY = keep(s.tangent_stiffness.createSubMatrix(idyn, iy))
        Jad = keep(s.tangent_stiffness.createSubMatrix(ialg, idyn))
        Jar = keep(s.tangent_stiffness.createSubMatrix(ialg, ireg))
        JaY = keep(s.tangent_stiffness.createSubMatrix(ialg, iy))
        if matmax(JaY, mpi) > 1e-12 or matmax(JdY, mpi) <= 1e-12:
            raise RuntimeError("Original rotating traction structure changed")
        if rowmax(s.tangent_mass, a, mpi) > 1e-12:
            raise RuntimeError("Original algebraic mass rows not zero")

        matnest = keep(PETSc.Mat().createNest(
            [[Mdd, None, JdY],
             [Jad, Jar, None]],
            comm=z.getComm(),
        ))
        matnest.assemble()
        if matnest.getType() != PETSc.Mat.Type.NEST:
            raise RuntimeError("Differentiated operator is not a PETSc MatNest")

        # PETSc Mat.convert("aij") WITHOUT an explicit out= converts
        # the input in place, destroying its native MatNest index sets.
        # Preserve the original nest AND its row/column maps for MPI
        # RHS packing and solution reconstruction.
        nest_row_is, nest_col_is = matnest.getNestISs()
        bordered = keep(matnest.convert("aij", out=PETSc.Mat()))
        if matnest.getType() != PETSc.Mat.Type.NEST:
            raise RuntimeError("Original MatNest lost during AIJ conversion")
        if bordered.getSize() != (n,n):
            raise RuntimeError("Unexpected differentiated operator dimensions")

        # CRITICAL MPI DETAIL:
        # MATNEST defaults to the *per-rank* block order, not to the
        # global concatenation [dynamic, regular algebraic, Y].
        # Its converted AIJ uses the original MATNEST IS permutations.
        # Writing a naive np.concatenate(forcing) or slicing a gathered
        # solution at global block sizes works in serial but silently
        # changes the physical RHS and variable meaning in MPI > 1.
        # Use the global index maps retrieved from the preserved Nest.
        if len(nest_row_is) != 2 or len(nest_col_is) != 3:
            raise RuntimeError("Unexpected native MatNest row/column blocks")

        # These prototype Vec ownership ranges describe the compact,
        # per-block PETSc ordering used by Mat.createSubMatrix.
        dyn_row_template = keep(Mdd.createVecLeft())
        alg_row_template = keep(Jar.createVecLeft())
        dyn_col_template = keep(Mdd.createVecRight())
        reg_col_template = keep(Jar.createVecRight())
        y_col_template = keep(JdY.createVecRight())
        row_blocks = (
            (nest_row_is[0], dyn_row_template, dynamic_rhs, "dynamic RHS"),
            (nest_row_is[1], alg_row_template, algebraic_time_rhs,
             "differentiated algebraic RHS"),
        )
        col_blocks = (
            (nest_col_is[0], dyn_col_template, len(d), "qdot"),
            (nest_col_is[1], reg_col_template, len(regular), "a_regular_dot"),
            (nest_col_is[2], y_col_template, len(y), "Y0"),
        )

        def local_native_nest_indices(iset, template, length, kind):
            # This is a mapping between two *different* parallel Vec
            # layouts. The returned IS indices are in the converted
            # AIJ coordinates, while template owns compact block
            # coordinates. Do not assume their global offsets coincide.
            mapped = np.asarray(iset.getIndices(), dtype=PETSc.IntType)
            first, stop = template.getOwnershipRange()
            if len(mapped) != stop-first or template.getSize() != length:
                raise RuntimeError(
                    f"Native MatNest {kind} ownership differs from subblock "
                    f"({len(mapped)} mapped vs {stop-first} owned)"
                )
            if np.any(mapped < 0) or np.any(mapped >= n):
                raise RuntimeError(f"Native MatNest {kind} mapping out of bounds")
            return mapped, first, stop

        def check_is_partition(blocks, expected):
            allids = []
            for iset, template, size, name in blocks:
                mapped, _, _ = local_native_nest_indices(
                    iset, template, size, name
                )
                allids.append(np.concatenate(
                    mpi.allgather(np.asarray(mapped, dtype=np.int64))
                ))
            union = np.concatenate(allids)
            if (len(union) != expected
                    or not np.array_equal(np.sort(union),
                                          np.arange(expected))):
                raise RuntimeError(
                    "MatNest global index sets do not partition converted AIJ"
                )
            return allids

        check_is_partition(
            ((iset, template, len(values), name)
             for iset, template, values, name in row_blocks), n
        )
        check_is_partition(col_blocks, n)

        rhs_true = keep(bordered.createVecLeft())
        rhs_true.set(0.)
        for iset, template, values, name in row_blocks:
            mapped, first, stop = local_native_nest_indices(
                iset, template, len(values), name
            )
            if len(mapped):
                rhs_true.setValues(mapped, values[first:stop],
                                   addv=PETSc.InsertMode.INSERT_VALUES)
        rhs_true.assemble()

        # Independently undo the IS permutation before factorization:
        # this also catches mistakes in RHS placement on TWO ranks.
        rhs_ai_full = full_vector(rhs_true, mpi)
        rhs_remapped = []
        for iset, template, values, name in row_blocks:
            mapped = np.asarray(iset.getIndices(), dtype=np.int64)
            rhs_remapped.append(np.concatenate(
                mpi.allgather(rhs_ai_full[mapped])
            ))
        rhs_roundtrip = np.concatenate(rhs_remapped)
        forcing = np.concatenate((dynamic_rhs, algebraic_time_rhs))
        rhs_layout_error = float(np.linalg.norm(rhs_roundtrip - forcing))
        if rhs_layout_error > 1e-13:
            raise RuntimeError(
                f"MatNest row IS mispackaged physical RHS: {rhs_layout_error}"
            )

        x = keep(bordered.createVecRight())
        x.set(0.)
        direct = keep(PETSc.KSP().create(comm=bordered.getComm()))
        direct.setOptionsPrefix("scrimp_issue28_true_rhs_")
        direct.setOperators(bordered)
        direct.setType(PETSc.KSP.Type.PREONLY)
        direct.getPC().setType(PETSc.PC.Type.LU)
        direct.getPC().setFactorSolverType("mumps")
        direct.solve(rhs_true, x)
        reason = int(direct.getConvergedReason())
        if reason <= 0:
            raise RuntimeError(f"True RHS bordered solve failed (reason={reason})")

        solution_ai_full = full_vector(x, mpi)
        recovered = []
        for iset, template, size, name in col_blocks:
            mapped, first, stop = local_native_nest_indices(
                iset, template, size, name
            )
            local_values = solution_ai_full[np.asarray(mapped, dtype=np.int64)]
            recovered.append(np.concatenate(mpi.allgather(local_values)))
        qdot, ardot, y0 = recovered
        correction = np.concatenate(recovered)
        if (len(qdot) != len(d) or len(ardot) != len(regular)
                or len(y0) != len(y) or not np.all(np.isfinite(correction))):
            raise RuntimeError("MatNest IS physical solution reconstruction failed")
        z_trial = keep(z.copy())
        localY = y[(y >= owned_lo) & (y < owned_hi)]
        if len(localY):
            j = np.searchsorted(y,localY)
            z_trial.setValues(
                np.asarray(localY,dtype=PETSc.IntType), y0[j],
                addv=PETSc.InsertMode.INSERT_VALUES,
            )
        z_trial.assemble()
        zdot = keep(z.duplicate())
        zdot.set(0.)
        localD = d[(d >= owned_lo) & (d < owned_hi)]
        if len(localD):
            j = np.searchsorted(d, localD)
            zdot.setValues(np.asarray(localD, dtype=PETSc.IntType),
                           qdot[j])
        zdot.assemble()
        residual = keep(z.duplicate())
        s.IFunction(None,0.,z_trial,zdot,residual)
        fd_verified = norm_sub(residual,idyn)
        fa_verified = norm_sub(residual,ialg)
        if not np.array_equal(full_vector(z_trial,mpi)[d],start_d):
            raise RuntimeError("True RHS trial moved the prescribed physical q0")
        if abs(fa_verified - fa_candidate)>1e-9:
            raise RuntimeError("Injecting traction Y0 modified algebraic Fa")
        # The *original* native PDE equations are the acceptance gate
        # for this diagnostic. Do not rely on the KSP residual alone:
        # the previous MPI bug gave a perfect KSP residual but
        # 0.3985857557 on the actual dynamic GetFEM IFunction.
        if fd_verified > 1e-8:
            raise RuntimeError(
                "Physical RHS MatNest candidate does not satisfy original "
                f"GetFEM dynamics: |Fd|={fd_verified:.9g}"
            )

        rec = keep(bordered.createVecLeft())
        bordered.mult(x,rec)
        rec.axpy(-1.,rhs_true)
        bordered_relative_residual = float(
            rec.norm()/max(1., rhs_true.norm())
        )

        # This report is diagnostic ONLY at *both* points. Even
        # a tiny candidate Fa raw coefficient norm cannot replace
        # mesh-aware physical boundary compatibility measurements.
        report = {
            "schema":"scrimp-v130-issue28-spinup-first-differentiated-v1",
            "status":"DIAGNOSTIC_SPINUP_T0_TRUE_RHS_SOLVED_NOT_INITIALIZED",
            "case":label,"nominal_mesh":mesh,"mpi_ranks":int(mpi.size),
            "original_q0_preserved_bitwise":True,
            "dofs":n,"dynamic_dofs":len(d),
            "regular_algebraic_dofs":len(regular),
            "traction_Y_dofs":len(y),
            "original_P2_geometry":True,
            "rotation_after_spinup":ROTATION,
            "initial_physical_p0_zero":bool(np.count_nonzero(initial[indices(("p",),iv)])==0),
            "physical_wall_facet_tangent":True,
            "spinup_duration":.02,
            "analytic_partial_t_Fa_at_t0_exact_zero":True,
            "native_wall_smoothstep_ratio":float(observed_ratio),
            "analytic_wall_smoothstep_ratio":float(analytic_ratio),
            "native_wall_time_change_U_delta_l2":float(u_time_change_one_l2),
            "native_wall_time_change_U_2delta_l2":float(u_time_change_two_l2),
            "native_time_change_other_rows_l2":float(other_time_change_l2),
            "initial_projection":projection,
            "native_Fd_at_zero_zdot_before_Y0_l2":fd_candidate,
            "native_algebraic_partial_t_Fa_l2":0.0,
            "dynamic_physical_rhs_l2":float(np.linalg.norm(dynamic_rhs)),
            "bordered_dynamic_algebraic_operator_shape":[n,n],
            "bordered_true_rhs_l2":float(rhs_true.norm()),
            "bordered_KSP_reason":reason,
            "bordered_relative_true_rhs_residual":bordered_relative_residual,
            "matnest_original_row_and_column_IS_used":True,
            "matnest_physical_rhs_roundtrip_l2":rhs_layout_error,
            "matnest_solution_reconstructed_by_block_IS":True,
            "candidate_qdot_coefficients":qdot.tolist(),
            "candidate_regular_algebraic_derivative_coefficients":ardot.tolist(),
            "candidate_traction_Y0_coefficients":y0.tolist(),
            "candidate_qdot_l2":float(np.linalg.norm(qdot)),
            "candidate_regular_algebraic_derivative_l2":float(np.linalg.norm(ardot)),
            "candidate_traction_Y0_l2":float(np.linalg.norm(y0)),
            "candidate_traction_Y0_maxabs":float(np.max(np.abs(y0))),
            "original_native_IFunction_Fd_with_candidate_derivative_Y0_l2":fd_verified,
            "original_native_IFunction_Fa_with_candidate_Y0_l2":fa_verified,
            "independent_spinup_t0_algebraic_compatibility_qualified":bool(fa_candidate<1e-11),
            "physical_rotating_Y0_accepted":False,
            "true_initialization_TS_executed":False,
            "true_physical_TS_executed":False,
            "mesh_independent_physical_constraint_check_qualified":False,
            "first_differentiated_rank_proven":False,
            "no_pseudo_or_physical_TS":True,
            "full_original_IFunction_checked":True,
            "fixed_initial_h0_p0_preserved_bitwise":True,
            "interpretation":(
                "At separate compatible physical p0=0 smoothstep spinup, "
                "true original Fd and exact analytically verified Fa,t=0 "
                "drive the native MatNest differentiated equations. "
                "Resulting diagnostic qdot0,Y0 are checked in the unchanged "
                "full GetFEM IFunction. No physical TS or accepted "
                "general high-index initializer."
            ),
        }
        if label == "spinup" and (
            fd_verified > 1e-10 or fa_verified > 1e-11
            or np.linalg.norm(y0) > 1e-8
            or np.linalg.norm(qdot) > 1e-8
            or np.linalg.norm(ardot) > 1e-8
            or bordered_relative_residual > 1e-10
        ):
            raise RuntimeError(
                "Compatible spinup t0 did not yield near-zero physical derivatives"
            )
        if mpi.rank == 0:
            print(
                f"ISSUE28 TRUE PHYSICAL RHS case={label} MPI={mpi.size}: "
                f"{report['status']}, projected |Fa|={fa_candidate:.9g}, "
                f"physical RHS={rhs_true.norm():.9g}, "
                f"|Y0|={np.linalg.norm(y0):.9g}, "
                f"native |Fd(qdot,Y0)|={fd_verified:.9g}, "
                f"native |Fa|={fa_verified:.9g}",
                flush=True,
            )
        return report
    finally:
        s.gf_model.set_time(0.)
        s.gf_model.to_variables(initial.copy())
        snap_global.copy(s.GlobalVec)
        s._explicit_rhs_ready = snapshot_ready
        s._explicit_rhs_time = snapshot_time
        if snap_rhs_copy is not None and snap_rhs is not None:
            snap_rhs_copy.copy(snap_rhs)
        if not mpi.allreduce(
            np.array_equal(
                np.asarray(s.gf_model.from_variables(),dtype=float),
                initial
            ), op=MPI.LAND
        ):
            raise RuntimeError("Read-only physical-RHS pilot did not roll back state")
        for obj in reversed(objects):
            obj.destroy()
        snap_global.destroy()
        if snap_rhs_copy is not None:
            snap_rhs_copy.destroy()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mesh",required=True,type=float,choices=(.2,.1))
    parser.add_argument("--out",required=True,type=Path)
    args=parser.parse_args()
    mpi=PETSc.COMM_WORLD.tompi4py()
    s=build_system(
        gamma=0.,mesh=args.mesh,geometry_order=2,height_order=2,
        final_time=.08,save_dt=.002,max_step=.002,
        ramp_start=.03,ramp_end=.06,
        init_step=True,ts_type="bdf",adapt_type="none",
        init_dt=1e-9,
        rotation=ROTATION,
        initial_profile="rigid",
        wall_velocity_profile="facet_tangent",
        algebraic_initial_guess="zero",start_from_rest_duration=.02,
        init_method="bdf",
    )
    result={}
    def intercept(self):
        result["report"]=audit(self,mpi,args.mesh)
        raise AuditComplete
    with patch.object(S.DPHS,"init_step",intercept):
        try:
            s.solve()
        except AuditComplete:
            pass
        else:
            raise RuntimeError("Native original solver init intercept missed")
    if "report" not in result:
        raise RuntimeError("No true physical RHS diagnostic result")
    if mpi.rank==0:
        args.out.mkdir(parents=True,exist_ok=True)
        (args.out/"spinup_first_derivative.json").write_text(
            json.dumps(result["report"],indent=2)+"\n",
            encoding="utf-8",
        )
    mpi.Barrier()


if __name__=="__main__":
    main()
