# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Issue #28: physical compatibility of the ORIGINAL rotating P2 SWE data.

Read-only native GetFEM integral diagnostics before SCRIMP's initializer.
The actual original h0,p0, curved disk and facet-tangent wall law stay
unchanged. Compare physical boundary L2 (normal/tangential mismatch) across
meshes independently of any LSQR. At the baseline mesh only, project
regular algebraic variables with distributed LSQR at two tolerances and
measure constitutive L2 errors, true IFunction Fa and normal equations.
Neither small L2 norms nor a converged LSQR alone certify exact DAE data.
No pseudo-BDF initialization, accepted state or physical TS step.
"""

import argparse
import json
from pathlib import Path
from unittest.mock import patch

import getfem as gf
import numpy as np
from mpi4py import MPI
from petsc4py import PETSc
import scrimp as S

from viscous_swe_2d_mpi import (
    build_system, H_INITIAL, RHO, GRAVITY, ROTATION,
)
from probe_v130_issue28_structure import (
    AuditComplete, indices, indexset, layout, norm_sub,
)

DOMAIN = 1
WALL = 10


def real_array(values, name):
    """Never silently discard an unexpected complex diagnostic quantity."""
    z = np.asarray(values)
    if not np.all(np.isfinite(z)):
        raise FloatingPointError(f"Nonfinite {name}")
    if np.iscomplexobj(z) and np.max(np.abs(z.imag), initial=0.) > (
        1e-12 * max(1., float(np.max(np.abs(z.real), initial=0.)))
    ):
        raise FloatingPointError(f"Unexpected imaginary component in {name}")
    return np.asarray(z.real, dtype=float)


def scalar_integral(system, mpi, expression, region):
    """Return both GetFEM raw per-rank values and their sum.

    Whether scalar assembly is local or replicated depends on the
    GetFEM MPI build. Ratios to the assembled measure are robust to
    either convention. The host comparator determines absolute L2
    normalization from the corresponding serial area/perimeter.
    """
    result = gf.asm_generic(
        system.domain._int_method[0], 0, expression, region, system.gf_model
    )
    values = real_array(result, expression).reshape(-1)
    if len(values) != 1:
        raise RuntimeError(f"Expected scalar integral: {expression}")
    raw = [float(x) for x in mpi.allgather(float(values[0]))]
    if not np.all(np.isfinite(raw)):
        raise FloatingPointError("Nonfinite MPI integral contributions")
    return {"per_rank": raw, "sum": float(np.sum(raw))}


def physical_integrals(system, mpi, rotation, projected):
    """Integrate real FE traces, not coefficient Euclidean norms."""
    omega = f"({rotation:.17g})"
    pnormal = "((p/rho).Normal)"
    ptangent = "((p/rho).Tangent)"
    mismatch = f"((p/rho)-{omega}*Tangent)"
    velocity_err = "(e_p-(p/rho))"
    wall_err = "(e_p-U)"
    wall_control = f"(U-U_n*Normal-{omega}*Tangent)"
    pressure_err = "(e_h-rho*g*h-0.5*(p.p)/rho)"
    expressions = {
        "area": ("1", DOMAIN),
        "perimeter": ("1", WALL),
        "initial_p_normal_squared": (f"{pnormal}*{pnormal}", WALL),
        "initial_p_tangent_mismatch_squared":
            (f"({ptangent}-{omega})*({ptangent}-{omega})", WALL),
        "initial_p_wall_vector_mismatch_squared":
            (f"{mismatch}.{mismatch}", WALL),
        "initial_p_tangent_squared": (f"{ptangent}*{ptangent}", WALL),
    }
    if projected:
        expressions.update({
            "interior_velocity_constitutive_squared":
                (f"{velocity_err}.{velocity_err}", DOMAIN),
            "interior_velocity_constitutive_h_weighted_squared":
                (f"h*({velocity_err}.{velocity_err})", DOMAIN),
            "interior_pressure_constitutive_squared":
                (f"{pressure_err}*{pressure_err}", DOMAIN),
            "boundary_ep_minus_U_squared":
                (f"{wall_err}.{wall_err}", WALL),
            "boundary_control_law_squared":
                (f"{wall_control}.{wall_control}", WALL),
            "boundary_normal_control_squared": ("U_n*U_n", WALL),
        })
    raw = {
        name: scalar_integral(system, mpi, expression, region)
        for name, (expression, region) in expressions.items()
    }
    area, perimeter = raw["area"]["sum"], raw["perimeter"]["sum"]
    if area <= 0 or perimeter <= 0:
        raise RuntimeError("Nonpositive native mesh measure")
    if not all(v["sum"] >= -1e-10 * max(1., area, perimeter)
               for k, v in raw.items() if k not in ("area", "perimeter")):
        raise RuntimeError("Negative physical squared integral")
    rms = {}
    for name, v in raw.items():
        if name in ("area", "perimeter"):
            continue
        measure = area if name.startswith("interior_") else perimeter
        rms[name.replace("_squared", "_rms")] = float(
            np.sqrt(max(0., v["sum"]) / measure)
        )
    return {"raw_integrals": raw, "rms": rms}


def native_fa(system, v, algebraic, iv, lo, hi, work):
    zero_dot = work(v.duplicate())
    zero_dot.set(0.)
    f = work(v.duplicate())
    system.IFunction(None, 0., v, zero_dot, f)
    ia = work(indexset(indices(algebraic, iv), lo, hi, v.getComm()))
    blocks = {}
    for name in algebraic:
        ix = work(indexset(indices((name,), iv), lo, hi, v.getComm()))
        blocks[name] = norm_sub(f, ix)
    return f, ia, float(norm_sub(f, ia)), blocks


def audit(system, mpi, mesh, rotation):
    differential, algebraic, iv, n = layout(system)
    if differential != ("h", "p") or "Y" not in algebraic:
        raise RuntimeError("Original SWE variable partition changed")
    original = real_array(system.gf_model.from_variables(), "initial GetFEM")
    if original.size != n:
        raise RuntimeError("Original initial state dimension mismatch")
    fixed = np.sort(indices(differential, iv))
    traction = np.sort(indices(("Y",), iv))
    regular = np.sort(indices(tuple(x for x in algebraic if x != "Y"), iv))
    lo, hi = system.tangent_stiffness.getOwnershipRange()
    created = []

    def keep(obj):
        created.append(obj)
        return obj

    snap_global = system.GlobalVec.copy()
    snap_ready = getattr(system, "_explicit_rhs_ready", False)
    snap_time = getattr(system, "_explicit_rhs_time", None)
    snap_rhs = getattr(system, "_explicit_rhs", None)
    snap_rhs_copy = snap_rhs.copy() if snap_rhs is not None else None

    try:
        z = keep(system.tangent_stiffness.createVecRight())
        if hi > lo:
            z.setValues(np.arange(lo, hi, dtype=PETSc.IntType),
                        original[lo:hi])
        z.assemble()
        system._sync_lumped_trial_state(0., z)
        boundary = physical_integrals(system, mpi, rotation, False)
        report = {
            "schema": "scrimp-v130-issue28-physical-compatibility-v1",
            "status": "DIAGNOSTIC_ONLY_NOT_INITIALIZED",
            "mpi_ranks": mpi.size, "mesh_nominal_h": mesh,
            "geometry_order": 2, "height_order": 2,
            "rotation": rotation, "dofs": n,
            "boundary_at_prescribed_q0": boundary,
            "projection": None,
            "original_q0_preserved_bitwise": False,
            "state_rollback_bitwise": False,
            "accepted_rotating_Y0": False,
            "exact_DAE_compatibility_qualified": False,
            "initialization_TS_executed": False,
            "physical_TS_executed": False,
        }

        # Reference solution at rest: exact e_h=0.5, Y_n=-0.5.
        # Rotating h=.2 is projected with the *same* 1716 regular
        # algebraic degrees as the successful physical-RHS pilot.
        if rotation == 0.:
            known = original.copy()
            for name, value in (("e_h", RHO*GRAVITY*H_INITIAL),
                                ("Y_n", -RHO*GRAVITY*H_INITIAL)):
                first, length = iv[name]
                known[first:first+length] = value
            candidate = keep(z.copy())
            if hi > lo:
                candidate.setValues(np.arange(lo, hi, dtype=PETSc.IntType),
                                    known[lo:hi])
            candidate.assemble()
            _, _, fa, blocks = native_fa(
                system, candidate, algebraic, iv, lo, hi, keep
            )
            system._sync_lumped_trial_state(0., candidate)
            integrals = physical_integrals(system, mpi, rotation, True)
            if fa > 1e-10:
                raise RuntimeError(f"Exact rest control has Fa={fa}")
            report["projection"] = {
                "type": "exact-analytical-rest",
                "native_Fa_l2": fa, "native_Fa_blocks_l2": blocks,
                "physical_integrals": integrals,
            }
        elif abs(mesh - .2) < 1e-14:
            zero_dot = keep(z.duplicate())
            zero_dot.set(0.)
            f = keep(z.duplicate())
            ia = keep(indexset(indices(algebraic, iv), lo, hi, z.getComm()))
            ir = keep(indexset(regular, lo, hi, z.getComm()))
            system.IFunction(None, 0., z, zero_dot, f)
            original_fa = norm_sub(f, ia)
            system._sync_lumped_trial_state(0., z)
            system.assemble_nl_stiffness()
            A = keep(system.tangent_stiffness.createSubMatrix(ia, ir))
            rhs = keep(A.createVecLeft())
            sub = f.getSubVector(ia)
            try:
                rhs.getArray()[:] = -real_array(
                    sub.getArray(readonly=True), "original algebraic RHS"
                )
            finally:
                f.restoreSubVector(ia, sub)
            rhs.assemble()
            trials = []
            standard_delta = None
            for kind, rtol, atol, max_it in (
                ("standard", 1e-12, 1e-13, 10000),
                ("tight", 1e-13, 1e-14, 20000),
            ):
                delta = keep(A.createVecRight())
                delta.set(0.)
                ksp = keep(PETSc.KSP().create(comm=A.getComm()))
                ksp.setOptionsPrefix(f"scrimp_issue28_phys_{kind}_")
                ksp.setOperators(A)
                ksp.setType("lsqr")
                ksp.getPC().setType("none")
                ksp.setTolerances(rtol=rtol, atol=atol, max_it=max_it)
                ksp.solve(rhs, delta)
                reason = int(ksp.getConvergedReason())
                correction_norm = float(delta.norm())
                if standard_delta is None:
                    standard_delta = keep(delta.copy())
                    coefficient_change = 0.
                else:
                    difference = keep(delta.copy())
                    difference.axpy(-1., standard_delta)
                    coefficient_change = float(
                        difference.norm() / max(1., standard_delta.norm())
                    )
                corr_local = real_array(
                    delta.getArray(readonly=True), "LSQR correction"
                )
                owned = regular[(regular >= lo) & (regular < hi)]
                if len(owned) != len(corr_local):
                    raise RuntimeError("PETSc reduced block ownership mismatch")
                candidate = keep(z.copy())
                if len(owned):
                    candidate.setValues(
                        np.asarray(owned, dtype=PETSc.IntType), corr_local,
                        addv=PETSc.InsertMode.ADD_VALUES
                    )
                candidate.assemble()
                system.IFunction(None, 0., candidate, zero_dot, f)
                fa = norm_sub(f, ia)
                blocks = {}
                for name in algebraic:
                    ix = keep(indexset(indices((name,), iv),
                                       lo, hi, z.getComm()))
                    blocks[name] = norm_sub(f, ix)
                fa_vec = f.getSubVector(ia)
                try:
                    gradient = keep(A.createVecRight())
                    A.multTranspose(fa_vec, gradient)
                    stationarity = float(gradient.norm())
                finally:
                    f.restoreSubVector(ia, fa_vec)
                system._sync_lumped_trial_state(0., candidate)
                trial_full = real_array(
                    system.gf_model.from_variables(), "projected trial"
                )
                fixed_preserved = np.array_equal(
                    trial_full[fixed], original[fixed]
                ) and np.array_equal(trial_full[traction], original[traction])
                if not mpi.allreduce(bool(fixed_preserved), op=MPI.LAND):
                    raise RuntimeError("Projection changed prescribed q0 or Y")
                integrals = physical_integrals(system, mpi, rotation, True)
                trials.append({
                    "kind": kind, "rtol": rtol, "atol": atol,
                    "max_iterations": max_it,
                    "converged_reason": reason,
                    "iterations": int(ksp.getIterationNumber()),
                    "reported_residual": float(ksp.getResidualNorm()),
                    "algebraic_correction_l2": correction_norm,
                    "relative_coefficient_change_from_standard": coefficient_change,
                    "native_Fa_l2": fa, "native_Fa_blocks_l2": blocks,
                    "normal_equation_AT_Fa_l2": stationarity,
                    "physical_integrals": integrals,
                })
            report["projection"] = {
                "type": "distributed-rectangular-LSQR-two-tolerances",
                "original_zero_guess_native_Fa_l2": original_fa,
                "trials": trials,
            }
        else:
            report["projection"] = {
                "type": "not-run-refinement-boundary-only",
                "reason": "Prescribed q0 trace norms do not require LSQR",
            }

        report["original_q0_preserved_bitwise"] = True
        report["state_rollback_bitwise"] = True
        if mpi.rank == 0:
            print(
                f"ISSUE28 PHYSICAL COMPATIBILITY h={mesh:g} "
                f"omega={rotation:g} MPI={mpi.size} "
                f"wall-normal-RMS={boundary['rms']['initial_p_normal_rms']:.9g} "
                f"wall-tangential-error-RMS="
                f"{boundary['rms']['initial_p_tangent_mismatch_rms']:.9g}",
                flush=True
            )
        return report
    finally:
        system.gf_model.set_time(0.)
        system.gf_model.to_variables(original.copy())
        snap_global.copy(system.GlobalVec)
        system._explicit_rhs_ready = snap_ready
        system._explicit_rhs_time = snap_time
        if snap_rhs_copy is not None and snap_rhs is not None:
            snap_rhs_copy.copy(snap_rhs)
        restored = real_array(
            system.gf_model.from_variables(), "restored GetFEM"
        )
        if not mpi.allreduce(bool(np.array_equal(restored, original)),
                             op=MPI.LAND):
            raise RuntimeError("Read-only compatibility audit failed rollback")
        for obj in reversed(created):
            obj.destroy()
        snap_global.destroy()
        if snap_rhs_copy is not None:
            snap_rhs_copy.destroy()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--mesh", type=float, choices=(.2, .1), required=True)
    p.add_argument("--case", choices=("rest", "rotating"), required=True)
    p.add_argument("--out", type=Path, required=True)
    args = p.parse_args()
    mpi = PETSc.COMM_WORLD.tompi4py()
    omega = 0. if args.case == "rest" else ROTATION
    system = build_system(
        gamma=0., mesh=args.mesh, geometry_order=2, height_order=2,
        final_time=.02, save_dt=.002, max_step=.002,
        ramp_start=.003, ramp_end=.009,
        init_step=True, ts_type="bdf", adapt_type="none",
        init_dt=1e-9, rotation=omega,
        initial_profile="tangent", wall_velocity_profile="facet_tangent",
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
            raise RuntimeError("Probe unexpectedly executed TS")
    if "report" not in result:
        raise RuntimeError("Native compatibility audit did not execute")
    if mpi.rank == 0:
        args.out.mkdir(parents=True, exist_ok=True)
        (args.out / "physical_compatibility.json").write_text(
            json.dumps(result["report"], indent=2)+"\n", encoding="utf-8"
        )
    mpi.Barrier()


if __name__ == "__main__":
    main()
