# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""#28: direct *weak* FE wall-mismatch moments at original rotating q0.

GetFEM order-one generic assembly evaluates the residual functional
  b(v) = integral_Gamma (p0/rho - omega*Rot*Normal).v
against the actual prescribed boundary U test functions. Distinguish
the normal and tangential components, check their reconstructed sum,
and compute a Riesz-scaled L2 norm via the native U-U wall mass block.

This is a READ-ONLY diagnostic at the original h0,p0, with curved P2
geometry, physical wall law, no LSQR, no Y0 and no TS time step.
A nonzero b is evidence of a nonzero *discrete weak* mismatch; alone,
it is not a numerical proof of full DAE incompatibility, which also
requires verifying exact algebraic elimination/trace-space relations.
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

from viscous_swe_2d_mpi import build_system, ROTATION
from probe_v130_issue28_structure import (
    AuditComplete, indices, indexset, layout,
)


def ensure_real(value, label):
    v = np.asarray(value)
    if not np.all(np.isfinite(v)):
        raise FloatingPointError(f"Nonfinite {label}")
    if np.iscomplexobj(v):
        real_max = max(1., float(np.max(abs(v.real), initial=0.)))
        imag_max = float(np.max(abs(v.imag), initial=0.))
        if imag_max > 1e-12 * real_max:
            raise FloatingPointError(f"Non-negligible imaginary {label}")
    return np.asarray(v.real, dtype=np.float64)


def integral(system, expression, region):
    result = ensure_real(gf.asm_generic(
        system.domain._int_method[0], 0, expression, region, system.gf_model
    ), expression).reshape(-1)
    if result.size != 1 or result[0] <= 0.:
        raise RuntimeError(f"Invalid native measure {expression}: {result}")
    return float(result[0])


def gather_scalar_measure(local, mpi):
    values = mpi.allgather(local)
    if mpi.size == 1:
        return local, "serial", values
    # The independently archived #28 physical-compatibility gate found
    # GetFEM scalar integrals to be replicated-global on this image;
    # retain raw contributions and independently check MPI1/2 later.
    if np.allclose(values, values[0], rtol=2e-10, atol=1e-12):
        return float(values[0]), "replicated-global", values
    return float(sum(values)), "summed-local", values


def assemble_wall_vector(system, mpi, iv, n, name, expression):
    raw = ensure_real(gf.asm_generic(
        system.domain._int_method[0], 1, expression, 10, system.gf_model
    ), name).reshape(-1)
    first, count = iv["U"]
    if raw.size == n:
        local = raw[first:first+count].copy()
        layout_name = "full-model-DOF-vector"
    elif raw.size == count:
        local = raw.copy()
        layout_name = "compact-U-vector"
    else:
        raise RuntimeError(
            f"GetFEM order-one {name} vector has size {raw.size}, "
            f"expected {n} or {count}"
        )
    all_vectors = mpi.allgather(local)
    if mpi.size == 1:
        native = local
        mpi_mode = "serial"
    else:
        same = all(
            np.allclose(v, all_vectors[0], rtol=2e-10, atol=1e-12)
            for v in all_vectors[1:]
        )
        # If the results are already global, summing doubles the load.
        # If they are rank-local contributions, they must be added.
        native = all_vectors[0] if same else np.sum(all_vectors, axis=0)
        mpi_mode = "replicated-global" if same else "summed-local"
    if not np.all(np.isfinite(native)):
        raise FloatingPointError(f"Nonfinite assembled {name}")
    return native, {
        "expression": expression,
        "gf_output_layout": layout_name,
        "mpi_vector_reduction_mode": mpi_mode,
        "per_rank_native_vector": [v.tolist() for v in all_vectors],
        "assembled_vector": native.tolist(),
        "coefficient_l2": float(np.linalg.norm(native)),
    }


def audit(s, mpi, mesh, rotation):
    differential, algebraic, iv, n = layout(s)
    if differential != ("h", "p") or "U" not in algebraic:
        raise RuntimeError("Unexpected original FE variable layout")
    original = ensure_real(s.gf_model.from_variables(), "initial state")
    if original.size != n:
        raise RuntimeError("Incorrect original physical state size")
    fixed = indices(differential, iv)
    traction = indices(("Y",), iv)
    owned_lo, owned_hi = s.tangent_stiffness.getOwnershipRange()
    objects = []

    def keep(obj):
        objects.append(obj)
        return obj

    saved_global = s.GlobalVec.copy()
    saved_ready = getattr(s, "_explicit_rhs_ready", False)
    saved_time = getattr(s, "_explicit_rhs_time", None)
    cached_rhs = getattr(s, "_explicit_rhs", None)
    saved_rhs = cached_rhs.copy() if cached_rhs is not None else None

    try:
        z = keep(s.tangent_stiffness.createVecRight())
        if owned_hi > owned_lo:
            z.setValues(np.arange(owned_lo, owned_hi, dtype=PETSc.IntType),
                        original[owned_lo:owned_hi])
        z.assemble()
        s._sync_lumped_trial_state(0., z)
        local_perimeter = integral(s, "1", 10)
        perimeter, perimeter_mode, raw_perimeters = (
            gather_scalar_measure(local_perimeter, mpi)
        )

        omega = f"({rotation:.17g})"
        strong_expression = (
            f"((p/rho)-{omega}*Tangent).((p/rho)-{omega}*Tangent)"
        )
        strong_raw = ensure_real(gf.asm_generic(
            s.domain._int_method[0], 0, strong_expression, 10, s.gf_model
        ), "strong wall mismatch").reshape(-1)
        if strong_raw.size != 1 or strong_raw[0] < -1e-12:
            raise RuntimeError("Invalid native squared strong wall mismatch")
        strong_total, strong_mode, strong_contributions = (
            gather_scalar_measure(float(max(0., strong_raw[0])), mpi)
        )
        strong_rms = float(np.sqrt(max(0., strong_total)/perimeter))
        exprs = {
            "full": f"((p/rho)-{omega}*Tangent).Test_U",
            "normal": "((p/rho).Normal)*(Normal.Test_U)",
            "tangential": (
                f"(((p/rho).Tangent)-{omega})*(Tangent.Test_U)"
            ),
        }
        vectors = {}
        moments = {}
        for name, expression in exprs.items():
            vectors[name], moments[name] = assemble_wall_vector(
                s, mpi, iv, n, name, expression
            )
        decomposition_error = float(np.linalg.norm(
            vectors["full"] - vectors["normal"] - vectors["tangential"]
        ))
        if decomposition_error > (
            1e-9 * max(1., float(np.linalg.norm(vectors["full"])))
        ):
            raise RuntimeError(
                "Native weak wall load does not decompose along n,t"
            )

        # Native derivative J_UU is just -integral_Gamma U.Test_U.
        # The Riesz mass is positive definite on the actual U trace
        # space, independently of whether the weak load is compatible.
        s._sync_lumped_trial_state(0., z)
        s.assemble_nl_stiffness()
        iu = keep(indexset(indices(("U",), iv), owned_lo, owned_hi,
                           z.getComm()))
        mass = keep(s.tangent_stiffness.createSubMatrix(iu, iu))
        if mass.getSize() != (iv["U"][1], iv["U"][1]):
            raise RuntimeError("Unexpected boundary U mass size")
        mass.scale(-1.)
        ksp = keep(PETSc.KSP().create(comm=mass.getComm()))
        ksp.setOptionsPrefix("scrimp_issue28_weak_wall_riesz_")
        ksp.setOperators(mass)
        ksp.setType(PETSc.KSP.Type.PREONLY)
        ksp.getPC().setType(PETSc.PC.Type.LU)
        ksp.getPC().setFactorSolverType("mumps")

        dual = {}
        for name, values in vectors.items():
            b = keep(mass.createVecLeft())
            lo, hi = b.getOwnershipRange()
            if hi > lo:
                b.setValues(np.arange(lo, hi, dtype=PETSc.IntType),
                            values[lo:hi])
            b.assemble()
            x = keep(mass.createVecRight())
            x.set(0.)
            if b.norm() > 0.:
                ksp.solve(b, x)
                if ksp.getConvergedReason() <= 0:
                    raise RuntimeError(f"Boundary Riesz KSP failed for {name}")
            residual = keep(mass.createVecLeft())
            mass.mult(x, residual)
            residual.axpy(-1., b)
            linear_residual = float(residual.norm())
            energy = float(np.real(b.dot(x)))
            if energy < -1e-12 or not np.isfinite(energy):
                raise RuntimeError("Boundary Riesz energy not nonnegative")
            if linear_residual > 1e-8 * max(1., b.norm()):
                raise RuntimeError(f"Boundary Riesz solve inaccurate for {name}")
            dual[name] = {
                "riesz_L2_projected_mismatch_l2": float(
                    np.sqrt(max(0., energy))
                ),
                "riesz_boundary_RMS": float(
                    np.sqrt(max(0., energy)/perimeter)
                ),
                "boundary_mass_solver_relative_residual": float(
                    linear_residual / max(1., b.norm())
                ),
                "boundary_mass_is_negative_native_J_UU": True,
            }

        if dual["full"]["riesz_boundary_RMS"] > strong_rms + 1e-8:
            raise RuntimeError(
                "Discrete boundary Riesz projection exceeds native strong L2"
            )
        now = ensure_real(s.gf_model.from_variables(), "post-measurement state")
        unchanged = np.array_equal(now[fixed], original[fixed]) and (
            np.array_equal(now[traction], original[traction])
        )
        if not mpi.allreduce(bool(unchanged), op=MPI.LAND):
            raise RuntimeError("Wall moment audit moved q0 or initial Y")

        result = {
            "schema": "scrimp-v130-issue28-weak-wall-moments-v1",
            "status": "NATIVE_WEAK_MOMENT_DIAGNOSTIC_ONLY",
            "mpi_ranks": mpi.size,
            "nominal_mesh": mesh, "geometry_order": 2,
            "rotation": rotation, "dofs": n,
            "boundary_U_dofs": iv["U"][1],
            "raw_per_rank_perimeter": raw_perimeters,
            "boundary_perimeter": perimeter,
            "perimeter_reduction_mode": perimeter_mode,
            "strong_wall_vector_mismatch_rms": strong_rms,
            "raw_per_rank_strong_wall_squared_integral": strong_contributions,
            "strong_squared_integral_reduction_mode": strong_mode,
            "moment_vectors": moments,
            "normal_tangential_decomposition_l2": decomposition_error,
            "boundary_mass_riesz": dual,
            "original_h0_p0_Y_preserved_bitwise": True,
            "entire_state_rollback_bitwise": True,
            "no_initialized_traction_Y0": True,
            "no_initialization_or_physical_TS": True,
            "algebraic_feasibility_proven": False,
            "interpretation": (
                "Direct GetFEM wall-mismatch test-function moments and "
                "Riesz mass dual norm at exactly original rotating q0. "
                "A nonzero weak moment needs corroboration by FE trace "
                "and algebraic elimination before declaring the original "
                "higher-index DAE initial state inconsistent."
            ),
        }
        if mpi.rank == 0:
            print(
                f"ISSUE28 WEAK WALL h={mesh:g} omega={rotation:g} "
                f"MPI={mpi.size} U_dofs={iv['U'][1]} "
                f"b_dual_RMS={dual['full']['riesz_boundary_RMS']:.10g} "
                f"normal={dual['normal']['riesz_boundary_RMS']:.10g} "
                f"tangent={dual['tangential']['riesz_boundary_RMS']:.10g}",
                flush=True
            )
        return result
    finally:
        s.gf_model.set_time(0.)
        s.gf_model.to_variables(original.copy())
        saved_global.copy(s.GlobalVec)
        s._explicit_rhs_ready = saved_ready
        s._explicit_rhs_time = saved_time
        if saved_rhs is not None and cached_rhs is not None:
            saved_rhs.copy(cached_rhs)
        restored = ensure_real(s.gf_model.from_variables(), "rollback state")
        if not mpi.allreduce(bool(np.array_equal(restored, original)),
                             op=MPI.LAND):
            raise RuntimeError("Wall moment diagnostic failed exact rollback")
        for obj in reversed(objects):
            obj.destroy()
        saved_global.destroy()
        if saved_rhs is not None:
            saved_rhs.destroy()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mesh", type=float, choices=(.2, .1), required=True)
    parser.add_argument("--case", choices=("rest", "rotating"), required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    mpi = PETSc.COMM_WORLD.tompi4py()
    rotation = 0. if args.case == "rest" else ROTATION
    system = build_system(
        gamma=0., mesh=args.mesh, geometry_order=2, height_order=2,
        final_time=.02, save_dt=.002, max_step=.002,
        ramp_start=.003, ramp_end=.009, init_step=True,
        ts_type="bdf", adapt_type="none", init_dt=1e-9,
        rotation=rotation, initial_profile="tangent",
        wall_velocity_profile="facet_tangent",
        algebraic_initial_guess="zero", init_method="bdf",
    )
    result = {}

    def intercept(self):
        result["report"] = audit(self, mpi, args.mesh, rotation)
        raise AuditComplete

    with patch.object(S.DPHS, "init_step", intercept):
        try:
            system.solve()
        except AuditComplete:
            pass
        else:
            raise RuntimeError("Unexpected physical time-stepping")
    if "report" not in result:
        raise RuntimeError("Original weak-wall audit was not called")
    if mpi.rank == 0:
        args.out.mkdir(parents=True, exist_ok=True)
        (args.out/"weak_wall.json").write_text(
            json.dumps(result["report"], indent=2)+"\n",
            encoding="utf-8"
        )
    mpi.Barrier()


if __name__ == "__main__":
    main()
