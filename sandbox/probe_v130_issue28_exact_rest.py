# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Issue #28: independent exact-rest control for nonlinear SWE geometry audit.

This tests the ORIGINAL, unchanged viscous weak forms, with h0=5, p0=0
and the physically prescribed wall at rest. Unlike the previous LSQR
feasibility test it starts from the exact known algebraic equilibrium:
  e_h=rho*g*h=0.5, e_p=0, Y_n=-e_h=-0.5,
  U_n=U=Y=0.
Negative Y_n follows SCRIMP's actual flow-brick sign convention:
  M_n[h] is registered as -h*Y_n*Test_Y_n, and
  -B_n[h]^T as -h*e_h*Test_Y_n.
All states, including dynamically prescribed h0,p0, are then restored.
No initialization TS, physical TS, or acceptance is performed.
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
    build_system, RHO, GRAVITY, H_INITIAL,
)
from probe_v130_issue28_structure import (
    AuditComplete, indices, indexset, layout, norm_sub,
)

TOL = 1e-11


def _norm_block(vec, variable, iv, lo, hi):
    iset = indexset(indices((variable,), iv), lo, hi, vec.getComm())
    try:
        return norm_sub(vec, iset)
    finally:
        iset.destroy()


def audit(system, mpi, order):
    dynamic, algebraic, iv, n = layout(system)
    if tuple(dynamic) != ("h", "p") or "Y_n" not in algebraic:
        raise ValueError("Unexpected original SWE structure")
    original = np.asarray(system.gf_model.from_variables(), dtype=float).copy()
    if original.shape != (n,) or not np.isfinite(original).all():
        raise ValueError("Invalid original GetFEM state")
    lo, hi = system.tangent_stiffness.getOwnershipRange()
    # The original h and p MUST be unchanged. The established
    # model has h0=5 and p0=0, as prescribed by the paper sandbox.
    h0 = indices(("h",), iv)
    p0 = indices(("p",), iv)
    if not np.allclose(original[h0], H_INITIAL, atol=1e-13, rtol=0):
        raise RuntimeError("Actual h0 differs from prescribed rest height")
    if not np.array_equal(original[p0], np.zeros(len(p0))):
        raise RuntimeError("Rest baseline p0 is not exactly zero")
    created = []
    def keep(obj):
        created.append(obj)
        return obj
    try:
        original_vec = keep(system.tangent_stiffness.createVecRight())
        owned_all = np.arange(lo, hi, dtype=PETSc.IntType)
        if len(owned_all):
            original_vec.setValues(owned_all, original[lo:hi])
        original_vec.assemble()
        reference = keep(original_vec.copy())
        fixed = indices(dynamic, iv)
        expected_pressure = float(RHO * GRAVITY * H_INITIAL)
        # Include a signed *negative* Y_n (dual normal-pressure port).
        # No PDE controls or prescribed state are changed.
        for name, value in (("e_h", expected_pressure),
                            ("Y_n", -expected_pressure)):
            var_indices = indices((name,), iv)
            owned = var_indices[(var_indices >= lo) & (var_indices < hi)]
            if len(owned):
                reference.setValues(
                    np.asarray(owned, dtype=PETSc.IntType),
                    np.full(len(owned), value, dtype=float),
                    addv=PETSc.InsertMode.INSERT_VALUES,
                )
        reference.assemble()
        deriv = keep(reference.duplicate())
        deriv.set(0)
        residual = keep(reference.duplicate())
        rows_d = keep(indexset(indices(dynamic, iv), lo, hi, reference.getComm()))
        rows_a = keep(indexset(indices(algebraic, iv), lo, hi, reference.getComm()))
        system.IFunction(None, 0., reference, deriv, residual)
        fa = norm_sub(residual, rows_a)
        fd = norm_sub(residual, rows_d)
        blocks = {
            name: _norm_block(residual, name, iv, lo, hi)
            for name in tuple(dynamic) + tuple(algebraic)
        }
        # A negative control protects the surprising minus sign
        # of Y_n. The wrong positive port sign must be detectable.
        wrong_sign = keep(reference.copy())
        yn = indices(("Y_n",), iv)
        owned_yn = yn[(yn >= lo) & (yn < hi)]
        if len(owned_yn):
            wrong_sign.setValues(
                np.asarray(owned_yn, dtype=PETSc.IntType),
                np.full(len(owned_yn), expected_pressure),
            )
        wrong_sign.assemble()
        system.IFunction(None, 0., wrong_sign, deriv, residual)
        wrong_sign_fa = norm_sub(residual, rows_a)
        system._sync_lumped_trial_state(0., reference)
        candidate = np.asarray(system.gf_model.from_variables(), dtype=float)
        same_d = bool(np.array_equal(candidate[fixed], original[fixed]))
        same_rest = bool(np.array_equal(candidate[p0], original[p0]))
        if not mpi.allreduce(same_d and same_rest, op=MPI.LAND):
            raise RuntimeError("Reference candidate modified prescribed h0,p0")
        if not all(np.isfinite(v) for v in (fa, fd, wrong_sign_fa)):
            raise FloatingPointError("Nonfinite reference residual")
        if fa > TOL or fd > TOL:
            raise AssertionError(
                f"Rest reference not a native equilibrium (Fa={fa:.9g}, Fd={fd:.9g})"
            )
        if wrong_sign_fa < 1e-4:
            raise AssertionError(
                "Negative control failed: wrong Y_n sign must create algebraic residual"
            )
        result = {
            "schema": "scrimp-v130-issue28-exact-rest-v1",
            "status": "EXACT_REST_NATIVE_PASS",
            "geometry_order": order,
            "mpi_ranks": int(mpi.size),
            "dofs": int(n),
            "h0": float(H_INITIAL),
            "p0": "identically zero",
            "e_h_reference": expected_pressure,
            "Y_n_reference": -expected_pressure,
            "F_a_reference_l2": float(fa),
            "F_d_reference_l2": float(fd),
            "F_reference_block_l2": blocks,
            "F_a_wrong_Yn_sign_l2": float(wrong_sign_fa),
            "tolerance_absolute": TOL,
            "dynamic_h0_p0_bitwise_preserved": True,
            "real_IFunction_evaluated": True,
            "least_squares_used": False,
            "initializer_TS_executed": False,
            "main_TS_executed": False,
            "interpretation": (
                "Known exact stationary solution calibrates a "
                "false LSQR incompatibility label; no rotating "
                "DAE initialization is qualified."
            ),
        }
        if mpi.rank == 0:
            print(
                f"ISSUE28 EXACT REST PASS P{order} MPI={mpi.size}: "
                f"|Fa|={fa:.9g}, |Fd|={fd:.9g}, "
                f"wrong-sign |Fa|={wrong_sign_fa:.9g}",
                flush=True,
            )
        return result
    finally:
        system.gf_model.set_time(0.)
        system.gf_model.to_variables(original.copy())
        if not mpi.allreduce(
            np.array_equal(
                np.asarray(system.gf_model.from_variables(), dtype=float),
                original
            ), op=MPI.LAND
        ):
            raise RuntimeError("Exact-rest diagnostic failed complete state rollback")
        for obj in reversed(created):
            obj.destroy()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--geometry-order", type=int, choices=(1, 2), required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    mpi = PETSc.COMM_WORLD.tompi4py()
    system = build_system(
        gamma=0., mesh=.2, geometry_order=args.geometry_order,
        final_time=.02, save_dt=.002, max_step=.002,
        ramp_start=.003, ramp_end=.009,
        init_step=True, ts_type="bdf", adapt_type="none",
        init_dt=1e-9, rotation=0.,
        initial_profile="tangent", wall_velocity_profile="facet_tangent",
        algebraic_initial_guess="zero", init_method="bdf", height_order=2,
    )
    result = {}
    def intercept(self):
        result["reference"] = audit(self, mpi, args.geometry_order)
        raise AuditComplete
    with patch.object(S.DPHS, "init_step", intercept):
        try:
            system.solve()
        except AuditComplete:
            pass
        else:
            raise RuntimeError("Unexpected TS initialization/physical solve")
    if "reference" not in result:
        raise RuntimeError("No exact stationary reference produced")
    if mpi.rank == 0:
        args.out.mkdir(parents=True, exist_ok=True)
        (args.out / "exact_rest.json").write_text(
            json.dumps(result["reference"], indent=2) + "\n",
            encoding="utf-8",
        )
    mpi.Barrier()


if __name__ == "__main__":
    main()
