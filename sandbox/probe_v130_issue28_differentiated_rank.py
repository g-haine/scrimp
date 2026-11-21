# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""#28: native first-differentiated-constraints rank diagnostic for Y0.

For the original open-loop viscous rotating SWE disk, write q=(h,p)
and a=(a_regular,Y), where Y is the original tangential wall traction.

The native algebraic rows have zero Y columns, while the dynamic
momentum rows couple to Y. Differentiating the algebraic equations
once leads to a *candidate* local bordered linear system:

  [ M_dd      0         J_dY  ] [qdot]     [dynamic RHS]
  [ J_ad    J_a,a_reg     0   ] [a_regdot] = [-partial_t F_a]
                                 [Y0]

with a matching 2x3 block layout. J blocks use the ACTUAL GetFEM
nonlinear tangent, M_dd the ACTUAL state-dependent mass matrix.

This pilot reports numerical invertibility using two synthetic
right-hand sides A*x_ref, not the physical RHS. It does NOT compute
a physical initial traction Y0 or determine the DAE index.

Two separate baselines:
 - exact stationary algebraic state (h0=5,p0=0), known native Fa=Fd=0;
 - unmodified ORIGINAL rotating h0,p0 at omega=0.2 and *zero*
   initial algebraic guess: its Fa is measured and may not vanish,
   so any linearized invertibility here is only STRUCTURAL.

No LSQR, hidden modification of q0, artificial forcing, pseudo
initializer or physical TS. Entire native GetFEM state is restored.
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
    build_system, H_INITIAL, RHO, GRAVITY, ROTATION,
)
from probe_v130_issue28_structure import (
    AuditComplete, indices, indexset, layout, matmax, norm_sub, rowmax,
)


REFERENCE_RESIDUAL_TOL = 1e-11
NUMERICAL_TEST_TOL = 1e-8


def audit(system, mpi, label):
    dnames, anames, iv, n = layout(system)
    if tuple(dnames) != ("h", "p") or "Y" not in anames:
        raise RuntimeError("Unexpected original open-loop SWE variables")

    d = np.sort(indices(dnames, iv))
    a = np.sort(indices(anames, iv))
    y = np.sort(indices(("Y",), iv))
    reg = np.setdiff1d(a, y, assume_unique=True)

    original = np.asarray(
        system.gf_model.from_variables(), dtype=float
    ).copy()
    if original.shape != (n,) or not np.isfinite(original).all():
        raise RuntimeError("Invalid initial GetFEM vector")
    fixed_q0 = original[d].copy()
    if label == "rest":
        trial = original.copy()
        for name, value in (
            ("e_h", float(RHO*GRAVITY*H_INITIAL)),
            ("Y_n", -float(RHO*GRAVITY*H_INITIAL)),
        ):
            lo_var, width = iv[name]
            trial[lo_var:lo_var+width] = value
        # The exact reference changes ONLY algebraic effort variables.
        if not np.array_equal(trial[d], fixed_q0):
            raise RuntimeError("Known resting equilibrium changed original q0")
        system.gf_model.set_time(0.)
        system.gf_model.to_variables(trial)
    elif label == "rotating":
        trial = original.copy()
    else:
        raise ValueError("Unknown original physical case")

    left, right = system.tangent_stiffness.getOwnershipRange()
    created = []

    def keep(obj):
        created.append(obj)
        return obj

    def set_from_full(vector, values):
        if right > left:
            vector.setValues(
                np.arange(left, right, dtype=PETSc.IntType),
                values[left:right],
            )
        vector.assemble()

    try:
        z = keep(system.tangent_stiffness.createVecRight())
        set_from_full(z, trial)
        zd = keep(z.duplicate())
        zd.set(0.)
        F = keep(z.duplicate())
        idof = keep(indexset(d, left, right, z.getComm()))
        ia = keep(indexset(a, left, right, z.getComm()))
        ir = keep(indexset(reg, left, right, z.getComm()))
        iy = keep(indexset(y, left, right, z.getComm()))

        system.IFunction(None, 0., z, zd, F)
        fa_norm = norm_sub(F, ia)
        fd_norm = norm_sub(F, idof)
        mass_alg_max = rowmax(system.tangent_mass, a, mpi)
        if mass_alg_max > 1e-12:
            raise RuntimeError(
                f"Original DAE unexpectedly has nonzero algebraic mass rows {mass_alg_max}"
            )
        if label == "rest" and (
            fa_norm > REFERENCE_RESIDUAL_TOL
            or fd_norm > REFERENCE_RESIDUAL_TOL
        ):
            raise RuntimeError(
                f"Known exact stationary solution failed: "
                f"Fa={fa_norm:.9g}, Fd={fd_norm:.9g}"
            )

        system._sync_lumped_trial_state(0., z)
        system.assemble_nl_stiffness()
        mass = system.tangent_mass
        jac = system.tangent_stiffness

        Mdd = keep(mass.createSubMatrix(idof, idof))
        JdY = keep(jac.createSubMatrix(idof, iy))
        Jad = keep(jac.createSubMatrix(ia, idof))
        Jar = keep(jac.createSubMatrix(ia, ir))
        JaY = keep(jac.createSubMatrix(ia, iy))
        m_ay = matmax(JaY, mpi)
        m_dy = matmax(JdY, mpi)
        if m_ay > 1e-12 or m_dy <= 1e-12:
            raise RuntimeError(
                "The original high-index traction null-column/dynamic coupling changed"
            )
        if (
            Mdd.getSize() != (len(d), len(d))
            or JdY.getSize() != (len(d), len(y))
            or Jad.getSize() != (len(a), len(d))
            or Jar.getSize() != (len(a), len(reg))
            or len(d)+len(a) != len(d)+len(reg)+len(y)
        ):
            raise RuntimeError("Unexpected differentiated block shape")

        # A genuine native PETSc sparse block operator; None denotes
        # zero blocks, not invented cross-couplings or forcing.
        bordered_nest = keep(PETSc.Mat().createNest(
            [[Mdd, None, JdY],
             [Jad, Jar, None]],
            comm=jac.getComm(),
        ))
        bordered_nest.assemble()
        bordered = keep(bordered_nest.convert("aij"))
        nrow, ncol = bordered.getSize()
        if (nrow, ncol) != (n, n):
            raise RuntimeError(
                f"Differentiated block operator size mismatch {nrow}x{ncol} vs {n}"
            )

        ksp = keep(PETSc.KSP().create(comm=bordered.getComm()))
        ksp.setOptionsPrefix("scrimp_issue28_index_rank_")
        ksp.setOperators(bordered)
        ksp.setType(PETSc.KSP.Type.PREONLY)
        ksp.getPC().setType(PETSc.PC.Type.LU)
        ksp.getPC().setFactorSolverType("mumps")
        random_trial_errors = []
        solver_reasons = []
        solver_failure = None
        manufactured_rhs = False
        try:
            for t in (0.011, 0.037):
                reference = keep(bordered.createVecRight())
                lo, hi = reference.getOwnershipRange()
                if hi > lo:
                    ix = np.arange(lo, hi, dtype=float)
                    reference.setValues(
                        np.arange(lo, hi, dtype=PETSc.IntType),
                        np.sin(1.+t*ix) + 0.2*np.cos(0.71+2*t*ix),
                    )
                reference.assemble()
                rhs = keep(bordered.createVecLeft())
                bordered.mult(reference, rhs)
                manufactured_rhs = True
                sol = keep(bordered.createVecRight())
                sol.set(0.)
                ksp.solve(rhs, sol)
                reason = int(ksp.getConvergedReason())
                solver_reasons.append(reason)
                delta = keep(sol.copy())
                delta.axpy(-1., reference)
                relative_x = float(delta.norm()/max(1., reference.norm()))
                reconstructed = keep(bordered.createVecLeft())
                bordered.mult(sol, reconstructed)
                reconstructed.axpy(-1., rhs)
                relative_r = float(reconstructed.norm()/max(1., rhs.norm()))
                random_trial_errors.append({
                    "relative_known_synthetic_solution_error": relative_x,
                    "relative_synthetic_linear_residual": relative_r,
                    "ksp_converged_reason": reason,
                })
        except (PETSc.Error, RuntimeError) as exc:
            solver_failure = f"{type(exc).__name__}: {exc}"

        checks = (
            solver_failure is None
            and len(random_trial_errors) == 2
            and all(
                row["ksp_converged_reason"] > 0
                and row["relative_known_synthetic_solution_error"] <= NUMERICAL_TEST_TOL
                and row["relative_synthetic_linear_residual"] <= NUMERICAL_TEST_TOL
                for row in random_trial_errors
            )
        )
        status = (
            "NUMERIC_REGULAR_AT_LINEARIZATION"
            if checks else "NOT_NUMERICALLY_CERTIFIED"
        )
        # Rank proxy is *not* an exact nullspace or DAE-index proof:
        # two synthetic RHS tests cannot replace singular value bounds.
        report = {
            "schema": "scrimp-v130-issue28-differentiated-rank-v1",
            "status": status,
            "initial_case": label,
            "mpi_ranks": int(mpi.size),
            "mesh": 0.2,
            "geometry_order": 2,
            "original_rotation": float(0.0 if label == "rest" else ROTATION),
            "dofs": int(n),
            "dynamic_dofs": int(len(d)),
            "algebraic_dofs": int(len(a)),
            "regular_algebraic_dofs": int(len(reg)),
            "traction_Y_dofs": int(len(y)),
            "differentiated_bordered_shape": [int(nrow), int(ncol)],
            "native_original_Fa_at_linearization_l2": fa_norm,
            "native_original_Fd_zdot_zero_at_linearization_l2": fd_norm,
            "mass_algebraic_row_max_abs": mass_alg_max,
            "J_aa_Y_column_max_abs": m_ay,
            "J_dynamic_Y_column_max_abs": m_dy,
            "synthetic_rhs_used_only_for_rank_diagnostic": manufactured_rhs,
            "synthetic_solve_tests": random_trial_errors,
            "synthetic_solve_failure": solver_failure,
            "original_q0_bitwise_preserved": True,
            "stationary_reference_exact_qualified": (
                label == "rest"
            ),
            "rotating_algebraic_compatibility_established": False,
            "physical_traction_Y0_calculated": False,
            "differentiation_index_proven": False,
            "real_initialization_TS_executed": False,
            "real_physical_TS_executed": False,
            "interpretation": (
                "The differentiated block operator consists of genuine "
                "GetFEM mass/tangent submatrices; numerical synthetic "
                "inversion is not a compatibility or index certificate. "
                "Rotating original q0 must satisfy Fa=0 independently."
            ),
        }
        if mpi.rank == 0:
            print(
                f"ISSUE28 DIFFERENTIATED RANK case={label} MPI={mpi.size}: "
                f"{status}, shape={nrow}x{ncol}, Y={len(y)}, "
                f"|Fa|={fa_norm:.9g}, synthetic={random_trial_errors}",
                flush=True,
            )
        return report
    finally:
        system.gf_model.set_time(0.)
        system.gf_model.to_variables(original.copy())
        restored = np.asarray(
            system.gf_model.from_variables(), dtype=float
        )
        for obj in reversed(created):
            obj.destroy()
        if not mpi.allreduce(
            bool(np.array_equal(restored, original)), op=MPI.LAND
        ):
            raise RuntimeError("Differentiated diagnostic did not roll back GetFEM")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", choices=("rest", "rotating"), required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    mpi = PETSc.COMM_WORLD.tompi4py()
    system = build_system(
        gamma=0., mesh=.2, geometry_order=2,
        height_order=2, final_time=.02,
        save_dt=.002, max_step=.002,
        ramp_start=.003, ramp_end=.009,
        init_step=True, ts_type="bdf", adapt_type="none",
        init_dt=1e-9,
        rotation=(0. if args.case == "rest" else ROTATION),
        initial_profile="tangent", wall_velocity_profile="facet_tangent",
        algebraic_initial_guess="zero", init_method="bdf",
    )
    report = {}

    def intercept(self):
        report["data"] = audit(self, mpi, args.case)
        raise AuditComplete

    with patch.object(S.DPHS, "init_step", intercept):
        try:
            system.solve()
        except AuditComplete:
            pass
        else:
            raise RuntimeError(
                "Missing native audit intercept or physical TS unexpectedly ran"
            )
    if "data" not in report:
        raise RuntimeError("No differentiated-block report generated")
    if mpi.rank == 0:
        args.out.mkdir(parents=True, exist_ok=True)
        (args.out/"differentiated_rank.json").write_text(
            json.dumps(report["data"], indent=2) + "\n",
            encoding="utf-8",
        )
    mpi.Barrier()


if __name__ == "__main__":
    main()
