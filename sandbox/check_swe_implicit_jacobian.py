# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Finite-difference audit of SCRIMP's full implicit SWE residual/Jacobian.

This is NOT a time simulation. No TS, SNES, KSP or MUMPS solver is created.

We check the PETSc contract
    J(z,zd,shift) @ v =
       d/d epsilon F(z+epsilon*v, zd+epsilon*shift*v)|epsilon=0
for the exact SCRIMP 2D nonlinear shallow-water model, at its initial
state and an optionally perturbed state. Thus F and J are tested in
the same GetFEM/MPI setting that exhibited divergence.

Examples:
    mpiexec -n 1 python -u sandbox/check_swe_implicit_jacobian.py
    mpiexec -n 2 python -u sandbox/check_swe_implicit_jacobian.py
    mpiexec -n 2 python -u sandbox/check_swe_implicit_jacobian.py \
        --state-perturb 0.05 --epsilon 1e-5

In addition to the directional error, the report checks:
- PETSc-to-GetFEM Newton-state synchronization on all MPI ranks;
- whether the GetFEM RHS directional derivative is identical on ranks;
- whether that derivative agrees with minus the nonlinear GetFEM
  tangent action (a separate check of GetFEM's RHS convention);
- which residual equation blocks account for the discrepancy.

The script reports the error and exits 0 by default, even if the
Jacobian is inconsistent: it is a diagnostic, not a green test.
Use --fail-above 1e-4 only when checking a candidate correction.
"""

import argparse
import hashlib
import json

from mpi4py import MPI

import numpy as np
from petsc4py import PETSc

from sandbox.viscous_swe_2d_mpi import build_system


def setup_model(mesh):
    # Reproduce DPHS.solve()'s initialization without creating a TS.
    swe = build_system(
        gamma=0.0, mesh=mesh, final_time=0.5, save_dt=0.1,
        max_step=0.01, init_step=False,
    )
    swe.set_linear_flags()
    swe.set_time_scheme()
    swe.time_scheme = PETSc.Options()
    swe.allocate_memory()
    swe.disable_all_bricks()
    swe.gf_model.set_time(0.0)
    swe.gf_model.first_iter()
    swe.assemble_mass()
    swe.assemble_nl_mass()
    swe.assemble_stiffness()
    swe.assemble_nl_stiffness()
    swe.assemble_rhs()
    swe.J = swe.tangent_stiffness.duplicate()
    swe.J.assemble()
    swe.F = swe.rhs.duplicate()
    swe.F.assemble()
    swe.buffer = swe.rhs.duplicate()
    swe.buffer.assemble()
    return swe


def from_global_array(template, values):
    result = template.duplicate()
    lo, hi = result.getOwnershipRange()
    result.getArray()[:] = values[lo:hi]
    return result


def report(swe, z, zd, direction, shift, epsilon):
    J = swe.J
    swe.IJacobian(None, 0.0, z, zd, shift, J, J)

    # Check the distributed PETSc Vec -> replicated GetFEM trial-state
    # scatter: errors here corrupt both nonlinear RHS and tangent assembly.
    mpi = MPI.COMM_WORLD
    global_trial = np.concatenate(
        mpi.allgather(z.getArray(readonly=True).real.copy())
    )
    getfem_trial = np.asarray(swe.gf_model.from_variables()).real.ravel()
    state_scatter_error = float(np.linalg.norm(getfem_trial - global_trial))
    state_scatter_errors = mpi.allgather(state_scatter_error)

    jv = direction.duplicate()
    J.mult(direction, jv)
    # Snapshot the *true* nonlinear GetFEM tangent action at the base
    # state before subsequent IFunction calls mutate the matrix.
    nl_action = direction.duplicate()
    swe.nl_stiffness.mult(direction, nl_action)

    zplus, zminus = z.copy(), z.copy()
    zdplus, zdminus = zd.copy(), zd.copy()
    fplus, fminus = swe.F.duplicate(), swe.F.duplicate()
    try:
        zplus.axpy(epsilon, direction)
        zminus.axpy(-epsilon, direction)
        zdplus.axpy(epsilon * shift, direction)
        zdminus.axpy(-epsilon * shift, direction)

        swe.IFunction(None, 0.0, zplus, zdplus, fplus)
        gf_rhs_plus = np.asarray(swe.gf_model.rhs()).copy()
        swe.IFunction(None, 0.0, zminus, zdminus, fminus)
        gf_rhs_minus = np.asarray(swe.gf_model.rhs()).copy()

        # GetFEM MPI-2 documents global, summed RHS vectors (unlike
        # distributed tangent matrices). Verify that claim at two points
        # in this *real* coupled problem, independently from PETSc's F.
        gf_rhs_fd = (gf_rhs_plus - gf_rhs_minus) / (2.0 * epsilon)
        mpi = MPI.COMM_WORLD
        rhs_by_rank = mpi.allgather(gf_rhs_fd)
        rhs_copy_disagreement = max(
            float(np.linalg.norm(rhs_k - rhs_by_rank[0]))
            for rhs_k in rhs_by_rank
        )
        if gf_rhs_fd.size != z.getSize():
            raise ValueError("GetFEM RHS dimension differs from PETSc")
        lo, hi = nl_action.getOwnershipRange()
        nl_local = np.asarray(nl_action.getArray(readonly=True))
        rhs_local = gf_rhs_fd[lo:hi]
        rhs_vs_tangent_sq = float(
            np.vdot(rhs_local + nl_local, rhs_local + nl_local).real
        )
        nl_sq = float(np.vdot(nl_local, nl_local).real)
        gf_rhs_derivative_vs_tangent = (
            mpi.allreduce(rhs_vs_tangent_sq, op=MPI.SUM) ** 0.5
            / max(mpi.allreduce(nl_sq, op=MPI.SUM) ** 0.5, 1e-14)
        )
        rhs_signatures = mpi.allgather({
            "rank": mpi.Get_rank(),
            "rhs_plus_norm": float(np.linalg.norm(gf_rhs_plus)),
            "rhs_minus_norm": float(np.linalg.norm(gf_rhs_minus)),
            "rhs_fd_norm": float(np.linalg.norm(gf_rhs_fd)),
            "rhs_fd_digest": hashlib.blake2b(
                np.ascontiguousarray(gf_rhs_fd).view(np.uint8),
                digest_size=8,
            ).hexdigest(),
        })

        fplus.axpy(-1.0, fminus)
        fplus.scale(0.5 / epsilon)
        difference = fplus.copy()
        difference.axpy(-1.0, jv)
        try:
            norm_fd = float(fplus.norm())
            norm_jv = float(jv.norm())
            norm_error = float(difference.norm())
            # Diagnose which variable equation blocks contribute to the
            # mismatch, using the *global* GetFEM variable numbering.
            groups = {}
            dlocal = difference.getArray(readonly=True)
            for name in ("h", "p", "e_h", "e_p",
                         "U_n", "Y_n", "U", "Y"):
                start, length = swe.gf_model.interval_of_variable(name)
                start, end = int(start), int(start + length)
                a, b = max(start, lo), min(end, hi)
                if b > a:
                    own = dlocal[a-lo:b-lo]
                    local_squared = float(np.vdot(own, own).real)
                else:
                    local_squared = 0.0
                groups[name] = (
                    mpi.allreduce(local_squared, op=MPI.SUM) ** 0.5
                )
        finally:
            difference.destroy()
        return {
            "zd_norm": float(zd.norm()),
            "shift": float(shift),
            "epsilon": float(epsilon),
            "finite_difference_norm": norm_fd,
            "jacobian_action_norm": norm_jv,
            "error_norm": norm_error,
            "relative_error": norm_error / max(norm_fd, norm_jv, 1e-14),
            "error_norm_by_equation": groups,
            "getfem_rhs_derivative_rank_difference": rhs_copy_disagreement,
            "getfem_rhs_derivative_vs_minus_nonlinear_tangent_relative":
                gf_rhs_derivative_vs_tangent,
            "getfem_rhs_rank_signatures": rhs_signatures,
            "petsc_to_getfem_state_scatter_error_by_rank": state_scatter_errors,
        }
    finally:
        for obj in (jv, nl_action, zplus, zminus, zdplus, zdminus, fplus, fminus):
            obj.destroy()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mesh", type=float, default=0.2)
    parser.add_argument("--epsilon", type=float, default=1e-5)
    parser.add_argument("--state-perturb", type=float, default=0.0,
                        help="Relative deterministic random state perturbation")
    parser.add_argument(
        "--residual-form", choices=("legacy", "direct_getfem"),
        default="legacy",
        help="DIAGNOSTIC ONLY: evaluate residual from GetFEM weak-form "
             "RHS and linear stiffness instead of K_tangent(z)*z+rhs. "
             "The production DPHS.IFunction is NOT modified.",
    )
    parser.add_argument(
        "--rhs-scale", type=float, default=1.0,
        help="DIAGNOSTIC ONLY: multiply PETSc RHS after GetFEM assembly; "
             "use 0.5 at 2 ranks to test whether the observed factor-two "
             "RHS bias fully explains the MPI finite-difference mismatch",
    )
    parser.add_argument("--fail-above", type=float, default=None,
                        help="Optional maximum relative directional FD error")
    args = parser.parse_args()
    if (args.mesh <= 0 or args.epsilon <= 0
            or args.state_perturb < 0 or args.rhs_scale <= 0):
        parser.error("mesh, epsilon, RHS scale positive; perturb >= 0")

    comm = PETSc.COMM_WORLD
    rank = comm.getRank()
    swe = setup_model(args.mesh)
    if args.residual_form == "direct_getfem":
        # Avoid touching the package while testing a suspected
        # mathematical error in the state-dependent-mass IFunction
        # branch. For a GetFEM nonlinear residual N(z), build_rhs
        # returns -N(z); its tangent matrix is N'(z), not N(z).
        # Thus the consistent implicit function is
        # M(z)*zd + K_linear*z - rhs_nonlinear(z).
        # The current IJacobian may still omit dM(z)/dz * zd.
        def direct_weak_form_residual(TS, t, z, zd, F):
            swe._sync_lumped_trial_state(t, z)
            if not swe.linear_mass:
                swe.assemble_nl_mass()
            swe.assemble_rhs()
            swe.tangent_mass.mult(zd, swe.buffer)
            swe.stiffness.multAdd(z, swe.buffer, F)
            F.axpy(-1, swe.rhs)

        swe.IFunction = direct_weak_form_residual

    if args.rhs_scale != 1.0:
        # This wrapper changes ONLY the isolated FD experiment, never
        # the SCRIMP package or the user's simulation. The raw GetFEM
        # RHS diagnostics in report() remain unscaled.
        real_assemble_rhs = swe.assemble_rhs

        def diagnostic_rhs():
            real_assemble_rhs()
            swe.rhs.scale(args.rhs_scale)

        swe.assemble_rhs = diagnostic_rhs

    global_z = np.asarray(swe.gf_model.from_variables()).real.ravel()
    n = len(global_z)

    rng = np.random.default_rng(513)
    direction_values = rng.standard_normal(n)
    direction_values /= np.linalg.norm(direction_values)

    if args.state_perturb:
        global_z = (
            global_z +
            args.state_perturb * max(np.linalg.norm(global_z), 1.0)
            * direction_values
        )

    vec_template = swe.J.createVecRight()
    z = from_global_array(vec_template, global_z)
    direction = from_global_array(vec_template, direction_values)
    nonzero_zd = from_global_array(
        vec_template, 0.1 * max(np.linalg.norm(global_z), 1.0)
        * direction_values,
    )
    zero_zd = vec_template.duplicate()
    zero_zd.set(0.0)

    try:
        results = []
        for zd, shift in (
            (zero_zd, 0.0),
            (zero_zd, 1.0),
            (nonzero_zd, 0.0),
            (nonzero_zd, 1.0),
        ):
            results.append(report(swe, z, zd, direction, shift, args.epsilon))

        # PETSc Vec.norm() is collective: ALL ranks must call it,
        # never only rank 0 in a branch used for reporting.
        state_norm = float(z.norm())
        if rank == 0:
            print(json.dumps({
                "dofs": n,
                "ranks": comm.getSize(),
                "state_norm": state_norm,
                "diagnostic_rhs_scale": args.rhs_scale,
                "diagnostic_residual_form": args.residual_form,
                "results": results,
            }, indent=2), flush=True)

        # No threshold by default: record failures as data, not as a
        # false-positive success or a deliberately red regression test.
        if args.fail_above is not None and any(
            result["relative_error"] > args.fail_above
            for result in results
        ):
            raise RuntimeError(
                "Implicit Jacobian does not match finite differences: "
                f"error exceeded {args.fail_above}"
            )
    finally:
        for obj in (z, direction, nonzero_zd, zero_zd, vec_template):
            obj.destroy()


if __name__ == "__main__":
    main()
