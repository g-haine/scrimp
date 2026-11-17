# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Opt-in PETSc matrix-memory helpers (no PETSc import at module load).

The DPHS owner retains every PETSc Vec/Mat across TS Newton callbacks.
Never destroy a matrix referenced by PETSc TS/KSP/PC; these helpers
only manage SCRIMP-owned tangent and Hessian scratch matrices.
"""


def alias_readonly_linear_tangent(linear_matrix, uninitialized_tangent):
    """Use the immutable assembled linear matrix as its own tangent.

    Only call this when there are NO state-dependent bricks for that
    operator, before PETSc TS creates its Jacobian or KSP. The tangent
    placeholder must not yet be initialized or referenced by PETSc.
    The independent TS Jacobian remains a separate Mat.
    """
    if linear_matrix is uninitialized_tangent:
        return linear_matrix
    uninitialized_tangent.destroy()
    return linear_matrix


def add_mass_hessian_fd_reused(owner, t, z, zd, jacobian, eta, structure):
    """Append the exact same central-FD Hessian using ONE matrix workspace.

    The first perturbation M(z+eta*zd) is copied into a persistent PETSc
    Mat. The second perturbation M(z-eta*zd) remains in the system's
    existing tangent_mass matrix. There is no need to duplicate it or
    allocate a third correction matrix.

    One persistent PETSc Vec is also enough for both perturbations,
    because z.copy(result=probe) restores the original state before
    each AXPY, avoiding accumulated perturbation roundoff.

    The GetFEM/PETSc base Newton trial state and tangent mass are
    restored unconditionally, matching the historical code.
    """
    probe = getattr(owner, "_mass_hessian_fd_probe", None)
    if probe is None:
        probe = z.duplicate()
        owner._mass_hessian_fd_probe = probe
        owner._mass_hessian_fd_probe_allocations = (
            getattr(owner, "_mass_hessian_fd_probe_allocations", 0) + 1
        )
    work = getattr(owner, "_mass_hessian_fd_matrix", None)
    try:
        z.copy(probe)
        probe.axpy(eta, zd)
        owner._sync_lumped_trial_state(t, probe)
        owner.assemble_nl_mass()

        if work is None:
            work = owner.tangent_mass.copy()
            owner._mass_hessian_fd_matrix = work
            owner._mass_hessian_fd_matrix_allocations = (
                getattr(owner, "_mass_hessian_fd_matrix_allocations", 0) + 1
            )
        else:
            # PETSc MatCopy(... DIFFERENT_NONZERO_PATTERN) accepts a
            # changing nonlinear tangent graph; do not assume SAME.
            owner.tangent_mass.copy(work, structure=structure)
            owner._mass_hessian_fd_matrix_reuses = (
                getattr(owner, "_mass_hessian_fd_matrix_reuses", 0) + 1
            )

        z.copy(probe)
        probe.axpy(-eta, zd)
        owner._sync_lumped_trial_state(t, probe)
        owner.assemble_nl_mass()

        # Overwrite the saved + perturbation with the correction.
        # The - perturbation is read directly from tangent_mass.
        work.axpy(-1.0, owner.tangent_mass, structure=structure)
        work.scale(0.5 / eta)
        jacobian.axpy(1.0, work, structure=structure)
        jacobian.assemble()
    finally:
        owner._sync_lumped_trial_state(t, z)
        owner.assemble_nl_mass()
