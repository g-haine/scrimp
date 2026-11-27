# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp


"""Native, opt-in 2D-wave qualification of SuperLU_DIST internal refinement.

Run the UNMODIFIED canonical examples.wave.main(), with the original
anisotropic mixed-boundary model, strict algebraic DAE initialisation,
10 fixed CN steps, four boundary powers and four genuine GetFEM
ParaView exports. Only change two PETSc factor-solver options:

    mat_superlu_dist_replacetinypivot = true
    mat_superlu_dist_iterrefine = true

These settings recovered the known 2230x2230 first-CN matrix on 1/2 MPI
ranks with original Px-b backward residual around 1e-16; this fact
does not yet qualify a full time-domain pH-DAE simulation. The driver
is never enabled by default and does not change SCRIMP or wave.py.

SCRIMP clears PETSc options at DPHS construction time. Thus install
the opt-in options AFTER the canonical DPHS.set_time_scheme, BEFORE
its TS/SNES/KSP setFromOptions. No shell-level PETSC_OPTIONS shortcut.
"""

import os

import scrimp as S


_OPTIONS = (
    "mat_superlu_dist_replacetinypivot",
    "mat_superlu_dist_iterrefine",
)


def main():
    if os.environ.get("SCRIMP_V130_WAVE2D_ENABLE_REFINED_DIST") != "1":
        raise SystemExit(
            "Research-only solver override requires explicit "
            "SCRIMP_V130_WAVE2D_ENABLE_REFINED_DIST=1"
        )
    original = S.DPHS.set_time_scheme
    invoked = [0]

    def enable_refined_superlu_dist(system, **kwargs):
        # The native canonical model still chooses its solver, mass,
        # initial DAE, CN/TS scheme and GMRES/PC LU just as before.
        result = original(system, **kwargs)
        actual_solver = system.time_scheme.getString(
            "pc_factor_mat_solver_type", ""
        ).strip().lower()
        if actual_solver != "superlu_dist":
            raise RuntimeError(
                "Refinement candidate MUST use real SuperLU_DIST, got "
                + repr(actual_solver)
            )
        if system.time_scheme.getString(
            "pc_type", ""
        ).strip().lower() != "lu":
            raise RuntimeError("Canonical candidate must preserve PETSc LU")
        for option in _OPTIONS:
            system.time_scheme[option] = "true"
            if not system.time_scheme.getBool(option, False):
                raise RuntimeError(
                    f"PETSc did not register required option {option!r}"
                )
        # DPHS._solve_impl() may call set_time_scheme() again with
        # zero kwargs; the numerical options must survive that call.
        # Recheck the two options each time, but announce the single
        # physical-system opt-in only once per DPHS instance/rank.
        if not getattr(system, "_v130_refined_dist_enabled", False):
            system._v130_refined_dist_enabled = True
            invoked[0] += 1
            print(
                "V130_WAVE2D_REFINED_SUPERLU_DIST_ENABLED "
                "solver=superlu_dist "
                "mat_superlu_dist_replacetinypivot=true "
                "mat_superlu_dist_iterrefine=true "
                f"instance={invoked[0]}",
                flush=True,
            )
        return result

    S.DPHS.set_time_scheme = enable_refined_superlu_dist
    try:
        from examples.wave import main as canonical_wave_main

        canonical_wave_main()
        if invoked[0] != 1:
            raise RuntimeError(
                "Refined solver options must be injected precisely once, "
                f"got {invoked[0]} calls"
            )
    finally:
        S.DPHS.set_time_scheme = original


if __name__ == "__main__":
    main()
