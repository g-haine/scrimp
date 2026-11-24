# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# See LICENSE for license information.
#
# github: https://github.com/g-haine/scrimp

"""Diagnostics of genuine SuperLU_DIST ordering on the unchanged 2D wave.

The PETSc native diagnostic identified KSP_DIVERGED_PC_FAILED=-11,
PC_FACTOR_NUMERIC_ZEROPIVOT=2, before any accepted TS step.

ONLY test non-regularizing *permutations* available in SuperLU_DIST,
without changing the wave equation, FE variables, DAE initialization,
time step, Jacobian, singularity threshold, or nonzero numerical
coefficients. This wrapper is research-only and is not the example
shown in Sphinx.

SCRIMP currently clears the PETSc options database on DPHS creation.
Consequently setting a PETSC_OPTIONS environment variable before
constructing DPHS would NOT be a reliable test of these options.
Instead we inject the selected genuine PETSc options **after**
the existing DPHS.set_time_scheme method and before TS.setFromOptions.
The public examples/wave.py file remains byte-for-byte unchanged.

The driver reuses examples.wave.main(), so the FULL 10-step original
Hamiltonian/physical four-port gates and true LU identity apply.
"""

import os

import scrimp as S


ORDERING_CASES = {
    "rowperm_awpm": {
        "mat_superlu_dist_rowperm": "LargeDiag_AWPM",
    },
    "colperm_mmd": {
        "mat_superlu_dist_colperm": "MMD_AT_PLUS_A",
    },
}


def main():
    selection = os.environ.get("SCRIMP_V130_SUPERLU_DIST_ORDERING", "")
    if selection not in ORDERING_CASES:
        raise SystemExit(
            "Set SCRIMP_V130_SUPERLU_DIST_ORDERING to "
            + " or ".join(sorted(ORDERING_CASES))
        )
    settings = ORDERING_CASES[selection]
    original = S.DPHS.set_time_scheme

    def with_ordering(system, **kwargs):
        result = original(system, **kwargs)
        actual_solver = system.time_scheme.getString(
            "pc_factor_mat_solver_type", ""
        ).strip().lower()
        if actual_solver != "superlu_dist":
            raise RuntimeError(
                "Ordering probe must use REAL SuperLU_DIST; found "
                + repr(actual_solver)
            )
        for key, value in settings.items():
            system.time_scheme[key] = value
            observed = system.time_scheme.getString(key, "")
            if observed != value:
                raise RuntimeError(
                    f"PETSc ordering option not installed: {key}={observed!r}"
                )
        if not getattr(system, "_wave2d_ordering_announced", False):
            print(
                "V130_SUPERLU_DIST_ORDERING_APPLIED "
                + selection + " " + repr(settings),
                flush=True,
            )
            system._wave2d_ordering_announced = True
        return result

    S.DPHS.set_time_scheme = with_ordering
    try:
        from examples.wave import main as original_wave_main

        original_wave_main()
    finally:
        S.DPHS.set_time_scheme = original


if __name__ == "__main__":
    main()
