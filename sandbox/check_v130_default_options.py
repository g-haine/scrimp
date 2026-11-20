# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Native PETSc.Options semantics for algebraic default without any solve."""

import scrimp as S


def main():
    # No initializer is requested; the implicit default must NOT make
    # previously valid init_step=False applications fail.
    noinit = S.DPHS("real")
    noinit.set_time_scheme(init_step=False)
    # PETSc.Options.hasName() means the option exists in the
    # database, not that the caller explicitly supplied it.
    # SCRIMP deliberately materializes its effective default so the
    # eventual _solve_impl() re-entry uses the same strategy.
    if noinit.time_scheme.getString(
        "init_step_ts_type", ""
    ).strip().lower() != "algebraic":
        raise RuntimeError(
            "No-init model did not register effective algebraic default"
        )
    if noinit.time_scheme.getString(
        "init_step", ""
    ).strip().lower() != "false":
        raise RuntimeError("init_step=False was not preserved")
    noinit.set_time_scheme()
    if noinit.time_scheme.getString(
        "init_step_ts_type", ""
    ).strip().lower() != "algebraic":
        raise RuntimeError(
            "Internal no-kwargs set_time_scheme changed the default"
        )
    print("PASS: default registered as algebraic; init_step=False "
          "remains legal through repeated set_time_scheme")

    # In contrast, explicitly requesting the projection while forbidding
    # initialization is contradictory and must raise an actionable error.
    invalid = S.DPHS("real")
    try:
        invalid.set_time_scheme(
            init_step=False, init_step_ts_type="algebraic"
        )
    except ValueError as exc:
        if "requires init_step=true" not in str(exc):
            raise
    else:
        raise RuntimeError(
            "Explicit algebraic + init_step=False unexpectedly accepted"
        )
    print("PASS: explicit conflicting algebraic + init_step=False rejected")

    # Preserve the original relaxation opt-in API.
    legacy = S.DPHS("real")
    legacy.set_time_scheme(
        init_step=True, init_step_ts_type="pseudo"
    )
    if legacy.time_scheme.getString(
        "init_step_ts_type", ""
    ) != "pseudo":
        raise RuntimeError("Legacy pseudo option was changed")
    print("PASS: explicit pseudo-time relaxation is still selectable")


if __name__ == "__main__":
    main()
