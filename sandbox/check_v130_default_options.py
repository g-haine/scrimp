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
    if noinit.time_scheme.hasName("init_step_ts_type"):
        raise RuntimeError(
            "No-init model unexpectedly acquired an explicit strategy"
        )
    if noinit.time_scheme.getString("init_step", "") != "false":
        raise RuntimeError("init_step=False was not preserved")
    print("PASS: implicit default + init_step=False is allowed")

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
