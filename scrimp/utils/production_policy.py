# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp


"""Production switches for conforming interfaces and MPI ParaView.

The two historical v1.3.0 pilot environment variables are retained
as backwards-compatible overrides for pinned native qualification
scripts; they are not required for ordinary production use.

These functions are importable without GetFEM, PETSc or MPI.
"""

import os


def _binary(name, value):
    if value not in ("0", "1"):
        raise ValueError(f"{name} must be '0' or '1', got {value!r}")
    return value == "1"


def auto_interface_enabled(environ=None):
    """Conforming internal interfaces are automatically compiled by default.

    No internal interface -> existing single-domain/multi-mesh model
    registration remains unchanged. Ambiguous or unsupported detected
    interfaces must fail closed (handled by the compiler, not here).

    SCRIMP_AUTO_INTERFACE=0 disables the production feature. The
    historical SCRIMP_V130_AUTO_INTERFACE_PILOT flag is a lower-priority
    backwards-compatible alias.
    """
    env = os.environ if environ is None else environ
    if "SCRIMP_AUTO_INTERFACE" in env:
        return _binary("SCRIMP_AUTO_INTERFACE", env["SCRIMP_AUTO_INTERFACE"])
    if "SCRIMP_V130_AUTO_INTERFACE_PILOT" in env:
        return _binary(
            "SCRIMP_V130_AUTO_INTERFACE_PILOT",
            env["SCRIMP_V130_AUTO_INTERFACE_PILOT"],
        )
    return True


def pv_export_mode(environ=None, mpi_size=1):
    """Return ('legacy' or 'root') for the current PETSc MPI communicator.

    Default MPI>=2 uses one writer instead of simultaneous writes to
    the same VTK/PVD files. Before root writes, the DPHS export method
    must independently check that every rank owns the same complete FE
    history; a partitioned state is an explicit unsupported case.

    SCRIMP_PV_EXPORT_MODE=legacy opts out and preserves historical
    multi-writer behavior. 'root' explicitly selects the one-writer
    policy. 'auto' is the default. The historical pilot flag remains
    a lower-priority alias for old Docker qualification scripts.
    """
    env = os.environ if environ is None else environ
    mode = env.get("SCRIMP_PV_EXPORT_MODE")
    if mode is None:
        old = env.get("SCRIMP_V130_PV_ROOT_ONLY_PILOT")
        if old is not None:
            mode = "root" if _binary(
                "SCRIMP_V130_PV_ROOT_ONLY_PILOT", old
            ) else "legacy"
        else:
            mode = "auto"
    if mode not in ("auto", "root", "legacy"):
        raise ValueError(
            "SCRIMP_PV_EXPORT_MODE must be 'auto', 'root' or 'legacy'"
        )
    return "root" if mode == "root" or (
        mode == "auto" and int(mpi_size) > 1
    ) else "legacy"
