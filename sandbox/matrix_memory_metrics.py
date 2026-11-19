# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Per-rank memory/ownership snapshots for short matrix-policy tests.

Linux RSS measures the whole process, including PETSc/MUMPS/GetFEM,
and is NOT the memory of any one Mat. Sample immediately after TS.solve
and BEFORE trajectory/Hamiltonian test-only materialization.
"""
from pathlib import Path
import resource


def snapshot(system, rank):
    values = {}
    for line in Path("/proc/self/status").read_text(
        encoding="ascii"
    ).splitlines():
        name, _, value = line.partition(":")
        if name in (
            "VmRSS", "VmHWM", "RssAnon", "RssFile", "VmSwap"
        ):
            words = value.strip().split()
            if len(words) == 2 and words[1] == "kB":
                values[name + "_kib"] = int(words[0])
    if "VmRSS_kib" not in values or "VmHWM_kib" not in values:
        raise RuntimeError("Missing Linux process memory counters")
    return {
        "mpi_rank": int(rank),
        "matrix_memory_mode":
            getattr(system, "_matrix_memory_mode", "legacy"),
        "linear_mass": bool(system.linear_mass),
        "linear_stiffness": bool(system.linear_stiffness),
        "tangent_mass_is_mass":
            system.tangent_mass is system.mass,
        "tangent_stiffness_is_stiffness":
            system.tangent_stiffness is system.stiffness,
        "hessian_fd_matrix_allocations": int(getattr(
            system, "_mass_hessian_fd_matrix_allocations", 0
        )),
        "hessian_fd_matrix_reuses": int(getattr(
            system, "_mass_hessian_fd_matrix_reuses", 0
        )),
        "hessian_fd_probe_allocations": int(getattr(
            system, "_mass_hessian_fd_probe_allocations", 0
        )),
        "ru_maxrss_kib": int(resource.getrusage(
            resource.RUSAGE_SELF
        ).ru_maxrss),
        **values,
    }
