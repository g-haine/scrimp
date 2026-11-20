# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Check MPI and SCRIMP imports without running a physics simulation.

Run from the project root (inside Dev Container / activated native conda):
    mpiexec -n 2 python simulations/mpi_check.py
"""

import scrimp  # noqa: F401: exercise SCRIMP initialization on each rank
from petsc4py import PETSc

comm = PETSc.COMM_WORLD
print(f"SCRIMP/MPI OK: rank {comm.getRank() + 1}/{comm.getSize()}", flush=True)
comm.barrier()
if comm.getRank() == 0:
    print("PASS: MPI startup and SCRIMP import")
