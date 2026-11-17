# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Fast native MPI smoke test: stream ONLY owned PETSc Vec segments.

The writer is started by the shell before mpiexec; no MPI process forks.
No physical PDE is simulated, no PETSc TS internals are changed.
"""
import argparse
import numpy as np
from mpi4py import MPI
from petsc4py import PETSc

from async_binary_trajectory import RankSender


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--socket", required=True)
    parser.add_argument("--n-global", type=int, default=17)
    args = parser.parse_args()

    comm = MPI.COMM_WORLD
    size, rank = comm.Get_size(), comm.Get_rank()
    vec = PETSc.Vec().createMPI(args.n_global, comm=PETSc.COMM_WORLD)
    lo, hi = vec.getOwnershipRange()
    sender = RankSender(
        args.socket, rank=rank, ranks=size, global_size=args.n_global,
        lo=lo, hi=hi, dtype=np.dtype(PETSc.ScalarType), max_queued=16,
    )
    try:
        for step in (0, 1, 4, 5):
            owned = vec.getArray()
            owned[:] = 100. * step + np.arange(
                lo, hi, dtype=PETSc.ScalarType
            )
            vec.assemble()
            sender.enqueue(
                step, 0.002 * step, vec.getArray(readonly=True)
            )
    finally:
        sender.close()
        vec.destroy()
    # Only a final barrier, outside any PETSc monitor/callback.
    comm.Barrier()
    if rank == 0:
        print(
            "PASS: native PETSc Vec local slices streamed from "
            f"{size} ranks without any all-gather",
            flush=True,
        )


if __name__ == "__main__":
    main()
