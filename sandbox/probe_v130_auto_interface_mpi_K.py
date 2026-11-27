# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp


"""MPI1/MPI2: capture truly owned PETSc rows of automatic Heat-Wave K.

Use the current outward-normal example in BOTH bulk permutations.
This opt-in pilot captures actual assembled stiffness just before DAE
initialization, avoiding unqualified parallel ParaView/PNG postprocess.
Each MPI rank saves ONLY its PETSc-owned CSR rows, without a dense gather.
The driver reconstructs a complete global matrix and checks parity.

This is NOT yet MPI2 DAE/TS/CN/energy qualification.
"""

import argparse
import json
import os
from pathlib import Path

import numpy as np
from mpi4py import MPI
from petsc4py import PETSc
from scipy.sparse import csr_matrix, save_npz

from scrimp.dphs import DPHS
from examples.heat_wave import heat_wave_eq


class CapturedNativeStiffness(Exception):
    """Expected synchronized stop before algebraic DAE initialization."""


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--heat", type=int, required=True, choices=(1, 2))
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    wave = 3 - args.heat
    comm = MPI.COMM_WORLD
    rank, size = comm.rank, comm.size
    if size not in (1, 2) or PETSc.COMM_WORLD.getSize() != size:
        raise RuntimeError("This probe accepts only actual MPI1/MPI2")
    if os.environ.get("SCRIMP_V130_AUTO_INTERFACE_PILOT") != "1":
        raise RuntimeError("Explicit automatic-interface pilot flag required")
    if not args.out.is_dir():
        raise RuntimeError("Output directory must exist and be shared by all ranks")

    original_init = DPHS.init_step
    calls = []

    def after_native_assembly(system):
        if calls:
            raise RuntimeError("Multiple assembly interceptions on one rank")
        compiler = system._internal_interface_compiler
        if compiler is None or not compiler.completed:
            raise RuntimeError("Automatic interface compilation missing")
        interfaces = compiler.report["interfaces"]
        if len(interfaces) != 1:
            raise RuntimeError("Expected exactly one logical interface")
        iface = interfaces[0]
        expected = {"Interface Heat": args.heat, "Interface Wave": wave}
        if (iface["logical_region"] != 10
                or iface["port_owners"] != expected
                or iface["facet_count"] != 78):
            raise RuntimeError("Wrong interface owners / oriented facets")
        ports = {}
        for name, region in expected.items():
            port = system.ports[name]
            if (port.get_region() != 10
                    or port.get_assembly_region() != iface["side_regions"][region]):
                raise RuntimeError("Wrong effective one-sided port region")
            ports[name] = port.get_assembly_region()
        if ports["Interface Heat"] == ports["Interface Wave"]:
            raise RuntimeError("Both interfaces unexpectedly use the same side")

        K = system.stiffness
        nrows, ncols = map(int, K.getSize())
        lo, hi = map(int, K.getOwnershipRange())
        n_model = int(system.gf_model.nbdof())
        expected_n = 7689 if args.heat == 1 else 8056
        if (nrows != ncols or nrows != n_model or nrows != expected_n
                or not (0 <= lo < hi <= nrows)):
            raise RuntimeError(
                f"Invalid PETSc owned range: {lo}:{hi}/{nrows}, model={n_model}"
            )
        indptr, indices, values = K.getValuesCSR()
        local = csr_matrix(
            (np.asarray(values).copy(), np.asarray(indices).copy(),
             np.asarray(indptr).copy()),
            shape=(hi - lo, ncols),
        )
        local.sum_duplicates()
        local.sort_indices()
        if not np.all(np.isfinite(local.data)):
            raise RuntimeError("Non-finite distributed assembled K coefficient")
        basename = f"rank{rank}"
        save_npz(args.out / (basename + ".npz"), local, compressed=True)
        local_manifest = {
            "rank": rank, "mpi_size": size, "heat_region": args.heat,
            "wave_region": wave, "global_shape": [nrows, ncols],
            "PETSc_owned_rows": [lo, hi],
            "local_csr_nnz": int(local.nnz),
            "local_shape": list(local.shape),
            "logical_region": 10,
            "oriented_port_regions": ports,
            "interface_compiler": compiler.report,
            "native_assembly_executed": True,
            "DAE_initialized": False,
            "CN_executed": False,
            "MPI2_energy_qualified": False,
        }
        (args.out / (basename + ".json")).write_text(
            json.dumps(local_manifest, indent=2) + "\n"
        )
        calls.append(local_manifest)
        # No rank may stop before all others have published its local
        # OWNED rows. Avoid rank-0-only callbacks across MPI collectives.
        all_ranges = comm.allgather((lo, hi, nrows, ncols))
        expected_ranges = sorted((a, b) for a, b, _, _ in all_ranges)
        if all_ranges and (
            expected_ranges[0][0] != 0
            or expected_ranges[-1][1] != nrows
            or any(a[1] != b[0]
                   for a, b in zip(expected_ranges, expected_ranges[1:]))
            or any((nr, nc) != (nrows, ncols)
                   for _, _, nr, nc in all_ranges)
        ):
            raise RuntimeError("PETSc rank-owned row ranges do not partition K")
        comm.Barrier()
        if rank == 0:
            print("MPI_AUTO_K_OWNED_ROWS_CAPTURED",
                  json.dumps({"heat": args.heat, "wave": wave,
                              "ranks": size, "shape": [nrows, ncols],
                              "owned_ranges": expected_ranges}), flush=True)
        raise CapturedNativeStiffness

    DPHS.init_step = after_native_assembly
    try:
        try:
            heat_wave_eq(
                heat_region=args.heat, wave_region=wave,
                factor_solver="mumps", final_time=.004,
                time_step=.001, save_step=.001,
            )
        except CapturedNativeStiffness:
            pass
        else:
            raise RuntimeError("Actual MPI GetFEM/PETSc assembly not intercepted")
    finally:
        DPHS.init_step = original_init
    if len(calls) != 1:
        raise RuntimeError("Wrong number of native assembly captures")


if __name__ == "__main__":
    main()
