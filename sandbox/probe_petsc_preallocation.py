# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Isolated PETSc matrix-assembly/preallocation probe on the real acoustic model.

This is a targeted diagnostic, NOT a new profiling infrastructure.
Use SCRIMP's existing TimingProfiler and convert_gmm_to_petsc() with the
existing examples/wave.py physical bricks. Never alters the PDE or the
production allocation path.

The distributed GetFEM bridge may contribute to any global row:
do not restrict transfer to PETSc-owned rows or discard stored zeros.

Run inside scrimp-mpi:local, with source/sandbox/examples mounted:
  mpiexec --bind-to none -n 2 python -u sandbox/probe_petsc_preallocation.py \
    --mesh 0.08 --prealloc 64 --output /artifacts/prealloc_2r.json

Compare the original diagonal-only PETSc allocation against PETSc AIJ
with 64 slots per local row in BOTH diagonal/off-diagonal blocks.
The two matrices must have the SAME stored PETSc CSR and coefficients.
The latter is deliberately an overestimate, not a production policy.
"""

import argparse
import json
import os
from pathlib import Path

import getfem as gf
from mpi4py import MPI
import numpy as np
from petsc4py import PETSc

from bench_acoustic_2d import _build_wave
from scrimp.utils.linalg import convert_gmm_to_petsc
from scrimp.utils.profiling import TimingProfiler


MPI_COMM = MPI.COMM_WORLD
PETSC_COMM = PETSc.COMM_WORLD


def _model_matrix(system, block):
    """Assemble one unmodified linear block from examples/wave.py."""
    system.disable_all_bricks()
    for brick in system.bricks.values():
        if brick.get_linear() and bool(brick.get_dt()) == (block == "mass"):
            brick.enable_id_bricks(system.gf_model)
    system.gf_model.assembly(option="build_matrix")
    # Keep assembled GetFEM brick state until PETSc transfer completes,
    # mirroring DPHS.assemble_mass()/assemble_stiffness().
    return system.gf_model.tangent_matrix()


def _raw_statistics(matrix, lo, hi):
    """Count local/nonlocal rows in GetFEM's globally numbered CSR."""
    copied = gf.Spmat("copy", matrix)
    copied.transpose()
    ptr, cols = (np.asarray(a) for a in copied.csc_ind())
    values = np.asarray(copied.csc_val())
    begin, end = int(ptr[lo]), int(ptr[hi])
    return {
        "stored_nnz": int(len(values)),
        "stored_zeros": int(len(values) - np.count_nonzero(values)),
        "stored_nnz_petsc_owned_rows": int(end - begin),
        "stored_nnz_other_rows": int(len(values) - (end - begin)),
        "nonzero_contributions_other_rows": int(
            np.count_nonzero(values[:begin]) + np.count_nonzero(values[end:])
        ),
    }


def _create_matrix(system, n, prealloc):
    if prealloc is None:
        # Production baseline, including initial diagonal allocation.
        mat = PETSc.Mat().create(comm=PETSC_COMM)
        system.init_matrix(mat)
        return mat

    # Experimental exact same PETSc AIJ type, but a larger initial
    # capacity per local row. The two NNZ values describe diagonal and
    # off-diagonal blocks in MPIAIJ, not GetFEM's global contribution.
    mat = PETSc.Mat().createAIJ(
        size=(n, n), nnz=(prealloc, prealloc), comm=PETSC_COMM
    )
    mat.setOption(PETSc.Mat.Option.FORCE_DIAGONAL_ENTRIES, True)
    zeros = mat.createVecLeft()
    try:
        zeros.set(0.0)
        mat.setDiagonal(zeros, addv=PETSc.InsertMode.ADD_VALUES)
        mat.assemble()
    finally:
        zeros.destroy()
    return mat


def _assemble(system, tangent, n, block, prealloc):
    mat = _create_matrix(system, n, prealloc)
    profiler = TimingProfiler(enabled=True)
    label = f"probe.bridge.{block}"
    try:
        convert_gmm_to_petsc(
            tangent, mat, profiler=profiler, profile_label=label,
        )
        csr = tuple(np.asarray(x).copy() for x in mat.getValuesCSR())
        # PETSc's default getInfo() does a GLOBAL_SUM collective.
        # LOCAL avoids an extra collective inside this diagnostic.
        info = {
            key: float(value) for key, value in
            mat.getInfo(PETSc.Mat.InfoType.LOCAL).items()
        }
        timing = profiler.report(MPI_COMM)  # collective, AFTER conversion
        return csr, info, timing
    finally:
        mat.destroy()


def _compare_csr(reference, candidate):
    """Strict structural check, numerical value check, collective result."""
    same_ptr = np.array_equal(reference[0], candidate[0])
    same_cols = np.array_equal(reference[1], candidate[1])
    same_size = len(reference[2]) == len(candidate[2])
    same_vals = (
        bool(np.allclose(reference[2], candidate[2], atol=1e-12,
                         rtol=1e-12))
        if same_size else False
    )
    diff = (
        float(np.max(np.abs(reference[2] - candidate[2])))
        if same_size and len(reference[2]) else
        (0.0 if same_size else float("inf"))
    )
    return {
        "stored_csr_indptr_equal": bool(same_ptr),
        "stored_csr_columns_equal": bool(same_cols),
        "stored_csr_values_close": bool(same_vals),
        "max_absolute_value_difference": diff,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mesh", type=float, default=0.08)
    parser.add_argument("--prealloc", type=int, default=64,
                        help="Slots per row for each MPIAIJ block")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.mesh <= 0 or args.prealloc <= 0:
        parser.error("--mesh and --prealloc must be positive")

    mode = os.environ.get(
        "SCRIMP_GETFEM_MATRIX_MODE", "replicated"
    ).strip().lower()
    if mode != "distributed":
        raise ValueError(
            "This MPI assembly probe requires "
            "SCRIMP_GETFEM_MATRIX_MODE=distributed"
        )

    system = _build_wave(args.mesh)
    system.gf_model.set_time(0.0)
    system.gf_model.first_iter()
    n = int(system.gf_model.nbdof())
    records = []

    for block in ("mass", "stiffness"):
        tangent = _model_matrix(system, block)
        # Matrix dimensions are global and PETSc row ownership is
        # independent of preallocation. Obtain range on a temporary Mat.
        owned_mat = PETSc.Mat().create(comm=PETSC_COMM)
        system.init_matrix(owned_mat)
        lo, hi = owned_mat.getOwnershipRange()
        owned_mat.destroy()

        raw = _raw_statistics(tangent, lo, hi)
        reference, default_info, default_timing = _assemble(
            system, tangent, n, block, None
        )
        candidate, prealloc_info, prealloc_timing = _assemble(
            system, tangent, n, block, args.prealloc
        )
        local_parity = _compare_csr(reference, candidate)
        gathered = MPI_COMM.allgather({
            "rank": MPI_COMM.rank,
            "owned_row_range": [int(lo), int(hi)],
            "raw_getfem": raw,
            "default_petsc_info": default_info,
            "prealloc_petsc_info": prealloc_info,
            "parity": local_parity,
        })
        records.append({
            "block": block,
            "by_rank": gathered,
            "baseline_profile": default_timing,
            "preallocated_profile": prealloc_timing,
        })
        if MPI_COMM.rank == 0:
            print(f"{block}: collected raw GetFEM/PETSc allocation data", flush=True)

        if not all(
            item["parity"]["stored_csr_indptr_equal"]
            and item["parity"]["stored_csr_columns_equal"]
            and item["parity"]["stored_csr_values_close"]
            for item in gathered
        ):
            raise RuntimeError(
                f"{block}: PETSc stored CSR differs after preallocation"
            )
        system.disable_all_bricks()

    report = {
        "model": "examples/wave.py wave_eq (unchanged linear bricks)",
        "matrix_mode": mode,
        "mesh": args.mesh,
        "degrees_of_freedom": n,
        "mpi_size": MPI_COMM.size,
        "prealloc_diagonal_nnz_per_row": args.prealloc,
        "prealloc_offdiagonal_nnz_per_row": args.prealloc,
        "blocks": records,
    }
    if MPI_COMM.rank == 0:
        rendered = json.dumps(report, indent=2)
        print(rendered, flush=True)
        if args.output:
            args.output.parent.mkdir(parents=True, exist_ok=True)
            args.output.write_text(rendered + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
