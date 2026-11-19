# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Isolated GetFEM distributed-matrix transfer: PETSc IJV vs native COO.

Uses the unchanged linear mass and stiffness bricks from examples/wave.py.
The reference path is SCRIMP's existing convert_gmm_to_petsc(), including
its initial diagonal PETSc allocation. The candidate path uses PETSc
setPreallocationCOO()/setValuesCOO() on the SAME GetFEM contributions.

PETSc COO accepts off-process rows and sums repeated coordinates. Owned
zero diagonal entries are provided explicitly to match SCRIMP.init_matrix
and preserve its stored CSR pattern (important for MUMPS symbolic reuse).

No production code is modified. Strict per-rank PETSc stored-CSR parity
is required before any candidate timing is considered valid. The test
measures cold and reversed-order trials to identify order effects.

Run only in the validated scrimp-mpi:local image:
  mpiexec --bind-to none -n 2 python -u sandbox/probe_petsc_coo.py \
      --mesh 0.08 --output /artifacts/coo_2r.json
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


def _linear_tangent(system, block):
    system.disable_all_bricks()
    for brick in system.bricks.values():
        if brick.get_linear() and bool(brick.get_dt()) == (block == "mass"):
            brick.enable_id_bricks(system.gf_model)
    system.gf_model.assembly(option="build_matrix")
    return system.gf_model.tangent_matrix()


def _csr_from_getfem(tangent, n):
    # GetFEM exposes CSC of the transpose, hence CSR of the original.
    copied = gf.Spmat("copy", tangent)
    copied.transpose()
    ptr, cols = copied.csc_ind()
    values = copied.csc_val()
    ptr = np.asarray(ptr, dtype=PETSc.IntType)
    cols = np.asarray(cols, dtype=PETSc.IntType)
    values = np.asarray(values, dtype=PETSc.ScalarType)
    if (len(ptr) != n + 1 or len(cols) != len(values)
            or int(ptr[0]) != 0 or int(ptr[-1]) != len(values)):
        raise ValueError("Unexpected GetFEM CSR shape/count")
    if len(cols) and (np.min(cols) < 0 or np.max(cols) >= n):
        raise ValueError("Invalid GetFEM global column indices")
    return ptr, cols, values


def _coo_from_csr(ptr, cols, values, n, lo, hi):
    counts = np.diff(ptr).astype(np.intp, copy=False)
    rows = np.repeat(np.arange(n, dtype=PETSc.IntType), counts)

    # The existing DPHS.init_matrix inserts structural zeros on ALL
    # diagonal entries. Keeping them avoids a false CSR pattern change.
    owned_diag = np.arange(lo, hi, dtype=PETSc.IntType)
    return (
        np.ascontiguousarray(np.concatenate((rows, owned_diag))),
        np.ascontiguousarray(np.concatenate((cols, owned_diag))),
        np.ascontiguousarray(np.concatenate((
            values,
            np.zeros(hi - lo, dtype=PETSc.ScalarType),
        ))),
    )


def _run_method(system, tangent, n, method):
    profiler = TimingProfiler(enabled=True)
    mat = None
    try:
        with profiler.measure(f"probe.{method}.total"):
            if method == "reference":
                with profiler.measure("probe.reference.initialize"):
                    mat = PETSc.Mat().create(comm=PETSC_COMM)
                    system.init_matrix(mat)
                with profiler.measure("probe.reference.bridge_total"):
                    convert_gmm_to_petsc(
                        tangent, mat, profiler=profiler,
                        profile_label="probe.reference.bridge",
                    )
            else:
                with profiler.measure("probe.coo.initialize"):
                    mat = PETSc.Mat().create(comm=PETSC_COMM)
                    mat.setSizes((n, n))
                    mat.setType(PETSc.Mat.Type.AIJ)
                with profiler.measure("probe.coo.prepare"):
                    # Matrix ownership is not necessarily initialized
                    # before PETSc COO preallocation. Use the default
                    # PETSc Vec layout on the same communicator.
                    layout = PETSc.Vec().createMPI(n, comm=PETSC_COMM)
                    try:
                        lo, hi = layout.getOwnershipRange()
                    finally:
                        layout.destroy()
                    ptr, cols, vals = _csr_from_getfem(tangent, n)
                    row, col, val = _coo_from_csr(
                        ptr, cols, vals, n, lo, hi
                    )
                with profiler.measure("probe.coo.preallocate"):
                    mat.setPreallocationCOO(row, col)
                if mat.getOwnershipRange() != (lo, hi):
                    raise RuntimeError(
                        "PETSc matrix and vector row layouts differ"
                    )
                with profiler.measure("probe.coo.set_values"):
                    # PETSc COO sums duplicate and off-process entries.
                    # INSERT replaces all values, like zeroEntries+ADD.
                    mat.setValuesCOO(
                        val, addv=PETSc.InsertMode.INSERT_VALUES
                    )
                # MatSetValuesCOO leaves the PETSc matrix assembled:
                # do not add an extra Mat.assemble() here.

        # Inspect only AFTER both assembly paths have completed.
        csr = tuple(np.array(x, copy=True) for x in mat.getValuesCSR())
        info = {
            k: float(v) for k, v in
            mat.getInfo(PETSc.Mat.InfoType.LOCAL).items()
        }
        timings = profiler.report(MPI_COMM)
        return csr, info, timings
    finally:
        if mat is not None:
            mat.destroy()


def _parity(reference, candidate):
    same_ptr = bool(np.array_equal(reference[0], candidate[0]))
    same_cols = bool(np.array_equal(reference[1], candidate[1]))
    shape_ok = reference[2].shape == candidate[2].shape
    values_close = (
        bool(np.allclose(reference[2], candidate[2],
                         atol=1e-12, rtol=1e-12))
        if shape_ok else False
    )
    max_diff = (
        float(np.max(np.abs(reference[2] - candidate[2])))
        if shape_ok and reference[2].size else
        (0.0 if shape_ok else None)
    )
    return {
        "stored_csr_indptr_equal": same_ptr,
        "stored_csr_columns_equal": same_cols,
        "stored_csr_values_close": values_close,
        "max_absolute_value_difference": max_diff,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mesh", type=float, default=0.08)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.mesh <= 0:
        parser.error("--mesh must be positive")
    if os.environ.get(
        "SCRIMP_GETFEM_MATRIX_MODE", "replicated"
    ).strip().lower() != "distributed":
        raise ValueError(
            "Probe requires SCRIMP_GETFEM_MATRIX_MODE=distributed"
        )

    system = _build_wave(args.mesh)
    system.gf_model.set_time(0.)
    system.gf_model.first_iter()
    n = int(system.gf_model.nbdof())
    records = []
    for block in ("mass", "stiffness"):
        tangent = _linear_tangent(system, block)
        trials = []
        for repeat in range(2):
            # Repeat in reversed order to distinguish real improvement
            # from PETSc caches and first-time MPI startup overhead.
            order = ("reference", "coo") if repeat == 0 else (
                "coo", "reference"
            )
            data = {}
            for method in order:
                csr, info, timings = _run_method(
                    system, tangent, n, method
                )
                data[method] = {
                    "csr": csr, "petsc_info_local": info,
                    "profile": timings,
                }
            parity = _parity(
                data["reference"]["csr"], data["coo"]["csr"]
            )
            parity_all = MPI_COMM.allgather(parity)
            if not all(
                p["stored_csr_indptr_equal"]
                and p["stored_csr_columns_equal"]
                and p["stored_csr_values_close"]
                for p in parity_all
            ):
                raise RuntimeError(
                    f"{block} trial {repeat}: PETSc CSR mismatch: "
                    f"{parity_all}"
                )
            trials.append({
                "trial": repeat,
                "order": list(order),
                "parity_by_rank": parity_all,
                "reference": {
                    "profile": data["reference"]["profile"],
                    "petsc_info_by_rank": MPI_COMM.allgather(
                        data["reference"]["petsc_info_local"]
                    ),
                },
                "coo": {
                    "profile": data["coo"]["profile"],
                    "petsc_info_by_rank": MPI_COMM.allgather(
                        data["coo"]["petsc_info_local"]
                    ),
                },
            })
            if MPI_COMM.rank == 0:
                print(
                    f"PASS {block} trial={repeat}, order={order}: "
                    "full PETSc stored-CSR parity",
                    flush=True,
                )

        system.disable_all_bricks()
        records.append({"block": block, "trials": trials})

    report = {
        "model": "unchanged examples/wave.py linear mass/stiffness",
        "mpi_size": MPI_COMM.size,
        "degrees_of_freedom": n,
        "mesh": args.mesh,
        "getfem_matrix_mode": "distributed",
        "candidate": "PETSc setPreallocationCOO + setValuesCOO",
        "blocks": records,
    }
    if MPI_COMM.rank == 0:
        payload = json.dumps(report, indent=2)
        print(payload, flush=True)
        if args.output:
            args.output.parent.mkdir(parents=True, exist_ok=True)
            args.output.write_text(payload + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
