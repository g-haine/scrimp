# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Measure native GetFEM matrix/RHS MPI conventions BEFORE SCRIMP scaling.

Independent of DPHS, PETSc matrix conversion, solver and MUMPS.
The exact 1D reference domain is [0,1], with u=1 and f=2:
  generic RHS for u*u*Test_u = -1,
  source-brick RHS for f=2    = +2,
  mixed generic + source     = +1,
  linear tangent of u*Test_u = 1 (one^T M one),
  nonlinear tangent of u^2   = 2 (one^T DN one).
These integrals must be unchanged by the process count.

The GetFEM MPI documentation says that global RHS vectors are copied
on all ranks whereas tangent matrices carry distributed contributions.
We test BOTH assertions independently: a complete RHS alone does not
justify division by MPI size!

Run with 1, 2, 4 MPI ranks under BOTH local GetFEM 5.4.x and 5.5:
  mpiexec -n 2 python sandbox/diagnose_getfem_mpi_assembly.py \
    --output sandbox/outputs/getfem_assembly_2.json
"""
import argparse
import hashlib
import json
from pathlib import Path

import getfem as gf
from mpi4py import MPI
import numpy as np

COMM = MPI.COMM_WORLD


def _prepare_model(kind):
    mesh = gf.Mesh("cartesian", np.linspace(0.0, 1.0, 17))
    fem = gf.MeshFem(mesh, 1)
    fem.set_fem(gf.Fem("FEM_PK(1,1)"))
    integration = gf.MeshIm(mesh, gf.Integ("IM_GAUSS1D(5)"))
    model = gf.Model("real")
    model.add_fem_variable("u", fem)
    if kind in ("generic", "mixed", "nonlinear_tangent"):
        model.add_nonlinear_term(integration, "u*u*Test_u")
    if kind in ("source", "mixed"):
        model.add_initialized_data("forcing", np.array([2.0]))
        model.add_source_term_brick(integration, "u", "forcing")
    if kind == "linear_tangent":
        model.add_linear_term(integration, "u*Test_u")
    model.set_variable("u", np.ones(model.nbdof()))
    # GetFEM's C++ context tracks dependencies by object identity.
    # Holding only Model is NOT enough: Mesh, MeshFem and MeshIm must
    # stay alive until all native assemblies and matrix probes finish.
    return model, mesh, fem, integration


def _matrix_fingerprint(model):
    raw = gf.Spmat("copy", model.tangent_matrix())
    raw.transpose()  # CSC of transpose => CSR of original matrix
    ptr, col = (np.asarray(v) for v in raw.csc_ind())
    val = np.asarray(raw.csc_val())
    h = hashlib.blake2b(digest_size=12)
    for part in (ptr, col, val):
        h.update(np.ascontiguousarray(part).tobytes())
    n = int(model.nbdof())
    if len(ptr) != n + 1:
        raise RuntimeError(
            f"Native tangent has {len(ptr)-1} rows instead of {n}"
        )
    return {
        "dofs": n,
        "raw_matrix_ones_integral": float(np.real(np.sum(val))),
        "raw_matrix_nnz_stored": int(len(val)),
        "raw_matrix_digest": h.hexdigest(),
    }


def _case(kind):
    model, mesh, fem, integration = _prepare_model(kind)
    # Keep the native GetFEM dependency graph alive for the entire case.
    # Local refs are intentional despite no explicit Python-side access.
    _dependencies = (mesh, fem, integration)
    n = int(model.nbdof())
    record = {"kind": kind, "dofs": n}
    if kind in ("generic", "source", "mixed"):
        model.assembly(option="build_rhs")
        native_rhs = np.asarray(model.rhs(), dtype=float).ravel().copy()
        record["raw_rhs_length"] = int(native_rhs.size)
        record["raw_rhs_integral"] = float(np.sum(native_rhs))
        # This is strictly a diagnostic, NOT an MPI reduction of the RHS:
        # GetFEM should already have performed the sum and replicated it.
        copies = COMM.allgather(native_rhs)
        record["rhs_copy_max_difference"] = max(
            float(np.linalg.norm(a - copies[0]))
            for a in copies
        ) if all(len(a) == len(copies[0]) for a in copies) else None
    if kind in ("linear_tangent", "nonlinear_tangent"):
        model.assembly(option="build_matrix")
        record.update(_matrix_fingerprint(model))
    per_rank = COMM.allgather(record)
    return {
        "kind": kind,
        "per_rank": per_rank,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args()
    results = [
        _case(kind) for kind in (
            "generic", "source", "mixed", "linear_tangent",
            "nonlinear_tangent"
        )
    ]
    report = {
        "mpi_size": COMM.size,
        "getfem_mpi_level": int(gf.util_mpi_parallelism_level()),
        "numpy_version": np.__version__,
        "cases": results,
    }
    if COMM.rank == 0:
        rendered = json.dumps(report, indent=2)
        print(rendered, flush=True)
        if args.output is not None:
            args.output.parent.mkdir(parents=True, exist_ok=True)
            args.output.write_text(rendered + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
