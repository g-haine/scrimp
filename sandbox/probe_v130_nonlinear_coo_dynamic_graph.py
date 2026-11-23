# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Native MPI1/2 COO reassembly with a deliberately changing matrix graph.

The SAME PETSc AIJ Mat receives GetFEM mass matrices alternating between
two disjoint scalar fields and a coupled-field operator. Its graph must
grow, shrink, grow, and shrink without stale off-diagonal blocks.

The production convert_gmm_to_petsc bridge is called with its default
arguments (no special COO opt-in). Checks stored CSR graph and values,
rank-local ownership, and matrix-entry parity against MPI1.

This is an isolated graph test; actual nonlinear pH Jacobians are
independently covered by the original shallow-water directional FD gate.
"""

import argparse
import json
import os
from pathlib import Path

import getfem as gf
from mpi4py import MPI
import numpy as np
from petsc4py import PETSc

from scrimp.utils.linalg import convert_gmm_to_petsc


def getfem_matrix(mesh, mim, *, coupled):
    model = gf.Model("real")
    fem = gf.MeshFem(mesh, 1)
    fem.set_fem(gf.Fem("FEM_PK(1,1)"))
    model.add_fem_variable("u", fem)
    model.add_fem_variable("v", fem)
    form = "u*Test_u+v*Test_v"
    if coupled:
        form += "+0.25*(u*Test_v+v*Test_u)"
    model.add_linear_term(mim, form)
    model.assembly(option="build_matrix")
    return model


def csr_snapshot(mat):
    n = mat.getSize()[0]
    lo, hi = mat.getOwnershipRange()
    ptr, cols, values = mat.getValuesCSR()
    ptr = np.asarray(ptr, dtype=np.int64)
    cols = np.asarray(cols, dtype=np.int64)
    values = np.asarray(values)
    if ptr.size != hi - lo + 1 or ptr[0] != 0 or ptr[-1] != values.size:
        raise RuntimeError("Invalid PETSc CSR row structure")
    if not np.all(np.isfinite(values)):
        raise RuntimeError("Nonfinite COO-updated matrix coefficients")
    if values.size and float(np.max(np.abs(values.imag))) > 1e-12:
        raise RuntimeError("Complex coefficients in real GetFEM mass")
    dense = np.zeros((hi-lo, n), dtype=np.float64)
    for row in range(hi-lo):
        dense[row, cols[ptr[row]:ptr[row+1]]] = values[ptr[row]:ptr[row+1]].real
    return {
        "ownership": [int(lo), int(hi)],
        "stored_nnz": int(values.size),
        "indptr": ptr.tolist(),
        "columns": cols.tolist(),
        "values_real": values.real.tolist(),
        "owned_dense_rows": dense.tolist(),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    mpi = MPI.COMM_WORLD
    if mpi.size not in (1,2) or PETSc.COMM_WORLD.getSize() != mpi.size:
        raise RuntimeError("MPI1/2 only")
    if int(gf.util_mpi_parallelism_level()) != 2:
        raise RuntimeError("GetFEM MPI level 2 required")
    if os.environ.get("SCRIMP_GETFEM_MATRIX_MODE") != "distributed":
        raise RuntimeError("This native graph test requires distributed GetFEM")
    if os.environ.get("SCRIMP_GETFEM_BRIDGE_DISTRIBUTED") != "coo":
        raise RuntimeError("This native graph test must use standard COO")
    args.out = args.out.with_name(
        f"dynamic_rank_{mpi.rank}.json"
    )
    if not args.out.parent.is_dir() or args.out.exists():
        raise RuntimeError("Host evidence path must be new and writable")

    mesh = gf.Mesh("regular simplices", np.linspace(0.,1.,6))
    mim = gf.MeshIm(mesh, gf.Integ("IM_GAUSS1D(5)"))
    model0 = getfem_matrix(mesh, mim, coupled=False)
    model1 = getfem_matrix(mesh, mim, coupled=True)
    n = int(model0.nbdof())
    if n != int(model1.nbdof()) or n < 6:
        raise RuntimeError("Unexpected GetFEM cross-coupled FE dimensions")
    A = PETSc.Mat().createAIJ(size=(n,n),comm=PETSc.COMM_WORLD)
    rows = []
    try:
        for label, source in (
            ("uncoupled", model0), ("coupled", model1),
            ("uncoupled_repeat", model0), ("coupled_repeat", model1),
        ):
            # Same allocated PETSc Mat for every NEW GetFEM tangent.
            convert_gmm_to_petsc(
                source.tangent_matrix(), A, comm=PETSc.COMM_WORLD
            )
            snapshot = csr_snapshot(A)
            row_nnz = mpi.allreduce(snapshot["stored_nnz"], op=MPI.SUM)
            snapshot.update({
                "label":label,
                "global_stored_nnz":int(row_nnz),
            })
            rows.append(snapshot)
            print(
                f"DYNAMIC COO MPI{mpi.size} rank={mpi.rank} "
                f"step={label} owned_nnz={snapshot['stored_nnz']} "
                f"global_nnz={row_nnz}", flush=True
            )
        for a,b in ((0,2),(1,3)):
            for key in ("indptr","columns"):
                if rows[a][key]!=rows[b][key]:
                    raise RuntimeError(
                        f"COO graph leaked between {rows[a]['label']} "
                        f"and {rows[b]['label']}"
                    )
            if not np.allclose(
                rows[a]["values_real"], rows[b]["values_real"],
                atol=1e-13,rtol=1e-13
            ):
                raise RuntimeError("Repeated COO matrix coefficient drift")
        if not rows[1]["global_stored_nnz"] > rows[0]["global_stored_nnz"]:
            raise RuntimeError(
                "Dynamic COO test did not actually change the graph"
            )
        if rows[0]["global_stored_nnz"] != rows[2]["global_stored_nnz"]:
            raise RuntimeError("COO failed to shrink sparsity graph")
        if rows[1]["global_stored_nnz"] != rows[3]["global_stored_nnz"]:
            raise RuntimeError("COO failed to restore expanded graph")

        report = {
            "schema":"scrimp-v130-nonlinear-coo-dynamic-pattern-v1",
            "status":"PASS_NATIVE_DYNAMIC_COO_REPREALLOCATION",
            "MPI":mpi.size,"rank":mpi.rank,"n":n,
            "same_petsc_mat_reused":True,
            "steps":rows,
            "GetFEM_matrix_mode":"distributed",
            "PETSc_transfer":"COO assembly, AIJ/CSR storage",
            "DAE_or_full_physics_qualified":False,
        }
        args.out.write_text(json.dumps(report,indent=2)+"\n",
                            encoding="utf-8")
        print(f"PASS DYNAMIC NONLINEAR COO GRAPH MPI{mpi.size} R{mpi.rank}",
              flush=True)
    finally:
        A.destroy()


if __name__ == "__main__":
    main()
