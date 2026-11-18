# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See LICENSE for license information.
#
# github: https://github.com/g-haine/scrimp

"""Read-only distributed variable-block coupling audit for assembled PETSc matrices.

Works with any disjoint, complete GetFEM variable layout; no variable names
are hardcoded. Reports directed block couplings, undirected connectivity,
and *candidates* for algebraic constraints. A zero diagonal block is
not by itself a proof of a Lagrange multiplier or singularity.

Uses Mat.getRow() on locally owned GLOBAL rows to avoid the MPIAIJ CSR
off-process column-numbering ambiguity. No global dense matrix or
PETSc vector/matrix gather. Only O(number_of_fields**2) counts and
norm summaries are reduced. This is a diagnostic, not a production PC.
"""
import argparse
import json
import os
from pathlib import Path

import numpy as np
from mpi4py import MPI
from petsc4py import PETSc

from bench_acoustic_2d import _build_wave
from scrimp.utils.ph_block_layout import discover_layout, classify_ph_formulation


def audit_blocks(A, fields, *, atol=0.0):
    """Return rank-invariant directed adjacency statistics for PETSc A."""
    comm = MPI.COMM_WORLD
    n = fields[-1]["stop"]
    nf = len(fields)
    if nf < 1 or A.getSize() != (n, n):
        raise ValueError("Matrix size and GetFEM variable layout disagree")
    if atol < 0 or not np.isfinite(atol):
        raise ValueError("Invalid nonzero threshold")

    starts = np.array([f["start"] for f in fields], dtype=np.int64)
    stops = np.array([f["stop"] for f in fields], dtype=np.int64)
    lo, hi = A.getOwnershipRange()
    stored = np.zeros((nf, nf), dtype=np.int64)
    nonzero = np.zeros((nf, nf), dtype=np.int64)
    fro2 = np.zeros((nf, nf), dtype=np.float64)
    diagonal_zero = np.zeros(nf, dtype=np.int64)
    diagonal_total = np.zeros(nf, dtype=np.int64)
    # getRow returns global column numbers on MPIAIJ, unlike raw local CSR.
    for row in range(lo, hi):
        source = int(np.searchsorted(starts, row, side="right") - 1)
        if source < 0 or row >= stops[source]:
            raise RuntimeError("Unclassified owned PETSc row")
        cols, values = A.getRow(row)
        try:
            cols = np.asarray(cols, dtype=np.int64)
            values = np.asarray(values)
            if cols.size and (
                np.any(cols < 0) or np.any(cols >= n)
                or not np.all(np.isfinite(values))
            ):
                raise RuntimeError("Invalid sparse row")
            if cols.size:
                targets = np.searchsorted(starts, cols, side="right") - 1
                if np.any(targets < 0) or np.any(cols >= stops[targets]):
                    raise RuntimeError("Unclassified PETSc column")
                stored[source] += np.bincount(targets, minlength=nf)
                active = np.abs(values) > atol
                nonzero[source] += np.bincount(
                    targets[active], minlength=nf
                )
                fro2[source] += np.bincount(
                    targets, weights=np.abs(values)**2, minlength=nf
                )
                diagonal = values[cols == row]
                is_zero = not np.any(np.abs(diagonal) > atol)
            else:
                is_zero = True
            diagonal_total[source] += 1
            diagonal_zero[source] += int(is_zero)
        finally:
            A.restoreRow(row, cols, values)

    # Reduce only small per-block statistics; never gather sparse entries.
    global_stored = np.empty_like(stored)
    global_nonzero = np.empty_like(nonzero)
    global_fro2 = np.empty_like(fro2)
    global_zero = np.empty_like(diagonal_zero)
    global_total = np.empty_like(diagonal_total)
    comm.Allreduce(stored, global_stored, op=MPI.SUM)
    comm.Allreduce(nonzero, global_nonzero, op=MPI.SUM)
    comm.Allreduce(fro2, global_fro2, op=MPI.SUM)
    comm.Allreduce(diagonal_zero, global_zero, op=MPI.SUM)
    comm.Allreduce(diagonal_total, global_total, op=MPI.SUM)
    if int(global_total.sum()) != n:
        raise RuntimeError("MPI ownership does not cover the full matrix")

    undirected = (global_nonzero > 0) | (global_nonzero.T > 0)
    visited = set()
    components = []
    for seed in range(nf):
        if seed in visited:
            continue
        todo = [seed]
        component = []
        visited.add(seed)
        while todo:
            i = todo.pop()
            component.append(i)
            for j in np.flatnonzero(undirected[i]):
                j = int(j)
                if j not in visited:
                    visited.add(j)
                    todo.append(j)
        components.append(sorted(fields[i]["name"] for i in component))

    edges = []
    for i in range(nf):
        for j in range(nf):
            if global_nonzero[i, j]:
                edges.append({
                    "row": fields[i]["name"],
                    "column": fields[j]["name"],
                    "entries": int(global_nonzero[i, j]),
                    "frobenius": float(np.sqrt(global_fro2[i, j])),
                })

    candidates = []
    for i, field in enumerate(fields):
        if global_nonzero[i, i] != 0:
            continue
        outgoing = [
            fields[j]["name"] for j in range(nf)
            if j != i and global_nonzero[i, j]
        ]
        incoming = [
            fields[j]["name"] for j in range(nf)
            if j != i and global_nonzero[j, i]
        ]
        candidates.append({
            "name": field["name"],
            "role_metadata": field["group"],
            "outgoing": outgoing,
            "incoming": incoming,
            "status": "zero_self_block_candidate_not_proven_constraint",
        })
    return {
        "stored_nnz": int(global_stored.sum()),
        "effective_nnz": int(global_nonzero.sum()),
        "fields": [
            {
                **f,
                "zero_diagonal_entries": int(global_zero[i]),
                "diagonal_entries": int(global_total[i]),
                "self_block_effective_nnz": int(global_nonzero[i, i]),
            }
            for i, f in enumerate(fields)
        ],
        "directed_edges": edges,
        "undirected_components": components,
        "zero_self_block_candidates": candidates,
    }


def run(args):
    mpi = MPI.COMM_WORLD
    if not (args.mesh > 0 and args.dt > 0 and
            args.shift1 > 0 and args.shift2 > 0 and
            args.shift1 != args.shift2):
        raise ValueError("Require positive mesh/dt and distinct positive shifts")
    if (os.environ.get("SCRIMP_GETFEM_MATRIX_MODE") != "distributed"
            or os.environ.get("SCRIMP_GETFEM_BRIDGE_DISTRIBUTED") != "coo"):
        raise RuntimeError("Distributed COO GetFEM-PETSc bridge required")

    system = _build_wave(args.mesh)
    system.set_time_scheme(
        ts_type="cn", ts_adapt_type="none",
        ts_exact_final_time="matchstep",
        dt=args.dt, dt_save=args.dt, t_f=2*args.dt,
        ksp_type="preonly", pc_type="lu",
        pc_factor_mat_solver_type="mumps",
    )
    system.solve()
    if (system.last_ts_converged_reason <= 0
            or system.last_ts_steps != 2
            or system.last_ts_rejections
            or system.last_ts_snes_failures):
        raise RuntimeError("Reference CN/MUMPS integration failed")
    if not system.linear_mass or not system.linear_stiffness:
        raise RuntimeError("Only linear time-invariant systems supported by this probe")

    fields = discover_layout(
        system.gf_model, states=system.states,
        costates=system.costates, ports=system.ports,
        controls=system.controls,
    )
    formulation = classify_ph_formulation(fields, costates=system.costates)
    audits = []
    for sigma in (args.shift1, args.shift2):
        A = system.stiffness.copy()
        A.axpy(
            sigma, system.mass,
            structure=PETSc.Mat.Structure.DIFFERENT_NONZERO_PATTERN,
        )
        A.assemble()
        try:
            result = audit_blocks(A, fields, atol=args.atol)
        finally:
            A.destroy()
        audits.append({"sigma": sigma, **result})
        if mpi.rank == 0:
            names = [v["name"] for v in result["zero_self_block_candidates"]]
            print(
                f"AUDIT sigma={sigma:g}: fields={len(fields)}, "
                f"edges={len(result['directed_edges'])}, "
                f"zero-self candidates={names}", flush=True,
            )
    if [f["name"] for f in audits[0]["fields"]] != [
        f["name"] for f in audits[1]["fields"]
    ]:
        raise RuntimeError("Variable layout changed across shifts")
    report = {
        "schema": "scrimp-ph-block-graph-audit-v1",
        "source_commit": os.getenv("SCRIMP_BENCH_GIT_SHA", ""),
        "mpi_ranks": mpi.size,
        "mesh_h": args.mesh,
        "dofs": int(system.gf_model.nbdof()),
        "formulation": formulation,
        "audits": audits,
        "limitations": [
            "Adjacency and diagonals are structural evidence, not a rank test.",
            "A zero self-block is not proof of a Lagrange multiplier.",
            "Only assembled A(sigma) for two shifts, not a nonlinear/time trajectory.",
            "No new preconditioner or change to production time integration.",
            "GetRow audit is O(matrix nnz) and diagnostic-only.",
        ],
    }
    if mpi.rank == 0:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2) + "\n")
        print("Saved:", args.output, flush=True)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--mesh", type=float, default=0.08)
    p.add_argument("--dt", type=float, default=0.002)
    p.add_argument("--shift1", type=float, default=1000)
    p.add_argument("--shift2", type=float, default=2000)
    p.add_argument("--atol", type=float, default=0.0)
    p.add_argument("--output", type=Path, required=True)
    run(p.parse_args())


if __name__ == "__main__":
    main()
