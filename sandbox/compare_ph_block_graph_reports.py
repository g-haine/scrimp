# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Validate distributed block-graph reports and propose *candidate* neighborhoods.

Pure standard-library postprocessor: safe to run on a host without PETSc.
Never infers the physical status of a variable from its name. A zero
diagonal/self block is only a structural constraint candidate; this
report does NOT authorize any solver change or Schur elimination.
"""
import argparse
import json
import math
from pathlib import Path


def edge_map(audit):
    edges = {}
    for edge in audit["directed_edges"]:
        key = (edge["row"], edge["column"])
        if key in edges:
            raise ValueError(f"Duplicate directed edge {key}")
        edges[key] = edge
    return edges


def validate_pair(left, right, *, rtol=1e-10, atol=1e-12):
    """Assert MPI-rank invariant graph statistics at matching shifts."""
    if left["dofs"] != right["dofs"]:
        raise AssertionError("Inconsistent global DoF count across MPI ranks")
    lseries, rseries = left["audits"], right["audits"]
    if len(lseries) != len(rseries):
        raise AssertionError("MPI reports contain different shift counts")
    if not lseries:
        raise AssertionError("No audited Jacobians")
    summaries = []
    for la, ra in zip(lseries, rseries):
        if la["sigma"] != ra["sigma"]:
            raise AssertionError("MPI reports use different shift ordering")
        for key in ("fields", "stored_nnz", "effective_nnz",
                    "undirected_components", "zero_self_block_candidates"):
            if la[key] != ra[key]:
                raise AssertionError(
                    f"Mismatch in {key} for sigma={la['sigma']}"
                )
        le, re = edge_map(la), edge_map(ra)
        if le.keys() != re.keys():
            raise AssertionError(f"Edge mismatch at sigma={la['sigma']}")
        max_relative = 0.0
        for key in le:
            if le[key]["entries"] != re[key]["entries"]:
                raise AssertionError(f"Nonzero-count mismatch on {key}")
            a = le[key]["frobenius"]
            b = re[key]["frobenius"]
            if not (math.isfinite(a) and math.isfinite(b)):
                raise AssertionError("Nonfinite Frobenius norm")
            if not math.isclose(a, b, rel_tol=rtol, abs_tol=atol):
                raise AssertionError(
                    f"Frobenius mismatch on {key}: {a} vs {b}"
                )
            max_relative = max(
                max_relative, abs(a-b)/max(abs(a), abs(b), atol)
            )
        summaries.append({
            "sigma": la["sigma"],
            "fields": len(la["fields"]),
            "directed_edges": len(le),
            "max_mpi_relative_frobenius_difference": max_relative,
        })
    # This does not require the graph to be invariant with sigma;
    # general nonlinear Jacobians may legitimately change sparsity.
    return summaries


def constraint_neighborhoods(report):
    """Return immediate directed coupling neighborhoods; no reordering."""
    result = []
    for audit in report["audits"]:
        edges = edge_map(audit)
        names = {f["name"] for f in audit["fields"]}
        for candidate in audit["zero_self_block_candidates"]:
            name = candidate["name"]
            if name not in names:
                raise AssertionError("Unknown variable in candidate list")
            outgoing = sorted(col for row, col in edges if row == name and col != name)
            incoming = sorted(row for row, col in edges if col == name and row != name)
            if incoming != sorted(candidate["incoming"]) or outgoing != sorted(candidate["outgoing"]):
                raise AssertionError("Candidate graph directions disagree")
            result.append({
                "sigma": audit["sigma"],
                "candidate": name,
                "role_metadata": candidate["role_metadata"],
                "row_couplings": outgoing,
                "column_couplings": incoming,
                "local_neighborhood": sorted(set([name] + incoming + outgoing)),
                "qualification": "structural_candidate_only",
            })
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("report_1r", type=Path)
    p.add_argument("report_2r", type=Path)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    a = json.loads(args.report_1r.read_text(encoding="utf-8"))
    b = json.loads(args.report_2r.read_text(encoding="utf-8"))
    if a["schema"] != "scrimp-ph-block-graph-audit-v1" or b["schema"] != a["schema"]:
        raise ValueError("Unexpected graph report schema")
    if a["mpi_ranks"] != 1 or b["mpi_ranks"] != 2:
        raise ValueError("Expected 1-rank and 2-rank input reports")
    if a["mesh_h"] != b["mesh_h"] or a["source_commit"] != b["source_commit"]:
        raise ValueError("Cannot compare different mesh or source commit")
    parity = validate_pair(a, b)
    neighborhoods = constraint_neighborhoods(a)
    result = {
        "schema": "scrimp-ph-block-graph-comparison-v1",
        "source_commit": a["source_commit"],
        "mesh_h": a["mesh_h"],
        "dofs": a["dofs"],
        "mpi_parity": "pass",
        "shifts": parity,
        "candidate_neighborhoods": neighborhoods,
        "limitations": [
            "Immediate adjacency does not define an invertible block.",
            "No automatic Lagrange-multiplier classification.",
            "No Schur-complement construction or solver change.",
            "Shift-invariant graph is not required for future nonlinear systems.",
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(f"PASS MPI graph parity: {len(parity)} shifts, {a['dofs']} DoFs")
    for item in neighborhoods:
        print(f"CANDIDATE sigma={item['sigma']:g}: {item['candidate']} "
              f"neighbors={item['local_neighborhood']}")
    print("Saved:", args.output)


if __name__ == "__main__":
    main()
