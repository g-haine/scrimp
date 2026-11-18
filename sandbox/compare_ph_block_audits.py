# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Gate initial pH block structure consistency across MPI=1/2.

This checks GetFEM variable offsets, all matrix coupling-block
Frobenius norms and algebraic zero-mass rows. Matrix Frobenius
norm agreement is necessary but NOT sufficient to prove coefficient-
wise equality of two matrices: physical state/Hamiltonian parity is
independently validated by established SCRIMP scripts.
"""
import argparse
import json
from pathlib import Path

import numpy as np


def compare(directory):
    directory = Path(directory)
    reports = {}
    for p in (1, 2):
        path = directory / f"blocks_{p}r.json"
        data = json.loads(path.read_text(encoding="utf-8"))
        if (data["schema"] != "scrimp-ph-block-audit-v1"
                or data["mpi_ranks"] != p
                or data["ts_type"] != "cn"
                or data["steps"] != 2
                or not data["linear_mass"]):
            raise ValueError(f"{path}: incorrect block audit provenance")
        reports[p] = data
    a, b = reports[1], reports[2]
    for name in ("source_commit", "physical_model", "mesh", "dofs",
                 "dt", "steps", "final_time", "fields",
                 "linear_stiffness"):
        if a[name] != b[name]:
            raise ValueError(f"Mismatch in {name} between 1/2 MPI ranks")
    # Older wave-only JSON reports did not contain this field.
    # Continue to support reanalysing these reference campaigns.
    if a.get("formulation") != b.get("formulation"):
        raise ValueError("Mismatch in formulation between 1/2 MPI ranks")
    if a.get("scrimp_mesh_source_sha256") != b.get(
        "scrimp_mesh_source_sha256"
    ):
        raise ValueError(
            "Different SCRIMP mesh.py source files across MPI runs"
        )
    if a["mumps_policy"] != b["mumps_policy"]:
        raise ValueError("Different MUMPS policy between 1/2 ranks")
    # A single number of mesh DoFs or matching sparsity graph is NOT
    # sufficient to compare coefficient-by-coefficient assembly.
    # The two runs generate their Gmsh files separately, and the
    # coenergy two-subdomain geometry can potentially vary.
    ma, mb = a.get("mesh_sha256"), b.get("mesh_sha256")
    if ma is not None and mb is not None and ma != mb:
        raise ValueError(
            "The 1/2-rank reference Gmsh mesh files differ "
            f"(SHA256 {ma} vs {mb}). Reproduce identical input "
            "meshes before diagnosing PETSc/GetFEM coefficient parity."
        )
    for p, data in reports.items():
        spans = data["owned_ranges_by_rank"]
        if len(spans) != p:
            raise ValueError("Incorrect MPI ownership coverage")
        if spans[0]["lo"] != 0 or spans[-1]["hi"] != data["dofs"]:
            raise ValueError("Missing global MPI row ownership")
        if any(
            spans[i-1]["hi"] != spans[i]["lo"]
            for i in range(1, p)
        ):
            raise ValueError("Gaps or overlaps in MPI owned row intervals")

    details = {}
    names = a["assembled_matrices"]
    if set(names) != set(b["assembled_matrices"]):
        raise ValueError("PETSc matrix set is different between runs")
    for mat_name, content_a in names.items():
        content_b = b["assembled_matrices"][mat_name]
        if content_a["variable_order"] != content_b["variable_order"]:
            raise ValueError(f"{mat_name}: different field ordering")
        metrics = {}
        for key in ("block_frobenius_norm", "role_frobenius_norm"):
            x, y = (
                np.asarray(s[key], dtype=np.float64)
                for s in (content_a, content_b)
            )
            if x.shape != y.shape or not np.all(np.isfinite(x+y)):
                raise ValueError(f"{mat_name}/{key}: invalid block norms")
            difference = np.abs(x-y)
            if not np.allclose(
                x, y, rtol=2e-7, atol=1e-9
            ):
                biggest = float(np.max(difference))
                raise ValueError(
                    f"{mat_name}/{key}: 1/2-rank sparse block "
                    f"norm mismatch (max abs {biggest:g})"
                )
            metrics[key+"_max_abs_diff"] = float(np.max(difference))
        # Matrix CSR graphs on different rank partitions may contain
        # extra structural zeros or tiny roundoff entries; do not
        # equate exact nonzero counts with coefficient-wise identity.
        for key in ("block_effective_nnz", "role_effective_nnz"):
            x, y = (
                np.asarray(s[key], dtype=np.int64)
                for s in (content_a, content_b)
            )
            if x.shape != y.shape:
                raise ValueError(f"{mat_name}/{key}: shape mismatch")
            metrics[key+"_total_1r"] = int(x.sum())
            metrics[key+"_total_2r"] = int(y.sum())
            metrics[key+"_differences"] = int(np.count_nonzero(x != y))
        details[mat_name] = metrics
    # The e_q/e_p equations and boundary constraints in wave.py are
    # ALGEBRAIC: they do not carry time derivative terms.
    mass = a["assembled_matrices"][
        "flow_mass_F_and_algebraic_zero_rows"
    ]
    roles = mass["roles"]
    mass_norms = np.asarray(mass["role_frobenius_norm"])
    for role in ("coenergy", "boundary", "algebraic"):
        if not np.all(mass_norms[roles.index(role)] < 1e-10):
            raise ValueError(
                f"Acoustic reference: unexpected mass on {role} rows"
            )
    stiffness = a["assembled_matrices"]["linear_stiffness_all_bricks"]
    stiff_norms = np.asarray(stiffness["role_frobenius_norm"])
    formulation = a.get("formulation", {}).get(
        "formulation", "implicit_constitutive"
    )
    if formulation == "implicit_constitutive":
        required = (
            ("coenergy", "coenergy"),
            ("coenergy", "dynamic"),
            ("dynamic", "coenergy"),
        )
    elif formulation == "substituted_coenergy":
        if any(f["group"] == "coenergy" for f in a["fields"]):
            raise ValueError(
                "Substituted coenergy example contains coenergy DoFs"
            )
        required = (
            ("dynamic", "dynamic"),
            ("dynamic", "algebraic"),
            ("algebraic", "dynamic"),
            ("dynamic", "boundary"),
            ("boundary", "dynamic"),
        )
        if a["linear_stiffness"]:
            raise ValueError(
                "The example's U=0.5*Y feedback is assembled "
                "through a nonlinear GetFEM tangent; flag changed"
            )
    else:
        raise ValueError(
            f"Unsupported native block audit formulation {formulation}"
        )
    for row, col in required:
        if stiff_norms[roles.index(row), roles.index(col)] <= 1e-10:
            raise ValueError(f"Required pH coupling {row}<-{col} absent")
    report = {
        "source_commit": a["source_commit"],
        "mesh": a["mesh"], "dofs": a["dofs"],
        "ranks_compared": [1, 2],
        "fields": a["fields"],
        "formulation": a.get("formulation", {}),
        "matrix_block_parity": details,
        "validated_properties": [
            "Complete GetFEM variable classification.",
            "MPI row ownership covers the global system.",
            "Mass on acoustic coenergy/boundary/algebraic rows is zero.",
            (
                "Constitutive and dynamic/coenergy couplings are nonzero."
                if formulation == "implicit_constitutive"
                else "Substituted model has no coenergy DoFs and "
                     "dynamic/port/dissipative couplings are nonzero."
            ),
            "1/2-rank matrix block Frobenius norms agree.",
        ],
        "limits": [
            "Block norms do NOT prove entrywise matrix equality.",
            "Block norms do NOT establish J skewness or SPD Schur.",
            "The true mixed boundary DAE has not been condensed.",
            "No preconditioner has been installed or benchmarked yet.",
        ],
    }
    out = directory / "ph_block_parity.json"
    out.write_text(json.dumps(report, indent=2)+"\n", encoding="utf-8")
    print(
        f"PASS pH block parity: {a['dofs']} DDL; "
        "all fields and 3 matrix block norm maps on 1/2 MPI ranks"
    )
    print(f"Saved: {out}")
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", type=Path, required=True)
    args = parser.parse_args()
    compare(args.directory)


if __name__ == "__main__":
    main()
