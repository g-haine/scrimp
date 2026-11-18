# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Pinpoint whether an MPI matrix mismatch arises before or after COO.

Read per-rank raw native GetFEM CSR and final PETSc owned-row CSR
recorded by audit_ph_blocks.py --capture-matrices. Reconstruct
full sparse matrices ONLY in this OFFLINE small-mesh diagnostic,
not inside the MPI production code. No factor-of-ranks correction.

A 1r-vs-2r raw GetFEM difference reveals a difference in native
assemblies (or geometry/numbering). An in-run difference between
summed GetFEM raw data and PETSc indicates a bridge or indexing error.
Compare the exact Gmsh mesh SHA256 before interpreting cross-run
differences. The output is a diagnostic, not an operator-parity pass.
"""
import argparse
import json
from pathlib import Path

import numpy as np
import scipy.sparse as sp


def load_raw(path):
    with np.load(path, allow_pickle=False) as z:
        nrows, ncols = map(int, z["shape"])
        matrix = sp.csr_matrix(
            (z["data"], z["indices"], z["indptr"]),
            shape=(nrows, ncols),
        )
    matrix.sum_duplicates()
    return matrix


def load_owned(folder, label, n):
    rows = []
    offset = 0
    for rank in range(n):
        path = folder / f"petsc_{label}_{n}r_rank{rank}.npz"
        with np.load(path, allow_pickle=False) as z:
            lo, hi = int(z["lo"]), int(z["hi"])
            shape = tuple(map(int, z["shape"]))
            if lo != offset or hi <= lo or shape[0] != shape[1]:
                raise ValueError(
                    f"{path}: invalid rank ownership {lo}:{hi}"
                )
            block = sp.csr_matrix(
                (z["data"], z["indices"], z["indptr"]),
                shape=(hi-lo, shape[1]),
            )
        block.sum_duplicates()
        rows.append(block)
        offset = hi
    matrix = sp.vstack(rows, format="csr")
    if matrix.shape != shape or matrix.shape[0] != offset:
        raise ValueError(f"PETSc rank ownership incomplete for {label}")
    return matrix


def metrics(test, reference):
    if test.shape != reference.shape:
        raise ValueError("Cannot compare different matrix sizes")
    d = (test - reference).tocsr()
    d.eliminate_zeros()
    ref_norm = float(sp.linalg.norm(reference))
    diff_norm = float(sp.linalg.norm(d))
    return {
        "reference_norm_fro": ref_norm,
        "test_norm_fro": float(sp.linalg.norm(test)),
        "difference_norm_fro": diff_norm,
        "relative_fro_error": diff_norm / max(ref_norm, 1e-300),
        "different_coefficient_locations": int(d.nnz),
        "reference_nnz": int(reference.nnz),
        "test_nnz": int(test.nnz),
    }


def feedback_block_metrics(stiffness, jacobian, fields):
    """Check dR_U/dY = +0.5 M_U from U=0.5*Y, coefficientwise.

    The linear SCRIMP control row is -M_U U. The GetFEM
    nonlinear tangent should therefore satisfy J_UY = -0.5*K_UU
    (and J_UU=K_UU) for this particular example, assuming the
    same boundary FE for U/Y. This is a model-specific audit,
    not a universal rule for user-defined feedback laws.
    """
    fields = {f["name"]: f for f in fields}
    if "U" not in fields or "Y" not in fields:
        return None
    u, y = fields["U"], fields["Y"]
    su = slice(int(u["start"]), int(u["stop"]))
    sy = slice(int(y["start"]), int(y["stop"]))
    if u["ndofs"] != y["ndofs"]:
        raise ValueError("Feedback U/Y must have equal size in this example")
    kuu = stiffness[su, su]
    juu = jacobian[su, su]
    juy = jacobian[su, sy]
    return {
        "feedback_full_J_UY_vs_minus_half_K_UU": metrics(
            juy, -0.5*kuu
        ),
        "control_full_J_UU_vs_linear_K_UU": metrics(
            juu, kuu
        ),
        "K_UU_frobenius_norm": float(sp.linalg.norm(kuu)),
        "J_UY_frobenius_norm": float(sp.linalg.norm(juy)),
    }


def analyze(folder):
    folder = Path(folder)
    reports = {
        n: json.loads(
            (folder / f"blocks_{n}r.json").read_text(encoding="utf-8")
        )
        for n in (1, 2)
    }
    a, b = reports[1], reports[2]
    for k in ("dofs", "fields", "physical_model", "mesh"):
        if a[k] != b[k]:
            raise ValueError(f"Different {k} in 1/2-rank assemblies")
    if not (a.get("capture_matrices") and b.get("capture_matrices")):
        raise ValueError(
            "Input reports lack --capture-matrices provenance"
        )
    same_mesh = (
        a.get("mesh_sha256") is not None
        and a["mesh_sha256"] == b.get("mesh_sha256")
    )
    if not same_mesh:
        print(
            "WARNING: mesh SHA256 differs/is missing; raw 1/2 "
            "rank comparisons cannot establish assembly parity",
            flush=True,
        )
    result = {
        "schema": "scrimp-getfem-bridge-localization-v1",
        "dofs": a["dofs"], "mesh": a["mesh"],
        "model": a["physical_model"],
        "same_mesh_sha256": same_mesh,
        "meshes": {
            "1rank": a.get("mesh_sha256"),
            "2rank": b.get("mesh_sha256"),
        },
        "getfem": {
            "1rank": {
                "module_path": a.get("getfem_module_path"),
                "mpi_level": a.get("getfem_mpi_level"),
            },
            "2rank": {
                "module_path": b.get("getfem_module_path"),
                "mpi_level": b.get("getfem_mpi_level"),
            },
        },
        "matrices": {},
        "interpretation": (
            "Compare raw native GetFEM contributions across rank "
            "counts and raw-sum vs PETSc within each rank count. "
            "No numerical scaling correction is applied."
        ),
    }
    fields = a["fields"]
    print(
        f"Native GetFEM vs PETSc COO MPI diagnostic: "
        f"{a['dofs']} DoFs, mesh SHA256 same={same_mesh}",
        flush=True,
    )
    # Additional *independent* check of the physical U=0.5 Y
    # feedback: compare the final TS Jacobian with its analytic
    # weak-form tangent; do not infer correctness from block norms.
    feedback = {}
    for ranks in (1, 2):
        jpath = folder / f"petsc_jacobian_{ranks}r_rank0.npz"
        if not jpath.exists():
            continue  # Backward compatibility with older archives.
        k = load_owned(folder, "stiffness", ranks)
        j = load_owned(folder, "jacobian", ranks)
        feedback[str(ranks)] = feedback_block_metrics(k, j, fields)
        if feedback[str(ranks)] is not None:
            error = feedback[str(ranks)][
                "feedback_full_J_UY_vs_minus_half_K_UU"
            ]["relative_fro_error"]
            print(
                f"  Feedback U=0.5*Y, {ranks}r: "
                f"exact coefficientwise tangent relative error "
                f"{error:.5g}",
                flush=True,
            )
    result["feedback_coefficientwise_by_ranks"] = feedback
    if len(feedback) == 2:
        j1 = load_owned(folder, "jacobian", 1)
        j2 = load_owned(folder, "jacobian", 2)
        u, y = (
            next(f for f in fields if f["name"] == name)
            for name in ("U", "Y")
        )
        sr = slice(int(u["start"]), int(u["stop"]))
        sc = slice(int(y["start"]), int(y["stop"]))
        result["feedback_2r_vs_1r"] = metrics(
            j2[sr, sc], j1[sr, sc]
        )
    for label in ("mass", "stiffness"):
        native_1 = load_raw(
            folder / f"raw_getfem_{label}_1r_rank0.npz"
        )
        native_rank0 = load_raw(
            folder / f"raw_getfem_{label}_2r_rank0.npz"
        )
        native_rank1 = load_raw(
            folder / f"raw_getfem_{label}_2r_rank1.npz"
        )
        native_2_sum = native_rank0 + native_rank1
        petsc_1 = load_owned(folder, label, 1)
        petsc_2 = load_owned(folder, label, 2)
        data = {
            "native_sum_2r_vs_native_1r":
                metrics(native_2_sum, native_1),
            "petsc_2r_vs_petsc_1r":
                metrics(petsc_2, petsc_1),
            "petsc_1r_vs_native_1r":
                metrics(petsc_1, native_1),
            "petsc_2r_vs_native_sum_2r":
                metrics(petsc_2, native_2_sum),
            "raw_2r_rank0_vs_raw_1r":
                metrics(native_rank0, native_1),
            "raw_2r_rank1_vs_raw_1r":
                metrics(native_rank1, native_1),
            "selected_fields": {},
        }
        for name in ("q", "p", "e_r", "U", "Y"):
            matches = [f for f in fields if f["name"] == name]
            if not matches:
                continue
            f = matches[0]
            sl = slice(int(f["start"]), int(f["stop"]))
            block1 = native_1[sl, sl]
            block2 = native_2_sum[sl, sl]
            block_pet2 = petsc_2[sl, sl]
            data["selected_fields"][name] = {
                "native_2r_vs_1r": metrics(block2, block1),
                "petsc_2r_vs_native_2r": metrics(block_pet2, block2),
            }
        result["matrices"][label] = data
        print(f"\n{label}:")
        for key in (
            "native_sum_2r_vs_native_1r",
            "petsc_2r_vs_petsc_1r",
            "petsc_1r_vs_native_1r",
            "petsc_2r_vs_native_sum_2r",
        ):
            m = data[key]
            print(
                f"  {key}: relative Fro error "
                f"{m['relative_fro_error']:.6g}, "
                f"nnz(diff)={m['different_coefficient_locations']}",
                flush=True,
            )
        for field, fld in data["selected_fields"].items():
            m = fld["native_2r_vs_1r"]
            if m["relative_fro_error"] > 1e-9:
                print(
                    f"  {field}<-{field}: native GetFEM 2r vs 1r "
                    f"relative Fro error {m['relative_fro_error']:.6g}",
                    flush=True,
                )
    path = folder / "native_getfem_vs_petsc_bridge.json"
    path.write_text(
        json.dumps(result, indent=2) + "\n", encoding="utf-8"
    )
    print(f"Saved: {path}", flush=True)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    args = parser.parse_args()
    analyze(args.directory)


if __name__ == "__main__":
    main()
