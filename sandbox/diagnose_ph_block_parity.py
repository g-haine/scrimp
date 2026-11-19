# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Explain *where* a pH-block MPI parity check failed, without rerunning MPI.

Reads the original 1/2-rank block-structure JSONs, including older
reports without mesh hashes. Norm disagreement cannot tell whether
GetFEM assembly, geometry, PETSc insertion or ordering caused it.
This tool therefore NEVER certifies physics or alters parity tolerances.

Example:
  python sandbox/diagnose_ph_block_parity.py \
      sandbox/outputs/ph-block-audit-wave_coenergy-h0p12-.../
"""
import argparse
import json
import math
from pathlib import Path


def _float(value, key):
    v = float(value)
    if not math.isfinite(v):
        raise ValueError(f"Nonfinite {key}: {value}")
    return v


def diagnose(folder, *, rtol=2e-7, atol=1e-9, top=12):
    folder = Path(folder)
    a = json.loads((folder / "blocks_1r.json").read_text(encoding="utf-8"))
    b = json.loads((folder / "blocks_2r.json").read_text(encoding="utf-8"))
    for key in ("dofs", "mesh", "fields", "physical_model"):
        if a[key] != b[key]:
            raise ValueError(f"Cannot compare different {key} across ranks")
    fields = [f["name"] for f in a["fields"]]
    if not fields:
        raise ValueError("The audit contains no physical fields")
    if set(a["assembled_matrices"]) != set(b["assembled_matrices"]):
        raise ValueError("The audit matrix sets differ")

    output = {
        "schema": "scrimp-ph-block-discrepancy-v1",
        "source_run": str(folder),
        "dofs": a["dofs"],
        "physical_model": a["physical_model"],
        "rtol": rtol, "atol": atol,
        "mesh_sha256_1r": a.get("mesh_sha256"),
        "mesh_sha256_2r": b.get("mesh_sha256"),
        "matrices": {},
    }
    print(f"pH block discrepancy, {a['dofs']} DoFs, "
          f"{a['physical_model']}", flush=True)
    if output["mesh_sha256_1r"] is None or output["mesh_sha256_2r"] is None:
        print("Mesh fingerprint: UNRECORDED in original reports; "
              "mesh identity not established.", flush=True)
    elif output["mesh_sha256_1r"] != output["mesh_sha256_2r"]:
        print("MESH FILES DIFFER: first test geometry/ordering before "
              "attributing matrix differences to MPI!", flush=True)
    else:
        print("Mesh file SHA256 matches across 1/2-rank runs.", flush=True)

    for label, m1 in a["assembled_matrices"].items():
        m2 = b["assembled_matrices"][label]
        nnz_a = int(m1["nnz_effective"])
        nnz_b = int(m2["nnz_effective"])
        entries = []
        for norm_key, nnz_key in (
            ("block_frobenius_norm", "block_effective_nnz"),
        ):
            matrix1, matrix2 = m1[norm_key], m2[norm_key]
            count1, count2 = m1[nnz_key], m2[nnz_key]
            n = len(fields)
            if (len(matrix1) != n or len(matrix2) != n
                    or len(count1) != n or len(count2) != n):
                raise ValueError(f"{label}: invalid matrix block dimensions")
            for i, row in enumerate(matrix1):
                if (len(row) != n or len(matrix2[i]) != n
                        or len(count1[i]) != n or len(count2[i]) != n):
                    raise ValueError(f"{label}: incorrect row dimensions")
                for j, v1 in enumerate(row):
                    v2 = _float(matrix2[i][j], label)
                    v1 = _float(v1, label)
                    abs_diff = abs(v1 - v2)
                    scale = max(abs(v1), abs(v2), 1e-300)
                    is_mismatch = abs_diff > atol + rtol*abs(v2)
                    entries.append({
                        "row": fields[i], "column": fields[j],
                        "frobenius_1r": v1, "frobenius_2r": v2,
                        "abs_diff": abs_diff,
                        "relative_diff": abs_diff/scale,
                        "ratio_2r_over_1r": (
                            v2/v1 if abs(v1) > 1e-300 else None
                        ),
                        "nonzero_count_1r": int(count1[i][j]),
                        "nonzero_count_2r": int(count2[i][j]),
                        "mismatch": is_mismatch,
                    })
        failures = [entry for entry in entries if entry["mismatch"]]
        by_abs = sorted(failures, key=lambda x:x["abs_diff"], reverse=True)
        output["matrices"][label] = {
            "total_nnz_1r": nnz_a,
            "total_nnz_2r": nnz_b,
            "mismatching_blocks": len(failures),
            "total_blocks": len(entries),
            "mismatches": by_abs,
            "all_block_norms": entries,
        }
        print(f"\n{label}", flush=True)
        print(f"  PETSc nonzeros: 1r={nnz_a:,} 2r={nnz_b:,}; "
              f"norm disagreements={len(failures)}/{len(entries)}", flush=True)
        for entry in by_abs[:top]:
            ratio = entry["ratio_2r_over_1r"]
            print(
                f"  {entry['row']:>5} <- {entry['column']:<5} : "
                f"{entry['frobenius_1r']:.9g} (1r) vs "
                f"{entry['frobenius_2r']:.9g} (2r), "
                f"abs Δ={entry['abs_diff']:.5g}, "
                f"rel Δ={entry['relative_diff']:.3g}, "
                f"ratio={ratio if ratio is not None else 'n/a'}; "
                f"nnz {entry['nonzero_count_1r']}/"
                f"{entry['nonzero_count_2r']}",
                flush=True,
            )

    output_path = folder / "ph_block_discrepancy.json"
    output_path.write_text(
        json.dumps(output, indent=2) + "\n", encoding="utf-8"
    )
    print(f"\nDetailed non-destructive report: {output_path}", flush=True)
    print("WARNING: equal block Frobenius norms/nnz do not prove "
          "entrywise equality; this report does not diagnose cause.",
          flush=True)
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("folder", type=Path)
    parser.add_argument("--top", type=int, default=12)
    args = parser.parse_args()
    if args.top < 1:
        parser.error("--top must be a positive integer")
    diagnose(args.folder, top=args.top)


if __name__ == "__main__":
    main()
