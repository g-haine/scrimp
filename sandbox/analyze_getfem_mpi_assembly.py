# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Classify unmodified GetFEM MPI assembly via analytical integrals.

Use JSON results from diagnose_getfem_mpi_assembly.py for P=1,2,4.
No corrections are applied. A complete vector on each MPI rank is
EXPECTED; its integral must not be multiplied by the process count.

Reports separately:
  generic nonlinear RHS factor, source-brick factor,
  mixed RHS consistency, tangent matrix distribution/replication.
Unknown / inconsistent conventions fail; known anomalous P-fold
generic contributions are reported, NOT hidden.
"""
import argparse
import json
import math
from pathlib import Path


def close(a, b):
    return math.isclose(a, b, rel_tol=2e-7, abs_tol=2e-7)


def classify_matrix(samples, exact, count):
    values = [float(row["raw_matrix_ones_integral"]) for row in samples]
    if all(close(v, exact) for v in values):
        return "replicated_full"
    if close(sum(values), exact):
        return "distributed_contributions"
    return "unexpected"


def classify_rhs(observed, exact, count):
    if close(observed, exact):
        return "physical"
    if close(observed, exact * count):
        return "multiplied_by_rank_count"
    return "unexpected"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reports", nargs="+", type=Path)
    parser.add_argument(
        "--matrix-mode", choices=("replicated", "distributed"),
        default="replicated",
        help="Configured SCRIMP_GETFEM_MATRIX_MODE for the tested build",
    )
    parser.add_argument(
        "--require-physical-rhs", action="store_true",
        help="Require native generic/source RHS to be unscaled (GetFEM 5.5)",
    )
    args = parser.parse_args()
    reports = [json.loads(p.read_text(encoding="utf-8")) for p in args.reports]
    sizes = [int(r["mpi_size"]) for r in reports]
    if len(set(sizes)) != len(sizes) or 1 not in sizes:
        parser.error("Require unique MPI sizes including P=1")
    if len(reports) < 2:
        parser.error("At least one P>1 observation is needed")

    incompatible = []
    print("| MPI | Generic RHS | Source RHS | Mixed RHS | Linear matrix | Nonlinear tangent |")
    print("|---:|---|---|---|---|---|")
    for r in sorted(reports, key=lambda x: x["mpi_size"]):
        p = int(r["mpi_size"])
        if int(r["getfem_mpi_level"]) != 2:
            incompatible.append(f"P={p}: GetFEM was not compiled at MPI level 2")
        cases = {
            c["kind"]: c["per_rank"] for c in r["cases"]
        }
        expected_cases = {
            "generic", "source", "mixed", "linear_tangent",
            "nonlinear_tangent"
        }
        if set(cases) != expected_cases:
            incompatible.append(f"P={p}: missing native assembly cases")
            continue

        scalar = {}
        for name, expected in (
            ("generic", -1.0),
            ("source", 2.0),
            ("mixed", 1.0),
        ):
            rows = cases[name]
            if len(rows) != p:
                incompatible.append(f"P={p}: {name} has {len(rows)} rank reports")
                continue
            if any(row["raw_rhs_length"] != row["dofs"] for row in rows):
                incompatible.append(
                    f"P={p}: {name} returned non-global RHS vectors"
                )
            if any(row["rhs_copy_max_difference"] is None
                   or row["rhs_copy_max_difference"] > 1e-8 for row in rows):
                incompatible.append(
                    f"P={p}: {name} RHS is not replicated consistently"
                )
            vals = [float(x["raw_rhs_integral"]) for x in rows]
            if not all(close(x, vals[0]) for x in vals):
                incompatible.append(
                    f"P={p}: {name} RHS integral differs across ranks: {vals}"
                )
            scalar[name] = vals[0]

        if len(scalar) != 3:
            continue
        g = classify_rhs(scalar["generic"], -1.0, p)
        s = classify_rhs(scalar["source"], 2.0, p)
        # Full model should be additive in raw GetFEM assembly, even
        # if generic/source contributions have different MPI scales.
        m = ("additive" if close(
            scalar["mixed"],
            scalar["generic"] + scalar["source"]
        ) else "NOT_ADDITIVE")
        lin = classify_matrix(cases["linear_tangent"], 1.0, p)
        nonlin = classify_matrix(cases["nonlinear_tangent"], 2.0, p)
        print(f"| {p} | {g} | {s} | {m} | {lin} | {nonlin} |")
        if (g == "unexpected" or s == "unexpected"
                or m != "additive"
                or lin == "unexpected" or nonlin == "unexpected"):
            incompatible.append(
                f"P={p}: unrecognized GetFEM assembly semantics"
            )
        if p > 1:
            observed = {lin, nonlin}
            required = (
                "distributed_contributions"
                if args.matrix_mode == "distributed" else "replicated_full"
            )
            if observed != {required}:
                incompatible.append(
                    f"P={p}: native GetFEM matrices={observed}, "
                    f"configured matrix mode={args.matrix_mode}; "
                    f"expected {required}. Select the correct PETSc bridge."
                )
        if args.require_physical_rhs and (g != "physical" or s != "physical"):
            incompatible.append(
                f"P={p}: native RHS not physical (generic={g}, "
                f"source={s}); do not assume GetFEM 5.5 is unscaled"
            )
        if g != s and p > 1:
            print(
                f"  P={p}: RHS generic/source scaling DIFFERS: "
                "do not divide the whole mixed RHS uniformly. "
                "Use separately calibrated source treatment."
            )

    if incompatible:
        print("\nFAIL: " + "; ".join(incompatible))
        raise SystemExit(1)
    print(
        "\nNative assembly conventions classified. "
        "P-fold factors, when present, are reported, not corrected. "
        "Confirm calibrated SCRIMP source handling and physical "
        "P=1/P>1 equivalence independently."
    )


if __name__ == "__main__":
    main()
