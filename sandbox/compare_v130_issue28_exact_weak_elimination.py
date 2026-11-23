# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Verify native exact-block algebraic elimination without accepting Y0."""

import argparse
import json
import math
from pathlib import Path

EXPECTED_IMAGE = "1d6d0c90115213298047b28237a7f38087563a35"
EXPECTED_WEAK_MOMENT_NORMS = {
    .2: 3.6797242326813234e-7,
    .1: 1.2763509037199882e-7,
}


def require(ok, message):
    if not ok:
        raise ValueError(message)


def real(value, name):
    require(type(value) in (int, float) and math.isfinite(value),
            f"Invalid {name}: {value!r}")
    return float(value)


def norm(values):
    return math.sqrt(math.fsum(real(x, "coefficient")**2 for x in values))


def relative_vector_error(a, b):
    require(len(a) == len(b), "Native algebraic vector dimensions differ")
    return norm([x-y for x, y in zip(a, b)])/max(1., norm(a), norm(b))


def validate(report, h, case, ranks):
    require(report.get("schema") ==
            "scrimp-v130-issue28-exact-weak-elimination-v1"
            and report.get("status") ==
            "EXACT_BLOCK_ELIMINATION_DIAGNOSTIC_ONLY",
            "Wrong native exact-block elimination schema")
    require(report.get("mpi_ranks") == ranks
            and abs(real(report.get("nominal_mesh"), "mesh") - h) < 1e-13
            and report.get("rotation") == (0. if case == "rest" else .2),
            "Native exact elimination run mismatched manifest")
    require(report.get("original_geometry_P2") is True
            and report.get("original_p_e_p_same_FE_dimension") is True
            and report.get("original_U_Y_same_FE_dimension") is True
            and report.get("original_fixed_h0_p0_preserved_bitwise") is True
            and report.get("original_traction_Y_preserved_bitwise") is True
            and report.get("entire_state_rollback_bitwise") is True
            and report.get("no_physical_TS") is True
            and report.get("no_pseudo_or_BDF_initializer") is True
            and report.get("consistent_rotating_initial_Y0_accepted") is False
            and report.get("generic_initialization_qualified") is False
            and report.get("differential_index_proven") is False,
            "Original physical data or diagnostic-only acceptance invalid")
    for field in (
        "native_Fep_after_exact_copy_l2",
        "native_Fep_after_wall_elimination_l2",
        "native_FU_after_wall_elimination_l2",
    ):
        require(0 <= real(report.get(field), field) < 1e-10,
                f"Native constitutive/control elimination failed: {field}")
    Fy = report["native_FY_coefficients"]
    b = report["independent_native_weak_U_moment_coefficients"]
    require(isinstance(Fy, list) and isinstance(b, list)
            and len(Fy) == len(b) and len(Fy) > 0,
            "Missing independent complete native weak moment coefficients")
    sign = report.get("native_FY_equals_signed_h0_wall_moment_sign")
    require(sign in (-1, 1), "Invalid actual GetFEM weak wall sign")
    discrepancy = norm([x - sign*5.*y for x, y in zip(Fy, b)])
    require(math.isclose(discrepancy,
                         real(report.get("native_FY_weak_moment_identity_error_l2"),
                              "weak identity error"),
                         rel_tol=1e-4, abs_tol=1e-12),
            "Archived native weak moment equation identity falsified")
    # Mirror the native probe's machine-roundoff-aware criterion;
    # use the archived old-run evidence even when the new optional
    # tolerance fields do not exist (safe resume of the six passed runs).
    identity_limit = 5e-15 + 2e-9 * (5.*norm(b))
    require(discrepancy <= identity_limit,
            "Native F_Y does not equal signed 5*original weak wall moment")
    if "native_FY_weak_moment_identity_tolerance_l2" in report:
        require(math.isclose(
            real(report["native_FY_weak_moment_identity_tolerance_l2"],
                 "recorded native identity tolerance"),
            identity_limit, rel_tol=1e-9, abs_tol=1e-18
        ), "Native probe and host comparator identity tolerances differ")
    if "native_FY_weak_moment_identity_relative_error" in report:
        require(math.isclose(
            real(report["native_FY_weak_moment_identity_relative_error"],
                 "recorded relative identity discrepancy"),
            discrepancy / max(5.*norm(b), 1e-30),
            rel_tol=1e-4, abs_tol=1e-12
        ), "Reported native relative identity discrepancy mismatches vectors")
    y_norm = norm(Fy)
    require(math.isclose(y_norm,
                         real(report.get("native_FY_after_both_exact_eliminations_l2"),
                              "native F_Y norm"),
                         rel_tol=1e-8, abs_tol=1e-12),
            "Native F_Y norm incompatible with actual full coefficient vector")
    expected_norm = real(
        report.get("expected_signed_h0_wall_moment_l2"), "expected moment norm"
    )
    require(math.isclose(expected_norm, 5.*norm(b),
                         rel_tol=1e-8, abs_tol=1e-12),
            "Expected signed original weak U moment norm falsified")
    if case == "rest":
        require(y_norm < 1e-10 and norm(b) < 1e-10,
                "Analytic exact rest must be algebraically compatible")
    else:
        require(norm(b) > 1e-10 and y_norm > 1e-10,
                "Rotating original weak obstruction unexpectedly absent")
        require(math.isclose(
            norm(b), EXPECTED_WEAK_MOMENT_NORMS[h],
            rel_tol=1e-4, abs_tol=1e-10),
            "Original prescribed physical wall moment changed since prior native gate")
    return {"Fy": Fy, "b": b, "Fy_l2": y_norm, "b_l2": norm(b),
            "ep_l2": report["native_Fep_after_wall_elimination_l2"],
            "U_l2": report["native_FU_after_wall_elimination_l2"],
            "weak_identity_error_l2": discrepancy,
            "sign": sign}


def compare(evidence):
    root = Path(evidence)
    manifest = json.loads((root/"manifest.json").read_text(encoding="utf-8"))
    require(manifest.get("schema") == "v130-issue28-exact-weak-elimination-v1"
            and manifest.get("image_revision") == EXPECTED_IMAGE
            and manifest.get("original_fixed_h0_p0") is True
            and manifest.get("no_physical_time_steps") is True
            and manifest.get("accepted_Y0") is False,
            "Wrong exact native algebraic elimination campaign")
    summary = {}
    for h in (.2, .1):
        for case in ("rest", "rotating"):
            per_rank = []
            for ranks in (1, 2):
                file = root/(
                    f"result_h{str(h).replace('.', 'p')}_{case}_{ranks}"
                )/"exact_weak_elimination.json"
                require(file.is_file(), f"Missing native elimination report {file}")
                per_rank.append(validate(
                    json.loads(file.read_text(encoding="utf-8")),
                    h, case, ranks
                ))
            first, second = per_rank
            require(first["sign"] == second["sign"],
                    "MPI changed native weak-flow/effort sign")
            require(relative_vector_error(first["Fy"], second["Fy"]) < 1e-9
                    and relative_vector_error(first["b"], second["b"]) < 1e-9,
                    "Original full weak residual/moment MPI vector parity failed")
            key = f"h{h:g}_{case}"
            summary[key] = {
                "native_F_Y_after_exact_e_p_U_elimination_l2":
                    first["Fy_l2"],
                "original_weak_U_moment_l2": first["b_l2"],
                "native_FY_equals_signed_h0_b": True,
                "original_F_ep_zero": first["ep_l2"] < 1e-10,
                "original_F_U_zero": first["U_l2"] < 1e-10,
                "mpi_full_vector_parity": True,
            }
            print(
                f"h={h:g} {case}: exact F_ep/F_U PASS, "
                f"F_Y={first['Fy_l2']:.9g} "
                f"independent 5*b={5*first['b_l2']:.9g}; "
                "MPI1/2 parity PASS", flush=True
            )
    result = {
        "schema": "scrimp-v130-issue28-exact-weak-elimination-comparison-v1",
        "status": "COMPLETE_BLOCK_ELIMINATION_DIAGNOSTIC_NOT_INITIALIZER",
        "cases": summary,
        "original_rotating_full_DAE_Y0_initialized": False,
        "differential_index_proven": False,
        "no_time_stepping": True,
        "interpretation": (
            "At original fixed rotating P2 q0 with h0=5, exact "
            "native F_ep and F_U block elimination leaves an "
            "independently corroborated nonzero native F_Y weak "
            "boundary constraint residual. This is conditional on "
            "nondegenerate same-basis FE mass operators; it does not "
            "initialize a physical rotating DAE or classify its index."
        ),
    }
    (root/"exact_weak_elimination_comparison.json").write_text(
        json.dumps(result, indent=2)+"\n", encoding="utf-8"
    )
    print("COMPLETE: exact algebraic block diagnosis, NO initial Y0 accepted.")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--evidence", type=Path, required=True)
    args = parser.parse_args()
    compare(args.evidence)


if __name__ == "__main__":
    main()
