# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Host comparison of #28/#6 native first-differentiated compatible spinup t0.

Insist on original physical C1 time law, correct MatNest MPI IS mapping,
true native IFunction, near-zero physical qdot0/Y0 and exact rollback.
This is not a production high-index DAE initializer or a physical TS.
"""

import argparse
import json
import math
from pathlib import Path

EXPECTED_IMAGE = "1d6d0c90115213298047b28237a7f38087563a35"
EXPECTED_DIMENSIONS = {
    .2: (3300, 1452, 1716, 132),
    .1: (10272, 4758, 5262, 252),
}


def require(ok, message):
    if not ok:
        raise ValueError(message)


def numeric(value, name):
    require(type(value) in (int, float) and math.isfinite(value),
            f"Missing or non-finite {name}: {value!r}")
    return float(value)


def norm(values):
    return math.sqrt(math.fsum(numeric(x, "vector")**2 for x in values))


def check(report, mesh, ranks):
    require(report.get("schema") ==
            "scrimp-v130-issue28-spinup-first-differentiated-v1"
            and report.get("status") ==
            "DIAGNOSTIC_SPINUP_T0_TRUE_RHS_SOLVED_NOT_INITIALIZED",
            "Wrong physical differentiated spinup report")
    require(report.get("case") == "spinup"
            and report.get("nominal_mesh") == mesh
            and report.get("mpi_ranks") == ranks
            and report.get("original_P2_geometry") is True
            and report.get("physical_wall_facet_tangent") is True
            and report.get("initial_physical_p0_zero") is True
            and report.get("spinup_duration") == .02
            and report.get("rotation_after_spinup") == .2,
            "Wrong physical spinup data, original geometry or MPI case")
    dim = EXPECTED_DIMENSIONS[mesh]
    require(tuple(report.get(x) for x in (
        "dofs", "dynamic_dofs",
        "regular_algebraic_dofs", "traction_Y_dofs",
    )) == dim,
            "Wrong original P2 finite element dimensions")
    projection = report.get("initial_projection")
    require(isinstance(projection, dict)
            and projection.get("type") == "known-exact-analytical-spinup-rest",
            "Initial compatibility is not the exact rest candidate")
    require(numeric(projection.get("native_final_Fa_l2"),
                    "initial Fa") < 1e-11,
            "Original spinup initial algebraic constraints are nonzero")
    require(report.get("independent_spinup_t0_algebraic_compatibility_qualified")
            is True
            and report.get("analytic_partial_t_Fa_at_t0_exact_zero") is True,
            "Original spinup state or exact C1 derivative not qualified")
    require(numeric(report.get("native_algebraic_partial_t_Fa_l2"),
                    "exact time derivative") == 0.,
            "The physical C1 wall starts with nonzero time derivative")

    tau = .02
    delta = 1e-4
    s = delta/tau
    g = lambda x: 3*x*x-2*x*x*x
    analytic = g(2*s)/g(s)
    measured = numeric(report.get("native_wall_smoothstep_ratio"),
                       "measured wall ratio")
    require(math.isclose(measured, analytic, rel_tol=2e-6),
            "Original native GetFEM wall does not have C1 smoothstep onset")
    require(math.isclose(
        numeric(report.get("analytic_wall_smoothstep_ratio"),
                "analytic wall ratio"),
        analytic, rel_tol=1e-12,
    ), "Wrong exact physical wall ramp definition")
    u1 = numeric(report.get("native_wall_time_change_U_delta_l2"), "U delta")
    u2 = numeric(report.get("native_wall_time_change_U_2delta_l2"),
                 "U 2delta")
    require(u1 > 1e-15 and u2 > u1
            and math.isclose(u2/u1, measured, rel_tol=1e-8)
            and numeric(report.get("native_time_change_other_rows_l2"),
                        "non-control temporal contamination") < 1e-11,
            "Physical native C1 control does not act only on port U")

    require(report.get("matnest_original_row_and_column_IS_used") is True
            and report.get("matnest_solution_reconstructed_by_block_IS") is True
            and report.get("no_pseudo_or_physical_TS") is True
            and report.get("original_q0_preserved_bitwise") is True
            and report.get("fixed_initial_h0_p0_preserved_bitwise") is True
            and report.get("full_original_IFunction_checked") is True
            and report.get("physical_rotating_Y0_accepted") is False
            and report.get("true_initialization_TS_executed") is False
            and report.get("true_physical_TS_executed") is False
            and report.get("first_differentiated_rank_proven") is False,
            "False physical/production initialization certification")
    require(report.get("bordered_dynamic_algebraic_operator_shape") ==
            [dim[0], dim[0]],
            "Malformed first-differentiated native bordered operator")
    require(type(report.get("bordered_KSP_reason")) is int
            and report["bordered_KSP_reason"] > 0,
            "Native bordered MUMPS direct solve failed")

    for field, limit in (
        ("native_Fd_at_zero_zdot_before_Y0_l2", 1e-10),
        ("native_algebraic_partial_t_Fa_l2", 1e-14),
        ("dynamic_physical_rhs_l2", 1e-10),
        ("bordered_true_rhs_l2", 1e-10),
        ("bordered_relative_true_rhs_residual", 1e-10),
        ("matnest_physical_rhs_roundtrip_l2", 1e-12),
        ("candidate_qdot_l2", 1e-8),
        ("candidate_regular_algebraic_derivative_l2", 1e-8),
        ("candidate_traction_Y0_l2", 1e-8),
        ("original_native_IFunction_Fd_with_candidate_derivative_Y0_l2", 1e-10),
        ("original_native_IFunction_Fa_with_candidate_Y0_l2", 1e-11),
    ):
        v = numeric(report.get(field), field)
        require(v >= 0 and v <= limit,
                f"Physical differentiated native gate fails: {field}={v:.9g}")

    vectors = []
    for label, field, expected in (
        ("qdot", "candidate_qdot_coefficients", dim[1]),
        ("a_regular_dot", "candidate_regular_algebraic_derivative_coefficients",
         dim[2]),
        ("Y0", "candidate_traction_Y0_coefficients", dim[3]),
    ):
        values = report.get(field)
        require(isinstance(values, list) and len(values) == expected,
                f"Incomplete native {label} coefficients")
        values = [numeric(x, label) for x in values]
        computed = norm(values)
        require(math.isclose(
            computed, numeric(report.get(
                {"qdot": "candidate_qdot_l2",
                 "a_regular_dot": "candidate_regular_algebraic_derivative_l2",
                 "Y0": "candidate_traction_Y0_l2"}[label]), label+" norm"),
            rel_tol=1e-8, abs_tol=1e-14,
        ), f"False native {label} coefficient norm")
        vectors.append(values)
    return {
        "Fa": projection["native_final_Fa_l2"],
        "physical_Fd_after_Y0":
            report["original_native_IFunction_Fd_with_candidate_derivative_Y0_l2"],
        "Y0": report["candidate_traction_Y0_l2"],
        "qdot": report["candidate_qdot_l2"],
        "wall_ratio": measured,
        "vectors": vectors,
    }


def compare(folder):
    root = Path(folder)
    manifest = json.loads((root/"manifest.json").read_text(encoding="utf-8"))
    require(manifest.get("schema") ==
            "v130-issue28-spinup-first-differentiated-v1"
            and manifest.get("image_revision") == EXPECTED_IMAGE
            and manifest.get("physical_spinup_only") is True
            and manifest.get("no_physical_time_steps") is True
            and manifest.get("no_accepted_Y0") is True,
            "Invalid original spinup differentiated gate provenance")
    summary = {}
    for mesh in (.2, .1):
        rank_reports = []
        for ranks in (1, 2):
            f = root/(
                f"result_h{str(mesh).replace('.', 'p')}_mpi{ranks}"
            )/"spinup_first_derivative.json"
            require(f.is_file(), f"Missing native spinup differentiated report {f}")
            rank_reports.append(
                check(json.loads(f.read_text(encoding="utf-8")), mesh, ranks)
            )
        r1, r2 = rank_reports
        for key in ("Fa", "physical_Fd_after_Y0", "wall_ratio", "qdot", "Y0"):
            require(math.isclose(
                numeric(r1[key], key), numeric(r2[key], key),
                rel_tol=1e-5, abs_tol=1e-11,
            ), f"Original native MPI first-differentiated mismatch in {key}")
        for xa, xb in zip(r1["vectors"], r2["vectors"]):
            require(norm([a-b for a, b in zip(xa, xb)]) <=
                    max(1e-11, 1e-5*norm(xa)),
                    "Native PETSc first-differentiated vector MPI parity failed")
        key = f"h{mesh:g}"
        summary[key] = {
            "native_MPI1_algebraic_Fa": r1["Fa"],
            "native_MPI1_physical_Fd_after_differentiation":
                r1["physical_Fd_after_Y0"],
            "candidate_qdot_l2": r1["qdot"],
            "candidate_Y0_l2": r1["Y0"],
            "native_C1_wall_ratio": r1["wall_ratio"],
            "mpi_parity": True,
        }
        print(
            f"mesh={mesh:g}: native differentiated spinup t0 MPI1/2 PASS "
            f"Fa={r1['Fa']:.9g} Fd={r1['physical_Fd_after_Y0']:.9g} "
            f"qdot={r1['qdot']:.9g} Y0={r1['Y0']:.9g}",
            flush=True
        )
    result = {
        "schema": "scrimp-v130-issue28-spinup-first-differentiated-comparison-v1",
        "status": "COMPLETE_DIFFERENTIATED_SPINUP_T0_DIAGNOSTIC_NO_TS",
        "cases": summary,
        "production_high_index_initializer_qualified": False,
        "physical_time_integration_qualified": False,
        "original_rotating_q0_accepted": False,
    }
    (root/"spinup_first_derivative_comparison.json").write_text(
        json.dumps(result, indent=2)+"\n", encoding="utf-8"
    )
    print("COMPLETE: differentiated t0 diagnostic; NO TS or Y0 acceptance.")
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--evidence", type=Path, required=True)
    args = p.parse_args()
    compare(args.evidence)


if __name__ == "__main__":
    main()
