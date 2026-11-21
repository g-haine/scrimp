# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See LICENSE for license information.

"""Host-only strict 8-run comparator for #28/#6 native startup gate.

Two separate claims: the CURRENT opt-in nonlinear_affine initializer
refuses high-index layouts without pseudo stepping, while a deliberately
distinct rest-to-spinup physical experiment has an exact stationary
t=0 DAE *candidate*. Neither claim establishes a physical trajectory
or a new production high-index initializer.
"""

import argparse
import json
import math
from pathlib import Path

OCI_REVISION = "1d6d0c90115213298047b28237a7f38087563a35"
CASES = ("compatible_rest_spinup", "incompatible_rotation")


def require(ok, message):
    if not ok:
        raise ValueError(message)


def real(value, label):
    require(type(value) in (int, float) and math.isfinite(value),
            f"Nonfinite or missing {label}: {value!r}")
    return float(value)


def check(report, case, mesh, ranks):
    require(report.get("schema") == "scrimp-v130-issue28-startup-gate-v1",
            "Wrong native startup schema")
    require(report.get("nominal_mesh") == mesh
            and report.get("mpi_ranks") == ranks
            and report.get("geometry_order") == 2
            and report.get("height_order") == 2,
            "Incorrect original P2 model or MPI configuration")
    require(report.get("original_differential_h0_p0_preserved_bitwise") is True
            and report.get("entire_GetFEM_and_GlobalVec_rollback_bitwise") is True
            and report.get("no_pseudo_or_physical_TS") is True
            and report.get("physical_time_evolution_qualified") is False
            and report.get("production_nonlinear_initializer_qualified") is False
            and report.get("physical_rotating_Y0_accepted") is False,
            "Unqualified TS/initial state acceptance, or physical-state mutation")
    if case == "compatible_rest_spinup":
        require(report.get("status") ==
                "COMPATIBLE_REST_SPINUP_T0_DIAGNOSTIC_ONLY",
                "Native positive control did not pass")
        require(report.get("rotation_after_spinup") == .2
                and report.get("spinup_duration") == .02
                and report.get("original_h0") == 5.
                and report.get("original_p0_zero") is True
                and report.get("original_physical_wall_facet_tangent") is True
                and report.get("initial_wall_rotation_exact_zero") is True
                and report.get("initial_wall_rotation_derivative_zero") is True,
                "Compatible physical positive control was changed")
        require(report.get("initial_algebraic_coenergy_e_h") == .5
                and report.get("initial_algebraic_normal_traction_Y_n") == -.5,
                "Initial analytical rest FE fields were changed")
        metrics = {}
        for key in (
            "native_IFunction_t0_Fd_l2", "native_IFunction_t0_Fa_l2",
            "native_IFunction_t0_full_l2",
        ):
            metrics[key] = real(report.get(key), key)
            require(0. <= metrics[key] < 1e-10,
                    f"Incompatible exact rest t0 physical DAE: {key}")
        blocks = report.get("native_IFunction_t0_blocks_l2")
        require(isinstance(blocks, dict)
                and {"h", "p", "e_h", "e_p", "U", "U_n", "Y", "Y_n"}
                <= set(blocks), "Missing original native IFunction row blocks")
        for name, value in blocks.items():
            require(real(value, name) < 1e-10,
                    f"Original t0 row {name} not compatible")
        u1 = real(report.get("native_F_U_at_delta_l2"), "first wall actuation")
        u2 = real(report.get("native_F_U_at_two_delta_l2"),
                  "second wall actuation")
        analytic = real(report.get("analytic_wall_smoothstep_ratio"),
                        "analytic smoothstep")
        ratio = real(report.get("native_wall_smoothstep_ratio"), "native ratio")
        a = .0001/.02
        target = (3*(2*a)**2-2*(2*a)**3)/(3*a*a-2*a*a*a)
        require(u1 > 1e-15 and u2 > u1 and
                math.isclose(ratio, u2/u1, rel_tol=2e-8)
                and math.isclose(ratio, analytic, rel_tol=2e-6)
                and math.isclose(analytic, target, rel_tol=1e-12),
                "Original physical facet-tangent smoothstep does not start O(t^2)")
        return {
            "status": report["status"], "Fd": metrics["native_IFunction_t0_Fd_l2"],
            "Fa": metrics["native_IFunction_t0_Fa_l2"],
            "u1": u1, "u2": u2, "ratio": ratio,
        }

    require(report.get("status") == "HIGH_INDEX_OPT_IN_REJECTED_BEFORE_TS",
            "Original high-index negative control was not rejected")
    require(report.get("original_rotating_q0_fixed") is True
            and report.get("original_wall_already_rotating_at_t0") is True
            and report.get("rotation") == .2
            and report.get("actual_production_init_step_dispatch") ==
            "nonlinear_affine"
            and report.get("expected_exception_type") ==
            "UnsupportedNonlinearAlgebraicDAE",
            "Original rotating negative experiment was altered")
    reason = report.get("rejection_reason")
    require(isinstance(reason, str) and "null columns" in reason,
            "Wrong production nonlinear-affine typed rejection reason")
    null_cols = report.get("native_algebraic_null_columns")
    require(type(null_cols) is int and null_cols > 0,
            "Missing algebraic multiplier null column obstruction")
    init = report.get("production_initialization_report")
    require(isinstance(init, dict)
            and init.get("strategy") == "nonlinear_affine"
            and init.get("status") == "REJECTED"
            and init.get("algebraic_null_columns") == null_cols
            and init.get("no_pseudo_ts_executed") is True
            and init.get("no_physical_ts_executed") is True
            and init.get("high_index_boundary_multipliers_supported") is False,
            "Production rejection report is not the original native initializer")
    return {"status": report["status"], "null_cols": null_cols,
            "reason": reason}


def compare(path):
    root = Path(path)
    manifest = json.loads((root/"manifest.json").read_text(encoding="utf-8"))
    require(manifest.get("schema") == "v130-issue28-startup-gate-v1"
            and manifest.get("image_revision") == OCI_REVISION
            and manifest.get("no_physical_time_steps") is True
            and manifest.get("accepted_Y0") is False
            and manifest.get("production_code_unchanged") is True,
            "Wrong native pinned-image startup-gate manifest")
    summary = {}
    for mesh in (.2, .1):
        for case in CASES:
            by_rank = []
            for ranks in (1, 2):
                file = root/(
                    f"result_h{str(mesh).replace('.', 'p')}_{case}_{ranks}"
                )/"startup_gate.json"
                require(file.is_file(), f"Missing native startup gate report {file}")
                by_rank.append(check(
                    json.loads(file.read_text(encoding="utf-8")),
                    case, mesh, ranks,
                ))
            serial, parallel = by_rank
            if case == "compatible_rest_spinup":
                for key in ("Fa", "Fd", "u1", "u2", "ratio"):
                    require(math.isclose(
                        serial[key], parallel[key], rel_tol=2e-7,
                        abs_tol=1e-12,
                    ), f"Native MPI1/2 startup positive mismatch: {key}")
            else:
                require(serial["null_cols"] == parallel["null_cols"],
                        "Native MPI1/2 null-column detection differs")
            summary[f"h{mesh:g}_{case}"] = {
                "mpi1": serial, "mpi2": parallel,
                "mpi_parity": True,
            }
            print(
                f"h={mesh:g} {case}: MPI1/2 PASS, "
                + (f"Fa={serial['Fa']:.9g}, Fd={serial['Fd']:.9g}, "
                   f"U_ratio={serial['ratio']:.9g}"
                   if case == "compatible_rest_spinup"
                   else f"typed-rejection=null-cols-{serial['null_cols']}"),
                flush=True
            )
    result = {
        "schema": "scrimp-v130-issue28-startup-gate-comparison-v1",
        "status": "NATIVE_STARTUP_SAFETY_AND_T0_FEASIBILITY_DIAGNOSTIC_ONLY",
        "cases": summary,
        "production_high_index_initializer_qualified": False,
        "physical_time_integration_qualified": False,
        "no_physical_time_steps": True,
    }
    (root/"startup_gate_comparison.json").write_text(
        json.dumps(result, indent=2)+"\n", encoding="utf-8",
    )
    print("COMPLETE: startup gate diagnostics only; no TS/Y0 acceptance.")
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--evidence", type=Path, required=True)
    args = p.parse_args()
    compare(args.evidence)


if __name__ == "__main__":
    main()
