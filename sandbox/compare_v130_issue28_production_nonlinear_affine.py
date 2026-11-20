# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Recheck saved #28 production nonlinear-affine native MPI1/2 evidence.

Runs on the HOST with stock Python, no PETSc, GetFEM, MPI or Docker.
Both JSON artifacts must come from actual prior native production
DPHS.init_step() runs; this comparator does NOT regenerate evidence.

Use the existing evidence folder (e.g. after a terminal script typo):
  python3 sandbox/compare_v130_issue28_production_nonlinear_affine.py \
      --evidence sandbox/outputs/v130-issue28-production-na-<commit>-<suffix>
"""

import argparse
import json
import math
from pathlib import Path
import re

EXPECTED_REVISION = "1d6d0c90115213298047b28237a7f38087563a35"
SCHEMA = "scrimp-v130-issue28-production-nonlinear-affine-shallow-v1"


def require(condition, message):
    if not condition:
        raise ValueError(message)


def finite(value, name):
    require(
        isinstance(value, (float, int)) and not isinstance(value, bool)
        and math.isfinite(value),
        f"Missing or nonfinite {name}: {value!r}",
    )
    return float(value)


def read_json(path):
    require(path.is_file(), f"Missing native evidence: {path}")
    with path.open("r", encoding="utf-8") as source:
        return json.load(source)


def check_evidence(folder):
    """Validate archived successful native results without rerunning them."""
    base = Path(folder).resolve()
    require(base.is_dir(), f"Evidence directory does not exist: {base}")

    manifest = read_json(base / "manifest.json")
    require(manifest.get("schema") == "issue28-production-nonlinear-affine-v1",
            "Unexpected evidence manifest schema")
    require(
        manifest.get("immutable_revision") == EXPECTED_REVISION,
        "Wrong Docker OCI revision; this evidence cannot qualify v1.3.0-rc3",
    )
    require(
        bool(re.fullmatch(r"[0-9a-f]{40}", manifest.get("source", ""))),
        "Missing or invalid source commit SHA in archived manifest",
    )
    require(manifest.get("opt_in_only") is True,
            "Manifest did not record explicit opt-in")
    require(manifest.get("physical_steps") == 0,
            "Manifest claims physical time stepping")

    results = []
    for ranks in (1, 2):
        path = base / f"result_{ranks}" / "production_nonlinear_affine.json"
        data = read_json(path)
        require(data.get("schema") == SCHEMA,
                f"MPI {ranks}: wrong production report schema")
        require(data.get("mpi_ranks") == ranks,
                f"MPI {ranks}: wrong rank count")
        require(data.get("model") == "original shipped examples/shallow_water.py",
                f"MPI {ranks}: model changed")
        require(finite(data.get("mesh"), "mesh") == 0.1,
                f"MPI {ranks}: original mesh changed")
        require(data.get("status") == "PASS_NATIVE_PRODUCTION_INITIALIZER",
                f"MPI {ranks}: production run was not successful")
        require(data.get("original_dynamic_h0_p0_bitwise") is True,
                f"MPI {ranks}: prescribed differential data changed")
        require(data.get("published_algebraic_dofs_modified") is True,
                f"MPI {ranks}: no algebraic coefficients were published")
        require(data.get("explicit_gyro_from_original_shallow_water") is True,
                f"MPI {ranks}: original explicit gyro not present")
        require(data.get("explicit_gyro_cache_ready") is True,
                f"MPI {ranks}: original explicit RHS not frozen")
        for field in (
            "main_physical_TS_executed",
            "pseudo_initial_TS_executed",
            "generic_nonlinear_algebraic_Newton_qualified",
            "high_index_rotating_Y_initialization_qualified",
        ):
            require(data.get(field) is False, f"MPI {ranks}: {field} must be False")

        report = data.get("real_production_init_step_report")
        require(isinstance(report, dict),
                f"MPI {ranks}: missing real production init report")
        require(report.get("strategy") == "nonlinear_affine",
                f"MPI {ranks}: wrong initializer strategy")
        require(report.get("status") == "PASS",
                f"MPI {ranks}: production initializer rejected trial")
        require(report.get("differential_initial_exact") is True,
                f"MPI {ranks}: differential data not preserved")
        require(report.get("accuracy_target_met") is True,
                f"MPI {ranks}: strict acceptance target missed")
        require(report.get("no_pseudo_ts_executed") is True,
                f"MPI {ranks}: pseudo init TS executed")
        require(report.get("no_physical_ts_executed") is True,
                f"MPI {ranks}: physical TS executed")
        require(report.get("algebraic_null_columns") == 0,
                f"MPI {ranks}: algebraic Jacobian structurally singular")
        require(report.get("linear_solver") == "PETSc PREONLY + MUMPS",
                f"MPI {ranks}: unexpected linear solver")
        require(report.get("linear_solver_converged_reason", 0) > 0,
                f"MPI {ranks}: PETSc direct solver did not converge")
        require(finite(report.get("max_abs_algebraic_mass_rows"), "mass rows") <= 1e-12,
                f"MPI {ranks}: unexpected algebraic mass rows")
        initial = finite(report.get("native_initial_algebraic_residual_l2"),
                         "initial Fa")
        require(initial > 1, f"MPI {ranks}: vacuous initial projection")
        final = finite(report.get("native_final_algebraic_residual_l2"),
                       "final Fa")
        published = finite(data.get("actual_published_IFunction_Fa_l2"),
                           "published Fa")
        accuracy = finite(report.get("accuracy_target"), "accuracy target")
        require(0 <= final <= min(1e-11, accuracy),
                f"MPI {ranks}: native initializer final residual too high")
        require(0 <= published <= min(1e-11, accuracy),
                f"MPI {ranks}: published native residual too high")
        require(math.isclose(published, final, rel_tol=1e-3, abs_tol=1e-11),
                f"MPI {ranks}: publication differs from init report")
        affine_error = finite(
            report.get("affine_algebraic_identity_error_l2"),
            "affine identity error",
        )
        affine_tol = finite(report.get("affine_identity_tolerance"),
                            "affine identity tolerance")
        require(0 <= affine_error <= min(1e-8, affine_tol),
                f"MPI {ranks}: affine-algebraic identity failed")
        require(data.get("dofs") == 4270
                and report.get("dofs") == 4270
                and report.get("differential_dofs") == 2035
                and report.get("algebraic_dofs") == 2235,
                f"MPI {ranks}: unexpected original model DOFs")
        results.append(data)

    serial, parallel = results
    for field in ("dofs", "mesh", "model"):
        require(serial[field] == parallel[field],
                f"MPI parity: {field} differs")

    # Important: Python's math.isclose uses rel_tol, not 'rtol'.
    for field in ("actual_published_IFunction_Fa_l2",):
        a = serial[field]
        b = parallel[field]
        require(math.isclose(a, b, rel_tol=1e-3, abs_tol=1e-11),
                f"MPI 1/2 parity failed for {field}: {a} vs {b}")

    for field in (
        "native_initial_algebraic_residual_l2",
        "native_final_algebraic_residual_l2",
        "affine_algebraic_identity_error_l2",
    ):
        a = serial["real_production_init_step_report"][field]
        b = parallel["real_production_init_step_report"][field]
        require(math.isclose(a, b, rel_tol=1e-3, abs_tol=1e-11),
                f"MPI 1/2 parity failed for {field}: {a} vs {b}")

    print(
        "PASS: archived *native* production nonlinear_affine initial projection, "
        f"MPI 1: |Fa|={serial['actual_published_IFunction_Fa_l2']:.9g}, "
        f"MPI 2: |Fa|={parallel['actual_published_IFunction_Fa_l2']:.9g}",
        flush=True,
    )
    print(
        "Original h0,p0 bitwise preserved; algebraic coefficients published; "
        "real IFunction acceptance met at both MPI ranks.",
        flush=True,
    )
    print(
        "No physical TS or nonlinear-in-algebraic Newton performed. "
        "Rotating high-index Y0 initialization remains unqualified.",
        flush=True,
    )
    print(f"Host evidence verified: {base}", flush=True)
    return manifest, serial, parallel


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--evidence", type=Path, required=True)
    args = parser.parse_args()
    check_evidence(args.evidence)


if __name__ == "__main__":
    main()
