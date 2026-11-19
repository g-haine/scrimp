# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Check 1/2-rank native MUMPS constitutive/Schur probe reports.

Acceptance is about algebraic identity and invertibility on tested
right-hand sides, NOT about scaling, power balance or spectral
conditioning of the future preconditioner. Native PETSc traces must
be checked separately for actual LU numeric factorization counts.
"""
import argparse
import json
from pathlib import Path

import numpy as np


def compare(directory):
    directory = Path(directory)
    records = {}
    for p in (1, 2):
        f = directory / f"schur_{p}r.json"
        result = json.loads(f.read_text(encoding="utf-8"))
        if (
            result["schema"] != "scrimp-constitutive-schur-probe-v1"
            or result["mpi_ranks"] != p
            or not result["reference_pde"]["mumps_full_system_unchanged"]
            or result["reference_pde"]["converged_reason"] <= 0
            or result["reference_pde"]["accepted_steps"] != 2
            or result["Q_solver"] != "PETSc KSPPREONLY + PCLU + MUMPS"
            or result["Q_factor_operator_configuration_count"] != 1
        ):
            raise ValueError(f"{f}: invalid native PETSc/MUMPS provenance")
        records[p] = result

    serial, parallel = records[1], records[2]
    for key in (
        "source_commit", "mesh_h", "dofs", "coenergy_dofs",
        "complement_dofs", "coenergy_fields", "complement_fields",
        "Q_solver",
    ):
        if serial[key] != parallel[key]:
            raise ValueError(f"1/2 rank constitutive mismatch: {key}")
    if serial["coenergy_dofs"] + serial["complement_dofs"] != serial["dofs"]:
        raise ValueError("Coenergy and complement do not cover all DoFs")
    if not serial["constitutive_Q_symmetric"] or not parallel[
        "constitutive_Q_symmetric"
    ]:
        raise ValueError("Constitutive Q is not symmetric")
    if min(
        serial["constitutive_Q_min_diagonal"],
        parallel["constitutive_Q_min_diagonal"],
    ) <= 0:
        raise ValueError("Constitutive Q has a nonpositive diagonal")
    norms = np.asarray([
        serial["constitutive_Q_frobenius_norm"],
        parallel["constitutive_Q_frobenius_norm"],
    ])
    if not np.all(np.isfinite(norms)) or not np.isclose(
        norms[0], norms[1], rtol=2e-9, atol=2e-11
    ):
        raise ValueError("1/2-rank constitutive Q norms differ")
    keys = (
        "relative_Aee_plus_Q",
        "Q_MUMPS_manufactured_relative_residual",
        "exact_schur_action_relative_error",
        "full_DAE_algebraic_relative_residual",
    )
    summary = []
    if len(serial["shift_results"]) != 2 or len(
        parallel["shift_results"]
    ) != 2:
        raise ValueError("Require precisely two distinct PETSc shifts")
    for i in range(2):
        shift = [record["shift_results"][i] for record in records.values()]
        if (
            shift[0]["sigma"] <= 0
            or shift[0]["sigma"] != shift[1]["sigma"]
        ):
            raise ValueError("PETSc shifts differ between MPI ranks")
        if i and shift[0]["sigma"] == serial["shift_results"][0]["sigma"]:
            raise ValueError("Shift values must differ")
        measurements = {}
        for key in keys:
            values = np.asarray([r[key] for r in shift], dtype=float)
            if not np.all(np.isfinite(values)) or (
                np.any(values < 0) or np.any(values > 5e-9)
            ):
                raise ValueError(
                    f"Schur {key} gate failed at shift "
                    f"{shift[0]['sigma']:g}: {values}"
                )
            measurements[key] = values.tolist()
        for result in shift:
            if (
                result["KSP_converged_reason"] <= 0
                or result["KSP_iterations"] < 0
                or not result["MUMPS_Q_operator_reused"]
            ):
                raise ValueError("Invalid constitutive MUMPS KSP solve")
        summary.append({"sigma": shift[0]["sigma"], **measurements})
    for result in records.values():
        error = float(
            result["constitutive_Q_manufactured_solution_relative_error"]
        )
        if not np.isfinite(error) or error > 5e-9 or error < 0:
            raise ValueError("Manufactured inverse-Q solve has drifted")
    out = {
        "schema": "scrimp-schur-MUMPS-parity-v1",
        "source_commit": serial["source_commit"],
        "mesh_h": serial["mesh_h"],
        "dofs": serial["dofs"],
        "coenergy_dofs": serial["coenergy_dofs"],
        "MPI_ranks": [1, 2],
        "Q_constitutive_norm_relative_difference": float(
            abs(norms[0]-norms[1]) / norms[0]
        ),
        "shift_parity": summary,
        "validated": [
            "Complete independent e / remaining x variable partition.",
            "PETSc Q=-K[e,e] independent of two distinct positive shifts.",
            "PETSc MUMPS solves the same Q for multiple RHS and shifts.",
            "Full uncondensed DAE and exact Schur actions agree.",
            "All boundary-variable DoFs remain in complement.",
        ],
        "not_yet_claimed": [
            "PETSc logging still needed to confirm LU factor reuse.",
            "No Krylov/Schur preconditioner or strong-scaling gain yet.",
            "No adaptive TS physical trajectory or energy-power gate yet.",
            "Q symmetry and positive diagonal are not an SPD proof.",
        ],
    }
    path = directory / "constitutive_schur_parity.json"
    path.write_text(json.dumps(out, indent=2)+"\n", encoding="utf-8")
    print(
        f"PASS MUMPS-Schur algebra: {serial['dofs']} DoFs, "
        f"{serial['coenergy_dofs']} coenergy, shifts "
        + ", ".join(f"{s['sigma']:g}" for s in summary)
        + "; 1/2 MPI ranks"
    )
    print(f"Saved: {path}")
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", type=Path, required=True)
    args = parser.parse_args()
    compare(args.directory)


if __name__ == "__main__":
    main()
