# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Compare complete 1-/2-rank coenergy CN states, powers, Hamiltonian.

This is an *independent* trajectory/physical comparison following
the demonstrated entrywise GetFEM/PETSc 1/2-rank matrix parity gate.
Do not claim adaptive-step correctness or full MUMPS scaling based
on a five-step fixed-CN microbenchmark.
"""
import argparse
import json
from pathlib import Path
import numpy as np


FIELDS = (
    "times", "states", "energy", "energy_getfem",
    "power_damping_mid", "power_boundary_mid",
    "cn_relative_residuals", "algebraic_relative_residuals",
    "relative_power_defects",
)
# Residual and defect arrays are nonnegative *diagnostics*: compare
# their absolute magnitudes, NOT their small relative differences.
PARITY = (
    "states", "energy", "energy_getfem",
    "power_damping_mid", "power_boundary_mid",
)


def compare(folder):
    folder = Path(folder)
    reports = {}
    arrays = {}
    for ranks in (1, 2):
        report = json.loads(
            (folder/f"coenergy_{ranks}r_physics.json").read_text(
                encoding="utf-8"
            )
        )
        if (
            report["schema"] != "scrimp-coenergy-physics-v1"
            or report["mpi_ranks"] != ranks
            or report["TS_converged_reason"] <= 0
            or report["formulation"] != "substituted_coenergy"
            or report["steps"] < 3
        ):
            raise ValueError("Invalid native coenergy physics provenance")
        reports[ranks] = report
        with np.load(
            folder/f"coenergy_{ranks}r_physics.npz",
            allow_pickle=False,
        ) as z:
            arrays[ranks] = {name: z[name] for name in FIELDS}

    a, b = reports[1], reports[2]
    for key in (
        "source_commit", "mesh", "dofs", "steps", "dt",
        "t_final", "full_state_snapshots", "formulation",
        "power_sign_convention",
    ):
        if a[key] != b[key]:
            raise ValueError(f"1/2 MPI physics mismatch: {key}")
    if not np.array_equal(
        arrays[1]["times"], arrays[2]["times"]
    ):
        raise ValueError("Mismatched saved accepted CN time grids")
    if (arrays[1]["states"].shape
            != (a["steps"]+1, a["dofs"])):
        raise ValueError("Incorrect full state snapshot shape")
    if arrays[2]["states"].shape != arrays[1]["states"].shape:
        raise ValueError("Incompatible MPI physical trajectory shapes")
    metrics = {}
    for key in PARITY:
        x, y = arrays[1][key], arrays[2][key]
        if x.shape != y.shape or not np.all(
            np.isfinite(x)
        ) or not np.all(np.isfinite(y)):
            raise ValueError(f"{key}: invalid MPI data")
        difference = np.linalg.norm(x-y)
        reference = max(
            np.linalg.norm(x), np.linalg.norm(y), 1e-300
        )
        relative = float(difference/reference)
        max_abs = float(np.max(np.abs(x-y)))
        if not np.allclose(x, y, rtol=5e-8, atol=2e-10):
            raise ValueError(
                f"1/2 MPI {key} physical drift: "
                f"rel={relative:.6g}, max_abs={max_abs:.6g}"
            )
        metrics[key] = {
            "relative_l2_difference": relative,
            "maximum_absolute_difference": max_abs,
        }
    for ranks, report in reports.items():
        for metric, measured in report["maxima"].items():
            if not np.isfinite(measured) or measured < 0 or measured > 5e-8:
                raise ValueError(
                    f"{ranks} MPI ranks: {metric} failed: {measured}"
                )
    out = {
        "schema": "scrimp-coenergy-physical-parity-v1",
        "source_commit": a["source_commit"],
        "mesh": a["mesh"],
        "dofs": a["dofs"],
        "dt": a["dt"], "steps": a["steps"],
        "mpi_ranks": [1, 2],
        "metrics": metrics,
        "per_rank_physical_gates": {
            str(r): reports[r]["maxima"] for r in reports
        },
        "validated": [
            "Full time-indexed dynamic, algebraic and boundary states.",
            "Finite-element Hamiltonian vs quadratic PETSc mass energy.",
            "Implicit CN DAE residual at all accepted steps.",
            "Damping and U=0.5Y boundary-feedback power midpoint integrals.",
            "Discrete energy balance after accounting for both power ports.",
        ],
        "limits": [
            "No variable-timestep comparison performed.",
            "No performance benchmark or preconditioner installed.",
            "State power-balance tests do not by themselves prove "
            "the PDE continuum spatial convergence.",
        ],
    }
    path = folder/"coenergy_physics_parity.json"
    path.write_text(json.dumps(out, indent=2)+"\n", encoding="utf-8")
    print(
        f"PASS coenergy physical parity: {a['dofs']} DoFs, "
        f"{a['steps']} steps, 1 vs 2 MPI; "
        f"state relative diff="
        f"{metrics['states']['relative_l2_difference']:.3g}, "
        f"Hamiltonian diff="
        f"{metrics['energy']['relative_l2_difference']:.3g}",
        flush=True,
    )
    print(f"Saved: {path}", flush=True)
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", type=Path, required=True)
    compare(parser.parse_args().directory)


if __name__ == "__main__":
    main()
