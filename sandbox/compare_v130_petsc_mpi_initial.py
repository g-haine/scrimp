# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Independent acceptance of the noncentralized PETSc DAE trial.

All four native runs must have identical physics per scenario at
MPI 1/2: a compatible nonzero polynomial state should start physical
CN exactly from the prescribed q/p (not pseudo-TS), preserve the
Hamiltonian + signed port-work identity, and yield MPI parity of
every full state. Incompatible published Gaussian must be explicitly
rejected by actual IFunction residual BEFORE physical CN starts.
"""

import argparse
import json
from pathlib import Path

import numpy as np


def _load(root, scenario, ranks):
    stem = f"distributed_{scenario}_{ranks}r"
    report = json.loads(
        (root / f"{stem}.json").read_text(encoding="utf-8")
    )
    if (
        report["scenario"] != scenario or report["mpi_ranks"] != ranks
        or report["pseudo_ts_executed"]
        or report["getfem_stiffness_gathered_to_rank_0"]
    ):
        raise RuntimeError(f"Invalid scientific provenance for {stem}")
    return report


def compare(root, *, max_scaled_state_discrepancy=5e-8):
    positive = {
        r: _load(root, "polynomial_compatible", r) for r in (1, 2)
    }
    negative = {r: _load(root, "published", r) for r in (1, 2)}
    basic = (
        "image_source_revision", "mesh_h", "dt", "steps",
    )
    ref = positive[1]
    for scenario, group in (
        ("polynomial_compatible", positive),
        ("published", negative),
    ):
        for ranks, x in group.items():
            for key in basic:
                if x[key] != ref[key]:
                    raise RuntimeError(
                        f"Different setup on {scenario}/{ranks}: {key}"
                    )
            if x["ksp_type"] != "lsqr":
                raise RuntimeError(
                    f"Unexpected KSP type in {scenario}/{ranks}"
                )
            if (
                x["algebraic_mass_row_max_abs"] > 1e-12
                or not x["differential_state_exactly_preserved"]
            ):
                raise RuntimeError(
                    f"Failed variable classification in {scenario}/{ranks}"
                )
            actual = x["native_algebraic_after_l2"]
            tol = x["native_residual_tolerance"]
            if scenario == "polynomial_compatible":
                if (
                    x["scientific_status"] != "PASS"
                    or not x["physical_cn_executed"]
                    or actual > tol
                ):
                    raise RuntimeError(
                        f"Compatible state not qualified at MPI={ranks}"
                    )
                power = x["power"]
                if (
                    power["absolute_port_work"] < 1e-9
                    or power["maximum_absolute_balance_defect"]
                    > x["power_balance_tolerance"]
                    or x["main_initial_vs_projected_l2"] > 5e-10
                ):
                    raise RuntimeError(
                        f"Non-vacuous port power or physical start "
                        f"failed MPI={ranks}"
                    )
            else:
                if (
                    x["scientific_status"] != "EXPECTED_REJECTION"
                    or not x["expected_incompatible_gaussian_rejected"]
                    or x["physical_cn_executed"]
                    or actual <= tol
                ):
                    raise RuntimeError(
                        "Published Gaussian incompatible-data guard "
                        f"not qualified at MPI={ranks}"
                    )

        a, b = group[1], group[2]
        for key in (
            "native_algebraic_before_l2",
            "native_algebraic_after_l2",
            "algebraic_state_correction_l2",
        ):
            if not np.isclose(
                a[key], b[key], atol=1e-8, rtol=5e-7,
            ):
                raise RuntimeError(
                    f"Different MPI algebraic solve {scenario}: {key}"
                )

    # Non-vacuous physics: ensure the new PETSc projection physically
    # reproduces the previously measured, independent centralized
    # diagnostic's signed boundary-work scale, not just parity
    # between two identical errors.
    for rank, r in positive.items():
        baseline_work = -5.39545892e-5
        measured = r["power"]["port_work_integral"][-1]
        if abs(measured-baseline_work) > 1e-8:
            raise RuntimeError(
                "New distributed initial state does not reproduce "
                f"independent rc3 wave boundary-work baseline: {measured}"
            )

    trajs = {}
    for ranks in (1, 2):
        with np.load(
            root / f"distributed_polynomial_compatible_{ranks}r.npz",
            allow_pickle=False,
        ) as f:
            trajs[ranks] = {
                k: np.array(f[k]) for k in (
                    "times", "states", "initial", "projected"
                )
            }
    a, b = trajs[1], trajs[2]
    if not np.array_equal(a["initial"], b["initial"]):
        raise RuntimeError("Prescribed initial GetFEM coefficients differ in MPI")
    if not np.allclose(
        a["times"], b["times"], atol=1e-12, rtol=0
    ):
        raise RuntimeError("MPI trajectory time grids differ")
    if a["states"].shape != b["states"].shape:
        raise RuntimeError("MPI trajectory snapshots differ in shape")

    errors = [
        float(np.linalg.norm(x-y)/max(
            1.0, np.linalg.norm(x)
        ))
        for x, y in zip(a["states"], b["states"])
    ]
    maximum = max(errors)
    if maximum > max_scaled_state_discrepancy:
        raise RuntimeError(
            "Distributed DAE full CN trajectory MPI disparity: "
            f"{maximum}"
        )
    for quantity in (
        "hamiltonian", "port_work_integral", "step_balance_defects",
    ):
        if not np.allclose(
            positive[1]["power"][quantity],
            positive[2]["power"][quantity],
            rtol=5e-8, atol=1e-9,
        ):
            raise RuntimeError(
                f"Distributed projection {quantity} differs in MPI"
            )
    for name in positive[1]["power"]["boundary_midpoint_powers"]:
        values1 = positive[1]["power"][
            "boundary_midpoint_powers"
        ][name]
        values2 = positive[2]["power"][
            "boundary_midpoint_powers"
        ][name]
        if not np.allclose(
            values1, values2, rtol=5e-8, atol=1e-9
        ):
            raise RuntimeError(
                f"Port {name} midpoint powers not MPI-equivalent"
            )

    summary = {
        "schema": "scrimp-v130-petsc-distributed-projection-verification-v1",
        "rc3_source_revision": ref["image_source_revision"],
        "distributed_kaa_no_rank0_gather": True,
        "original_pseudo_time_ts_bypassed": True,
        "compatible_polynomial": "PASS",
        "incompatible_gaussian": "EXPECTED_REJECTION",
        "power_balance_with_nonzero_port_work": "PASS",
        "full_CN_trajectory_MPI_1_2": "PASS",
        "max_scaled_full_trajectory_MPI_discrepancy": maximum,
        "reproduces_independent_scientific_signed_port_work": True,
    }
    (root / "distributed_dae_comparison.json").write_text(
        json.dumps(summary, indent=2)+"\n",
        encoding="utf-8",
    )
    print(
        "SCRIMP DAE projected initial state: PETSc-distributed Kaa "
        "without pseudo-time TS"
    )
    print("| MPI | Compatible residual | Gaussian residual | "
          "Gaussian rejected before CN | CN max power defect | "
          "Integrated port work |")
    print("| ---: | ---: | ---: | --- | ---: | ---: |")
    for ranks in (1, 2):
        p, n = positive[ranks], negative[ranks]
        print(
            f"| {ranks} | {p['native_algebraic_after_l2']:.8g} | "
            f"{n['native_algebraic_after_l2']:.8g} | "
            f"{n['scientific_status']} | "
            f"{p['power']['maximum_absolute_balance_defect']:.8g} | "
            f"{p['power']['port_work_integral'][-1]:.8g} |"
        )
    print(
        "PASS: actual PETSc-distributed algebraic block, "
        "prescribed q0/p0, full CN MPI parity, and non-vacuous "
        "Hamiltonian + signed boundary work"
    )
    print(
        "PASS: incompatible Gaussian explicitly rejected at t0 "
        "without pseudo-time evolution or CN"
    )
    return summary


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("output_dir", type=Path)
    args = p.parse_args()
    compare(args.output_dir)


if __name__ == "__main__":
    main()
