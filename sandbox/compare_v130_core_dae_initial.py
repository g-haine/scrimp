# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Scientific acceptance of the PRODUCTION opt-in DAE initializer.

Uses the immutable *independent* rc3 prototype .npz/JSON as a real
trajectory reference, not merely a hard-coded terminal-state scalar.
The default pseudo and explicit negative-time tiny BDF2 legacy paths
are separately observed; only algebraic is expected to preserve q0/p0.
"""

import argparse
import json
from pathlib import Path

import numpy as np


def _read(folder, mode, scenario, ranks):
    stem = f"core_{mode}_{scenario}_{ranks}r"
    result = json.loads(
        (folder / f"{stem}.json").read_text(encoding="utf-8")
    )
    if (
        result["mode"] != mode
        or result["scenario"] != scenario
        or result["mpi_ranks"] != ranks
    ):
        raise RuntimeError(f"Wrong result identity: {stem}")
    return result


def _npz(folder, name):
    with np.load(folder / name, allow_pickle=False) as src:
        return {key: np.array(src[key]) for key in src.files}


def _relative(a, b):
    a = np.asarray(a)
    b = np.asarray(b)
    if a.shape != b.shape:
        raise RuntimeError(
            f"Cannot compare different state layouts: {a.shape}, {b.shape}"
        )
    return float(
        np.linalg.norm(a-b) / max(1., np.linalg.norm(b))
    )


def compare(out, reference, legacy_reference):
    positive = {
        ranks: _read(out, "algebraic", "polynomial_compatible", ranks)
        for ranks in (1, 2)
    }
    negative = {
        ranks: _read(out, "algebraic", "published", ranks)
        for ranks in (1, 2)
    }
    summary = {
        "schema": "scrimp-v130-core-dae-qualification-v1",
        "producer": "production DPHS.init_step opt-in, no monkeypatch",
        "qualified_prototype_reference": str(reference),
        "qualified_pseudo_reference": str(legacy_reference),
        "modes": {},
        "prototype_full_state_parity": {},
    }
    projected_states = {}
    for ranks in (1, 2):
        p = positive[ranks]
        n = negative[ranks]
        if (
            p["status"] != "PASS"
            or not p["physical_cn_executed"]
            or p["init_report"]["status"] != "PASS"
            or not p["init_report"]["differential_initial_exact"]
            or p["init_report"]["pseudo_ts_executed"]
            or p["init_report"]["native_final_residual_l2"]
               > p["init_report"]["feasibility_tolerance"]
            or not p["init_report"]["accuracy_target_met"]
            or p["init_report"]["mass_algebraic_max_abs"] > 1e-12
        ):
            raise RuntimeError(
                f"Production algebraic physical/DAE gate failed: MPI={ranks}"
            )
        if (
            n["status"] != "EXPECTED_REJECTION"
            or n["physical_cn_executed"]
            or not n["initial_data_preserved"]
            or n["init_report"]["status"] != "UNRESOLVED"
            or n["init_report"]["native_final_residual_l2"]
               <= n["init_report"]["feasibility_tolerance"]
        ):
            raise RuntimeError(
                f"Production incompatible Gaussian not atomically "
                f"rejected: MPI={ranks}"
            )

        actual = _npz(
            out, f"core_algebraic_polynomial_compatible_{ranks}r.npz"
        )
        frozen = _npz(
            reference, f"distributed_polynomial_compatible_{ranks}r.npz"
        )
        if not np.array_equal(actual["initial"], frozen["initial"]):
            raise RuntimeError(
                "The production run does not use the frozen "
                "prototype's exact prescribed GetFEM initial state"
            )
        if not np.allclose(
            actual["times"], frozen["times"], atol=1e-12, rtol=0
        ):
            raise RuntimeError(
                "Production CN time grid differs from frozen prototype"
            )
        # Independent complete 11-snapshot physical trajectory baseline;
        # do not loosen the pre-existing 5e-8 acceptance gate.
        full_error = max(
            _relative(a,b) for a,b in zip(
                actual["states"], frozen["states"]
            )
        )
        init_error = _relative(
            actual["projected"], frozen["projected"]
        )
        if full_error > 5e-8 or init_error > 5e-8:
            raise RuntimeError(
                "Production DAE differs from frozen distributed "
                f"prototype MPI={ranks}: t0={init_error:.9g}, "
                f"trajectory={full_error:.9g}"
            )
        power = p["physical_power"]
        if (
            power["absolute_port_work"] <= 1e-9
            or power["maximum_absolute_balance_defect"]
               > p["power_balance_tolerance"]
        ):
            raise RuntimeError(
                f"Production DAE nonvacuous power gate failed MPI={ranks}"
            )
        reference_power = json.loads(
            (reference / f"distributed_polynomial_compatible_{ranks}r.json"
             ).read_text(encoding="utf-8")
        )["power"]
        for key in (
            "hamiltonian", "port_work_integral",
            "step_balance_defects",
        ):
            if not np.allclose(
                power[key], reference_power[key],
                rtol=5e-8, atol=1e-9
            ):
                raise RuntimeError(
                    f"Production {key} differs from frozen physical "
                    f"reference MPI={ranks}"
                )
        for name, values in reference_power[
            "boundary_midpoint_powers"
        ].items():
            if name not in power["boundary_midpoint_powers"] or (
                not np.allclose(
                    power["boundary_midpoint_powers"][name],
                    values, rtol=5e-8, atol=1e-9,
                )
            ):
                raise RuntimeError(
                    f"Production signed {name} boundary powers "
                    f"changed versus frozen prototype MPI={ranks}"
                )
        summary["prototype_full_state_parity"][f"{ranks}r"] = {
            "initial_scaled_l2": init_error,
            "complete_trajectory_max_scaled_l2": full_error,
            "native_algebraic_residual": p["init_report"][
                "native_final_residual_l2"
            ],
            "energy_balance_max_defect": power[
                "maximum_absolute_balance_defect"
            ],
            "signed_integrated_port_work": power[
                "port_work_integral"
            ][-1],
            "negative_gaussian_rejected": True,
        }
        projected_states[ranks] = actual

    cross_mpi_error = max(
        _relative(a,b) for a,b in zip(
            projected_states[1]["states"],
            projected_states[2]["states"],
        )
    )
    if cross_mpi_error > 5e-8:
        raise RuntimeError(
            "Production full CN 1/2 MPI trajectory parity fails: "
            f"{cross_mpi_error}"
        )
    summary["new_core_MPI_1_2_parity_max_scaled_l2"] = cross_mpi_error

    # Legacy paths are deliberately not required to preserve q0/p0:
    # they integrate a tiny pseudo-time or backward BDF2 interval.
    for mode in ("pseudo", "bdf2"):
        legacy = _read(out, mode, "polynomial_compatible", 1)
        if (
            legacy["status"] != "PASS"
            or not legacy["physical_cn_executed"]
            or legacy["physical_power"]["absolute_port_work"] <= 1e-9
        ):
            raise RuntimeError(
                f"Historical initializer {mode} no longer runs in core"
            )
        if mode == "bdf2":
            if (
                legacy["initialization_options"][
                    "init_step_from_negative_time"
                ].lower() != "true"
                or legacy["initialization_options"][
                    "init_step_bdf_order"
                ] != "2"
            ):
                raise RuntimeError("BDF2 negative-time smoke not configured")
        summary["modes"][mode] = {
            "differential_initial_drift": legacy[
                "initial_dynamic_drift_l2"
            ],
            "power_balance_max_defect": legacy[
                "physical_power"
            ]["maximum_absolute_balance_defect"],
            "scope": (
                "legacy small time integration, no exact fidelity guarantee"
            ),
        }

    # At least the pseudo mode must remain numerically identical
    # to the measured, independent, unchanged pre-integration path.
    legacy = _npz(
        out, "core_pseudo_polynomial_compatible_1r.npz"
    )
    frozen_pseudo = _npz(
        legacy_reference, "wave_pseudo_1r.npz"
    )
    if not np.allclose(
        legacy["times"], frozen_pseudo["times"],
        rtol=0, atol=1e-12,
    ):
        raise RuntimeError("Historical pseudo time grid changed")
    old_pseudo_error = max(
        _relative(a,b) for a,b in zip(
            legacy["states"], frozen_pseudo["trajectory"]
        )
    )
    if old_pseudo_error > 5e-8:
        raise RuntimeError(
            "Historical pseudo TS main CN changed versus frozen "
            f"physical baseline: {old_pseudo_error}"
        )
    summary["legacy_pseudo_reference_max_scaled_l2"] = old_pseudo_error
    summary["qualified_opt_in_core_DAE"] = "PASS"
    summary["promote_to_default"] = "NOT YET — separate gate required"
    (out / "core_dae_qualification.json").write_text(
        json.dumps(summary, indent=2)+"\n", encoding="utf-8"
    )
    print(
        "SCRIMP core opt-in algebraic initialization qualification"
    )
    print("| MPI | F_a(t0) | Frozen CN max scaled L2 | "
          "CN power defect | Gaussian rejected |")
    print("|---:|---:|---:|---:|---|")
    for ranks in (1,2):
        p = summary["prototype_full_state_parity"][f"{ranks}r"]
        print(
            f"| {ranks} | {p['native_algebraic_residual']:.9g} | "
            f"{p['complete_trajectory_max_scaled_l2']:.9g} | "
            f"{p['energy_balance_max_defect']:.9g} | YES |"
        )
    print(
        f"PASS: legacy pseudo TS matches independent reference "
        f"to {old_pseudo_error:.9g}; legacy BDF2 runs with negative "
        f"t0-dt and has documented dynamic drift"
    )
    print(
        "PASS: core opt-in DAE native 1/2 MPI full physics; "
        "default-promotion gate remains separate"
    )
    return summary


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--reference", type=Path, required=True)
    p.add_argument("--legacy-reference", type=Path, required=True)
    args = p.parse_args()
    compare(args.out, args.reference, args.legacy_reference)


if __name__ == "__main__":
    main()
