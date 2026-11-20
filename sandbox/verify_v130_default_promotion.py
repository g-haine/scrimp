# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Strict final native gate for changing SCRIMP's DAE default.

The existing independent physical comparators must PASS unchanged:
conservative 2D wave versus frozen diagnostic and pseudo/bdf
controls, and dissipative 2D wave with signed damping work. This
wrapper additionally proves *no init_step_ts_type option* was
specified for either the successful or rejected new-default case.
"""

import argparse
import json
from pathlib import Path


def verify(folder):
    wave = json.loads(
        (folder / "core_dae_qualification.json").read_text(
            encoding="utf-8"
        )
    )
    dissipative = json.loads(
        (folder / "second_linear_dae_comparison.json").read_text(
            encoding="utf-8"
        )
    )
    if (
        wave["qualified_opt_in_core_DAE"] != "PASS"
        or dissipative["scientific_status"] != "PASS"
    ):
        raise RuntimeError(
            "One or both independent default physics comparators failed"
        )
    cases = {}
    for ranks in (1, 2):
        for scenario in ("polynomial_compatible", "published"):
            stem = f"core_algebraic_{scenario}_{ranks}r"
            report = json.loads(
                (folder / f"{stem}.json").read_text(encoding="utf-8")
            )
            if (
                not report.get("implicit_default_requested")
                or report.get("explicit_init_step_ts_type_present")
                or "init_step_ts_type" in report["initialization_options"]
                or report["mode"] != "algebraic"
            ):
                raise RuntimeError(
                    "A wave test supplied an explicit initializer "
                    f"instead of exercising the actual default: {stem}"
                )
            if scenario == "published":
                if (
                    report["status"] != "EXPECTED_REJECTION"
                    or report["physical_cn_executed"]
                    or not report["initial_data_preserved"]
                ):
                    raise RuntimeError(
                        "Default algebraic initializer silently modified "
                        f"or accepted the incompatible Gaussian MPI={ranks}"
                    )
            elif (
                report["status"] != "PASS"
                or report["init_report"]["strategy"] != "algebraic"
                or report["init_report"]["pseudo_ts_executed"]
                or not report["init_report"]["differential_initial_exact"]
            ):
                raise RuntimeError(
                    "Implicit default did not select the native "
                    f"algebraic projection MPI={ranks}"
                )
            cases[stem] = report["status"]
        stem = f"dissipative_core_{ranks}r"
        report = json.loads(
            (folder / f"{stem}.json").read_text(encoding="utf-8")
        )
        if (
            not report.get("implicit_default_requested")
            or report.get("explicit_init_step_ts_type_present")
            or report["dae_projection"]["strategy"] != "algebraic"
            or report["dae_projection"]["pseudo_ts_executed"]
        ):
            raise RuntimeError(
                "Dissipative wave supplied an explicit strategy or "
                f"used pseudo in the new default: {stem}"
            )
        cases[stem] = "PASS"
    for mode in ("pseudo", "bdf2"):
        stem = f"core_{mode}_polynomial_compatible_1r"
        report = json.loads(
            (folder / f"{stem}.json").read_text(encoding="utf-8")
        )
        if (
            report["status"] != "PASS"
            or report.get("implicit_default_requested")
            or not report.get("explicit_init_step_ts_type_present")
            or report["initialization_options"]["init_step_ts_type"]
               != ("pseudo" if mode == "pseudo" else "bdf")
        ):
            raise RuntimeError(
                f"Historical explicit {mode} option no longer runs"
            )
        cases[stem] = "PASS"
    summary = {
        "schema": "scrimp-v130-default-algebraic-migration-v1",
        "status": "PASS",
        "actual_default_option_omitted": True,
        "qualified_linear_conservative_and_dissipative_models": True,
        "historical_pseudo_and_BDF2_explicitly_usable": True,
        "incompatible_gaussian_rejected_before_CN_without_modification":
            True,
        "passed_independent_full_CN_and_signed_power_comparators": True,
        "new_default_is_released": False,
        "cases": cases,
        "wave_reference_max_scaled_1r": wave[
            "prototype_full_state_parity"
        ]["1r"]["complete_trajectory_max_scaled_l2"],
        "wave_reference_max_scaled_2r": wave[
            "prototype_full_state_parity"
        ]["2r"]["complete_trajectory_max_scaled_l2"],
        "dissipative_MPI_max_scaled_full_CN":
            dissipative["max_scaled_full_CN_MPI_discrepancy"],
    }
    (folder / "default_promotion_qualification.json").write_text(
        json.dumps(summary, indent=2)+"\n",encoding="utf-8",
    )
    print("SCRIMP v1.3.0 algebraic DAE default native migration:")
    print("PASS: no init_step_ts_type option on 2 linear models, MPI 1/2")
    print("PASS: published Gaussian rejected without change to q0/p0")
    print("PASS: independent complete CN, signed boundary and damping power")
    print("PASS: historical pseudo and tiny negative-time BDF2 still explicit")
    print("PASS: strict native gate of staged default candidate")
    print("Release readiness and large-MPI memory remain separate gates")
    return summary


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("output_dir", type=Path)
    args = p.parse_args()
    verify(args.output_dir)


if __name__ == "__main__":
    main()
