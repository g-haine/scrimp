# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Host-only comparator for original #28 physical FE trace compatibility.

Do not infer admissibility of the rotating DAE state from small residuals
or from an apparent mesh trend. Full native GetFEM MPI evidence required.
"""

import argparse
import json
import math
from pathlib import Path

EXPECTED_IMAGE = "1d6d0c90115213298047b28237a7f38087563a35"
FIELDS = (
    "initial_p_normal_squared",
    "initial_p_tangent_mismatch_squared",
    "initial_p_wall_vector_mismatch_squared",
    "initial_p_tangent_squared",
)


def require(ok, message):
    if not ok:
        raise ValueError(message)


def number(x, name):
    require(type(x) in (float, int) and math.isfinite(x),
            f"Nonfinite {name}: {x}")
    return float(x)


def close(a, b, rel=1e-3, absolute=1e-8):
    return math.isclose(a, b, rel_tol=rel, abs_tol=absolute)


def checked_integrals(sample, ranks):
    raw = sample["raw_integrals"]
    rms = sample["rms"]
    measures = {}
    for name, obj in raw.items():
        local = obj.get("per_rank")
        require(isinstance(local, list) and len(local) == ranks,
                f"Wrong MPI scalar integral ownership for {name}")
        arr = [number(v, name) for v in local]
        total = number(obj.get("sum"), f"sum {name}")
        require(close(total, math.fsum(arr), rel=1e-12, absolute=1e-10),
                f"Raw GetFEM MPI integral sum mismatch {name}")
        require(total >= -1e-10,
                f"Negative squared GetFEM integral {name}")
        measures[name] = total
    area, length = measures["area"], measures["perimeter"]
    require(area > 0 and length > 0, "Nonpositive physical domain or boundary")
    for name in measures:
        if name in ("area", "perimeter"):
            continue
        key = name.replace("_squared", "_rms")
        ratio = math.sqrt(max(0., measures[name])
                          / (area if name.startswith("interior_") else length))
        require(close(number(rms.get(key), key), ratio,
                      rel=2e-8, absolute=1e-11),
                f"Incorrect physically scaled RMS for {name}")
    for field in FIELDS:
        require(field in measures, f"Missing independent original trace {field}")
    # Rot*n is an orthogonal tangent: the vector defect decomposes into
    # actual normal and tangential contributions on the curved boundary.
    error = abs(
        measures["initial_p_wall_vector_mismatch_squared"]
        - measures["initial_p_normal_squared"]
        - measures["initial_p_tangent_mismatch_squared"]
    )
    require(error < 1e-7 * max(1., measures["initial_p_wall_vector_mismatch_squared"]),
            "Native boundary normal/tangent decomposition inconsistent")
    return measures


def mode_against_serial(serial, parallel):
    """Distinguish rank-local assembly from globally replicated scalars.

    Calibration uses TWO independent known measures, area and perimeter.
    It is unsafe to assume MPI scalar assembly returns local contributions.
    """
    local = all(close(parallel[k], serial[k], rel=2e-3)
                for k in ("area", "perimeter"))
    replicated = all(close(parallel[k] / 2, serial[k], rel=2e-3)
                     for k in ("area", "perimeter"))
    require(local != replicated,
            "MPI GetFEM scalar integral reduction is ambiguous or inconsistent")
    return "summed-local" if local else "replicated-global"


def compare(evidence):
    folder = Path(evidence)
    manifest = json.loads((folder/"manifest.json").read_text(encoding="utf-8"))
    require(manifest.get("schema") == "v130-issue28-physical-compatibility-v1",
            "Unexpected physical compatibility campaign")
    require(manifest.get("image_revision") == EXPECTED_IMAGE,
            "Unqualified image revision")
    require(manifest.get("no_physical_time_steps") is True
            and manifest.get("physical_initial_data_accepted") is False,
            "Unqualified physical initial data were accepted")

    overview = {}
    for mesh in (.2, .1):
        for case in ("rest", "rotating"):
            data, measured = [], []
            for ranks in (1, 2):
                file = folder/f"result_h{str(mesh).replace('.', 'p')}_{case}_{ranks}"/"physical_compatibility.json"
                require(file.is_file(), f"Missing native evidence: {file}")
                r = json.loads(file.read_text(encoding="utf-8"))
                require(r.get("schema") == "scrimp-v130-issue28-physical-compatibility-v1"
                        and r.get("status") == "DIAGNOSTIC_ONLY_NOT_INITIALIZED",
                        "Unexpected compatibility probe status")
                require(r.get("mpi_ranks") == ranks
                        and close(r.get("mesh_nominal_h"), mesh, rel=1e-12)
                        and r.get("rotation") == (0. if case=="rest" else .2),
                        f"Wrong physical experiment {mesh}/{case}/{ranks}")
                require(r.get("geometry_order") == 2
                        and r.get("height_order") == 2,
                        "Original curved P2 mesh or P2 state changed")
                require(r.get("original_q0_preserved_bitwise") is True
                        and r.get("state_rollback_bitwise") is True,
                        "Original physical h0,p0 changed or state not restored")
                require(r.get("accepted_rotating_Y0") is False
                        and r.get("exact_DAE_compatibility_qualified") is False
                        and r.get("initialization_TS_executed") is False
                        and r.get("physical_TS_executed") is False,
                        "Audit falsely claims initialized DAE")
                measured.append(checked_integrals(
                    r["boundary_at_prescribed_q0"], ranks
                ))
                projection = r.get("projection")
                require(isinstance(projection, dict), "No projection metadata")
                if case == "rest":
                    require(projection.get("type") == "exact-analytical-rest",
                            "Missing exact rest positive reference")
                    require(number(projection.get("native_Fa_l2"), "rest Fa")
                            < 1e-10, "Exact rest algebraic reference failed")
                    checked_integrals(projection["physical_integrals"], ranks)
                    for name in r["boundary_at_prescribed_q0"]["rms"]:
                        require(close(
                            projection["physical_integrals"]["rms"][name],
                            r["boundary_at_prescribed_q0"]["rms"][name],
                            rel=2e-8, absolute=1e-11),
                            "Rest projection altered physical initial trace")
                    for name in ("initial_p_normal_rms",
                                 "initial_p_tangent_mismatch_rms"):
                        require(r["boundary_at_prescribed_q0"]["rms"][name] < 1e-9,
                                "Nonzero velocity at exact physical rest")
                elif mesh == .2:
                    require(projection.get("type") ==
                            "distributed-rectangular-LSQR-two-tolerances",
                            "Missing independent rotating tolerance test")
                    trials = projection.get("trials", [])
                    require([v.get("kind") for v in trials] == ["standard", "tight"],
                            "Missing paired LSQR tolerances")
                    for trial in trials:
                        number(trial.get("native_Fa_l2"), "LSQR Fa")
                        number(trial.get("normal_equation_AT_Fa_l2"), "LSQR normal")
                        require(number(trial.get("algebraic_correction_l2"),
                                       "LSQR correction") >= 0.
                                and number(trial.get(
                                    "relative_coefficient_change_from_standard"
                                ), "LSQR sensitivity") >= 0.,
                                "Invalid algebraic LSQR coefficient sensitivity")
                        require(isinstance(trial.get("iterations"), int)
                                and trial["iterations"] >= 0,
                                "Invalid LSQR iteration count")
                        checked_integrals(trial["physical_integrals"], ranks)
                        for name in r["boundary_at_prescribed_q0"]["rms"]:
                            require(close(
                                trial["physical_integrals"]["rms"][name],
                                r["boundary_at_prescribed_q0"]["rms"][name],
                                rel=2e-8, absolute=1e-11),
                                "LSQR trial altered prescribed physical trace")
                else:
                    require(projection.get("type") ==
                            "not-run-refinement-boundary-only",
                            "Refinement must not silently project h0,p0")
                data.append(r)

            mode = mode_against_serial(measured[0], measured[1])
            for name in measured[0]:
                if mode == "replicated-global":
                    p = measured[1][name] / 2
                else:
                    p = measured[1][name]
                require(close(measured[0][name], p, rel=2e-3,
                              absolute=1e-9),
                        f"GetFEM MPI physical integral mismatch {name}")
            rms1 = data[0]["boundary_at_prescribed_q0"]["rms"]
            rms2 = data[1]["boundary_at_prescribed_q0"]["rms"]
            for name in rms1:
                require(close(rms1[name], rms2[name], rel=2e-3,
                              absolute=1e-9),
                        f"Physical RMS MPI 1/2 parity failure for {name}")
            key = f"h{mesh:g}_{case}"
            overview[key] = {
                "mpi2_scalar_assembly_mode": mode,
                "nominal_h": mesh, "case": case,
                "wall_normal_rms": rms1["initial_p_normal_rms"],
                "wall_tangential_mismatch_rms":
                    rms1["initial_p_tangent_mismatch_rms"],
                "wall_vector_mismatch_rms":
                    rms1["initial_p_wall_vector_mismatch_rms"],
                "mpi_physical_rms_parity_passed": True,
                "projection": "none" if mesh != .2 and case == "rotating"
                              else "exact" if case == "rest" else "two-LSQR",
            }
            print(
                f"h={mesh:g} {case}: normal RMS="
                f"{overview[key]['wall_normal_rms']:.9g}, "
                f"tangent mismatch RMS="
                f"{overview[key]['wall_tangential_mismatch_rms']:.9g}, "
                f"GetFEM MPI2 scalar mode={mode}", flush=True
            )
    coarse = overview["h0.2_rotating"]
    fine = overview["h0.1_rotating"]
    result = {
        "schema": "scrimp-v130-issue28-physical-compatibility-summary-v1",
        "status": "COMPLETE_PHYSICAL_TRACE_DIAGNOSTIC_NOT_INITIALIZER",
        "cases": overview,
        "normal_RMS_fine_over_coarse":
            fine["wall_normal_rms"] / max(coarse["wall_normal_rms"], 1e-100),
        "tangential_RMS_fine_over_coarse":
            fine["wall_tangential_mismatch_rms"]
            / max(coarse["wall_tangential_mismatch_rms"], 1e-100),
        "exact_rotating_DAE_initial_compatibility_qualified": False,
        "rotating_traction_Y0_accepted": False,
        "physical_time_integration_executed": False,
        "interpretation": (
            "Measured independent physically scaled P2 boundary traces "
            "with MPI parity, two nominal meshes and fixed original q0. "
            "Mesh trends and LSQR sensitivity are diagnostics, not a "
            "proof of exact algebraic feasibility or an initialized DAE."
        ),
    }
    (folder/"physical_compatibility_comparison.json").write_text(
        json.dumps(result, indent=2)+"\n", encoding="utf-8"
    )
    print("COMPLETE: physical trace diagnostic only; rotating DAE/Y0 NOT qualified.")
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--evidence", type=Path, required=True)
    args = p.parse_args()
    compare(args.evidence)


if __name__ == "__main__":
    main()
