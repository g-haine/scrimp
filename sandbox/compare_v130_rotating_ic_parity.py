# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Native DATA comparator for SCRIMP checkpoint equivalence, not energy accuracy.

Continuous PETSc/GetFEM TSCN trajectory is the oracle. Each simulation
keeps its own (possibly imperfect) Hamiltonian balance. The sole gate is
whether pausing/resuming produces the SAME states and diagnostics.
"""

import argparse
import json
from pathlib import Path

import numpy as np


TIME_TOL = 1e-12
FULL_STATE_L2_TOL = 1e-8
OBS_ABS_TOL = 1e-8
FA_TOL = 1e-8
CN_TOL = 1e-7

TRACE_FIELDS = (
    "H", "volume", "p_normal", "p_viscous", "d_grad", "d_div",
    "boundary_work", "viscous_dissipation",
    "power_balance", "volume_balance",
)


def load_segment(directory, expected_phase):
    directory = Path(directory)
    report = json.loads(
        (directory / "initial_state_report.json").read_text(encoding="utf-8")
    )
    if (report.get("schema") != "scrimp-v130-rotating-controlled-ic-native-v1"
            or report.get("phase") != expected_phase
            or report.get("qualification_mode") != "restart_parity"
            or report.get("physical_relative_power_gate_applied") is not False
            or report.get("physical_energy_defect_retained_not_modified") is not True
            or report.get("status") !=
               "ACCEPTED_COMPATIBLE_ROTATING_DAE_INITIAL_STATE_NOT_STEADY"):
        raise ValueError(f"Not a qualified accepted native {expected_phase} report")
    if (report["PETSc_reason"] <= 0 or report["PETSc_rejected"]
            or report["SNES_failures"]
            or report["original_Fa_max_l2"] >= FA_TOL
            or report["original_TSCN_max_l2"] >= CN_TOL
            or report["minimum_positive_height_FE_coefficient"] <= 0
            or report["original_volume_balance_defect"] >= 1e-8
            or report["stationary_rotating_state_certified"] is not False):
        raise ValueError(f"Original PETSc/DAE segment is not accepted: {directory}")
    with np.load(directory / "initial_state_trace.npz", allow_pickle=False) as a:
        required = ("times", "states", "full_state_start", "full_state_final",
                    *TRACE_FIELDS)
        if any(k not in a.files for k in required):
            raise ValueError("Native segment lacks complete trajectory/power traces")
        trace = {key: a[key].copy() for key in required}
    t, z = trace["times"], trace["states"]
    if (t.ndim != 1 or z.ndim != 2 or len(t) < 2
            or z.shape != (len(t), int(report["ndof"]))
            or int(report["accepted_steps"]) != len(t)-1
            or any(not np.all(np.isfinite(a)) for a in trace.values())
            or any(trace[key].shape != t.shape for key in TRACE_FIELDS)
            or not np.array_equal(z[0], trace["full_state_start"])
            or not np.array_equal(z[-1], trace["full_state_final"])
            or not np.all(np.diff(t) > 0)):
        raise ValueError(f"Malformed original trajectory: {directory}")
    if (abs(float(t[0])-float(report["physical_time_initial"])) > TIME_TOL
            or abs(float(t[-1])-float(report["physical_time_final"])) > TIME_TOL):
        raise ValueError("Physical clock/report mismatch")
    checkpoint_path = directory / "accepted_state.npz"
    with np.load(checkpoint_path, allow_pickle=False) as cp:
        if not np.array_equal(cp["state"], z[-1]):
            raise ValueError("Production DPHS checkpoint differs from final FE state")
    if expected_phase == "resume":
        initialization = report["initialization"]
        if (initialization.get("no_pseudo_BDF_or_projection") is not True
                or initialization.get(
                    "all_GetFEM_and_GlobalVec_coefficients_restored_bitwise"
                ) is not True
                or initialization.get("PETSc_TS_history_restored") is not False):
            raise ValueError("Cold resume did not preserve exact public-API contract")
    return report, trace


def _measure(name, reference, obtained, tolerance, out):
    a, b = np.asarray(reference), np.asarray(obtained)
    if a.shape != b.shape or not np.all(np.isfinite(a)) or not np.all(np.isfinite(b)):
        raise ValueError(f"{name}: shape or finite check failed")
    delta = float(np.max(np.abs(a-b)))
    out[name + "_maximum_absolute_difference"] = delta
    out[name + "_tolerance"] = tolerance
    if delta > tolerance:
        raise ValueError(
            f"{name}: continuous/split mismatch {delta:.9g} > {tolerance:.9g}"
        )


def compare(seed_dir, continuous_dir, resumed_dir):
    seed, s = load_segment(seed_dir, "seed")
    continuous, c = load_segment(continuous_dir, "seed")
    resumed, r = load_segment(resumed_dir, "resume")
    if not all(doc["gamma"] == .25 for doc in (seed, continuous, resumed)):
        raise ValueError("Not the original gamma=0.25 controlled model")
    if (len({doc["mesh"] for doc in (seed, continuous, resumed)}) != 1
            or len({doc["original_FE_fingerprint_sha256"]
                    for doc in (seed, continuous, resumed)}) != 1
            or len({doc["dt"] for doc in (seed, continuous, resumed)}) != 1):
        raise ValueError("Different model/mesh/FE/control/timestep fingerprint")
    if (abs(float(s["times"][0])) > TIME_TOL
            or abs(float(c["times"][0])) > TIME_TOL
            or abs(float(r["times"][0]-s["times"][-1])) > TIME_TOL
            or abs(float(r["times"][-1]-c["times"][-1])) > TIME_TOL):
        raise ValueError("Not the same physical time range with a pause")
    stitched_times = np.concatenate((s["times"], r["times"][1:]))
    if (len(stitched_times) != len(c["times"])
            or not np.allclose(stitched_times, c["times"],
                               rtol=0., atol=TIME_TOL)):
        raise ValueError("Actual PETSc accepted timestamps disagree")
    if not np.array_equal(s["states"][-1], r["states"][0]):
        raise ValueError("Cold MPI checkpoint not restored bitwise")
    with np.load(Path(seed_dir) / "accepted_state.npz", allow_pickle=False) as cp:
        if not np.array_equal(cp["state"], r["states"][0]):
            raise ValueError("Cold restart not equal to source checkpoint bytes")

    # Compare the complete *trajectory*, not only z(t_final).
    stitched_states = np.concatenate((s["states"],r["states"][1:]), axis=0)
    err_by_time = np.linalg.norm(c["states"]-stitched_states, axis=1)
    err = float(np.max(err_by_time))
    if not np.isfinite(err) or err > FULL_STATE_L2_TOL:
        raise ValueError(f"Full FE trajectory diverged: max L2 = {err:.9g}")

    result = {
        "schema": "scrimp-v130-rotating-swe-checkpoint-parity-v1",
        "status": "PASS_RESTART_EQUIVALENT_TO_UNINTERRUPTED",
        "mesh": seed["mesh"],
        "source_ranks": seed["MPI"],
        "reference_ranks": continuous["MPI"],
        "resume_ranks": resumed["MPI"],
        "time_of_pause": float(s["times"][-1]),
        "final_time": float(c["times"][-1]),
        "total_steps": len(c["times"])-1,
        "all_original_DAE_CN_steps_accepted": True,
        "all_FE_variables_restored_bitwise_at_pause": True,
        "maximum_full_FE_trajectory_L2_difference": err,
        "final_full_FE_L2_difference": float(err_by_time[-1]),
        "original_FE_fingerprint_sha256":
            seed["original_FE_fingerprint_sha256"],
        "checkpoint_API_history_restored": False,
        "physical_energy_defect_observed_not_gated": True,
        "stationary_equilibrium_claimed": False,
    }
    for key in ("H", "volume", "p_normal", "p_viscous", "d_grad", "d_div"):
        stitched = np.concatenate((s[key],r[key][1:]))
        _measure(key, c[key], stitched, OBS_ABS_TOL, result)
    # Restarted integrated work/dissipation start at zero at the pause.
    # Recover global accumulated values while preserving native
    # midpoint physical port and volume samples from each segment.
    for key in ("boundary_work", "viscous_dissipation", "volume_balance"):
        stitched = np.concatenate((s[key],s[key][-1]+r[key][1:]))
        _measure(key, c[key], stitched, OBS_ABS_TOL, result)
    stitched_H = np.concatenate((s["H"],r["H"][1:]))
    stitched_W = np.concatenate((
        s["boundary_work"], s["boundary_work"][-1]+r["boundary_work"][1:]
    ))
    stitched_D = np.concatenate((
        s["viscous_dissipation"],
        s["viscous_dissipation"][-1]+r["viscous_dissipation"][1:]
    ))
    stitched_defect = stitched_H-stitched_H[0]+stitched_D-stitched_W
    _measure("H_plus_D_minus_W", c["power_balance"],
             stitched_defect, OBS_ABS_TOL, result)
    # The defect's MAGNITUDE is not an acceptance criterion. Only its
    # equality in uninterrupted and resumed runs matters for SCRIMP.
    result["continuous_final_energy_defect"] = float(c["power_balance"][-1])
    result["resumed_final_energy_defect"] = float(stitched_defect[-1])
    result["continuous_segment_power_defect_relative_to_exchange"] = float(
        continuous["original_power_defect_over_exchange"]
    )
    result["resumed_segment_power_defect_relative_to_exchange"] = float(
        resumed["original_power_defect_over_exchange"]
    )
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--seed", type=Path, required=True)
    p.add_argument("--continuous", type=Path, required=True)
    p.add_argument("--resume", type=Path, required=True)
    p.add_argument("--report", type=Path, required=True)
    args=p.parse_args()
    report=compare(args.seed,args.continuous,args.resume)
    if args.report.exists():
        raise FileExistsError("Never overwrite original checkpoint-parity evidence")
    args.report.parent.mkdir(parents=True,exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2)+"\n", encoding="utf-8")
    print("PASS ORIGINAL NATIVE ROTATING SWE CHECKPOINT PARITY "
          f"mesh={report['mesh']} "
          f"MPI={report['source_ranks']}->{report['resume_ranks']} "
          f"max_full_state_L2={report['maximum_full_FE_trajectory_L2_difference']:.9g} "
          f"power_defect_preserved={report['H_plus_D_minus_W_maximum_absolute_difference']:.9g}",
          flush=True)


if __name__=="__main__":
    main()
