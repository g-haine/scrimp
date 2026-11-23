# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Independent HOST-ONLY true native continuous-vs-split diagnostic.

Compare a successfully accepted t=0..0.02 seed and a t=0..0.024
continuous control to a t=0.02..0.024 restart diagnostic (which may
have failed a strict local energy-exchange relative gate).

This script DOES NOT certify an energy-inaccurate checkpoint; it
distinguishes a physical-step/integrator defect from a cross-MPI
restart defect. Only uses recorded original FE vectors and original
GetFEM physical energy traces; never solves a synthetic PDE.
"""

import argparse
import json
import math
from pathlib import Path

import numpy as np


def segment(directory):
    path = Path(directory)
    precheck = path / "physical_precheck_trace.npz"
    if precheck.is_file():
        report = json.loads(
            (path / "physical_precheck_report.json").read_text(encoding="utf-8")
        )
        assert report["status"] == "NATIVE_PHYSICAL_DIAGNOSTIC_NOT_A_CHECKPOINT"
        metadata = {
            "phase": report["phase"], "gamma": report["gamma"],
            "mesh": report["mesh"], "mpi": report["ranks"],
            "dt": None,
            "fingerprint": report["source_model_fingerprint"],
            "Fa": report["maximum_actual_original_Fa_residual"],
            "CN": report["maximum_actual_original_CN_residual"],
        }
        filename = precheck
    else:
        report = json.loads(
            (path / "initial_state_report.json").read_text(encoding="utf-8")
        )
        if report["status"] != "ACCEPTED_COMPATIBLE_ROTATING_DAE_INITIAL_STATE_NOT_STEADY":
            raise ValueError("Expected genuinely accepted native seed or continuous")
        metadata = {
            "phase": report["phase"], "gamma": report["gamma"],
            "mesh": report["mesh"], "mpi": report["MPI"],
            "dt": report["dt"],
            "fingerprint": report["original_FE_fingerprint_sha256"],
            "Fa": report["original_Fa_max_l2"],
            "CN": report["original_TSCN_max_l2"],
        }
        filename = path / "initial_state_trace.npz"
    with np.load(filename, allow_pickle=False) as a:
        fields = {key: a[key].copy() for key in (
            "times", "full_state_start", "full_state_final",
            "H", "boundary_work", "viscous_dissipation",
            "power_balance", "volume", "volume_balance"
        )}
    t = fields["times"]
    if (t.ndim != 1 or len(t) < 2 or not np.all(np.isfinite(t))
            or not np.all(np.diff(t) > 0)):
        raise ValueError("Malformed native physical time trace")
    if (any(v.shape != t.shape or not np.all(np.isfinite(v))
            for k, v in fields.items()
            if k not in ("times", "full_state_start", "full_state_final"))
            or fields["full_state_start"].ndim != 1
            or fields["full_state_final"].shape !=
               fields["full_state_start"].shape):
        raise ValueError("Malformed physical states/power arrays")
    return metadata, fields


def physical_energy(fields):
    return (
        float(fields["H"][-1] - fields["H"][0]),
        float(fields["boundary_work"][-1]),
        float(fields["viscous_dissipation"][-1]),
    )


def compare(seed, continuous, resumed):
    sm, s = segment(seed)
    cm, c = segment(continuous)
    rm, r = segment(resumed)
    if [sm["phase"],cm["phase"],rm["phase"]] != ["seed","seed","resume"]:
        raise ValueError("Required original seed, longer continuous seed, restart")
    if (not (sm["gamma"] == cm["gamma"] == rm["gamma"] == .25)
            or not (sm["mesh"] == cm["mesh"] == rm["mesh"])
            or not (sm["fingerprint"] == cm["fingerprint"]
                    == rm["fingerprint"])):
        raise ValueError("Not identical original physical GetFEM model/control")
    dt = np.diff(c["times"])
    if (not np.allclose(dt, dt[0], atol=1e-12, rtol=0.)
            or not np.allclose(np.diff(s["times"]), dt[0],
                               atol=1e-12, rtol=0.)
            or not np.allclose(np.diff(r["times"]), dt[0],
                               atol=1e-12, rtol=0.)
            or abs(float(s["times"][0])) > 1e-12
            or abs(float(c["times"][0])) > 1e-12
            or abs(float(s["times"][-1])-float(r["times"][0])) > 1e-12
            or abs(float(c["times"][-1])-float(r["times"][-1])) > 1e-12):
        raise ValueError("Physical time grids are not continuous and identical")
    if not np.array_equal(s["full_state_final"], r["full_state_start"]):
        raise ValueError("Cold restart did not begin with BITWISE original seed")
    if not np.array_equal(s["full_state_start"], c["full_state_start"]):
        raise ValueError("Continuous and seed initial FE states differ")
    zerr = float(np.linalg.norm(c["full_state_final"] - r["full_state_final"]))
    hjump = float(r["H"][0] - s["H"][-1])
    verr = float(r["volume"][0] - s["volume"][-1])
    hdiff = float(c["H"][-1] - r["H"][-1])
    cw, sw, rw = physical_energy(c), physical_energy(s), physical_energy(r)
    split_h_change = float(r["H"][-1] - s["H"][0])
    split_work = sw[1] + rw[1]
    split_dissipation = sw[2] + rw[2]
    split_defect = split_h_change + split_dissipation - split_work
    continuous_defect = cw[0] + cw[2] - cw[1]
    work_diff = split_work - cw[1]
    diss_diff = split_dissipation - cw[2]
    defect_diff = split_defect - continuous_defect
    result = {
        "status": "NATIVE_RESTART_PARITY_DIAGNOSTIC_NOT_ENERGY_ACCEPTANCE",
        "seed_mpi": sm["mpi"], "continuous_mpi": cm["mpi"],
        "resume_mpi": rm["mpi"], "nominal_mesh": sm["mesh"],
        "source_full_FE_end_state_bitwise_restored": True,
        "continuous_vs_restarted_final_full_FE_l2": zerr,
        "restart_H_jump": hjump, "restart_volume_jump": verr,
        "continuous_vs_restarted_final_H": hdiff,
        "continuous_vs_split_original_boundary_work": work_diff,
        "continuous_vs_split_original_dissipation": diss_diff,
        "continuous_H_plus_D_minus_W": continuous_defect,
        "split_H_plus_D_minus_W": split_defect,
        "continuous_vs_split_energy_defect_difference": defect_diff,
        "resume_original_CN_max": rm["CN"],
        "resume_original_Fa_max": rm["Fa"],
        "physical_segment_relative_energy_gate_not_overridden": True,
    }
    # These thresholds qualify parity only. They do NOT reclassify
    # the failed original physical energy defect as acceptable.
    if (zerr >= 1e-8 or abs(hjump) >= 1e-8
            or abs(verr) >= 1e-8 or abs(hdiff) >= 1e-8
            or abs(work_diff) >= 1e-8 or abs(diss_diff) >= 1e-8
            or abs(defect_diff) >= 1e-8
            or rm["CN"] >= 1e-7 or rm["Fa"] >= 1e-8):
        raise ValueError(
            "Native continuous/split physical parity REJECTED; "
            + json.dumps(result, indent=2)
        )
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--seed", type=Path, required=True)
    p.add_argument("--continuous", type=Path, required=True)
    p.add_argument("--resume", type=Path, required=True)
    args = p.parse_args()
    result = compare(args.seed, args.continuous, args.resume)
    print(json.dumps(result, indent=2))
    print(
        "RESTART PARITY PASS (only). Physical relative energy gate "
        "remains INDEPENDENT; no new checkpoint has been accepted."
    )


if __name__ == "__main__":
    main()
