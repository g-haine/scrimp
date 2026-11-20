# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Strict native MPI parity and energy+boundary+damping audit.

The second DAE is a different dissipative port-Hamiltonian model:
verify its original GetFEM weak constitutive terms, nonzero dissipated
work, signed boundary powers, actual CN evolution, and 1/2-rank parity.
Never treat source-only tests as native scientific evidence.
"""

import argparse
import json
from pathlib import Path

import numpy as np


def _load(path, ranks):
    stem = f"dissipative_core_{ranks}r"
    record = json.loads(
        (path / f"{stem}.json").read_text(encoding="utf-8")
    )
    with np.load(path / f"{stem}.npz", allow_pickle=False) as f:
        traj = {k: np.asarray(f[k]).copy() for k in f.files}
    return record, traj


def _scaled(a, b):
    if a.shape != b.shape:
        raise RuntimeError(f"Incompatible vector shapes: {a.shape} {b.shape}")
    return float(np.linalg.norm(a-b)/max(1., np.linalg.norm(a)))


def compare(path):
    records = {}
    data = {}
    for ranks in (1, 2):
        info, traj = _load(path, ranks)
        if info["rank_count"] != ranks or info["steps"] != 10:
            raise RuntimeError("Wrong MPI/step identity for dissipative DAE")
        report = info["dae_projection"]
        if (
            report["strategy"] != "algebraic"
            or report["status"] != "PASS"
            or not report["differential_initial_exact"]
            or report["pseudo_ts_executed"]
            or not report["accuracy_target_met"]
            or report["native_final_residual_l2"] >
               report["feasibility_tolerance"]
        ):
            raise RuntimeError(f"Algebraic DAE FAIL MPI={ranks}")
        if any(info["initial_dynamic_drift"].values()):
            raise RuntimeError(
                f"Prescribed dissipative initial states changed MPI={ranks}"
            )
        t = traj["times"]
        expected = np.arange(11)*.002
        if not np.allclose(t, expected, rtol=0, atol=1e-12):
            raise RuntimeError("Dissipative time grid not full fixed CN")
        physics = info["physical"]
        if (
            physics["absolute_dissipation_work"] < 1.e-9
            or physics["integrated_dissipation"] < 0
            or physics["absolute_signed_boundary_work"] < 1.e-9
            or physics["max_absolute_balance_defect"]
               > info["power_tolerance"]
        ):
            raise RuntimeError(
                f"Non-vacuous dissipative pH power gate FAIL MPI={ranks}"
            )
        records[ranks], data[ranks] = info, traj
    a,b=data[1],data[2]
    if not np.array_equal(a["initial"],b["initial"]):
        raise RuntimeError("Prescribed dissipative initial data not MPI identical")
    if not np.allclose(
        a["times"],b["times"], rtol=0,atol=1.e-12
    ):
        raise RuntimeError("Different dissipative MPI time grids")
    diffs = [
        _scaled(x,y) for x,y in zip(a["states"],b["states"])
    ]
    maximum=max(diffs)
    if maximum>5e-8:
        raise RuntimeError(
            "Dissipative complete CN 1/2-rank trajectory mismatch "
            f"{maximum:.9g} > 5e-8"
        )
    p1,p2=records[1]["physical"],records[2]["physical"]
    for key in ("hamiltonian", "balance_defects", "step_defects"):
        if not np.allclose(
            p1[key],p2[key],rtol=5e-8,atol=1e-9
        ):
            raise RuntimeError(f"Independent dissipative {key} MPI mismatch")
    for name in p1["midpoint_port_powers"]:
        if name not in p2["midpoint_port_powers"] or not np.allclose(
            p1["midpoint_port_powers"][name],
            p2["midpoint_port_powers"][name],
            rtol=5e-8,atol=1e-9
        ):
            raise RuntimeError(
                f"Dissipative port {name} MPI power mismatch"
            )
    summary={
        "schema":"scrimp-v130-independent-dissipative-linear-dae-v1",
        "scientific_status":"PASS",
        "model":"original examples/wave_dissipative.py with compatible p0",
        "qualified_scope":"linear dissipative DAE on h=0.2; CN .002 x 10; MPI 1/2",
        "native_residual_1r":records[1]["dae_projection"][
            "native_final_residual_l2"
        ],
        "native_residual_2r":records[2]["dae_projection"][
            "native_final_residual_l2"
        ],
        "prescribed_dynamic_state_preserved":True,
        "max_scaled_full_CN_MPI_discrepancy":maximum,
        "max_absolute_signed_power_defect_1r":p1[
            "max_absolute_balance_defect"
        ],
        "max_absolute_signed_power_defect_2r":p2[
            "max_absolute_balance_defect"
        ],
        "signed_boundary_work_1r":p1[
            "integrated_signed_boundary_work"
        ],
        "signed_dissipation_1r":p1["integrated_dissipation"],
        "source_only_substitute":False,
    }
    (path/"second_linear_dae_comparison.json").write_text(
        json.dumps(summary,indent=2)+"\n",encoding="utf-8"
    )
    print("SCRIMP second independent dissipative linear pH-DAE native audit")
    print("| MPI | true F_a(t0) | dissipation work | signed boundary work "
          "| max |H-H0+Wdiss+Wbnd| |")
    print("|---:|---:|---:|---:|---:|")
    for ranks in (1,2):
        m=records[ranks]
        phys=m["physical"]
        print(f"| {ranks} | "
              f"{m['dae_projection']['native_final_residual_l2']:.9g} | "
              f"{phys['integrated_dissipation']:.9g} | "
              f"{phys['integrated_signed_boundary_work']:.9g} | "
              f"{phys['max_absolute_balance_defect']:.9g} |")
    print(f"PASS: second linear DAE, positive true dissipation, "
          f"nonzero boundary work; full CN MPI parity {maximum:.9g}")
    return summary


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("folder",type=Path)
    args=p.parse_args()
    compare(args.folder)


if __name__=="__main__":
    main()
