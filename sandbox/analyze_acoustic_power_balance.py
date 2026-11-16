# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Compare the SAME driven wave's Hamiltonian and port power at 1/2/4 MPI ranks.

This is a numerical-invariance test, not a performance benchmark. The
per-run discrete balance is checked by check_acoustic_energy_mpi.py.
Require pointwise Hamiltonian and signed boundary powers to agree
across all ranks, not just the final state norms.
"""
import argparse
import json
from pathlib import Path

import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reports", nargs="+", type=Path)
    parser.add_argument("--rtol", type=float, default=1e-5)
    parser.add_argument("--atol", type=float, default=1e-9)
    args = parser.parse_args()

    records = {
        int(report["mpi_size"]): report
        for report in (
            json.loads(p.read_text(encoding="utf-8")) for p in args.reports
        )
    }
    if set(records) != {1, 2, 4} or len(args.reports) != 3:
        parser.error("Need exactly one report each for MPI 1, 2 and 4")
    ref = records[1]
    config = (
        "degrees_of_freedom", "rectangle_mesh_h", "dt",
        "requested_final_time", "left_boundary_drive",
    )
    for ranks, data in sorted(records.items()):
        for key in config:
            if data[key] != ref[key]:
                raise SystemExit(
                    f"P={ranks}: changed {key}: {data[key]} != {ref[key]}"
                )
        if data["ts_converged_reason"] <= 0 or data["ts_snes_failures"]:
            raise SystemExit(f"P={ranks}: time integration did not converge")
        if data["maximum_absolute_balance_defect"] > data["balance_tolerance"]:
            raise SystemExit(f"P={ranks}: wave power identity violated")
        keys = ("times", "hamiltonian",
                "port_work_integral", "step_balance_defects")
        for key in keys:
            if not np.allclose(data[key], ref[key],
                               rtol=args.rtol, atol=args.atol):
                raise SystemExit(
                    f"P={ranks}: {key} differs from P=1. "
                    f"max error={np.max(np.abs(np.asarray(data[key]) - np.asarray(ref[key]))):.3g}"
                )
        if data["boundary_midpoint_powers"].keys() != ref["boundary_midpoint_powers"].keys():
            raise SystemExit(f"P={ranks}: boundary port set changed")
        for name, series in data["boundary_midpoint_powers"].items():
            if not np.allclose(
                series, ref["boundary_midpoint_powers"][name],
                rtol=args.rtol, atol=args.atol,
            ):
                raise SystemExit(
                    f"P={ranks}: boundary port {name} power differs from P=1"
                )

    print(
        f"Validated driven acoustic CN Hamiltonian + boundary power at 1/2/4 "
        f"ranks, DOFs={ref['degrees_of_freedom']}, "
        f"dt={ref['dt']}, t_f={ref['requested_final_time']}"
    )
    print("| MPI ranks | ΔH | Integrated signed port power | Max |H-H0+W| | Relative defect |")
    print("|---:|---:|---:|---:|---:|")
    for p, data in sorted(records.items()):
        print(
            f"| {p} | {data['energy_change']:.8e} | "
            f"{data['port_work_integral'][-1]:.8e} | "
            f"{data['maximum_absolute_balance_defect']:.3e} | "
            f"{data['relative_balance_defect']:.3e} |"
        )


if __name__ == "__main__":
    main()
