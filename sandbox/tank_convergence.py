# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3

"""Time-step and mesh refinement study for sandbox/sloshing_tank.py.

Run from the repository root:
    python sandbox/tank_convergence.py
    python sandbox/tank_convergence.py --mesh 0.05 0.025 --dt 0.008 0.004 0.002

The table reports absolute defects (not exact solution errors), on the
same interval [0, t_f]. Temporal convergence should be examined at
fixed mesh, and spatial convergence separately at fixed time step.
"""

import argparse
import csv
import math
from pathlib import Path

from sloshing_tank import sloshing_tank


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mesh", nargs="+", type=float, default=[0.05, 0.025])
    parser.add_argument("--dt", nargs="+", type=float,
                        default=[0.008, 0.004, 0.002])
    parser.add_argument("--final-time", type=float, default=0.2)
    parser.add_argument("--csv", type=Path,
                        default=Path("outputs/tank_convergence.csv"))
    args = parser.parse_args()
    if any(h <= 0 for h in args.mesh) or any(dt <= 0 for dt in args.dt):
        parser.error("Mesh sizes and time steps must be positive")

    rows = []
    for h in args.mesh:
        for dt in args.dt:
            print(f"\n== Coupled tank: h={h:g}, dt={dt:g}, T={args.final_time:g} ==")
            tank = sloshing_tank(
                mesh_size=h,
                final_time=args.final_time,
                time_step=dt,
                save_step=dt,
                adaptive=False,
                plot=False,
            )
            metrics = dict(tank.sloshing_metrics)
            if not math.isclose(
                metrics["last_saved_time"], args.final_time, abs_tol=1e-10
            ):
                raise RuntimeError(
                    f"Last saved time {metrics['last_saved_time']} "
                    f"differs from {args.final_time}; "
                    "check TS final-time option and monitor."
                )
            rows.append(metrics)

    args.csv.parent.mkdir(parents=True, exist_ok=True)
    with args.csv.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    print("\nFixed-step time/space refinement (max absolute defects):")
    print(
        f"{'mesh':>8} {'dt':>9} {'mass':>12} "
        f"{'momentum':>13} {'energy':>13} {'constraint':>13}"
    )
    for row in rows:
        print(
            f"{row['mesh_size']:8.4f} {row['time_step']:9.4g} "
            f"{row['mass_error']:12.3e} "
            f"{row['impulse_exact_error']:13.3e} "
            f"{row['energy_balance_error']:13.3e} "
            f"{row['constraint_error']:13.3e}"
        )

    print("\nConvergence ratios at fixed mesh (energy-balance defect):")
    for h in args.mesh:
        subset = sorted(
            (row for row in rows if row["mesh_size"] == h),
            key=lambda row: row["time_step"],
            reverse=True,
        )
        for coarser, finer in zip(subset, subset[1:]):
            d1, d2 = coarser["time_step"], finer["time_step"]
            e1, e2 = coarser["energy_balance_error"], finer["energy_balance_error"]
            if e2 > 0 and e1 > 0:
                observed = math.log(e1/e2)/math.log(d1/d2)
                print(f"  h={h:g}: dt {d1:g} -> {d2:g}, "
                      f"error ratio {e1/e2:.3g}, observed order {observed:.2f}")
            else:
                print(f"  h={h:g}: zero / invalid defect at dt={d2:g}")
    print(f"\nResults written to {args.csv}")


if __name__ == "__main__":
    main()
