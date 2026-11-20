# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Issue #6: inspect true MPI 1/2/4 viscous SWE state, power and output evidence.

This is a deliberately bounded short smoke, not a proof of long-time
paper reproduction, pH energy preservation or absence of race conditions
in every other output path. Refuses missing/nonfinite evidence and
rank-dependent generic port powers. Physical balance uses the h-weighted
boundary powers and both viscous dissipations written by the driver.
"""

import argparse
import json
from pathlib import Path

import numpy as np

RANKS = (1, 2, 4)
ENERGY_KEYS = (
    "t", "hamiltonian", "kinetic", "potential",
    "diss_grad", "diss_div", "power_normal", "power_viscous",
    "supplied", "dissipated", "balance", "relative_balance_exchange",
    "volume", "cumulative_outflow", "volume_balance",
)


def _finite_array(v, label):
    x = np.asarray(v, dtype=float)
    if not x.size or not np.all(np.isfinite(x)):
        raise ValueError(f"{label}: missing/nonfinite samples")
    return x


def _assert_close(a, b, label, rtol=2e-5, atol=1e-8):
    if a.shape != b.shape or not np.allclose(a, b, rtol=rtol, atol=atol):
        delta = (float(np.max(np.abs(a-b))) if a.shape == b.shape else float("inf"))
        raise ValueError(f"{label}: serial/MPI or inter-rank discrepancy, max |delta|={delta:.9g}")


def _case(output, name, ranks, tf, exchange_limit, volume_limit, activity):
    root = output / f"{name}_{ranks}r"
    metrics = json.loads((root / "metrics.json").read_text(encoding="utf-8"))
    with np.load(root / "energies.npz") as archive:
        energy = {key: _finite_array(archive[key], f"{name}/{ranks}r/{key}")
                  for key in ENERGY_KEYS}
    with np.load(root / "state_endpoints.npz") as archive:
        states = {key: _finite_array(archive[key], f"{name}/{ranks}r/{key}")
                  for key in ("initial", "final")}
    t = energy["t"]
    if t.ndim != 1 or len(t) < 2 or np.any(np.diff(t) <= 0):
        raise ValueError(f"{name}/{ranks}r: invalid time grid")
    if abs(t[0]) > 1e-12 or abs(t[-1] - tf) > 1e-9:
        raise ValueError(f"{name}/{ranks}r: incorrect final time {t[-1]}")
    if any(v.shape != t.shape for v in energy.values()):
        raise ValueError(f"{name}/{ranks}r: incompatible physical time series")
    if int(metrics["mpi_ranks"]) != ranks:
        raise ValueError(f"{name}/{ranks}r: incorrect world-size provenance")
    if (int(metrics["petsc_ts_converged_reason"]) <= 0
            or int(metrics["petsc_ts_accepted_steps"]) <= 0
            or int(metrics["petsc_ts_rejected_steps"]) != 0
            or int(metrics["petsc_ts_snes_failures"]) != 0):
        raise ValueError(f"{name}/{ranks}r: incomplete or failed PETSc TS")
    if float(metrics["minimum_height_coefficient"]) <= 0:
        raise ValueError(f"{name}/{ranks}r: nonpositive sampled h DOF")
    for key in ("diss_grad", "diss_div"):
        if float(np.min(energy[key])) < -1e-10:
            raise ValueError(f"{name}/{ranks}r: negative viscous dissipation ({key})")
    h, work, diss = (energy[k] for k in ("hamiltonian", "supplied", "dissipated"))
    reconstructed = h - h[0] + diss - work
    _assert_close(energy["balance"], reconstructed,
                  f"{name}/{ranks}r: signed physical power identity",
                  rtol=1e-11, atol=1e-8)
    eb = float(metrics["maximum_absolute_balance_error_over_actual_exchange"])
    vb = float(metrics["maximum_absolute_relative_volume_balance_error"])
    if not np.isfinite(eb) or eb > exchange_limit:
        raise ValueError(f"{name}/{ranks}r: physical power defect {eb:.9g} > {exchange_limit}")
    if not np.isfinite(vb) or vb > volume_limit:
        raise ValueError(f"{name}/{ranks}r: volume defect {vb:.9g} > {volume_limit}")
    if activity:
        magnitude = max(float(np.max(np.abs(energy[key])))
                        for key in ("diss_grad", "diss_div", "power_normal",
                                    "power_viscous"))
        if magnitude <= 1e-10:
            raise ValueError(f"{name}/{ranks}r: all physical powers are vacuous")
    powers = []
    for rank in range(ranks):
        path = root / f"mpi_port_powers_rank{rank}.json"
        doc = json.loads(path.read_text(encoding="utf-8"))
        if int(doc["rank"]) != rank or int(doc["world_size"]) != ranks:
            raise ValueError(f"{name}/{ranks}r: rank provenance mismatch")
        time_rank = _finite_array(doc["times"], f"{name}/{ranks}r/rank{rank}/t")
        _assert_close(time_rank, t, f"{name}/{ranks}r/rank{rank} times",
                      rtol=0, atol=1e-12)
        port = {key: _finite_array(v, f"{name}/{ranks}r/{rank}/{key}")
                for key, v in doc["powers"].items()}
        if not port:
            raise ValueError(f"{name}/{ranks}r: no generic algebraic ports")
        if any(v.shape != t.shape for v in port.values()):
            raise ValueError(f"{name}/{ranks}r: incomplete generic port traces")
        powers.append(port)
    for rank, power in enumerate(powers[1:], start=1):
        if power.keys() != powers[0].keys():
            raise ValueError(f"{name}/{ranks}r: rank-{rank} port names differ")
        for key in powers[0]:
            _assert_close(power[key], powers[0][key],
                          f"{name}/{ranks}r: rank-{rank} generic port {key}")
    return {"energy": energy, "states": states, "metrics": metrics,
            "ports": powers[0], "max_exchange_defect": eb,
            "max_volume_defect": vb}


def compare(root, cases, final_time, exchange_limit, volume_limit, require_activity):
    results = {}
    for case in cases:
        for ranks in RANKS:
            results[case, ranks] = _case(
                root, case, ranks, final_time, exchange_limit,
                volume_limit, require_activity,
            )
    summary = {"schema": "scrimp-v130-issue6-viscous-mpi-smoke-v1",
               "status": "PASS", "source": "native outputs",
               "cases": cases, "mpi_ranks": list(RANKS),
               "requires_further_long_time_qualification": True,
               "exact_discrete_power_conservation_claimed": False,
               "comparisons": {}}
    for case in cases:
        ref = results[case, 1]
        for ranks in (2, 4):
            cur = results[case, ranks]
            for key in ENERGY_KEYS:
                _assert_close(cur["energy"][key], ref["energy"][key],
                              f"{case} {ranks}r vs 1r physical {key}")
            for key in ("initial", "final"):
                # Compare full coefficients, not only energetics.
                _assert_close(cur["states"][key], ref["states"][key],
                              f"{case} {ranks}r vs 1r full {key}",
                              rtol=5e-7, atol=1e-8)
            if cur["ports"].keys() != ref["ports"].keys():
                raise ValueError(f"{case}: MPI {ranks}r port names differ")
            for port in ref["ports"]:
                _assert_close(cur["ports"][port], ref["ports"][port],
                              f"{case} {ranks}r vs 1r generic port {port}")
        summary["comparisons"][case] = {
            "max_exchange_scaled_power_balance_defect":
                max(results[case, n]["max_exchange_defect"] for n in RANKS),
            "max_relative_volume_balance_defect":
                max(results[case, n]["max_volume_defect"] for n in RANKS),
            "mpi_full_state_initial_final_parity": True,
            "generic_port_power_rank_parity": True,
            "h_dof_positive_at_saved_times": True,
        }
    path = root / "issue6_mpi_comparison.json"
    path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print("PASS: issue #6 bounded native MPI 1/2/4 viscous SWE smoke")
    print(json.dumps(summary["comparisons"], indent=2))
    print("NOT a long-time physical release qualification.")
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--cases", nargs="+", choices=("open", "closed"),
                        default=["open", "closed"])
    parser.add_argument("--final-time", type=float, required=True)
    parser.add_argument("--max-exchange-defect", type=float, default=0.05)
    parser.add_argument("--max-relative-volume-defect", type=float, default=1e-4)
    parser.add_argument("--require-activity", action="store_true")
    args = parser.parse_args()
    if (args.final_time <= 0 or args.max_exchange_defect <= 0
            or args.max_relative_volume_defect <= 0):
        parser.error("time and defect limits must be positive")
    compare(args.output_dir, args.cases, args.final_time,
            args.max_exchange_defect, args.max_relative_volume_defect,
            args.require_activity)


if __name__ == "__main__":
    main()
