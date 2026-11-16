# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Audit the discrete Hamiltonian/port-power identity of the real 2D wave.

The wave example is linear, non-dissipative, and uses Crank-Nicolson.
Its quadratic Hamiltonian and algebraic control ports should satisfy,
at each consecutive saved step, the midpoint discrete identity

    H[n] - H[n-1] + dt * sum_ports (flow*effort)[n-1/2] = 0.

The PLUS sign comes directly from DPHS.plot_powers(): all four wave
control ports are flagged dissipative=True and contribute +integrated
flow*effort to the Hamiltonian balance. This is NOT internal damping.

Hamiltonian terms are the very same ones in examples/wave.py.
All GetFEM scalar assemblies are evaluated on every rank and their
rank-to-rank equality is checked rather than silently reduced or /P'd.
The optional left-boundary drive overrides only a prescribed control,
without changing the wave PDE or variational bricks.

Use identical mesh/dt/final-time/drive for MPI 1,2,4; compare the
JSON records with analyze_acoustic_power_balance.py.
"""
import argparse
import json
from pathlib import Path
from unittest.mock import patch

import getfem as gf
from mpi4py import MPI
import numpy as np
import scrimp as S

from sandbox.bench_acoustic_2d import _build_wave


def _build_with_left_drive(mesh, left_drive):
    if left_drive == "published":
        return _build_wave(mesh)

    original = S.DPHS.set_control

    def overridden_set_control(system, name, expression):
        if name == "Boundary control (left)":
            expression = left_drive
        return original(system, name, expression)

    # Only the published example's prescribed boundary signal is changed.
    # The model, FE choices, port conventions and variational forms remain.
    with patch.object(S.DPHS, "set_control", overridden_set_control):
        return _build_wave(mesh)


def _global_scalar(value, mpi, description):
    """Reject rank-local or duplicated scalar assembly without guessing."""
    value = complex(np.asarray(value).item())
    if not np.isfinite(value.real) or not np.isfinite(value.imag):
        raise RuntimeError(f"{description}: non-finite GetFEM integral {value}")
    vals = mpi.allgather(value)
    if any(not np.isclose(x, vals[0], atol=1e-10, rtol=1e-10)
           for x in vals[1:]):
        raise RuntimeError(
            f"{description}: GetFEM scalar integral differs per rank: {vals}. "
            "Must classify native assembly before applying any MPI reduction."
        )
    if abs(value.imag) > 1e-10 * max(1.0, abs(value.real)):
        raise RuntimeError(f"{description}: unexpectedly complex value {value}")
    return float(value.real)


def _audit(system, mpi):
    t = np.asarray(system.solution["t"], dtype=float)
    z = np.asarray(system.solution["z"])
    if t.ndim != 1 or t.size < 2 or len(z) != len(t):
        raise RuntimeError("Need at least two matching saved time/state snapshots")
    if not np.all(np.diff(t) > 0):
        raise RuntimeError("Saved time grid must be strictly increasing")

    # Copy the actual wave.py Hamiltonian terms (created there after solve).
    # Since we intercepted wave.solve(), they are registered here instead.
    system.hamiltonian.set_name("Mechanical energy")
    system.hamiltonian.add_term(
        S.Term("Potential energy", "0.5*q.T.q", [1])
    )
    system.hamiltonian.add_term(
        S.Term("Kinetic energy", "0.5*p*p/rho", [1])
    )
    H_native = system.get_Hamiltonian()
    H = np.array([
        _global_scalar(h, mpi, f"Hamiltonian at t={tk:g}")
        for tk, h in zip(t, H_native)
    ])

    # Following DPHS.plot_powers(), these algebraic, non-substituted
    # control ports all have dissipative=True and contribute with the
    # PLUS sign to H + int(flow*effort). No other dissipative mechanism.
    ports = list(system.controls.values())
    if len(ports) != 4:
        raise RuntimeError(f"Expected four acoustic boundary ports, got {len(ports)}")
    for port in ports:
        if (not port.get_algebraic() or port.get_substituted()
                or not port.get_dissipative()):
            raise RuntimeError(
                f"Unexpected port power convention for {port.get_name()}"
            )

    mim = system.domain._int_method
    model = system.gf_model
    powers = {p.get_name(): [] for p in ports}
    balance = [0.0]
    cumulative_port_work = [0.0]
    step_defects = []
    for k in range(1, t.size):
        tm = 0.5 * (t[k] + t[k - 1])
        zm = 0.5 * (z[k] + z[k - 1])
        model.set_time(float(tm))
        model.to_variables(zm)
        p_sum = 0.0
        for port in ports:
            expression = port.get_flow() + "*" + port.get_effort()
            value = gf.asm(
                "generic", mim[port.get_mesh_id()], 0,
                expression, port.get_region(), model
            )
            scalar = _global_scalar(
                value, mpi, f"{port.get_name()} power at t={tm:g}"
            )
            powers[port.get_name()].append(scalar)
            p_sum += scalar

        dt = float(t[k] - t[k - 1])
        w = dt * p_sum
        cumulative_port_work.append(cumulative_port_work[-1] + w)
        defect = float(H[k] - H[k - 1] + w)
        step_defects.append(defect)
        balance.append(float(H[k] - H[0] + cumulative_port_work[-1]))

    H_scale = max(float(np.max(np.abs(H))),
                  float(sum(abs(np.diff(cumulative_port_work)))),
                  1e-12)
    result = {
        "times": t.tolist(),
        "hamiltonian": H.tolist(),
        "hamiltonian_terms": {
            term.get_description(): [
                float(x) for x in term.get_values()
            ] for term in system.hamiltonian.get_terms()
        },
        "boundary_midpoint_powers": powers,
        "port_work_integral": cumulative_port_work,
        "balance_minus_initial": balance,
        "step_balance_defects": step_defects,
        "maximum_absolute_balance_defect": float(np.max(np.abs(balance))),
        "maximum_absolute_step_defect": float(
            np.max(np.abs(step_defects))
        ),
        "balance_reference_scale": H_scale,
        "relative_balance_defect": float(
            np.max(np.abs(balance)) / H_scale
        ),
        "absolute_port_work": float(
            sum(abs(np.diff(cumulative_port_work)))
        ),
        "energy_change": float(H[-1] - H[0]),
    }
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mesh", type=float, default=0.2)
    parser.add_argument("--dt", type=float, default=0.002)
    parser.add_argument("--final-time", type=float, default=0.02)
    parser.add_argument(
        "--drive-left", default="0.08*sin(20*t)",
        help="'published' for the original near-zero pulse, or a "
             "GetFEM expression for the left boundary U_L",
    )
    parser.add_argument("--output", type=Path)
    parser.add_argument("--rtol", type=float, default=1e-5)
    parser.add_argument("--atol", type=float, default=1e-9)
    parser.add_argument(
        "--diagnose-only", action="store_true",
        help="Record raw energy/power evidence without balance gate; "
             "for first-time audits on an unvalidated GetFEM build",
    )
    args = parser.parse_args()
    if min(args.mesh, args.dt, args.final_time) <= 0:
        parser.error("mesh, dt and final time must be positive")
    if args.rtol < 0 or args.atol < 0:
        parser.error("rtol/atol must be non-negative")

    mpi = MPI.COMM_WORLD
    system = _build_with_left_drive(args.mesh, args.drive_left)
    system.set_time_scheme(
        ts_type="cn", ts_adapt_type="none",
        ts_exact_final_time="matchstep",
        dt=args.dt, dt_save=args.dt, t_f=args.final_time,
        ksp_type="preonly", pc_type="lu",
        pc_factor_mat_solver_type="mumps",
    )
    system.solve()
    audit = _audit(system, mpi)
    audit.update({
        "mpi_size": mpi.Get_size(),
        "degrees_of_freedom": int(system.gf_model.nbdof()),
        "rectangle_mesh_h": args.mesh,
        "dt": args.dt,
        "requested_final_time": args.final_time,
        "left_boundary_drive": args.drive_left,
        "ts_converged_reason": int(system.last_ts_converged_reason),
        "ts_steps": int(system.last_ts_steps),
        "ts_snes_failures": int(system.last_ts_snes_failures),
        "balance_tolerance": (
            args.atol + args.rtol * audit["balance_reference_scale"]
        ),
    })
    if mpi.Get_rank() == 0:
        output = json.dumps(audit, indent=2)
        print(output, flush=True)
        if args.output:
            args.output.parent.mkdir(parents=True, exist_ok=True)
            args.output.write_text(output + "\n", encoding="utf-8")

    if args.diagnose_only:
        return
    if int(system.last_ts_converged_reason) <= 0:
        raise RuntimeError("Acoustic TS failed to converge")
    if audit["maximum_absolute_balance_defect"] > audit["balance_tolerance"]:
        raise RuntimeError(
            "Discrete wave Hamiltonian + boundary-work balance FAILED: "
            f"defect={audit['maximum_absolute_balance_defect']:.6g}, "
            f"tolerance={audit['balance_tolerance']:.6g}, "
            f"H change={audit['energy_change']:.6g}, "
            f"port work={audit['port_work_integral'][-1]:.6g}"
        )
    if args.drive_left != "published" and audit["absolute_port_work"] < 1e-9:
        raise RuntimeError(
            "Prescribed left-boundary drive did not generate a measurable "
            "boundary-power transfer; this would be a vacuous energy check"
        )


if __name__ == "__main__":
    main()
