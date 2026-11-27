# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

# SCRIMP v1.3.0 — canonical 2D heat physical diagnostic, not a model rewrite.
# Reuse the unmodified examples/heat.py entry point; collect real GetFEM data.
# No exact BDF2 energy identity is assumed.
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import getfem as gf
from petsc4py import PETSc

from examples.heat import heat_eq


def real_series(seq, name, n):
    values = np.asarray(seq)
    if values.shape != (n,):
        raise RuntimeError(f"{name}: expected {(n,)}, got {values.shape}")
    if np.iscomplexobj(values):
        if np.max(np.abs(values.imag)) > 1e-10:
            raise RuntimeError(f"{name}: significant imaginary component")
        values = values.real
    values = np.asarray(values, dtype=np.float64)
    if not np.all(np.isfinite(values)):
        raise RuntimeError(f"{name}: nonfinite GetFEM values")
    return values


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    comm = PETSc.COMM_WORLD.tompi4py()
    rank, size = comm.rank, comm.size

    # This CALL preserves all original FEMs, controls, time scheme,
    # initial-data expression, and canonical Hamiltonian figure.
    system = heat_eq(time_stepper="bdf", initializer="pseudo")
    time = np.asarray(system.solution["t"], dtype=np.float64)
    states = np.asarray(system.solution["z"])
    n = len(time)
    if (n != 501 or time.shape != (501,)
            or states.shape != (501, 5830)
            or not np.all(np.isfinite(states))
            or not np.all(np.diff(time) > 0)
            or not np.allclose(time, np.arange(501)*0.01, atol=1e-9, rtol=0)):
        raise RuntimeError(f"Unexpected original 500-step heat timeline/FE layout: "
                           f"{time.shape}, {states.shape}")
    if (system.last_ts_converged_reason <= 0
            or system.last_ts_steps != 500
            or system.last_ts_rejections != 0
            or system.last_ts_snes_failures != 0):
        raise RuntimeError("Canonical PETSc BDF2 convergence gate failed")
    if system.time_scheme.getString("ts_type", "").lower() != "bdf":
        raise RuntimeError("Original heat example no longer uses BDF")

    H = real_series(system.get_Hamiltonian(), "H", n)
    # Independently evaluate original scalar energy and physical quadratic
    # Fourier dissipation on every accepted native FE state.
    H_integral = real_series(
        0.5*np.asarray(system.get_quantity("T*rho*T", region=1)),
        "independent FE Hamiltonian", n
    )
    dissipation = real_series(
        system.get_quantity("f_Q.Lambda.f_Q", region=1),
        "physical Fourier dissipation", n
    )
    flux_dissipation = real_series(
        system.get_quantity("f_Q.J_Q", region=1),
        "actual constitutive flux-product", n
    )

    # Physical outward boundary flux for the ENERGY balance is
    # B_r = int_{Gamma_r} T*(J_Q dot n) ds. Signs are uniform here.
    # It is distinct from the signed u*y convention of a control port,
    # notably on the Neumann boundary (Gamma_Left).
    bflux = {
        region: real_series(
            system.get_quantity("T*(J_Q.Normal)", region=region),
            f"physical outward heat boundary flux region {region}", n
        ) for region in (10, 11, 12, 13)
    }
    # The default plot computed generic control-port powers on all ranks.
    raw_port_powers = {}
    for name, port in system.ports.items():
        if port.get_algebraic() and not port.get_substituted():
            if port.get_power_sampling() != "instantaneous":
                raise RuntimeError(f"Unexpected port-power sampling on {name}")
            raw_port_powers[name] = real_series(
                port.get_power(), "port power "+name, n
            )

    # Independent diagnostic: trapezoidal quadrature on accepted BDF
    # endpoint data is APPROXIMATE, never an exact BDF2 identity.
    boundary_total = np.sum(np.stack(list(bflux.values())), axis=0)
    q_total = dissipation + boundary_total
    increments = 0.5*np.diff(time)*(q_total[1:] + q_total[:-1])
    accum = np.r_[0., np.cumsum(increments)]
    defect = H - H[0] + accum
    energy_scale = max(
        1e-12, abs(float(H[-1]-H[0])),
        abs(float(np.sum(0.5*np.diff(time)*(dissipation[1:]+dissipation[:-1])))),
        abs(float(accum[-1]))
    )

    # Full accepted states support independent cross-MPI comparison, including
    # initial physical vs pseudo-initialized values; do not force equality
    # of the user's symbolic initial function after pseudo consistency.
    if np.iscomplexobj(states):
        if np.max(np.abs(states.imag)) > 1e-10:
            raise RuntimeError("Complex FE values unexpectedly significant")
        states = states.real
    states = np.ascontiguousarray(states, dtype=np.float64)
    if np.min(H) <= 0 or np.min(dissipation) < -1e-10:
        raise RuntimeError("Nonphysical Hamiltonian or negative Fourier dissipation")

    # Correct discrete pH power: the flow/effort products on the control
    # ports, not the pointwise T*J_Q.normal traces. DG T need not equal U_D
    # on the Dirichlet facets; U_L*Y_L is the Neumann contribution.
    if ("Heat flux" not in raw_port_powers or
            len(raw_port_powers) != 5):
        raise RuntimeError("Unexpected original heat power-port topology")
    control_power = np.sum(np.stack([
        val for key, val in raw_port_powers.items()
        if key != "Heat flux"
    ]), axis=0)
    total_discrete_power = dissipation + control_power
    port_work = np.r_[0., np.cumsum(
        0.5*np.diff(time)*(total_discrete_power[1:]
                          + total_discrete_power[:-1])
    )]
    port_balance_approx = H - H[0] + port_work

    # Stronger exact BDF2 G-stability identity for a *quadratic* H:
    # (a, 3a-4b+c)_M = H(a)-H(b)
    #   + H(2a-b)-H(2b-c)+H(a-2b+c).
    # Hence the original BDF2 weak equations imply
    #   [H(a)-H(b)+H(2a-b)-H(2b-c)+H(a-2b+c)]/(2dt)
    #   + D(a) + SUM P_control(a) = 0
    # for each step n>=2. Every H below is a REAL GetFEM mass integral
    # at an authentic combination of the original accepted FE vectors.
    # It does not reuse the approximate trapezoidal time quadrature.
    def quadratic_hamiltonian(state, current_t):
        system.gf_model.set_time(float(current_t))
        system.gf_model.to_variables(state)
        return float(gf.asm_generic(
            system.domain._int_method[0], 0, "0.5*T*rho*T", 1,
            system.gf_model
        ))

    bdf2_power_balance_discrete = np.full(n, np.nan)
    extrap_prev = quadratic_hamiltonian(
        2*states[1]-states[0], time[1]
    )
    # PETSc BDF2 starts with its one-step implicit formula.
    bdf1_first_step_balance = (
        (H[1]-H[0] +
         quadratic_hamiltonian(states[1]-states[0], time[1]))
        /(time[1]-time[0]) + total_discrete_power[1]
    )
    for k in range(2, n):
        extrap_current = quadratic_hamiltonian(
            2*states[k]-states[k-1], time[k]
        )
        second_diff = quadratic_hamiltonian(
            states[k]-2*states[k-1]+states[k-2], time[k]
        )
        bdf2_power_balance_discrete[k] = (
            (H[k]-H[k-1]+extrap_current-extrap_prev+second_diff)
            /(2*(time[k]-time[k-1]))
            + total_discrete_power[k]
        )
        extrap_prev = extrap_current
    if not np.all(np.isfinite(bdf2_power_balance_discrete[2:])):
        raise RuntimeError("Nonfinite native BDF2 quadratic power identity")

    args.out.mkdir(parents=True, exist_ok=True)
    trace = args.out/f"rank{rank}.npz"
    report = args.out/f"rank{rank}.json"
    if trace.exists() or report.exists():
        raise FileExistsError("Refuse to overwrite original native heat evidence")

    np.savez_compressed(
        trace, time=time, states=states, H=H,
        H_independent=H_integral, dissipation=dissipation,
        flux_dissipation=flux_dissipation,
        boundary_10=bflux[10], boundary_11=bflux[11],
        boundary_12=bflux[12], boundary_13=bflux[13],
        energy_defect_approximate=defect,
        port_balance_approximate=port_balance_approx,
        control_power=control_power,
        discrete_bdf2_power_defect=bdf2_power_balance_discrete,
        **{"port_"+str(i): v for i, v in enumerate(raw_port_powers.values())}
    )
    record = {
        "schema": "scrimp-v130-heat-original-physical-diagnostic-v1",
        "status": "NATIVE_HEAT_PHYSICS_DIAGNOSTIC_ONLY",
        "original_model": "examples/heat.py",
        "original_FEM_bricks_controls_scheme_unchanged": True,
        "rank": rank, "MPI": size, "accepted_steps": int(system.last_ts_steps),
        "SNES_failures": int(system.last_ts_snes_failures),
        "rejected_steps": int(system.last_ts_rejections),
        "total_DOFs": int(states.shape[1]),
        "final_time": float(time[-1]),
        "H_initial": float(H[0]), "H_final": float(H[-1]),
        "max_H_independent_difference": float(np.max(np.abs(H-H_integral))),
        "min_Fourier_dissipation": float(np.min(dissipation)),
        "max_constitutive_dissipation_difference": float(
            np.max(np.abs(dissipation-flux_dissipation))),
        "outward_boundary_flux_integral_final": {
            str(k):float(v[-1]) for k, v in bflux.items()
        },
        "generic_port_order": list(raw_port_powers),
        "max_BDF_trapezoid_energy_defect": float(np.max(np.abs(defect))),
        "final_BDF_trapezoid_energy_defect": float(defect[-1]),
        "relative_BDF_trapezoid_defect": float(np.max(np.abs(defect))/energy_scale),
        "max_correct_port_trapezoid_energy_defect": float(
            np.max(np.abs(port_balance_approx))),
        "final_correct_port_trapezoid_energy_defect": float(
            port_balance_approx[-1]),
        "max_exact_BDF2_quadratic_power_defect": float(
            np.max(np.abs(bdf2_power_balance_discrete[2:]))),
        "rms_exact_BDF2_quadratic_power_defect": float(
            np.sqrt(np.mean(bdf2_power_balance_discrete[2:]**2))),
        "BDF1_initial_step_power_defect": float(
            bdf1_first_step_balance),
        "BDF2_uses_exact_quadratic_mass_identity": True,
        "BDF_trapezoid_is_exact_energy_identity": False,
        "energy_balance_physically_accepted": False,
        "states_sha256": hashlib.sha256(states.tobytes()).hexdigest()
    }
    report.write_text(json.dumps(record, indent=2)+"\n", encoding="utf-8")
    print("HEAT PHYSICAL DIAGNOSTIC",json.dumps({
        k:record[k] for k in (
            "rank","MPI","accepted_steps","H_initial","H_final",
            "max_H_independent_difference",
            "max_constitutive_dissipation_difference",
            "max_BDF_trapezoid_energy_defect",
            "relative_BDF_trapezoid_defect",
            "max_correct_port_trapezoid_energy_defect",
            "max_exact_BDF2_quadratic_power_defect",
            "BDF1_initial_step_power_defect")
    }), flush=True)


if __name__ == "__main__":
    main()
