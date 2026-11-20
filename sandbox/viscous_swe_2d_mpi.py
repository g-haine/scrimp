# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Experimental MPI benchmark: 2D viscous rotational Saint-Venant equations.

Circular basin; open-loop and normal-output-feedback control based on
Cardoso-Ribeiro, Lefevre, Haine and Matignon, "Passivity-Based Boundary
Control for the 2D Shallow Water Equations with Viscous Dissipation".

This is a reproducibility experiment, NOT yet a validated SCRIMP example.
Run the smoke tests first and compare the power-balance defect against
the published results before using it as a release benchmark.
"""

import argparse
import json
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from petsc4py import PETSc

import scrimp as S

comm = PETSc.COMM_WORLD
rank = comm.getRank()

RHO = 1.0
GRAVITY = 0.1
VISCOSITY = 1e-3
H_INITIAL = 5.0
H_TARGET = 4.0
ROTATION = 0.2
REGION_DOMAIN = 1
REGION_WALL = 10


def ramp_expression(start=75.0, end=200.0):
    """C1 smooth decreasing boundary velocity, 1 until 75 s, 0 after 200 s.

    The paper shows a smooth decay in this interval but does not provide
    its explicit analytic expression; this is an approximation, NOT
    necessarily the exact function used in the publication.
    """
    if start < 0 or end <= start:
        raise ValueError("ramp start must be nonnegative and smaller than ramp end")
    s = f"(min(1., max(0., (t-{start})/{end-start})))"
    return f"(1. - 3.*({s})*({s}) + 2.*({s})*({s})*({s}))"


def _set_diagnostic_initial_guess(swe, name, expression):
    """Set an ALGEBRAIC GetFEM variable's starting value, not the PDE state.

    SCRIMP.set_initial_value() intentionally accepts only dynamic
    states. For this diagnostic, use the existing port FEM to evaluate
    a given field, restrict boundary DOFs exactly as set_initial_value
    does, and verify the GetFEM variable interval. This is not a
    generally DAE-consistent initialization or a production option.
    """
    candidates = [
        port for port in swe.ports.values()
        if name in (port.get_flow(), port.get_effort())
    ]
    if len(candidates) != 1:
        raise RuntimeError(
            f"Expected one SCRIMP port for algebraic variable {name}"
        )
    port = candidates[0]
    fem = port.get_fem()
    values = np.asarray(fem.eval(expression, globals(), locals())).reshape(-1)
    if port.get_region() is not None:
        dofs = fem.basic_dof_on_region(port.get_region())
        values = values[dofs]
    start, count = swe.gf_model.interval_of_variable(name)
    if values.size != int(count):
        raise RuntimeError(
            f"Initial guess for {name}: got {values.size} coefficients; "
            f"GetFEM expects {count}"
        )
    if not np.all(np.isfinite(values)):
        raise RuntimeError(f"Nonfinite initial guess for {name}")
    swe.gf_model.set_variable(name, values)
    actual = np.asarray(swe.gf_model.from_variables())[
        int(start):int(start + count)
    ]
    if not np.allclose(actual, values, rtol=1e-12, atol=1e-14):
        raise RuntimeError(f"GetFEM did not retain the {name} initial guess")
    if rank == 0:
        print(
            f"[SCRIMP algebraic guess] {name}: "
            f"norm={np.linalg.norm(actual):.9g}, "
            f"maxabs={np.max(np.abs(actual)):.9g}",
            flush=True,
        )


def build_system(gamma, mesh, final_time, save_dt, max_step,
                 ramp_start=75.0, ramp_end=200.0, init_step=True,
                 ts_type="arkimex", adapt_type="basic",
                 init_dt=1e-9, rotation=ROTATION,
                 initial_profile="tangent",
                 wall_velocity_profile="facet_tangent",
                 algebraic_initial_guess="zero",
                 start_from_rest_duration=0.0,
                 init_method="bdf", init_pseudo_steps=1,
                 height_order=2, geometry_order=1):
    if type(geometry_order) is not int or geometry_order not in (1, 2):
        raise ValueError("geometry_order must be integer 1 or 2")
    if height_order not in (2, 3):
        raise ValueError("height_order must be 2 or 3")
    if gamma not in (0.0, 0.25):
        raise ValueError("gamma must be 0.0 (open) or 0.25 (closed)")
    if not np.isfinite(rotation):
        raise ValueError("rotation must be finite")
    if initial_profile not in ("tangent", "rigid", "legacy_opposite"):
        raise ValueError("unknown initial momentum profile")
    if wall_velocity_profile not in ("facet_tangent", "rigid_trace"):
        raise ValueError("unknown wall velocity profile")
    if not np.isfinite(start_from_rest_duration) or start_from_rest_duration < 0:
        raise ValueError("start_from_rest_duration must be finite and nonnegative")
    if start_from_rest_duration > 0 and (
        gamma != 0.0 or initial_profile != "rigid"
        or wall_velocity_profile != "facet_tangent"
        or algebraic_initial_guess != "zero"
    ):
        raise ValueError(
            "start_from_rest_duration is diagnostic-only and requires "
            "the open case, rigid initial profile, physical facet_tangent "
            "wall and zero algebraic initial guesses"
        )
    if algebraic_initial_guess not in ("zero", "ep", "ep_wall"):
        raise ValueError("unknown algebraic initial guess")
    if init_method not in ("bdf", "pseudo"):
        raise ValueError("unknown initialization method")
    if init_pseudo_steps < 1:
        raise ValueError("init_pseudo_steps must be >= 1")
    if algebraic_initial_guess != "zero" and initial_profile != "rigid":
        raise ValueError("algebraic initial guess requires rigid profile")
    if algebraic_initial_guess == "ep_wall" and (
        gamma != 0.0 or wall_velocity_profile != "rigid_trace"
    ):
        raise ValueError(
            "ep_wall initial guess requires open case and rigid_trace "
            "diagnostic wall velocity to make both traces identical"
        )
    if wall_velocity_profile == "rigid_trace" and (
        gamma != 0.0 or initial_profile != "rigid"
    ):
        raise ValueError(
            "rigid_trace is a diagnostic restricted to the open case "
            "with rigid initial momentum"
        )
    if (mesh <= 0 or final_time <= 0 or save_dt <= 0
            or max_step <= 0 or not np.isfinite(init_dt) or init_dt <= 0):
        raise ValueError("mesh and time parameters must be positive")

    swe = S.DPHS("real")
    swe.gf_model.add_macro("div(v)", "Trace(Grad(v))")
    swe.gf_model.add_macro("D(v)", "Sym(Grad(v))")
    # IMPORTANT: GetFEM GWFL nested matrices list COLUMNS, not ROWS.
    # [[0,-1],[1,0]] means the actual matrix [[0,1],[-1,0]],
    # clockwise Rot*n=(n_y,-n_x) for an outward boundary normal.
    # The prior [[0,1],[-1,0]] acted counterclockwise and
    # imposed the OPPOSITE boundary velocity to the prescribed
    # clockwise initial momentum p0, as proven natively in #28.
    swe.gf_model.add_macro("Rot", "[[0,-1],[1,0]]")
    swe.gf_model.add_macro("Tangent", "(Rot*Normal)")
    swe.gf_model.add_macro("Curl2D(v)", "div(Rot*v)")
    swe.gf_model.add_macro("Gyro(v)", "Curl2D(v)*Rot")
    swe.set_domain(S.Domain(
        "Disk", {"R": 1.0, "h": mesh, "geometry_order": geometry_order}
    ))

    h = S.State("h", "Fluid height", "scalar-field")
    p = S.State("p", "Linear momentum", "vector-field")
    swe.add_state(h)
    swe.add_state(p)
    swe.add_costate(S.CoState("e_h", "Pressure", h))
    swe.add_costate(S.CoState("e_p", "Velocity", p))

    control_normal = S.Control_Port(
        "Normal boundary velocity", "U_n", "Normal velocity",
        "Y_n", "Normal pressure", "scalar-field",
        region=REGION_WALL, position="effort",
    )
    control_vector = S.Control_Port(
        "Wall velocity", "U", "Boundary velocity",
        "Y", "Viscous boundary traction", "vector-field",
        region=REGION_WALL, position="effort",
    )
    swe.add_control_port(control_normal)
    swe.add_control_port(control_vector)
    # The published PFEM experiments use k=2 for p/boundary fields
    # and k+1=3 for h/e_h. The original sandbox used P2 throughout.
    # Keep that historical baseline as the default; select 3 for
    # fidelity with the published mixed-order discretization.
    swe.add_FEM(S.FEM("h", height_order, FEM="CG"))
    for name in ("p", control_normal.get_name(), control_vector.get_name()):
        swe.add_FEM(S.FEM(name, 2, FEM="CG"))
    for name, value in (("rho", RHO), ("g", GRAVITY), ("mu", VISCOSITY)):
        swe.add_parameter(S.Parameter(name, name, "scalar-field", str(value), "h"))

    bricks = [
        # Left-hand-side mass matrices (including the h-dependent one).
        S.Brick("M_h", "h * Test_h", [REGION_DOMAIN], dt=True, position="flow"),
        S.Brick("M_p[h]", "h * p . Test_p", [REGION_DOMAIN],
                dt=True, linear=False, position="flow"),
        # Rotational interconnection and nonlinear constitutive relations.
        S.Brick("G[h,p]", "h * (Gyro(p) * e_p) . Test_p", [REGION_DOMAIN],
                linear=False, position="effort"),
        S.Brick("-M_e_h", "-e_h * Test_e_h", [REGION_DOMAIN]),
        S.Brick("Q_h", "rho*g*h*Test_e_h", [REGION_DOMAIN]),
        S.Brick("K_h", "0.5*(p.p)/rho*Test_e_h", [REGION_DOMAIN],
                linear=False),
        S.Brick("-M_e_p[h]", "-h*e_p.Test_e_p", [REGION_DOMAIN],
                linear=False),
        S.Brick("Q_p[h]", "h*(p/rho).Test_e_p", [REGION_DOMAIN],
                linear=False),
        # PFEM grad formulation; normal and viscous boundary ports.
        S.Brick("M_n[h]", "h*Y_n*Test_Y_n", [REGION_WALL],
                linear=False, position="flow"),
        S.Brick("D[h]", "h*e_p.Grad(Test_h)", [REGION_DOMAIN],
                linear=False, position="effort"),
        S.Brick("-B_n[h]", "-h*U_n*Test_h", [REGION_WALL],
                linear=False, position="effort"),
        S.Brick("-D[h]^T", "-h*Grad(e_h).Test_p", [REGION_DOMAIN],
                linear=False, position="effort"),
        S.Brick("-B_n[h]^T", "-h*e_h*Test_Y_n", [REGION_WALL],
                linear=False, position="effort"),
        # Viscous contributions of the symmetric gradient and divergence.
        S.Brick("-2mu R_Grad[h]", "-2*mu*h*D(e_p):D(Test_p)",
                [REGION_DOMAIN], linear=False, position="effort"),
        S.Brick("-2mu R_div[h]", "-2*mu*h*div(e_p)*div(Test_p)",
                [REGION_DOMAIN], linear=False, position="effort"),
        S.Brick("B_mu[h]", "h*Y.Test_p", [REGION_WALL],
                linear=False, position="effort"),
        S.Brick("M_b[h]", "h*U.Test_Y", [REGION_WALL],
                linear=False, position="flow"),
        S.Brick("B_mu[h]^T", "h*e_p.Test_Y", [REGION_WALL],
                linear=False, position="effort"),
    ]
    for brick in bricks:
        swe.add_brick(brick)

    f = ramp_expression(ramp_start, ramp_end)
    if start_from_rest_duration > 0:
        # DIAGNOSTIC ONLY. Instead of initializing a nonzero rotating
        # dynamic state and an already moving wall, start from the
        # validated stationary state. Then smoothly drive the *physical*
        # facet-tangential wall to the requested rotation. The decreasing
        # long-time ramp is retained by multiplying both envelopes.
        # At t=0 the rising envelope and its derivative both vanish.
        # This changes the initial-value/forcing experiment but NOT the
        # SCRIMP formulations, KSP/MUMPS options or BDF2 consistency gate.
        tau = start_from_rest_duration
        q = f"(min(1., max(0., t/{tau})))"
        rising = f"(3.*({q})*({q}) - 2.*({q})*({q})*({q}))"
        f = f"({f})*({rising})"
    if gamma == 0.0:
        swe.set_control(control_normal.get_name(), "0.")
    else:
        # u_n = -(1-f(t))*gamma*(y^0 + y^mu.n - y*),
        # where y^0 = Y_n, y^mu = Y and y* = -rho*g*h*.
        # An implicit algebraic brick is used: feedback is state-dependent,
        # so it must NOT be implemented as an external source term.
        normal_law = (
            f"(-U_n - ({gamma})*(1.-({f}))"
            f"*(Y_n + Y.Normal + {RHO*GRAVITY*H_TARGET}))*Test_U_n"
        )
        swe.add_brick(S.Brick("Normal output feedback", normal_law,
                              [REGION_WALL], linear=False))

    # Full boundary velocity: U = U_n*n + rotation*f(t)*Rot*n.
    # The ramp is DECREASING: f(0)=1 for the default start=75 s.
    # The intended clockwise Rot is encoded column-first in GetFEM:
    # Rot=[[0,-1],[1,0]]. Match the initial p sign,
    # but beware: on a faceted disk, the imposed Rot*Normal is not
    # exactly the trace of rigid-body p/rho=rotation*(y,-x).
    # Curved P2 geometry reduces the remaining normal-trace mismatch.
    # Even after correcting the Rot sign, an initial DAE consistency
    # solve must still check feasibility; do not infer from geometry.
    #
    # DIAGNOSTIC-ONLY 'rigid_trace': match the initial rigid p/rho
    # exactly on the polygonal mesh, using omega*Rot*X instead of
    # omega*Rot*Normal at t=0 (the decreasing ramp f(0)=1).
    # This is deliberately NOT a physically admissible no-penetration
    # condition on a faceted circle: Rot*X has a normal component.
    # It isolates the proposed geometric incompatibility without
    # changing GetFEM assembly, PETSc, or the BDF2 initializer.
    wall_vector = (
        "Tangent" if wall_velocity_profile == "facet_tangent"
        else "(Rot*X)"
    )
    vector_law = (
        f"(-U + U_n*Normal + ({rotation})*({f})*{wall_vector}).Test_U"
    )
    swe.add_brick(S.Brick("Tangential wall actuation", vector_law,
                          [REGION_WALL], linear=False))

    swe.set_initial_value("h", str(H_INITIAL))
    # Tangent=Rot*Normal, Rot=[[0,-1],[1,0]] in GWFL column syntax.
    # On the unit circle
    # Tangent=(y,-x). With e_p=p/rho the boundary condition needs
    # p/rho to point in this SAME direction; the old (-y,+x)
    # initialization was opposite to the imposed wall velocity.
    # Even with this corrected sign, the prescribed tangent based on
    # polygonal face normals is not exactly a smooth rigid-rotation trace.
    radius = "np.sqrt(x*x+y*y)"
    eps = 1e-6
    amplitude = RHO * (0.0 if start_from_rest_duration > 0 else rotation)
    if initial_profile == "tangent":
        # Approximately constant tangential speed; regularized
        # center, and approximate normals on the polygonal boundary.
        p_x = f"({amplitude})*y/({radius}+{eps})"
        p_y = f"-({amplitude})*x/({radius}+{eps})"
    elif initial_profile == "rigid":
        # Smooth rigid-body rotation, speed |rotation| at r=1,
        # zero at the center. This matches an ideal circular boundary,
        # not the finite-element disk's facetwise Rot*Normal exactly.
        p_x = f"({amplitude})*y"
        p_y = f"-({amplitude})*x"
    else:  # legacy_opposite preserves the earlier, opposite p0 sign
        p_x = f"-({amplitude})*y/({radius}+{eps})"
        p_y = f"({amplitude})*x/({radius}+{eps})"
    swe.set_initial_value("p", f"[{p_x}, {p_y}]")

    if algebraic_initial_guess != "zero":
        # The dynamic state p is already initialized but the algebraic
        # e_p and U default to zero. Test whether warm-starting these
        # variables removes the 1/init_dt startup impulse in Y.
        # This is a targeted initial-GUESS experiment, not a computed
        # consistent DAE initial condition.
        initial_velocity = f"[({rotation})*y, -({rotation})*x]"
        _set_diagnostic_initial_guess(swe, "e_p", initial_velocity)
        if algebraic_initial_guess == "ep_wall":
            _set_diagnostic_initial_guess(swe, "U", initial_velocity)

    # Keep the paper-inspired ARKIMEX-A2 setting as the baseline.
    # Distinguish truncation-error rejections from a truly failed
    # nonlinear solve with TSAdapt 'none', i.e. diagnostic fixed steps.
    # A fixed-step run is NOT a physical/accuracy validation.
    time_scheme = dict(
        ts_type=ts_type, ts_adapt_type=adapt_type,
        ksp_type="preonly", pc_type="lu", pc_factor_mat_solver_type="mumps",
        t_0=0.0, t_f=final_time, dt=max_step,
        dt_save=save_dt, ts_adapt_dt_max=max_step,
        ts_adapt_dt_min=1e-8,
        # State-dependent mass: use one tiny, fixed BDF2 step
        # from -1e-9 to t=0 to initialize the algebraic variables.
        # PETSc BDF2 bootstraps with an internal half-step.
        # This replaces the historical pseudo-TS *for this model*.
        init_step=init_step,
        init_step_nb_iter=init_pseudo_steps if init_method == "pseudo" else 1,
        init_step_ts_type=init_method, init_step_bdf_order=2,
        init_step_dt=init_dt,
        init_step_from_negative_time=(init_method == "bdf"),
        # This guard is now checked for both BDF2 and pseudo by
        # DPHS.init_step() when explicitly configured, preventing
        # a misleading "success" from earlier ~1e24 pseudo-TS blowups.
        init_step_max_norm_amplification=100.0,
    )
    if ts_type == "arkimex":
        time_scheme["ts_arkimex_type"] = "a2"
    elif ts_type == "bdf":
        time_scheme["ts_bdf_order"] = 2
    swe.set_time_scheme(**time_scheme)
    return swe


def integrate_midpoint(t, quantity):
    """Piecewise integral of midpoint-state power samples, as in MCSS_2024.

    With get_quantity(CN=True), sample k > 0 uses 0.5*(z[k-1]+z[k]),
    so multiplying sample k by dt is a midpoint-rule approximation.
    """
    energy = np.zeros_like(t, dtype=float)
    energy[1:] = np.cumsum(np.diff(t) * np.asarray(quantity)[1:])
    return energy


def diagnostics(swe, gamma):
    swe.hamiltonian.set_name("Mechanical Hamiltonian")
    swe.hamiltonian.add_term(
        S.Term("Kinetic energy", "0.5*h*(p.p)/rho", [REGION_DOMAIN])
    )
    swe.hamiltonian.add_term(
        S.Term("Potential energy", "0.5*rho*g*h*h", [REGION_DOMAIN])
    )
    swe.compute_Hamiltonian()  # GetFEM assembly: called on ALL MPI ranks.
    t = np.asarray(swe.solution["t"], dtype=float)
    kinetic = np.asarray(swe.hamiltonian.get_terms()[0].get_values())
    potential = np.asarray(swe.hamiltonian.get_terms()[1].get_values())
    H = kinetic + potential

    # Separate boundary inflows and nonnegative volume dissipation;
    # see equation (13) of the accompanying article.
    d_grad = np.asarray(
        swe.get_quantity("2*mu*h*D(e_p):D(e_p)", CN=True)
    )
    d_div = np.asarray(
        swe.get_quantity("2*mu*h*div(e_p)*div(e_p)", CN=True)
    )
    p_normal = np.asarray(
        swe.get_quantity("h*U_n*Y_n", CN=True, region=REGION_WALL)
    )
    p_viscous = np.asarray(
        swe.get_quantity("h*U.Y", CN=True, region=REGION_WALL)
    )
    dissipated = integrate_midpoint(t, d_grad + d_div)
    supplied = integrate_midpoint(t, p_normal + p_viscous)
    # dH/dt = P_normal + P_viscous - D_grad - D_div.
    balance = H - H[0] + dissipated - supplied
    relative_balance = balance / max(abs(float(H[0])), 1e-14)
    # A tiny error relative to the large background potential energy
    # can hide a large error relative to the actual energy exchange.
    exchanged_scale = max(
        float(np.max(np.abs(H - H[0]))),
        float(np.max(np.abs(supplied))),
        float(np.max(np.abs(dissipated))),
        1e-14,
    )
    relative_balance_exchange = balance / exchanged_scale
    # An independent conservation law for the impermeable open case:
    # d/dt int_Omega h = -int_Gamma h U_n.
    volume = np.asarray(swe.get_quantity("h"), dtype=float)
    normal_volume_outflow = np.asarray(
        swe.get_quantity("h*U_n", CN=True, region=REGION_WALL),
        dtype=float,
    )
    cumulative_outflow = integrate_midpoint(t, normal_volume_outflow)
    volume_balance = volume - volume[0] + cumulative_outflow
    relative_volume_balance = volume_balance / max(
        abs(float(volume[0])), 1e-14
    )
    lyapunov = np.asarray(
        swe.get_quantity(
            f"0.5*h*(p.p)/rho + 0.5*rho*g*(h-{H_TARGET})*(h-{H_TARGET})"
        )
    )
    return locals()


def save_outputs(swe, result, out):
    if rank != 0:
        return
    out.mkdir(parents=True, exist_ok=True)
    d = result
    t = d["t"]
    np.savez_compressed(
        out / "energies.npz", t=t, kinetic=d["kinetic"],
        potential=d["potential"], hamiltonian=d["H"],
        diss_grad=d["d_grad"], diss_div=d["d_div"],
        power_normal=d["p_normal"], power_viscous=d["p_viscous"],
        supplied=d["supplied"], dissipated=d["dissipated"],
        balance=d["balance"], relative_balance=d["relative_balance"],
        relative_balance_exchange=d["relative_balance_exchange"],
        volume=d["volume"], normal_volume_outflow=d["normal_volume_outflow"],
        cumulative_outflow=d["cumulative_outflow"],
        volume_balance=d["volume_balance"],
        relative_volume_balance=d["relative_volume_balance"],
        lyapunov=d["lyapunov"],
    )
    # Positivity of the height FE *coefficients* is a necessary,
    # not sufficient, condition for pointwise positive P2 height.
    # Preserve this weaker meaning explicitly in the metric name.
    h_start, h_count = swe.gf_model.interval_of_variable("h")
    h_start, h_count = int(h_start), int(h_count)
    height_coefficients = [
        np.asarray(z, dtype=float)[h_start:h_start + h_count]
        for z in swe.solution["z"]
    ]
    if not height_coefficients or any(v.size != h_count for v in height_coefficients):
        raise RuntimeError("Cannot verify stored Saint-Venant height coefficients")
    if not all(np.all(np.isfinite(v)) for v in height_coefficients):
        raise RuntimeError("Nonfinite height coefficients in saved trajectory")
    np.savez_compressed(
        out / "state_endpoints.npz",
        initial=np.asarray(swe.solution["z"][0], dtype=float),
        final=np.asarray(swe.solution["z"][-1], dtype=float),
    )
    metrics = {
        "minimum_height_coefficient": min(float(np.min(v)) for v in height_coefficients),
        "maximum_height_coefficient": max(float(np.max(v)) for v in height_coefficients),
        "petsc_ts_converged_reason": int(swe.last_ts_converged_reason),
        "petsc_ts_accepted_steps": int(swe.last_ts_steps),
        "petsc_ts_rejected_steps": int(swe.last_ts_rejections),
        "petsc_ts_snes_failures": int(swe.last_ts_snes_failures),
        "gamma": float(d["gamma"]),
        "mpi_ranks": comm.getSize(),
        "steps_saved": int(t.size),
        "final_time_reached": float(t[-1]),
        "initial_hamiltonian": float(d["H"][0]),
        "final_hamiltonian": float(d["H"][-1]),
        "maximum_absolute_relative_power_balance_error": float(
            np.max(np.abs(d["relative_balance"]))
        ),
        "final_absolute_relative_power_balance_error": float(
            abs(d["relative_balance"][-1])
        ),
        "maximum_absolute_balance_error_over_actual_exchange": float(
            np.max(np.abs(d["relative_balance_exchange"]))
        ),
        "maximum_absolute_volume_balance_error": float(
            np.max(np.abs(d["volume_balance"]))
        ),
        "maximum_absolute_relative_volume_balance_error": float(
            np.max(np.abs(d["relative_volume_balance"]))
        ),
        "initial_volume": float(d["volume"][0]),
        "final_volume": float(d["volume"][-1]),
        "final_hamiltonian_change": float(d["H"][-1]-d["H"][0]),
        "integrated_boundary_work": float(d["supplied"][-1]),
        "integrated_viscous_dissipation": float(d["dissipated"][-1]),
    }
    (out / "metrics.json").write_text(json.dumps(metrics, indent=2) + "\n")
    fig, axs = plt.subplots(3, 1, figsize=(10, 9), sharex=True)
    axs[0].plot(t, d["H"], label="Hamiltonian")
    axs[0].plot(t, d["kinetic"], label="Kinetic")
    axs[0].plot(t, d["potential"], label="Potential")
    axs[0].legend()
    axs[0].set_ylabel("Energy")
    axs[1].plot(t, d["supplied"], label="Boundary supply")
    axs[1].plot(t, d["dissipated"], label="Viscous dissipation")
    axs[1].legend()
    axs[1].set_ylabel("Integrated power")
    axs[2].plot(t, d["relative_balance"], label="Relative balance defect")
    axs[2].axhline(0, lw=0.6, color="black")
    axs[2].legend()
    axs[2].set_ylabel("(H - H0 + D - S)/|H0|")
    axs[2].set_xlabel("Time (s)")
    for ax in axs:
        ax.grid(True)
    fig.tight_layout()
    fig.savefig(out / "power_balance.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(9, 4))
    ax.plot(t, d["lyapunov"])
    ax.set(xlabel="Time (s)", ylabel="Lyapunov functional")
    ax.grid(True)
    fig.tight_layout()
    fig.savefig(out / "lyapunov.png", dpi=180)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", choices=("open", "closed"), required=True)
    parser.add_argument("--mesh", type=float, default=0.05)
    parser.add_argument(
        "--height-order", type=int, choices=(2, 3), default=2,
        help="CG polynomial order for h and e_h. The existing P2/P2 "
             "sandbox baseline is the default; select 3 to reproduce "
             "the paper's P3(h,e_h)/P2(p,e_p,boundary) choice.",
    )
    parser.add_argument(
        "--geometry-order", type=int, choices=(1, 2), default=1,
        help="Gmsh geometric order (1=straight-sided historical disk, "
             "2=experimental curved quadratic circle), independent of "
             "P2/P3 field FEM. Must verify GetFEM import and constraints.",
    )
    parser.add_argument("--final-time", type=float, default=250.0)
    parser.add_argument("--save-dt", type=float, default=1.0)
    parser.add_argument("--max-step", type=float, default=0.01)
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument("--no-init-step", action="store_true",
                        help="Diagnostic only: bypass BDF2 consistency initialization")
    parser.add_argument(
        "--rotation", type=float, default=ROTATION,
        help="Initial swirl and imposed tangential wall speed; "
             "zero gives a stationary-water diagnostic",
    )
    parser.add_argument(
        "--initial-profile",
        choices=("tangent", "rigid", "legacy_opposite"),
        default="tangent",
        help="Initial momentum profile: tangent approximates circular "
             "wall velocity (default); rigid matches an ideal circle but "
             "not polygonal facet normals exactly; legacy_opposite "
             "reproduces the historical sign mistake",
    )
    parser.add_argument(
        "--wall-velocity-profile",
        choices=("facet_tangent", "rigid_trace"),
        default="facet_tangent",
        help="Boundary tangential actuator: facet_tangent=Rot*Normal "
             "(physical default); rigid_trace=Rot*X is a diagnostic "
             "matching rigid initial p/rho at t=0, but has nonzero "
             "normal velocity on a polygonal boundary and MUST NOT "
             "be treated as a physically valid no-penetration model",
    )
    parser.add_argument(
        "--algebraic-initial-guess",
        choices=("zero", "ep", "ep_wall"),
        default="zero",
        help="DIAGNOSTIC ONLY: zero keeps historical algebraic initial "
             "guesses; ep initializes e_p to p/rho; ep_wall additionally "
             "initializes U to the matching rigid trace. ep_wall "
             "requires open case, rigid profile and rigid_trace wall. "
             "No production DAE initializer is changed.",
    )
    parser.add_argument(
        "--init-dt", type=float, default=1e-9,
        help="Tiny BDF2 initialization step, from -init-dt to t=0; "
             "default 1e-9",
    )
    parser.add_argument(
        "--newton-block-trace", action="store_true",
        help="Log per-GetFEM-variable Newton state norms and fluid "
             "height range to identify unstable initialization",
    )
    parser.add_argument(
        "--ts-type", choices=("arkimex", "cn", "bdf"),
        default="arkimex",
        help="PETSc time integrator; default ARKIMEX-A2",
    )
    parser.add_argument(
        "--adapt-type", choices=("basic", "none"),
        default="basic",
        help="PETSc TSAdapt algorithm; 'none' tests fixed-step behavior",
    )
    parser.add_argument(
        "--petsc-monitor", action="store_true",
        help="Enable PETSc TSAdapt/SNES/KSP convergence logging",
    )
    parser.add_argument("--ramp-start", type=float, default=75.0)
    parser.add_argument("--ramp-end", type=float, default=200.0)
    parser.add_argument(
        "--init-method", choices=("bdf", "pseudo"), default="bdf",
        help="DIAGNOSTIC initialization integrator, bdf (default, "
             "one step from t=-init_dt to zero) or PETSc pseudo (steady-"
             "state pseudo time, NOT a DAE-consistency projection). "
             "The physical simulation time integrator is unchanged.",
    )
    parser.add_argument(
        "--init-pseudo-steps", type=int, default=1,
        help="Maximum PETSc pseudo initialization steps (default 1). "
             "Only used with --init-method pseudo.",
    )
    parser.add_argument(
        "--start-from-rest-duration", type=float, default=0.0,
        help="DIAGNOSTIC ONLY: initialize rigid momentum at zero and "
             "multiply the physical facet-tangential wall rotation by "
             "a C1 smoothstep rising from 0 at t=0 to 1 at the stated "
             "positive duration. The original long-time decaying ramp "
             "is retained. Requires --case open --initial-profile rigid "
             "--wall-velocity-profile facet_tangent and default zero "
             "algebraic guesses. Disabled by default.",
    )
    args = parser.parse_args()
    if args.init_pseudo_steps < 1:
        parser.error("--init-pseudo-steps must be >= 1")
    if args.start_from_rest_duration < 0 or not np.isfinite(
        args.start_from_rest_duration
    ):
        parser.error("--start-from-rest-duration must be finite and >= 0")
    if args.start_from_rest_duration > 0 and (
        args.case != "open" or args.initial_profile != "rigid"
        or args.wall_velocity_profile != "facet_tangent"
        or args.algebraic_initial_guess != "zero"
    ):
        parser.error(
            "--start-from-rest-duration requires --case open "
            "--initial-profile rigid --wall-velocity-profile "
            "facet_tangent --algebraic-initial-guess zero"
        )
    if args.algebraic_initial_guess != "zero" and (
        args.initial_profile != "rigid"
    ):
        parser.error("nonzero algebraic guess requires --initial-profile rigid")
    if args.algebraic_initial_guess == "ep_wall" and (
        args.case != "open" or args.wall_velocity_profile != "rigid_trace"
    ):
        parser.error(
            "ep_wall guess requires --case open --wall-velocity-profile "
            "rigid_trace"
        )
    if args.wall_velocity_profile == "rigid_trace" and (
        args.case != "open" or args.initial_profile != "rigid"
    ):
        parser.error(
            "--wall-velocity-profile rigid_trace is diagnostic-only "
            "and requires --case open --initial-profile rigid"
        )
    gamma = 0.0 if args.case == "open" else 0.25
    out = args.output_dir or (
        Path(__file__).resolve().parent / "outputs" / "viscous_swe_2d"
        / args.case
    )
    if rank == 0:
        out.mkdir(parents=True, exist_ok=True)
    comm.barrier()

    if args.newton_block_trace:
        os.environ["SCRIMP_NEWTON_BLOCK_TRACE"] = "1"
    if rank == 0:
        print(
            "Saint-Venant diagnostic configuration: "
            f"TS={args.ts_type}, TSAdapt={args.adapt_type}, "
            f"pseudo_init={not args.no_init_step}, "
            f"dt={args.max_step}, t_f={args.final_time}, "
            f"height_order={args.height_order}, "
            f"geometry_order={args.geometry_order}, "
            f"init_dt={args.init_dt}, init_method={args.init_method}, "
            f"init_pseudo_steps={args.init_pseudo_steps}, "
            f"rotation={args.rotation}, "
            f"initial_profile={args.initial_profile}, "
            f"wall_velocity_profile={args.wall_velocity_profile}, "
            f"algebraic_initial_guess={args.algebraic_initial_guess}, "
            f"start_from_rest_duration={args.start_from_rest_duration}, "
            f"ramp_start={args.ramp_start}, ramp_end={args.ramp_end}",
            flush=True,
        )
    swe = build_system(gamma, args.mesh, args.final_time,
                       args.save_dt, args.max_step,
                       args.ramp_start, args.ramp_end,
                       init_step=not args.no_init_step,
                       ts_type=args.ts_type,
                       adapt_type=args.adapt_type,
                       init_dt=args.init_dt,
                       rotation=args.rotation,
                       initial_profile=args.initial_profile,
                       wall_velocity_profile=args.wall_velocity_profile,
                       algebraic_initial_guess=args.algebraic_initial_guess,
                       start_from_rest_duration=args.start_from_rest_duration,
                       init_method=args.init_method,
                       init_pseudo_steps=args.init_pseudo_steps,
                       height_order=args.height_order,
                       geometry_order=args.geometry_order)
    if args.petsc_monitor:
        # PETSc's -ts_monitor/-ts_adapt_monitor are VIEWER-valued
        # command-line flags, not Boolean-valued TS settings.
        # set_time_scheme(ts_monitor=True) stores the string "true":
        # TSSetFromOptions() interprets that as a viewer name and fails
        # with PETSC_ERR_SUP (56): "Unsupported viewer true".
        #
        # Store the flags with *no* argument as PETSc expects.
        # They apply to both the optional pseudo-TS and the main TS.
        options = PETSc.Options()
        for key in (
            "ts_monitor", "ts_adapt_monitor",
            "snes_converged_reason", "ksp_converged_reason",
        ):
            options.setValue(key, None)
    swe.solve()
    # SCRIMP currently allows its TS event handler to stop on minimum
    # timestep without raising.  Exit code 0 must not be interpreted as a
    # successful PDE simulation if t_f was never reached (observed: only
    # 1e-8 of the requested 0.5 seconds, with huge power-balance defect).
    # A monitor snapshot is NOT the authoritative TS completion time.
    # The published MCSS model reached 0.0003 while its last saved
    # monitor point was 0.0001. Check the actual PETSc convergence
    # reason and end time before trusting any power-balance output.
    reached = float(getattr(swe, "last_ts_time", float("-inf")))
    reason = int(getattr(swe, "last_ts_converged_reason", 0))
    if reason <= 0 or not np.isfinite(reached) or (
        reached < args.final_time - 1e-10 * max(1.0, abs(args.final_time))
    ):
        raise RuntimeError(
            "Saint-Venant integration did not converge to requested "
            f"t_f={args.final_time:.9g}: PETSc t={reached:.9g}, "
            f"reason={reason}. "
            "Do not interpret energy figures as validated."
        )
    result = diagnostics(swe, gamma)
    save_outputs(swe, result, out)
    comm.barrier()

    # Exercise issue #6 code path on every MPI rank. Note: SCRIMP's generic
    # port-power plot does not include the h weights used in the physical
    # boundary fluxes above, so power_balance.png is the reference diagnostic.
    swe.plot_Hamiltonian(
        with_powers=True, save_figure=True,
        filename=f"viscous_swe_{args.case}_Hamiltonian.png",
    )
    # Issue #6 MPI-native evidence: independently record the generic
    # (UNWEIGHTED) SCRIMP port powers actually assembled on EACH rank.
    # Never infer rank-local power equality from rank-0 diagnostics.
    # Names are per rank, so there are no competing output writers.
    port_powers = {}
    for port in swe.ports.values():
        if port.get_algebraic() and not port.get_substituted():
            values = np.asarray(port.get_power(), dtype=float)
            if values.shape != (len(swe.solution["t"]),) or not np.all(np.isfinite(values)):
                raise RuntimeError(f"Invalid rank-{rank} power trace for {port.get_name()}")
            port_powers[port.get_name()] = values.tolist()
    (out / f"mpi_port_powers_rank{rank}.json").write_text(
        json.dumps({
            "rank": rank, "world_size": comm.getSize(),
            "times": [float(t) for t in swe.solution["t"]],
            "powers": port_powers,
            "note": "generic unweighted port power; NOT physical h-weighted SWE balance",
        }, indent=2) + "\n",
        encoding="utf-8",
    )
    comm.barrier()

    if rank == 0:
        print(json.dumps({
            "case": args.case,
            "output_dir": str(out),
            "maximum_relative_balance_error":
                float(np.max(np.abs(result["relative_balance"]))),
            "maximum_balance_error_over_actual_exchange":
                float(np.max(np.abs(
                    result["relative_balance_exchange"]
                ))),
            "maximum_relative_volume_balance_error":
                float(np.max(np.abs(
                    result["relative_volume_balance"]
                ))),
        }, indent=2))


if __name__ == "__main__":
    main()
