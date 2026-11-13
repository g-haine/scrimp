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


def build_system(gamma, mesh, final_time, save_dt, max_step,
                 ramp_start=75.0, ramp_end=200.0, init_step=True):
    if gamma not in (0.0, 0.25):
        raise ValueError("gamma must be 0.0 (open) or 0.25 (closed)")
    if mesh <= 0 or final_time <= 0 or save_dt <= 0 or max_step <= 0:
        raise ValueError("mesh and time parameters must be positive")

    swe = S.DPHS("real")
    swe.gf_model.add_macro("div(v)", "Trace(Grad(v))")
    swe.gf_model.add_macro("D(v)", "Sym(Grad(v))")
    swe.gf_model.add_macro("Rot", "[[0,1],[-1,0]]")
    swe.gf_model.add_macro("Tangent", "(Rot*Normal)")
    swe.gf_model.add_macro("Curl2D(v)", "div(Rot*v)")
    swe.gf_model.add_macro("Gyro(v)", "Curl2D(v)*Rot")
    swe.set_domain(S.Domain("Disk", {"R": 1.0, "h": mesh}))

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
    for name in ("h", "p", control_normal.get_name(), control_vector.get_name()):
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

    # Full boundary velocity: u = u_n n + (ROTATION*f(t)) t.
    # Keep U_n as an unknown so the feedback equation is solved implicitly.
    vector_law = f"(-U + U_n*Normal + ({ROTATION})*({f})*Tangent).Test_U"
    swe.add_brick(S.Brick("Tangential wall actuation", vector_law,
                          [REGION_WALL], linear=False))

    swe.set_initial_value("h", str(H_INITIAL))
    # Approximate the radial-constant tangential speed 0.2 from the paper:
    # regularize the singular direction at the disk center.
    radius = "np.sqrt(x*x+y*y)"
    eps = 1e-6
    swe.set_initial_value(
        "p",
        f"[-{RHO*ROTATION}*y/({radius}+{eps}),"
        f" {RHO*ROTATION}*x/({radius}+{eps})]",
    )

    # PETSc ARKIMEX A2 is the TR-BDF2-type second-order DIRK scheme
    # used in the published computation; MUMPS supports MPI factorization.
    swe.set_time_scheme(
        ts_type="arkimex", ts_arkimex_type="a2",
        ksp_type="preonly", pc_type="lu", pc_factor_mat_solver_type="mumps",
        t_0=0.0, t_f=final_time, dt=max_step,
        dt_save=save_dt, ts_adapt_dt_max=max_step,
        ts_adapt_dt_min=1e-8,
        init_step=init_step, init_step_nb_iter=2,
        init_step_ts_type="pseudo", init_step_dt=1e-9,
    )
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
        lyapunov=d["lyapunov"],
    )
    metrics = {
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
    parser.add_argument("--final-time", type=float, default=250.0)
    parser.add_argument("--save-dt", type=float, default=1.0)
    parser.add_argument("--max-step", type=float, default=0.01)
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument("--no-init-step", action="store_true",
                        help="Diagnostic only: bypass pseudo-TS initial consistency solve")
    parser.add_argument("--ramp-start", type=float, default=75.0)
    parser.add_argument("--ramp-end", type=float, default=200.0)
    args = parser.parse_args()
    gamma = 0.0 if args.case == "open" else 0.25
    out = args.output_dir or (
        Path(__file__).resolve().parent / "outputs" / "viscous_swe_2d"
        / args.case
    )
    if rank == 0:
        out.mkdir(parents=True, exist_ok=True)
    comm.barrier()

    swe = build_system(gamma, args.mesh, args.final_time,
                       args.save_dt, args.max_step,
                       args.ramp_start, args.ramp_end,
                       init_step=not args.no_init_step)
    swe.solve()
    # SCRIMP currently allows its TS event handler to stop on minimum
    # timestep without raising.  Exit code 0 must not be interpreted as a
    # successful PDE simulation if t_f was never reached (observed: only
    # 1e-8 of the requested 0.5 seconds, with huge power-balance defect).
    recorded_times = np.asarray(swe.solution["t"], dtype=float)
    reached = float(recorded_times[-1]) if recorded_times.size else float("-inf")
    if not np.isfinite(reached) or (
        reached < args.final_time - 1e-10 * max(1.0, abs(args.final_time))
    ):
        raise RuntimeError(
            f"Saint-Venant integration stopped at t={reached:.9g}, "
            f"before requested t_f={args.final_time:.9g}. "
            "Do not interpret any resulting energy figures as validated."
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
    if rank == 0:
        print(json.dumps({
            "case": args.case,
            "output_dir": str(out),
            "maximum_relative_balance_error":
                float(np.max(np.abs(result["relative_balance"]))),
        }, indent=2))


if __name__ == "__main__":
    main()
