# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""The first SCRIMP tutorial: a vibrating string with Neumann controls.

The physical strain q=dw/dx and momentum p=rho*dw/dt satisfy
q_t=partial_x(e_p), p_t=partial_x(e_q), e_q=T*q, e_p=p/rho.
The Hamiltonian is H=integral_0^1 (T*q*q+p*p/rho)/2 dx.

Both boundary controls are OUTWARD tractions:
U_L=-e_q(0), U_R=e_q(1), and Y_L=e_p(0), Y_R=e_p(1).
The boundary power balance reads dH/dt = U_L*Y_L+U_R*Y_R.
The negative flow-side algebraic Y rows encode SCRIMP's
minus-observation convention in the mixed Dirac equations.

Scientific background (not an exact reproduction of paper figures):
Haine, Matignon & Serhani (2023), DOI 10.4208/ijnam2023-1005;
Ferraro, Fournie & Haine (2024), DOI 10.1016/j.ifacol.2024.08.267.
"""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import scrimp as S


def wave_1d(*, final_time=0.04, time_step=0.004, mesh_size=0.05,
            factor_solver=None, plot=False):
    """Build and solve the genuine PFEM port-Hamiltonian string DAE.

    The quick defaults are suitable for a first installation smoke,
    not for a spatial/temporal convergence study. Omitting the solver
    selects the existing MUMPS default. SuperLU is SERIAL ONLY.
    """
    final_time, time_step, mesh_size = (
        float(final_time), float(time_step), float(mesh_size)
    )
    if not 0 < time_step <= final_time or mesh_size <= 0:
        raise ValueError("Require 0 < dt <= final_time and mesh_size > 0")
    steps = final_time/time_step
    if abs(steps-round(steps)) > 1e-10:
        raise ValueError("The fixed-CN tutorial requires final_time/dt integer")
    if factor_solver not in (None, "mumps", "superlu", "superlu_dist"):
        raise ValueError("Unsupported PETSc direct LU solver")

    system = S.DPHS("real")
    system.set_domain(S.Domain("Interval", {"L": 1., "h": mesh_size}))

    q = S.State("q", "Strain", "scalar-field")
    p = S.State("p", "Linear momentum", "scalar-field")
    system.add_state(q)
    system.add_state(p)
    system.add_costate(S.CoState("e_q", "Stress", q))
    system.add_costate(S.CoState("e_p", "Velocity", p))

    for side, region, suffix in (("left", 10, "L"), ("right", 11, "R")):
        system.add_control_port(S.Control_Port(
            f"Boundary control ({side})", f"U_{suffix}",
            "Outward traction", f"Y_{suffix}", "Velocity",
            "scalar-field", region=region,
        ))

    # Mixed FE spaces on one spatial interval, not lumped scalar ports.
    for finite_element in (
        S.FEM("q", 2), S.FEM("p", 1),
        S.FEM("Boundary control (left)", 1),
        S.FEM("Boundary control (right)", 1),
    ):
        system.add_FEM(finite_element)

    system.add_parameter(S.Parameter(
        "T", "String tension", "scalar-field", "1", "q",
    ))
    system.add_parameter(S.Parameter(
        "rho", "Linear mass density", "scalar-field", "1+x*(1-x)", "p",
    ))

    bricks = [
        # Time-derivative mass matrix M.
        S.Brick("M_q", "q*Test_q", [1], dt=True, position="flow"),
        S.Brick("M_p", "p*Test_p", [1], dt=True, position="flow"),
        # Flow-side boundary entries are -y, not +y.
        S.Brick("-M_Y_L", "-Y_L*Test_Y_L", [10], position="flow"),
        S.Brick("-M_Y_R", "-Y_R*Test_Y_R", [11], position="flow"),
        # Antisymmetric Dirac interconnection and outward traction.
        S.Brick("D", "Grad(e_p)*Test_q", [1], position="effort"),
        S.Brick("-D^T", "-e_q*Grad(Test_p)", [1], position="effort"),
        S.Brick("B_L", "U_L*Test_p", [10], position="effort"),
        S.Brick("B_R", "U_R*Test_p", [11], position="effort"),
        S.Brick("-B_L^T", "-e_p*Test_Y_L", [10], position="effort"),
        S.Brick("-B_R^T", "-e_p*Test_Y_R", [11], position="effort"),
        # Linear constitutive laws, implicit zero-residual form.
        S.Brick("-M_e_q", "-e_q*Test_e_q", [1]),
        S.Brick("CR_q", "q*T*Test_e_q", [1]),
        S.Brick("-M_e_p", "-e_p*Test_e_p", [1]),
        S.Brick("CR_p", "p/rho*Test_e_p", [1]),
    ]
    for brick in bricks:
        system.add_brick(brick)

    # Compatible t=0 data: the tapered Gaussian strain vanishes
    # at each boundary, where the two imposed tractions are zero.
    system.set_control("Boundary control (left)", "sin(2*pi*t)")
    system.set_control("Boundary control (right)", "0.")
    system.set_initial_value(
        "q", "2*x*(1-x)*np.exp(-50*(x-0.5)*(x-0.5))",
    )
    system.set_initial_value("p", "0.")

    options = {}
    if factor_solver is not None:
        options["pc_factor_mat_solver_type"] = factor_solver
    system.set_time_scheme(
        init_step_ts_type="algebraic", init_step=True,
        ts_type="cn", ts_adapt_type="none",
        ts_exact_final_time="matchstep",
        ksp_type="gmres", pc_type="lu",
        t_0=0., t_f=final_time, dt=time_step, dt_save=time_step,
        **options,
    )
    system.solve()

    system.hamiltonian.set_name("Mechanical energy of a vibrating string")
    for term in (
        S.Term("Potential energy", "0.5*q*T*q", [1]),
        S.Term("Kinetic energy", "0.5*p*p/rho", [1]),
    ):
        system.hamiltonian.add_term(term)

    if plot:
        # Called collectively in MPI; only rank 0 plots.
        system.plot_Hamiltonian(
            with_powers=True, save_figure=True,
            filename="wave_1d_energy_and_power.png",
        )
    return system


def observe_wave_1d(system):
    """Measure independent genuine GetFEM Hamiltonian and port powers.

    CN step powers P[k>=1] are evaluated on mean FE states. They
    must be integrated as sum(dt[k]*P[k]); P[0] is instantaneous,
    and must NOT be included in a trapezoidal stage-power integral.
    """
    times = np.asarray(system.solution["t"], dtype=float)
    states = np.asarray(system.solution["z"], dtype=float)
    H = np.asarray(system.get_Hamiltonian(), dtype=float)
    if not system.powers_computed:
        system.compute_powers()
    sides = ("Boundary control (left)", "Boundary control (right)")
    powers = {
        side: np.asarray(system.ports[side].get_power(), dtype=float)
        for side in sides
    }
    if any(system.ports[side].get_power_sampling() != "cn_step"
           for side in sides):
        raise RuntimeError("This physical check requires actual CN step powers")
    if (times.ndim != 1 or len(times) < 2
            or states.ndim != 2 or states.shape[0] != len(times)
            or H.shape != times.shape
            or any(v.shape != times.shape for v in powers.values())
            or not np.all(np.isfinite(times))
            or not np.all(np.diff(times) > 0)
            or not np.all(np.isfinite(states))
            or not np.all(np.isfinite(H))
            or any(not np.all(np.isfinite(v)) for v in powers.values())):
        raise RuntimeError("Nonfinite or malformed accepted FE/port history")

    P = sum(powers.values())
    step_work = np.diff(times)*P[1:]
    delta_H = np.diff(H)
    defect = delta_H-step_work
    cumulative_work = np.concatenate(([0.], np.cumsum(step_work)))
    max_defect = float(np.max(np.abs(defect)))
    max_boundary_power = float(np.max(np.abs(P[1:])))
    if (system.last_ts_converged_reason <= 0
            or system.last_ts_steps != len(times)-1
            or system.last_ts_rejections != 0
            or system.last_ts_snes_failures != 0
            or not np.all(H >= -1e-12)):
        raise RuntimeError("Real PETSc CN convergence or energy positivity failed")
    scale = max(1., float(np.max(np.abs(H))))
    init = system.initialization_report
    gates = {
        "PETSc_accepted_fixed_CN": True,
        "strict_DAE_init_without_changing_dynamic_data":
            isinstance(init, dict)
            and init.get("status") == "PASS"
            and bool(init.get("differential_initial_exact")),
        "nonzero_boundary_work": max_boundary_power > 1e-9,
        "stepwise_signed_energy_balance": max_defect <= 1e-8*scale,
        "cumulative_signed_energy_balance":
            float(np.max(np.abs(H-H[0]-cumulative_work))) <= 1e-8*scale,
    }
    return {
        "times": times, "states": states, "H": H,
        "power_left": powers[sides[0]], "power_right": powers[sides[1]],
        "signed_work": cumulative_work, "energy_defect": defect,
        "max_abs_step_energy_defect": max_defect,
        "max_abs_boundary_power": max_boundary_power,
        "gates": gates,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--solver", choices=("mumps", "superlu", "superlu_dist"),
        default="mumps", help="PETSc direct LU backend",
    )
    parser.add_argument("--mesh-size", type=float, default=0.05)
    parser.add_argument("--final-time", type=float, default=0.04)
    parser.add_argument("--time-step", type=float, default=0.004)
    parser.add_argument("--output", type=Path, default=None,
                        help="Optional persistent result directory")
    parser.add_argument("--plot", action="store_true")
    parser.add_argument("--check", action="store_true",
                        help="Fail on genuine DAE/CN/energy defects")
    args = parser.parse_args()

    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    if args.solver == "superlu" and comm.size != 1:
        parser.error("SuperLU is sequential; use SuperLU_DIST or MUMPS with MPI")
    system = wave_1d(
        final_time=args.final_time, time_step=args.time_step,
        mesh_size=args.mesh_size, factor_solver=args.solver,
        plot=args.plot,
    )
    report = observe_wave_1d(system)
    if args.output is not None:
        args.output.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(
            args.output/f"rank{comm.rank}.npz",
            **{k:v for k,v in report.items() if isinstance(v, np.ndarray)},
        )
        metadata = {
            "schema": "scrimp-v130-pedagogical-wave-1d-native-v1",
            "rank": comm.rank,
            "MPI": comm.size,
            "requested_solver": args.solver,
            "selected_PETSc_option": system.time_scheme.getString(
                "pc_factor_mat_solver_type", ""
            ).lower(),
            "full_state_DOFs": int(report["states"].shape[1]),
            "accepted_steps": int(system.last_ts_steps),
            "initialization_report": system.initialization_report,
            "physical_gates": report["gates"],
            "max_abs_step_energy_defect":
                report["max_abs_step_energy_defect"],
            "max_abs_boundary_power": report["max_abs_boundary_power"],
            "state_SHA256": hashlib.sha256(
                np.ascontiguousarray(report["states"]).tobytes()
            ).hexdigest(),
        }
        (args.output/f"rank{comm.rank}.json").write_text(
            json.dumps(metadata, indent=2)+"\n", encoding="utf-8",
        )
    passed = all(report["gates"].values())
    if comm.rank == 0:
        print(json.dumps({
            "status": (
                "PASS_WAVE_1D_NATIVE_ENERGY_AND_PORT_WORK"
                if passed else "FAIL_WAVE_1D_NATIVE_PHYSICS"
            ),
            "solver": args.solver, "mpi_ranks": comm.size,
            "accepted_steps": int(system.last_ts_steps),
            "max_energy_defect": report["max_abs_step_energy_defect"],
            "max_boundary_power": report["max_abs_boundary_power"],
            "physical_gates": report["gates"],
        }, indent=2), flush=True)
    if args.check and not passed:
        raise RuntimeError("Genuine original PFEM energy/DAE gate failed")


if __name__ == "__main__":
    main()
