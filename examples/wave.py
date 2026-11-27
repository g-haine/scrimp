# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp


"""Conservative 2D anisotropic membrane with mixed boundary ports.

2D mixed-boundary anisotropic membrane, energy variables q=grad w,
p=rho*w_t, coenergy e_q=T q and e_p=p/rho. This preserves the
original model's four boundary ports and its actual GetFEM GWFL
Hamiltonian 0.5*q.T.q, where T is the full tensor parameter.

Physical outward-normal traction: u_N=n.e_q at bottom/right/top;
observation y_N=e_p. Left side prescribes u_D=e_p (velocity) and
observes y_D=n.e_q (traction). The pH energy balance is

    H(t)-H(0) = integral (u_N*y_N + u_D*y_D) dt.

The negative flow-side algebraic rows impose Y_N=e_p and U_D=e_p:
this is the same physical sign convention as qualified wave_1d.py.
The *physical* sign of the measured work is audited independently
against the original GetFEM Hamiltonian; no success is inferred
solely from these formulae.

This tutorial starts from EXACTLY zero, compatible with zero
initial Dirichlet velocity and initial normal stress, and uses
a left-boundary velocity that is zero at t=0 but nontrivial soon
after. The original model, before promotion, passed real MUMPS MPI1/2
10-step TSCN physical and shared-ParaView tests at cf42089d;
other LU backends require independent native acceptance.

Scientific context: Haine, Matignon, Serhani (2023), anisotropic
heterogeneous wave PFEM, DOI 10.4208/ijnam2023-1005; and
Brugnoli, Haine, Matignon (2022) on mixed boundary control.
"""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import scrimp as S


# TUTORIAL 2A-00: boundary names and GetFEM region identifiers
BOUNDARIES = (
    ("bottom", 10, "B"),
    ("right", 11, "R"),
    ("top", 12, "T"),
    ("left", 13, "L"),
)


def wave_eq(
    *, mesh_size=0.2, final_time=0.08, time_step=0.008,
    factor_solver=None,
):
    """Build and integrate the conservative mixed-boundary 2D pH wave.

    The time/mesh defaults are intentionally for native installation
    acceptance, not for a convergence order or published benchmark.
    """
    mesh_size = float(mesh_size)
    final_time = float(final_time)
    time_step = float(time_step)
    if not (np.isfinite(mesh_size) and mesh_size > 0
            and np.isfinite(final_time) and final_time > 0
            and np.isfinite(time_step) and 0 < time_step <= final_time):
        raise ValueError("Finite positive mesh/time values required")
    ratio = final_time / time_step
    if abs(ratio-round(ratio)) > 1e-10:
        raise ValueError("Physical fixed CN requires final_time/dt integer")
    if factor_solver not in (None, "mumps", "superlu", "superlu_dist"):
        raise ValueError("Unknown PETSc LU backend")

    # TUTORIAL 2A-01: DPHS, domain, energy and co-energy variables
    system = S.DPHS("real")
    system.set_output_name("wave_2d")
    system.set_domain(S.Domain("Rectangle", {
        "L": 2., "l": 1., "h": mesh_size,
    }, mesh_output_dir=system.get_output_path("mesh",create=True)))

    q = S.State("q", "Strain", "vector-field")
    p = S.State("p", "Linear momentum", "scalar-field")
    system.add_state(q)
    system.add_state(p)
    system.add_costate(S.CoState("e_q", "Stress", q))
    system.add_costate(S.CoState("e_p", "Velocity", p))

    # TUTORIAL 2A-02: Neumann and velocity-controlled boundary ports
    for side, region, suffix in BOUNDARIES:
        # Only the LEFT boundary is velocity-controlled (flow port).
        is_left = side == "left"
        system.add_control_port(S.Control_Port(
            f"Boundary control ({side})",
            f"U_{suffix}",
            "Velocity" if is_left else "Outward normal traction",
            f"Y_{suffix}",
            "Outward normal traction" if is_left else "Velocity",
            "scalar-field",
            region=region,
            position="flow" if is_left else "effort",
        ))

    # TUTORIAL 2A-03: partitioned finite elements, anisotropic material
    for fem in (
        S.FEM("q", 1, "DG"),
        S.FEM("p", 2, "CG"),
        *(S.FEM(f"Boundary control ({side})", 1, "DG")
          for side, _, _ in BOUNDARIES),
    ):
        system.add_FEM(fem)

    # This T is symmetric positive definite throughout [0,2] x [0,1].
    # GWFL nested matrices list COLUMNS; it is symmetric here.
    system.add_parameter(S.Parameter(
        "T", "Anisotropic membrane tension", "tensor-field",
        "[[5+x,x*y],[x*y,2+y]]", "q",
    ))
    system.add_parameter(S.Parameter(
        "rho", "Areal mass density", "scalar-field",
        "3-x", "p",
    ))

    # TUTORIAL 2A-04: mixed Dirac and constitutive bricks
    bricks = [
        S.Brick("M_q", "q.Test_q", [1], dt=True, position="flow"),
        S.Brick("M_p", "p*Test_p", [1], dt=True, position="flow"),
        # Signed algebraic observation in *flow*:
        #   -Y_N = -e_p -> Y_N=e_p (Neumann outputs).
        S.Brick("-M_Y_B", "-Y_B*Test_Y_B", [10], position="flow"),
        S.Brick("-M_Y_R", "-Y_R*Test_Y_R", [11], position="flow"),
        S.Brick("-M_Y_T", "-Y_T*Test_Y_T", [12], position="flow"),
        # Mixed left velocity input: -U_L=-e_p -> U_L=e_p.
        S.Brick("-M_U_L", "-U_L*Test_Y_L", [13], position="flow"),
        S.Brick("D", "Grad(e_p).Test_q", [1], position="effort"),
        S.Brick("-D^T", "-e_q.Grad(Test_p)", [1], position="effort"),
        S.Brick("B_B", "U_B*Test_p", [10], position="effort"),
        S.Brick("B_R", "U_R*Test_p", [11], position="effort"),
        S.Brick("B_T", "U_T*Test_p", [12], position="effort"),
        S.Brick("B_L", "Y_L*Test_p", [13], position="effort"),
        S.Brick("-B_B^T", "-e_p*Test_Y_B", [10], position="effort"),
        S.Brick("-B_R^T", "-e_p*Test_Y_R", [11], position="effort"),
        S.Brick("-B_T^T", "-e_p*Test_Y_T", [12], position="effort"),
        S.Brick("-B_L^T", "-e_p*Test_Y_L", [13], position="effort"),
        S.Brick("-M_e_q", "-e_q.Test_e_q", [1]),
        S.Brick("CR_q", "q.T.Test_e_q", [1]),
        S.Brick("-M_e_p", "-e_p*Test_e_p", [1]),
        S.Brick("CR_p", "p/rho*Test_e_p", [1]),
    ]
    for brick in bricks:
        system.add_brick(brick)

    # TUTORIAL 2A-05: compatible physical initial and boundary inputs
    # Exactly compatible at t=0 (zero stress and zero velocity).
    # The left boundary injects energy immediately after t=0.
    for side, _, _ in BOUNDARIES:
        value = (
            "0.2*sin(4*pi*t)*sin(pi*y)"
            if side == "left" else "0."
        )
        system.set_control(f"Boundary control ({side})", value)
    system.set_initial_value("q", "[0., 0.]")
    system.set_initial_value("p", "0.")

    # TUTORIAL 2A-06: strict DAE initialization, fixed Crank-Nicolson
    options = {}
    if factor_solver is not None:
        options["pc_factor_mat_solver_type"] = factor_solver
    system.set_time_scheme(
        init_step_ts_type="algebraic", init_step=True,
        ts_type="cn", ts_adapt_type="none",
        ts_exact_final_time="matchstep",
        ksp_type="gmres", pc_type="lu",
        t_0=0., t_f=final_time,
        dt=time_step, dt_save=time_step,
        **options,
    )
    system.solve()

    # TUTORIAL 2A-07: exact tensor Hamiltonian from the original GetFEM field
    system.hamiltonian.set_name(
        "Anisotropic membrane mechanical energy (mixed boundary control)"
    )
    system.hamiltonian.add_term(S.Term(
        "Potential energy (q dot T q)", "0.5*q.T.q", [1],
    ))
    system.hamiltonian.add_term(S.Term(
        "Kinetic energy", "0.5*p*p/rho", [1],
    ))
    # TUTORIAL 2A-08: end of code shown progressively in Sphinx
    return system


def measure_wave_2d(system):
    """Measure ORIGINAL GetFEM energy and all four signed port powers."""
    t = np.asarray(system.solution["t"], dtype=float)
    z = np.asarray(system.solution["z"], dtype=float)
    H = np.asarray(system.get_Hamiltonian(), dtype=float)
    if not system.powers_computed:
        system.compute_powers()
    port_power = {
        side: np.asarray(
            system.ports[f"Boundary control ({side})"].get_power(),
            dtype=float,
        ) for side, _, _ in BOUNDARIES
    }
    sampling = {
        side: system.ports[f"Boundary control ({side})"].get_power_sampling()
        for side, _, _ in BOUNDARIES
    }
    if any(v != "cn_step" for v in sampling.values()):
        raise RuntimeError(f"Expected REAL step-averaged CN powers: {sampling}")
    if (t.ndim != 1 or len(t) < 2
            or z.ndim != 2 or z.shape[0] != len(t)
            or H.shape != t.shape
            or any(p.shape != t.shape for p in port_power.values())
            or not np.all(np.isfinite(t))
            or not np.all(np.diff(t) > 0)
            or not np.all(np.isfinite(z))
            or not np.all(np.isfinite(H))
            or any(not np.all(np.isfinite(p)) for p in port_power.values())):
        raise RuntimeError("Malformed native accepted FE/Hamiltonian/port history")

    total = sum(port_power.values())
    step_work = np.diff(t) * total[1:]
    step_deltaH = np.diff(H)
    defects = step_deltaH - step_work
    W = np.concatenate(([0.], np.cumsum(step_work)))
    init = system.initialization_report
    init_ok = (
        isinstance(init, dict)
        and init.get("status") == "PASS"
        and bool(init.get("differential_initial_exact"))
    )
    cv = (
        system.last_ts_converged_reason > 0
        and system.last_ts_steps == len(t) - 1
        and system.last_ts_rejections == 0
        and system.last_ts_snes_failures == 0
    )
    magnitude = max(1., float(np.max(np.abs(H))))
    gates = {
        "native_TS_fixed_CN_complete_accepted_samples": bool(cv),
        "original_differential_initial_state_preserved": init_ok,
        "strictly_positive_anisotropic_Hamiltonian_after_excitation":
            bool(np.min(H) >= -1e-10 and np.max(H) > 1e-9),
        "nontrivial_signed_left_boundary_work":
            bool(np.max(np.abs(port_power["left"][1:])) > 1e-8),
        "zero_power_from_three_unforced_Neumann_boundaries":
            all(float(np.max(np.abs(port_power[side]))) < 1e-10
                for side in ("bottom", "right", "top")),
        "each_CN_step_mechanical_energy_balance":
            bool(np.max(np.abs(defects)) <= 1e-8*magnitude),
        "cumulative_CN_mechanical_energy_balance":
            bool(np.max(np.abs(H-H[0]-W)) <= 1e-8*magnitude),
    }
    return dict(
        t=t, z=z, H=H,
        power_bottom=port_power["bottom"],
        power_right=port_power["right"],
        power_top=port_power["top"],
        power_left=port_power["left"],
        signed_work=W,
        energy_step_defects=defects,
        max_abs_step_defect=float(np.max(np.abs(defects))),
        max_abs_cumulative_defect=float(np.max(np.abs(H-H[0]-W))),
        max_abs_left_power=float(np.max(np.abs(port_power["left"]))),
        gates=gates,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--solver", default="mumps",
                        choices=("mumps", "superlu", "superlu_dist"))
    parser.add_argument("--mesh-size", type=float, default=.2)
    parser.add_argument("--final-time", type=float, default=.08)
    parser.add_argument("--time-step", type=float, default=.008)
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--export", action="store_true")
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args()

    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    if comm.size > 1 and args.solver == "superlu":
        parser.error("PETSc SuperLU is sequential; use SuperLU_DIST/MUMPS")
    if (args.export or args.check) and args.output is None:
        parser.error("--check/--export require --output for raw native evidence")

    system = wave_eq(
        mesh_size=args.mesh_size,
        final_time=args.final_time,
        time_step=args.time_step,
        factor_solver=args.solver,
    )
    try:
        observations = measure_wave_2d(system)
    except Exception as error:
        if args.output:
            args.output.mkdir(parents=True, exist_ok=True)
            np.savez_compressed(
                args.output / f"rank{comm.rank}.failed.npz",
                times=np.asarray(system.solution["t"], dtype=float),
                states=np.asarray(system.solution["z"], dtype=float),
            )
            (args.output / f"rank{comm.rank}.failed.json").write_text(
                json.dumps({
                    "status": "UNQUALIFIED_NATIVE_FE_TRAJECTORY",
                    "error": f"{type(error).__name__}: {error}",
                    "accepted_steps": system.last_ts_steps,
                }, indent=2)+"\n",
                encoding="utf-8",
            )
        raise

    if args.output is not None:
        args.output.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(
            args.output / f"rank{comm.rank}.npz",
            **{key:value for key,value in observations.items()
               if isinstance(value, np.ndarray)},
        )
        metadata = {
            "status": (
                "PASS_ORIGINAL_FE_PHYSICS"
                if all(observations["gates"].values())
                else "FAIL_ORIGINAL_FE_PHYSICS"
            ),
            "rank": comm.rank,
            "MPI": comm.size,
            "requested_LU": args.solver,
            "selected_LU": system.time_scheme.getString(
                "pc_factor_mat_solver_type", ""
            ).lower(),
            "DOFs": int(observations["z"].shape[1]),
            "steps": int(system.last_ts_steps),
            "TS_reason": int(system.last_ts_converged_reason),
            "DAE_initialization": system.initialization_report,
            "max_abs_step_defect": observations["max_abs_step_defect"],
            "max_abs_cumulative_defect":
                observations["max_abs_cumulative_defect"],
            "max_abs_left_power": observations["max_abs_left_power"],
            "physical_gates": observations["gates"],
            "state_SHA256": hashlib.sha256(
                np.ascontiguousarray(observations["z"]).tobytes()
            ).hexdigest(),
        }
        (args.output/f"rank{comm.rank}.json").write_text(
            json.dumps(metadata, indent=2, default=str)+"\n",
            encoding="utf-8",
        )
    if args.check and not all(observations["gates"].values()):
        if comm.rank == 0:
            print(json.dumps({
                "status": "FAIL_WAVE2D_MIXED_NATIVE_PHYSICS",
                "physical_gates": observations["gates"],
                "max_step_defect": observations["max_abs_step_defect"],
                "max_cumulative_defect":
                    observations["max_abs_cumulative_defect"],
            }, indent=2), flush=True)
        raise RuntimeError(
            "Original GetFEM tensor Hamiltonian/mixed boundary work "
            "failed scientific acceptance; no PASS assumed"
        )
    # Regular end-user invocation writes original GetFEM figures/data
    # under ./outputs/wave_2d/, never rank-named probe folders.
    # Explicit --output retains historical per-rank scientific
    # evidence semantics and existing runner paths unchanged.
    if args.output is None:
        system.plot_Hamiltonian(
            save_figure=True,show=False,
            filename="wave2d_hamiltonian.png"
        )
        for variable in ("q", "p", "e_q", "e_p"):
            system.export_to_pv(variable,t="All")
        if comm.rank==0:
            print("SCRIMP wave_2d user outputs:",system._output_root,flush=True)

    if args.export:
        # Required collective: root-only writer first validates rank
        # replication before producing shared native VTK/PVD.
        for name in ("q", "p", "e_q", "e_p"):
            system.export_to_pv(name, path=str(args.output), t="All")
    if comm.rank == 0:
        print(json.dumps({
            "status": "PASS_WAVE2D_MIXED_NATIVE_PHYSICS",
            "solver": args.solver,
            "MPI": comm.size,
            "steps": int(system.last_ts_steps),
            "DOFs": observations["z"].shape[1],
            "gates": observations["gates"],
            "max_abs_step_defect": observations["max_abs_step_defect"],
            "max_abs_cumulative_defect":
                observations["max_abs_cumulative_defect"],
        }, indent=2), flush=True)


if __name__ == "__main__":
    main()
