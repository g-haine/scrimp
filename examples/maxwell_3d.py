# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
"""Small genuine 3D Maxwell cavity using patched intrinsic GetFEM elements.

Normalized lossless Maxwell: eps=mu=1, no currents or port drive.
E (electric field) in Nedelec H(curl), B (magnetic induction) in RT0
H(div), on a genuine tetrahedral 3D ball.

The weak equations are
    (dE/dt, v) = (B, curl v),
    (dB/dt, w) = -(curl E, w).
The first equation has natural zero tangential magnetic boundary work;
this is a weak homogeneous closure, not a separately enforced PEC
tangential E condition. The discrete coupling is skew-adjoint, hence
the *exact semidiscrete and CN identity* H(t)=H(0), with
H=0.5 integral(|E|^2+|B|^2).

Run only inside the pinned SCRIMP GetFEM 5.5 MPI image:
    python -m examples.maxwell_3d --check --export --output /artifacts
This example is not an electromagnetic convergence or divergence
benchmark. It deliberately tests the Nedelec and RT GetFEM patches.
"""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import scrimp as S

# TUTORIAL M3D-01: explicit three-dimensional curl (GWFL has no native Curl)
CURL3 = (
    "[Grad(v)(3,2)-Grad(v)(2,3);"
    "Grad(v)(1,3)-Grad(v)(3,1);"
    "Grad(v)(2,1)-Grad(v)(1,2)]"
)


def maxwell_3d(*, mesh_size=0.7, final_time=0.2, time_step=0.02,
               factor_solver="mumps"):
    """Build and solve the two-field Maxwell system on a tetrahedral ball."""
    values = (mesh_size, final_time, time_step)
    if any(not np.isfinite(float(v)) or float(v) <= 0 for v in values):
        raise ValueError("Mesh and time parameters must be finite and positive")
    if abs(final_time/time_step-round(final_time/time_step)) > 1e-10:
        raise ValueError("Choose an integer number of fixed CN steps")
    if factor_solver not in ("mumps", "superlu", "superlu_dist"):
        raise ValueError("Unknown PETSc LU backend")

    model = S.DPHS("real")
    model.gf_model.add_macro("Curl3(v)", CURL3)

    # TUTORIAL M3D-02: real 3D tetrahedral mesh, provided by SCRIMP
    model.set_domain(S.Domain("Ball", {"R": 1., "h": float(mesh_size)}))
    if model.domain.get_dim()[0] != 3:
        raise RuntimeError("Maxwell requires a genuine 3D volume mesh")

    # TUTORIAL M3D-03: coenergy = state (eps=mu=1)
    E = S.State("E", "Electric field", "vector-field")
    B = S.State("B", "Magnetic induction", "vector-field")
    model.add_state(E)
    model.add_state(B)
    model.add_costate(S.CoState("e_E", "Electric field", E, substituted=True))
    model.add_costate(S.CoState("e_B", "Magnetic induction", B, substituted=True))

    # TUTORIAL M3D-04: true 3D intrinsic vector elements, NOT CG vectors
    # S.FEM uses a raw GetFEM descriptor when order<0.
    # Both these GetFEM families use the SCRIMP-v1.3.0 patched image.
    model.add_FEM(S.FEM("E", -1, "FEM_NEDELEC(3)"))
    model.add_FEM(S.FEM("B", 0, "RT"))  # FEM_RTK(3,0)

    # S.FEM has a legacy silent CG fallback. Fail CLOSED if that happens.
    for name, requested in (("E", "FEM_NEDELEC(3)"),
                            ("B", "FEM_RTK(3,0)")):
        import getfem as gf
        actual = {f.char() for f in model.ports[name].get_fem().fem()}
        reference = gf.Fem(requested).char()
        if actual != {reference}:
            raise RuntimeError(
                f"Required patched {requested} not active for {name}: {actual}"
            )
        fem = model.ports[name].get_fem()
        if fem.qdim() != 3 or fem.nbdof() == 0:
            raise RuntimeError(f"Invalid intrinsic 3D field for {name}")

    # TUTORIAL M3D-05: two positive mass blocks, two skew curl blocks
    for brick in (
        S.Brick("M_E", "E.Test_E", [1], dt=True, position="flow"),
        S.Brick("M_B", "B.Test_B", [1], dt=True, position="flow"),
        S.Brick("C^T", "B.Curl3(Test_E)", [1], position="effort"),
        S.Brick("-C", "-Curl3(E).Test_B", [1], position="effort"),
    ):
        model.add_brick(brick)

    # TUTORIAL M3D-06: smooth solenoidal rotation (exact Nedelec mode),
    # zero magnetic field; curl E=[0,0,2] triggers nontrivial dynamics.
    model.set_initial_value("E", "[-y, x, 0.]")
    model.set_initial_value("B", "[0., 0., 0.]")

    model.set_time_scheme(
        # Both states carry a positive mass: this is an ODE, NOT a DAE.
        # The algebraic initializer requires algebraic unknowns;
        # a pseudo-TS or artificial DAE spin-up would be wrong here.
        init_step=False,
        ts_type="cn", ts_adapt_type="none",
        ts_exact_final_time="matchstep",
        t_0=0., t_f=float(final_time),
        dt=float(time_step), dt_save=float(time_step),
        ksp_type="gmres", pc_type="lu",
        pc_factor_mat_solver_type=factor_solver,
    )
    model.solve()

    # TUTORIAL M3D-07: true GetFEM 3D electromagnetic Hamiltonian
    model.hamiltonian.set_name("3D electromagnetic energy")
    model.hamiltonian.add_term(S.Term("Electric", "0.5*E.E", [1]))
    model.hamiltonian.add_term(S.Term("Magnetic", "0.5*B.B", [1]))
    return model


def measure_maxwell(system):
    """Audit real FE energy, nontrivial B, and discrete Gauss diagnostic."""
    import getfem as gf

    t = np.asarray(system.solution["t"], dtype=float)
    z = np.asarray(system.solution["z"], dtype=float)
    H = np.asarray(system.get_Hamiltonian(), dtype=float)
    if t.ndim != 1 or len(t) < 2 or z.ndim != 2 or len(z) != len(t):
        raise RuntimeError("Invalid original FE Maxwell time/state arrays")
    if H.shape != t.shape or not np.all(np.diff(t) > 0):
        raise RuntimeError("Invalid Maxwell energy sampling")
    if not all(np.all(np.isfinite(x)) for x in (t, z, H)):
        raise RuntimeError("Nonfinite Maxwell FE data")
    mim = system.domain._int_method[0]
    electric, magnetic, gauss_B = [], [], []
    for tk, zk in zip(t, z):
        system.gf_model.set_time(float(tk))
        system.gf_model.to_variables(zk)
        electric.append(float(gf.asm(
            "generic", mim, 0, "0.5*E.E", 1, system.gf_model)))
        magnetic.append(float(gf.asm(
            "generic", mim, 0, "0.5*B.B", 1, system.gf_model)))
        gauss_B.append(float(gf.asm(
            "generic", mim, 0,
            "Trace(Grad(B))*Trace(Grad(B))", 1, system.gf_model)))
    He = np.asarray(electric)
    Hb = np.asarray(magnetic)
    div_B_squared = np.asarray(gauss_B)
    steps = np.diff(H)
    defect = H-H[0]
    scale = max(1.0, float(H[0]))
    gates = {
        "real_3d_mesh": system.domain.get_dim()[0] == 3,
        "authentic_fixed_CN_steps": bool(
            system.last_ts_converged_reason > 0
            and system.last_ts_steps == len(t)-1
            and system.last_ts_rejections == 0
            and system.last_ts_snes_failures == 0),
        "pure_differential_Maxwell_no_artificial_DAE": bool(
            len(system.states) == 2 and
            all(not p.get_algebraic() for p in system.ports.values())),
        "positive_initial_em_energy": bool(
            H[0] > 1e-7 and min(np.min(He), np.min(Hb)) >= -1e-10),
        "electric_and_magnetic_exchange": bool(
            Hb[-1] > 1e-6*H[0] and He[-1] < He[0]),
        "independent_fe_energy_matches": bool(
            np.max(np.abs(He+Hb-H)) < 1e-9*scale),
        "every_CN_step_conserves_energy": bool(
            np.max(np.abs(steps)) < 1e-9*scale),
        "cumulative_em_energy_conserved": bool(
            np.max(np.abs(defect)) < 1e-9*scale),
        # Gauss(B): curl E belongs to RT0 on the same tetrahedral mesh;
        # B(0)=0, hence the compatible discrete div B must stay zero.
        # This guards the *second* Maxwell constraint independently
        # of the skew/energy check, not only as a plotted diagnostic.
        "magnetic_gauss_law": bool(
            np.min(div_B_squared) >= -1e-12*scale and
            np.max(div_B_squared) <= 1e-8*scale),
    }
    return dict(t=t, z=z, H=H, electric_energy=He, magnetic_energy=Hb,
                div_B_squared=div_B_squared, step_energy_defect=steps,
                cumulative_energy_defect=defect, gates=gates)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mesh-size", type=float, default=0.7)
    parser.add_argument("--time-step", type=float, default=0.02)
    parser.add_argument("--final-time", type=float, default=0.2)
    parser.add_argument("--solver", choices=("mumps", "superlu", "superlu_dist"),
                        default="mumps")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--export", action="store_true")
    args = parser.parse_args()
    if (args.check or args.export) and args.output is None:
        parser.error("--check/--export require --output")
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    if comm.size > 1 and args.solver == "superlu":
        parser.error("Sequential SuperLU cannot run in MPI2")

    system = maxwell_3d(mesh_size=args.mesh_size,
                        final_time=args.final_time,
                        time_step=args.time_step, factor_solver=args.solver)
    if args.output is None:
        return system
    args.output.mkdir(parents=True, exist_ok=True)
    obs = measure_maxwell(system)
    np.savez_compressed(
        args.output/f"rank{comm.rank}.npz",
        **{k: v for k, v in obs.items() if isinstance(v, np.ndarray)}
    )
    metadata = {
        "status": ("PASS_NATIVE_MAXWELL3D" if all(obs["gates"].values())
                   else "FAIL_NATIVE_MAXWELL3D"),
        "rank": comm.rank, "MPI": comm.size,
        "solver": args.solver,
        "DOFs": int(obs["z"].shape[1]),
        "steps": int(system.last_ts_steps),
        "TS_reason": int(system.last_ts_converged_reason),
        "DAE_initialization": system.initialization_report,
        "physical_gates": obs["gates"],
        "initial_energy": float(obs["H"][0]),
        "final_energy": float(obs["H"][-1]),
        "final_magnetic_energy": float(obs["magnetic_energy"][-1]),
        "max_step_energy_defect": float(np.max(np.abs(
            obs["step_energy_defect"]))),
        "max_cumulative_energy_defect": float(np.max(np.abs(
            obs["cumulative_energy_defect"]))),
        "max_divB_squared": float(np.max(obs["div_B_squared"])),
        "state_SHA256": hashlib.sha256(
            np.ascontiguousarray(obs["z"]).tobytes()).hexdigest(),
    }
    (args.output/f"rank{comm.rank}.json").write_text(
        json.dumps(metadata, indent=2)+"\n", encoding="utf-8"
    )
    if args.check and not all(obs["gates"].values()):
        raise RuntimeError("Native Maxwell 3D physical acceptance FAILED: "
                           + json.dumps(metadata))
    if args.export:
        for name in ("E", "B"):
            system.export_to_pv(name, path=str(args.output), t="All")
    if comm.rank == 0:
        print(json.dumps(metadata, indent=2), flush=True)
    return system


if __name__ == "__main__":
    main()
