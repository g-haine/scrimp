# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Short 1-/2-rank native physics gate for examples/wave_coenergy.py.

The full coenergy DAE includes independently assembled q,p, a
substituted dissipative port e_r and boundary feedback U=0.5*Y.
MUMPS and production SCRIMP equations are unchanged. With constant
coefficients and fixed-step Crank-Nicolson, exact implicit residual is

    M (z[n+1]-z[n])/dt + K_full (z[n+1]+z[n])/2 = 0,

where K_full = system.tangent_stiffness includes the linear GetFEM
feedback tangent, unlike system.stiffness. Every equality is checked
against owned PETSc CSR rows; no new nonlinear residual assumption is
introduced into production code.

The energy balance uses the *midpoint* powers

    H = 0.5 z^T M z,
    P_r = e_r^T K[e_r,e_r] e_r,
    P_boundary = U^T (-K[U,U]) Y,

    H[n+1]-H[n] + dt*(P_r+P_boundary) == 0.

The full state is already gathered by SCRIMP's normal memory trajectory
on every rank. This tiny diagnostic validates against independent
GetFEM Hamiltonian terms from the published example, including their
original coefficient and subdomain definitions. 1-/2-rank trajectory
and power parity is checked by a second, offline script.

Not an adaptive-step test: the goal is FIRST to establish that the
canonical CN power balance and dissipative/feedback constraints hold
after the Concentric geometry fix.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path

from mpi4py import MPI
import numpy as np
from scipy.sparse import csr_matrix
import scrimp as S

from audit_ph_blocks import _build_coenergy_wave
from scrimp.utils.ph_block_layout import (
    classify_ph_formulation, discover_layout,
)


def _local_operator(mat, n):
    lo, hi = map(int, mat.getOwnershipRange())
    ptr, cols, values = mat.getValuesCSR()
    if tuple(mat.getSize()) != (n, n):
        raise ValueError("Unexpected PETSc operator dimensions")
    return lo, hi, csr_matrix(
        (values, cols, ptr), shape=(hi-lo, n)
    )


def _global_norm(local_vec, mpi):
    x = np.asarray(local_vec)
    value = float(np.vdot(x, x).real)
    return float(np.sqrt(mpi.allreduce(value, op=MPI.SUM)))


def _quadratic(m, lo, hi, x, mpi):
    local = float(np.vdot(x[lo:hi], m.dot(x)).real)
    return float(mpi.allreduce(local, op=MPI.SUM))


def _midpoint_port_powers(k, lo, hi, mid, indices, mpi):
    """Distributed positive (dissipation + feedback) port powers."""
    e = indices["e_r"]
    u = indices["U"]
    y = indices["Y"]
    if (u.stop-u.start) != (y.stop-y.start):
        raise ValueError("U and Y must share the boundary FE size")
    eri, erj = max(lo,e.start), min(hi,e.stop)
    uri, urj = max(lo,u.start), min(hi,u.stop)
    local_r = 0.
    local_b = 0.
    if eri < erj:
        r = slice(eri-lo, erj-lo)
        local_r = float(np.vdot(
            mid[eri:erj],
            k[r, e].dot(mid[e]),
        ).real)
    if uri < urj:
        r = slice(uri-lo, urj-lo)
        local_b = -float(np.vdot(
            mid[uri:urj],
            k[r, u].dot(mid[y]),
        ).real)
    return (
        float(mpi.allreduce(local_r, op=MPI.SUM)),
        float(mpi.allreduce(local_b, op=MPI.SUM)),
    )


def verify(args):
    mpi = MPI.COMM_WORLD
    rank, ranks = mpi.Get_rank(), mpi.Get_size()
    if args.mesh <= 0 or args.dt <= 0 or args.steps < 3:
        raise ValueError("Need positive mesh/dt and >=3 CN steps")
    if os.environ.get("SCRIMP_GETFEM_MATRIX_MODE") != "distributed":
        raise ValueError("Require GetFEM distributed matrix mode")
    if os.environ.get("SCRIMP_GETFEM_BRIDGE_DISTRIBUTED") != "coo":
        raise ValueError("Require COO GetFEM-to-PETSc bridge")
    system = _build_coenergy_wave(args.mesh)
    system.set_time_scheme(
        ts_type="cn", ts_adapt_type="none",
        ts_exact_final_time="matchstep",
        dt=args.dt, dt_save=args.dt,
        t_f=args.steps * args.dt,
        ksp_type="preonly", pc_type="lu",
        pc_factor_mat_solver_type="mumps",
    )
    system.solve()
    if (
        system.last_ts_converged_reason <= 0
        or system.last_ts_steps != args.steps
        or system.last_ts_rejections
        or system.last_ts_snes_failures
    ):
        raise RuntimeError("Full MUMPS/CN coenergy time solve failed")
    times = np.asarray(system.solution["t"], dtype=float)
    states = np.stack(
        [np.asarray(z, dtype=float) for z in system.solution["z"]]
    )
    n = int(system.gf_model.nbdof())
    if times.shape != (args.steps+1,) or states.shape != (
        args.steps+1, n
    ):
        raise RuntimeError(
            "All accepted CN steps and their complete state vectors "
            "must be saved (including t=0)"
        )
    expected = np.arange(args.steps+1)*args.dt
    if not np.allclose(times, expected, rtol=0, atol=2e-12):
        raise RuntimeError("Unexpected accepted CN time grid")
    fingerprint = hashlib.sha256(states.tobytes()).hexdigest()
    all_digest = mpi.allgather(fingerprint)
    if len(set(all_digest)) != 1:
        raise RuntimeError("SCRIMP full saved states differ across ranks")

    fields = discover_layout(
        system.gf_model, states=system.states,
        costates=system.costates, ports=system.ports,
        controls=system.controls,
    )
    formulation = classify_ph_formulation(
        fields, costates=system.costates
    )
    if formulation["formulation"] != "substituted_coenergy":
        raise ValueError("The coenergy example is no longer substituted")
    indices = {
        f["name"]: slice(f["start"], f["stop"])
        for f in fields
    }
    if not {"q", "p", "e_r", "U", "Y"} <= set(indices):
        raise ValueError("A dynamic, damping or boundary DoF is missing")
    lm, hm, m = _local_operator(system.mass, n)
    lk, hk, k = _local_operator(system.tangent_stiffness, n)
    if (lm, hm) != (lk, hk):
        raise ValueError("Mass/K row ownership differs")
    lo, hi = lm, hm

    energy = np.array([
        0.5*_quadratic(m, lo, hi, state, mpi)
        for state in states
    ], dtype=float)
    if not np.all(np.isfinite(energy)) or np.min(energy) <= 0:
        raise ValueError("Invalid/nonpositive acoustic Hamiltonian")

    # These terms are defined AFTER solve() in the published example.
    # They are the independent GetFEM finite-element H reference,
    # not reconstructed from PETSc's mass matrix.
    system.hamiltonian.add_term(
        S.Term("Potential energy", "0.5*q.Tinv.q", [1, 2])
    )
    system.hamiltonian.add_term(
        S.Term("Kinetic energy", "0.5*p*p*rho", [1, 2])
    )
    system.compute_Hamiltonian()
    fe_energy = np.asarray(system.get_Hamiltonian(), dtype=float)
    fe_energy_rel = float(
        np.linalg.norm(fe_energy-energy)
        / max(np.linalg.norm(energy), 1e-300)
    )
    if fe_energy.shape != energy.shape or fe_energy_rel > 5e-9:
        raise RuntimeError(
            "GetFEM physical H does not equal 0.5 z^T M z: "
            f"{fe_energy_rel:.6g}"
        )

    power_d = []
    power_b = []
    cn_residuals = []
    algebraic_residuals = []
    defects = []
    for j in range(args.steps):
        dt = float(times[j+1] - times[j])
        mid = .5*(states[j]+states[j+1])
        dz = (states[j+1]-states[j])/dt
        mdot = m.dot(dz)
        kmid = k.dot(mid)
        residual = mdot + kmid
        resid_rel = (
            _global_norm(residual, mpi)
            / max(
                _global_norm(mdot, mpi)+_global_norm(kmid, mpi),
                1e-300,
            )
        )
        algebraic_rows = np.zeros(hi-lo, dtype=bool)
        for fld in ("U", "Y", "e_r"):
            interval = indices[fld]
            a, b = max(lo, interval.start), min(hi, interval.stop)
            if a < b:
                algebraic_rows[a-lo:b-lo] = True
        constraint_rel = (
            _global_norm(kmid[algebraic_rows], mpi)
            / max(_global_norm(kmid, mpi), 1e-300)
        )
        pr, pb = _midpoint_port_powers(
            k, lo, hi, mid, indices, mpi
        )
        balance = (energy[j+1]-energy[j]) + dt*(pr+pb)
        scale = max(
            abs(energy[j]), abs(energy[j+1]),
            abs(dt*(pr+pb)), 1e-300,
        )
        defects.append(float(abs(balance)/scale))
        power_d.append(pr)
        power_b.append(pb)
        cn_residuals.append(float(resid_rel))
        algebraic_residuals.append(float(constraint_rel))

    # Conservativity of the lossless Dirac terms and correct signs:
    # the only allowed positive powers here are nu damping in
    # the disk and the dissipative U=0.5Y boundary feedback.
    if min(power_d) < -1e-11 or min(power_b) < -1e-11:
        raise RuntimeError("Negative physical dissipated power")
    maxima = {
        "cn_relative_residual": max(cn_residuals),
        "algebraic_relative_residual": max(algebraic_residuals),
        "relative_discrete_power_defect": max(defects),
        "independent_GetFEM_H_relative_error": fe_energy_rel,
    }
    if max(maxima.values()) > 5e-8:
        raise RuntimeError(
            "Coenergy physics gate failure: " + str(maxima)
        )
    if rank == 0:
        args.output_dir.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(
            args.output_dir / f"coenergy_{ranks}r_physics.npz",
            times=times, states=states, energy=energy,
            energy_getfem=fe_energy,
            power_damping_mid=np.asarray(power_d),
            power_boundary_mid=np.asarray(power_b),
            cn_relative_residuals=np.asarray(cn_residuals),
            algebraic_relative_residuals=np.asarray(
                algebraic_residuals
            ),
            relative_power_defects=np.asarray(defects),
        )
        report = {
            "schema": "scrimp-coenergy-physics-v1",
            "source_commit": os.environ.get("SCRIMP_BENCH_GIT_SHA", ""),
            "mpi_ranks": ranks, "mesh": args.mesh,
            "dofs": n, "steps": args.steps,
            "t_final": float(times[-1]), "dt": args.dt,
            "complete_state_sha256": fingerprint,
            "full_state_snapshots": len(times),
            "TS_converged_reason": int(
                system.last_ts_converged_reason
            ),
            "effective_mumps_policy": getattr(
                system, "_mumps_reset_policy", None
            ),
            "formulation": formulation["formulation"],
            "power_sign_convention": (
                "P_damping=+er^T K_erer er; "
                "P_boundary=U^T (-K_UU)Y; "
                "H[n+1]-H[n]+dt*(P_damping+P_boundary)=0"
            ),
            "maxima": maxima,
        }
        (args.output_dir/f"coenergy_{ranks}r_physics.json").write_text(
            json.dumps(report, indent=2)+"\n", encoding="utf-8"
        )
        print(
            f"PASS coenergy physics {ranks}r: {n} DoFs, "
            f"{args.steps} CN steps, "
            f"max CN residual={maxima['cn_relative_residual']:.3g}, "
            "max discrete power defect="
            f"{maxima['relative_discrete_power_defect']:.3g}, "
            f"GetFEM H error={fe_energy_rel:.3g}",
            flush=True,
        )
    return maxima


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mesh", type=float, default=.12)
    parser.add_argument("--dt", type=float, default=.002)
    parser.add_argument("--steps", type=int, default=5)
    parser.add_argument("--output-dir", type=Path, required=True)
    verify(parser.parse_args())


if __name__ == "__main__":
    main()
