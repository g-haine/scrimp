# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Issue #28.0: native algebraic-structure audit of original viscous SWE.

This intercepts production DPHS initialization AFTER real GetFEM/PETSc
mass, stiffness, nonlinear tangent and RHS assembly, BEFORE any
initialization or physical TS step. It holds prescribed h0 and p0 fixed.
The experiment is not a nonlinear DAE initializer or an index proof.
"""

import argparse
import json
from pathlib import Path
from unittest.mock import patch

import numpy as np
from mpi4py import MPI
from petsc4py import PETSc
import scrimp as S
from viscous_swe_2d_mpi import build_system


class AuditComplete(Exception):
    """Expected interruption before TS time stepping."""


def layout(s):
    names = [x.get_name() for x in s.states.values()]
    dynamic = tuple(names)
    for p in s.ports.values():
        if p.get_substituted():
            raise ValueError("Substituted port not supported")
        choices = (p.get_flow(), p.get_effort()) if p.get_algebraic() else (p.get_effort(),)
        names.extend(x for x in choices if x)
    if len(names) != len(set(names)) or set(dynamic) != {"h", "p"} or "Y" not in names:
        raise ValueError("Unexpected SWE unknown partition")
    intervals = {n: tuple(map(int, s.gf_model.interval_of_variable(n))) for n in names}
    n = int(s.gf_model.nbdof())
    cover = np.zeros(n, dtype=np.int8)
    for start, width in intervals.values():
        if start < 0 or width < 1 or start + width > n:
            raise ValueError("Invalid variable interval")
        cover[start:start + width] += 1
    if not np.all(cover == 1):
        raise ValueError("Variables fail to cover GetFEM DOFs exactly")
    return dynamic, tuple(x for x in names if x not in dynamic), intervals, n


def indices(names, intervals):
    return np.concatenate([
        np.arange(intervals[n][0], intervals[n][0] + intervals[n][1], dtype=np.int64)
        for n in names
    ])


def indexset(rows, lo, hi, comm):
    subset = rows[(rows >= lo) & (rows < hi)]
    return PETSc.IS().createGeneral(np.asarray(subset, dtype=PETSc.IntType), comm=comm)


def matmax(mat, mpi):
    _, _, values = mat.getValuesCSR()
    maximum = float(np.max(np.abs(values))) if len(values) else 0.0
    return float(mpi.allreduce(maximum, op=MPI.MAX))


def rowmax(mat, rows, mpi):
    ptr, _, vals = mat.getValuesCSR()
    lo, hi = mat.getOwnershipRange()
    maximum = 0.0
    for row in rows[(rows >= lo) & (rows < hi)]:
        values = vals[ptr[row - lo]:ptr[row - lo + 1]]
        if len(values):
            maximum = max(maximum, float(np.max(np.abs(values))))
    return float(mpi.allreduce(maximum, op=MPI.MAX))


def norm_sub(v, idx):
    sub = v.getSubVector(idx)
    try:
        return float(sub.norm())
    finally:
        v.restoreSubVector(idx, sub)


def inspect(s, mpi):
    differential, algebraic, intervals, n = layout(s)
    d = indices(differential, intervals)
    a = indices(algebraic, intervals)
    y = indices(("Y",), intervals)
    original = np.asarray(s.gf_model.from_variables(), dtype=float).copy()
    lo, hi = s.tangent_stiffness.getOwnershipRange()
    if s.tangent_mass.getOwnershipRange() != (lo, hi):
        raise RuntimeError("PETSc tangent mass/stiffness partition mismatch")
    owned_y = y[(y >= lo) & (y < hi)]
    owned = np.arange(lo, hi, dtype=PETSc.IntType)
    made = []
    def keep(obj):
        made.append(obj)
        return obj
    try:
        z = keep(s.tangent_stiffness.createVecRight())
        if hi > lo:
            z.setValues(owned, original[lo:hi])
        z.assemble()
        zd = keep(z.duplicate())
        zd.set(0)
        f = keep(z.duplicate())
        isd = keep(indexset(d, lo, hi, z.getComm()))
        isa = keep(indexset(a, lo, hi, z.getComm()))
        isy = keep(indexset(y, lo, hi, z.getComm()))
        s.IFunction(None, 0.0, z, zd, f)
        fa_initial = norm_sub(f, isa)
        blocks = {}
        for name in algebraic:
            sub = keep(indexset(indices((name,), intervals), lo, hi, z.getComm()))
            blocks[name] = norm_sub(f, sub)
        s._sync_lumped_trial_state(0.0, z)
        s.assemble_nl_stiffness()
        k = s.tangent_stiffness
        kaa = keep(k.createSubMatrix(isa, isa))
        kay = keep(k.createSubMatrix(isa, isy))
        kdy = keep(k.createSubMatrix(isd, isy))
        cutoff = 1e-12 * max(1.0, matmax(kaa, mpi))
        mass_algebraic_max = rowmax(s.tangent_mass, a, mpi)
        kay_max = matmax(kay, mpi)
        kdy_max = matmax(kdy, mpi)
        # Independently perturb Y only and evaluate the actual IFunction.
        v = keep(z.duplicate())
        v.set(0)
        if len(owned_y):
            v.setValues(np.asarray(owned_y, dtype=PETSc.IntType),
                        np.ones(len(owned_y)))
        v.assemble()
        if v.norm() <= 0:
            raise RuntimeError("No traction Y DOFs available")
        v.scale(1.0 / v.norm())
        zp = keep(z.copy())
        zm = keep(z.copy())
        zp.axpy(1e-6, v)
        zm.axpy(-1e-6, v)
        fp = keep(z.duplicate())
        fm = keep(z.duplicate())
        s.IFunction(None, 0.0, zp, zd, fp)
        s.IFunction(None, 0.0, zm, zd, fm)
        fp.axpy(-1, fm)
        fp.scale(5e5)
        jv = keep(z.duplicate())
        k.mult(v, jv)
        fd_a, fd_d = norm_sub(fp, isa), norm_sub(fp, isd)
        j_a, j_d = norm_sub(jv, isa), norm_sub(jv, isd)
        fp.axpy(-1, jv)
        err_a, err_d = norm_sub(fp, isa), norm_sub(fp, isd)
        s._sync_lumped_trial_state(0.0, z)
        restored = np.asarray(s.gf_model.from_variables(), dtype=float)
        unchanged = all(np.array_equal(
            original[intervals[name][0]:sum(intervals[name])],
            restored[intervals[name][0]:sum(intervals[name])]
        ) for name in differential)
        if not mpi.allreduce(unchanged, op=MPI.LAND):
            raise RuntimeError("Differential initial state changed in audit")
        if not all(np.isfinite(q) for q in (
            mass_algebraic_max, kay_max, kdy_max, fd_a, fd_d, err_a, err_d)):
            raise FloatingPointError("Nonfinite native algebraic evidence")
        if mass_algebraic_max > cutoff:
            raise RuntimeError("Nonzero algebraic rows in nonlinear tangent mass")
        if kay_max > cutoff or kdy_max <= cutoff:
            raise RuntimeError("Expected open-loop Y algebraic null column/dynamic coupling absent")
        if err_a > 1e-5 * max(1.0, fd_a, j_a) or err_d > 1e-5 * max(1.0, fd_d, j_d):
            raise RuntimeError("Native tangent and Y-directional residual FD disagree")
        report = {
            "schema": "scrimp-v130-issue28-structure-v1",
            "status": "PASS", "mpi_ranks": mpi.size,
            "dofs": n, "dynamic_dofs": len(d),
            "algebraic_dofs": len(a), "traction_Y_dofs": len(y),
            "intervals": {name: list(item) for name, item in intervals.items()},
            "initial_Fa_l2": fa_initial, "initial_Fa_block_l2": blocks,
            "max_abs_algebraic_mass_rows": mass_algebraic_max,
            "max_abs_Jaa_Y_columns": kay_max,
            "max_abs_Jda_Y_columns": kdy_max,
            "Y_direction_Fa_FD_l2": fd_a, "Y_direction_Fd_FD_l2": fd_d,
            "Y_direction_Fa_J_l2": j_a, "Y_direction_Fd_J_l2": j_d,
            "Y_direction_Fa_FD_error": err_a,
            "Y_direction_Fd_FD_error": err_d,
            "differential_initial_bitwise_preserved": True,
            "main_TS_executed": False,
            "initializer_TS_executed": False,
            "nonlinear_initializer_qualified": False,
            "differentiation_index_proven": False,
            "conclusion": "Original open-loop Y has a null algebraic Jacobian column but nonzero dynamic coupling.",
        }
        if mpi.rank == 0:
            print("ISSUE28 STRUCTURE PASS: Y absent from Fa, present in Fd", flush=True)
            print(json.dumps({
                "Fa_initial_l2": fa_initial,
                "mass_a_max": mass_algebraic_max,
                "Jaa_Y_max": kay_max, "Jda_Y_max": kdy_max,
                "Y_fd_Fa_l2": fd_a, "Y_fd_Fd_l2": fd_d,
            }, indent=2), flush=True)
        return report
    finally:
        for obj in reversed(made):
            obj.destroy()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--mesh", type=float, default=0.2)
    args = parser.parse_args()
    if not np.isfinite(args.mesh) or args.mesh <= 0:
        parser.error("Invalid mesh resolution")
    mpi = PETSc.COMM_WORLD.tompi4py()
    system = build_system(
        gamma=0., mesh=args.mesh, final_time=.02, save_dt=.002,
        max_step=.002, ramp_start=.003, ramp_end=.009,
        init_step=True, ts_type="bdf", adapt_type="none",
        init_dt=1e-9, rotation=.2,
        initial_profile="tangent", wall_velocity_profile="facet_tangent",
        algebraic_initial_guess="zero", init_method="bdf", height_order=2,
    )
    evidence = {}
    def intercept(self):
        evidence["report"] = inspect(self, mpi)
        raise AuditComplete()
    with patch.object(S.DPHS, "init_step", intercept):
        try:
            system.solve()
        except AuditComplete:
            pass
        else:
            raise RuntimeError("No physical TS must execute during structural audit")
    if "report" not in evidence:
        raise RuntimeError("Native structural audit did not execute")
    if mpi.rank == 0:
        args.out.mkdir(parents=True, exist_ok=True)
        (args.out / "structure.json").write_text(
            json.dumps(evidence["report"], indent=2) + "\n", encoding="utf-8"
        )
    mpi.Barrier()


if __name__ == "__main__":
    main()
