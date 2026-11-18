# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Native MPI read-only audit of the FULL pH-DAE matrix structure.

Two short fixed-CN steps are used ONLY to materialize the normal
GetFEM/PETSc matrices. This is not an optimization tied to fixed dt:
the eventual Schur preconditioner will use TS's varying shift.

Inspect GetFEM's actual unknown ranges and every coupling block in
PETSc mass, stiffness and last implicit Jacobian, including control/
boundary variables. No changes to the operator, no global vector
gather and no dense Mat conversion. This first audit does not infer
J, R, M_c or Q_c by guessing from variable names: it records the
numerical blocks that must be identified and signed before designing
a field-split preconditioner.
"""
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
from unittest.mock import patch

from mpi4py import MPI
import numpy as np

from bench_acoustic_2d import _build_wave
import scrimp as S
from scrimp.utils.ph_block_layout import (
    GROUPS, discover_layout, classify_ph_formulation,
    owned_group_ranges,
    classify_owned_csr, aggregate_role_pairs,
)


class _StopBeforeExampleSolve(Exception):
    def __init__(self, system):
        self.system = system


def _build_coenergy_wave(mesh):
    """Build examples/wave_coenergy.py without its published long run."""
    source = (
        Path(__file__).resolve().parents[1]
        / "examples" / "wave_coenergy.py"
    )
    spec = importlib.util.spec_from_file_location(
        "_scrimp_audit_wave_coenergy", source
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    original_domain = S.Domain

    def short_mesh(kind, params, *a, **kw):
        updated = dict(params)
        if kind == "Concentric":
            updated["h"] = mesh
        return original_domain(kind, updated, *a, **kw)

    def stop_before_solve(system):
        raise _StopBeforeExampleSolve(system)

    with patch.object(S, "Domain", short_mesh), patch.object(
        S.DPHS, "solve", stop_before_solve
    ):
        try:
            module.wave_coenergy_eq()
        except _StopBeforeExampleSolve as captured:
            return captured.system
    raise RuntimeError(
        "wave_coenergy.py did not call DPHS.solve"
    )


def audit(args):
    mpi = MPI.COMM_WORLD
    rank, ranks = mpi.Get_rank(), mpi.Get_size()
    if not (args.mesh > 0 and args.dt > 0 and args.steps >= 2):
        raise ValueError("Require mesh>0, dt>0 and at least two CN steps")
    if os.environ.get("SCRIMP_GETFEM_MATRIX_MODE") != "distributed":
        raise ValueError("Native block audit requires distributed GetFEM")
    if os.environ.get("SCRIMP_GETFEM_BRIDGE_DISTRIBUTED") != "coo":
        raise ValueError("Native block audit requires PETSc COO bridge")

    system = (
        _build_coenergy_wave(args.mesh)
        if args.model == "wave_coenergy"
        else _build_wave(args.mesh)
    )
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
        or system.last_ts_rejections or system.last_ts_snes_failures
    ):
        raise RuntimeError("Reference MUMPS/CN acoustic solve did not converge")

    fields = discover_layout(
        system.gf_model,
        states=system.states,
        costates=system.costates,
        ports=system.ports,
        controls=system.controls,
    )
    formulation = classify_ph_formulation(
        fields, costates=system.costates
    )
    if args.model == "wave_coenergy":
        if formulation["formulation"] != "substituted_coenergy":
            raise RuntimeError(
                "Coenergy example unexpectedly has an explicit -Qc block"
            )
    elif formulation["formulation"] != "implicit_constitutive":
        raise RuntimeError(
            "Mixed wave benchmark unexpectedly lacks independent coenergies"
        )

    # Compare the EXACT Gmsh .msh bytes between independent 1- and
    # 2-rank runs. Same 5790 DoFs do NOT prove identical geometry,
    # cell connectivity, numbering or physical-region assignments.
    # Domain("Concentric") regenerates the file for every run.
    mesh_file = (
        Path.cwd() / "outputs" / "mesh" /
        ("Concentric.msh" if args.model == "wave_coenergy"
         else "Rectangle.msh")
    )
    if rank == 0:
        if not mesh_file.is_file():
            raise FileNotFoundError(
                f"Gmsh reference mesh missing at {mesh_file}"
            )
        mesh_digest = hashlib.sha256(
            mesh_file.read_bytes()
        ).hexdigest()
        mesh_size = mesh_file.stat().st_size
    else:
        mesh_digest = None
        mesh_size = None
    mesh_digest = mpi.bcast(mesh_digest, root=0)
    mesh_size = mpi.bcast(mesh_size, root=0)

    n = int(system.gf_model.nbdof())
    matrices = {
        "flow_mass_F_and_algebraic_zero_rows": system.mass,
        "linear_stiffness_all_bricks": system.stiffness,
        "implicit_TS_Jacobian_last_shift": system.J,
    }
    by_matrix = {}
    ownership = []
    for label, mat in matrices.items():
        if tuple(mat.getSize()) != (n, n):
            raise ValueError(f"{label}: unexpected PETSc matrix shape")
        lo, hi = map(int, mat.getOwnershipRange())
        if mpi.allreduce(hi-lo, op=MPI.SUM) != n:
            raise ValueError(f"{label}: incomplete distributed matrix rows")
        if label == "flow_mass_F_and_algebraic_zero_rows":
            ranges = owned_group_ranges(fields, lo, hi)
            ownership = mpi.gather(
                {"rank": rank, "lo": lo, "hi": hi,
                 "role_owned_ranges": ranges},
                root=0,
            )

        stored, effective, fro2 = classify_owned_csr(
            fields, lo, hi, mat.getValuesCSR()
        )
        global_stored = np.zeros_like(stored)
        global_effective = np.zeros_like(effective)
        global_fro2 = np.zeros_like(fro2)
        mpi.Allreduce(stored, global_stored, op=MPI.SUM)
        mpi.Allreduce(effective, global_effective, op=MPI.SUM)
        mpi.Allreduce(fro2, global_fro2, op=MPI.SUM)
        global_role_nnz = aggregate_role_pairs(
            fields, global_effective
        )
        global_role_fro2 = aggregate_role_pairs(
            fields, global_fro2
        )
        by_matrix[label] = {
            "nnz_stored": int(np.sum(global_stored)),
            "nnz_effective": int(np.sum(global_effective)),
            "variable_order": [f["name"] for f in fields],
            "block_stored_nnz": global_stored.tolist(),
            "block_effective_nnz": global_effective.tolist(),
            "block_frobenius_norm": np.sqrt(global_fro2).tolist(),
            "roles": list(GROUPS),
            "role_effective_nnz": global_role_nnz.astype(
                np.int64
            ).tolist(),
            "role_frobenius_norm": np.sqrt(
                global_role_fro2
            ).tolist(),
        }

    if rank == 0:
        outcome = {
            "schema": "scrimp-ph-block-audit-v1",
            "source_commit": os.environ.get("SCRIMP_BENCH_GIT_SHA", ""),
            "mesh": args.mesh,
            "mpi_ranks": ranks,
            "dofs": n,
            "physical_model": (
                "examples/wave_coenergy.py:wave_coenergy_eq"
                if args.model == "wave_coenergy"
                else "examples/wave.py:wave_eq"
            ),
            "mesh_sha256": mesh_digest,
            "mesh_file_size_bytes": mesh_size,
            "formulation": formulation,
            "ts_type": "cn", "dt": args.dt, "steps": args.steps,
            "final_time": float(system.last_ts_time),
            "linear_mass": bool(system.linear_mass),
            "linear_stiffness": bool(system.linear_stiffness),
            "mumps_policy": getattr(
                system, "_mumps_reset_policy", None
            ),
            "fields": fields,
            "owned_ranges_by_rank": ownership,
            "candidate_preconditioner_partition": (
                {
                    "first": "coenergy",
                    "second": "dynamic + boundary + algebraic",
                    "scope": (
                        "Exploratory ONLY. Boundary saddle-point "
                        "rows require further analysis."
                    ),
                }
                if formulation["formulation"] == "implicit_constitutive"
                else {
                    "first": "dynamic",
                    "second": "boundary + algebraic",
                    "scope": (
                        "NO independent constitutive -Qc exists. "
                        "Explore mass/Dirac/port split; MUMPS global "
                        "remains the only enabled solver."
                    ),
                }
            ),
            "assembled_matrices": by_matrix,
            "caveats": [
                "The final TS Jacobian contains the actual last shift.",
                "The stiffness includes Dirac, constitutive and "
                "boundary/control bricks: it is not pure J-R.",
                "Variable roles do not alone identify mathematical "
                "J,R,M_c,Q_c, and no block reduction is applied.",
                "The explicit zero-mass rows are algebraic equations "
                "and MUST be retained by any preconditioner.",
                "Nonzero patterns and norms do not by themselves "
                "prove discrete adjointness or energy positivity.",
                "A feedback law may use a GetFEM nonlinear brick "
                "even if it is mathematically linear; a zero "
                "linear_stiffness flag does not prove a nonlinear PDE.",
            ],
        }
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(
            json.dumps(outcome, indent=2) + "\n", encoding="utf-8"
        )
        print(
            f"Gmsh input mesh: SHA256={mesh_digest}, "
            f"size={mesh_size} bytes"
        )
        print(
            f"PASS pH blocks: {ranks} MPI ranks, {n} DoFs, "
            f"{len(fields)} distinct GetFEM variable fields"
        )
        print("GetFEM unknown ordering:")
        for field in fields:
            print(
                f"  {field['name']:7s} {field['group']:10s} "
                f"[{field['start']:>8d}, {field['stop']:>8d})"
            )
        print("MPI role couplings (nonzero entries) in K:")
        data = by_matrix["linear_stiffness_all_bricks"]
        for role, row in zip(GROUPS, data["role_effective_nnz"]):
            print(f"  {role:10s} " + " ".join(
                f"{v:10d}" for v in row
            ))
        print(f"Saved: {args.output}")
    return by_matrix


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--model", choices=("wave", "wave_coenergy"),
        default="wave",
        help="Original full wave DAE or substituted-coenergy wave",
    )
    parser.add_argument("--mesh", type=float, default=0.08)
    parser.add_argument("--dt", type=float, default=0.002)
    parser.add_argument("--steps", type=int, default=2)
    parser.add_argument("--output", type=Path, required=True)
    audit(parser.parse_args())


if __name__ == "__main__":
    main()
