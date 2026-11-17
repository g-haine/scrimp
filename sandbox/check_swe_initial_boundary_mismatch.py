# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Measure the initial rigid-rotation / polygonal-wall trace mismatch.

No SCRIMP DPHS, TS, SNES, KSP, MUMPS, initialization or PDE solve is
performed. We generate exactly the same SCRIMP Disk meshes as the
viscous Saint-Venant sandbox and use GetFEM boundary quadrature to
measure the *continuous* initial velocity trace vs the prescribed
wall velocity on the DISCRETE mesh.

For the default decreasing boundary ramp f(0)=1, unit disk R=1,
rho=1, and clockwise Rot(x,y)=(y,-x):

    u_initial = omega * Rot * X
    u_wall    = omega * Rot * Normal       (with U_n=0)

In particular, u_initial.Normal is generally nonzero on a polygon,
even though it vanishes on the exact circle. The squared discrepancy
decomposes into its normal and tangential parts:

    |u_initial-u_wall|^2 = omega^2 * |X-Normal|^2
    (u_initial.Normal)^2 = omega^2 * (y*n_x-x*n_y)^2
    ((u_initial-u_wall).Rot*Normal)^2
                          = omega^2 * (X.Normal-1)^2

This diagnostic measures continuous geometric mismatch; it does NOT
measure finite-element projection error, the multiplier traction Y,
the stability/inf-sup constants, or a discrete physical solution.

Example (one MPI rank only):
  python -u sandbox/check_swe_initial_boundary_mismatch.py \\
      --meshes 0.2 0.15 0.1 --rotation 0.02 \\
      --output sandbox/outputs/symbolic-reuse/boundary_trace.json
"""

import argparse
import json
import math
from pathlib import Path

import getfem as gf
from petsc4py import PETSc
import scrimp as S

REGION_WALL = 10


def assembly(mim, expr):
    """Evaluate one scalar boundary integral with GetFEM GWFL."""
    value = gf.asm_generic(mim, 0, expr, REGION_WALL)
    return float(value)


def mismatch_for_mesh(mesh_size, rotation):
    # Domain uses the same SCRIMP Disk(R=1,h=mesh_size) / Gmsh mesh
    # generation and MeshIm as viscous_swe_2d_mpi.build_system().
    domain = S.Domain("Disk", {"R": 1.0, "h": mesh_size})
    mesh = domain.get_mesh()[0]
    mim = domain._int_method[0]

    # The rotation omega and the initial-momentum rho cancel from the
    # geometric factors, which are useful to compare meshes independently.
    tangent_sq = (
        "((X(1)*Normal(1)+X(2)*Normal(2)-1)"
        "*(X(1)*Normal(1)+X(2)*Normal(2)-1))"
    )
    normal_sq = (
        "((X(2)*Normal(1)-X(1)*Normal(2))"
        "*(X(2)*Normal(1)-X(1)*Normal(2)))"
    )
    full_sq = (
        "((X(1)-Normal(1))*(X(1)-Normal(1))"
        "+(X(2)-Normal(2))*(X(2)-Normal(2)))"
    )
    boundary_measure = assembly(mim, "1")
    if not (boundary_measure > 0 and math.isfinite(boundary_measure)):
        raise RuntimeError("GetFEM boundary integration returned invalid length")

    geometric_full_sq = assembly(mim, full_sq)
    geometric_normal_sq = assembly(mim, normal_sq)
    geometric_tangent_sq = assembly(mim, tangent_sq)
    for name, value in (
        ("full", geometric_full_sq),
        ("normal", geometric_normal_sq),
        ("tangential", geometric_tangent_sq),
    ):
        if value < -1e-10 or not math.isfinite(value):
            raise RuntimeError(f"Invalid GetFEM {name} integral: {value}")

    # Orthogonal decomposition validates the three boundary expressions
    # and the signed normal convention, independently of solver behavior.
    discrepancy = abs(
        geometric_full_sq - geometric_normal_sq - geometric_tangent_sq
    )
    if discrepancy > 1e-8 * max(
        geometric_full_sq, geometric_normal_sq, geometric_tangent_sq, 1e-12
    ):
        raise RuntimeError(
            "Normal/tangential decomposition failed: "
            f"full={geometric_full_sq}, normal={geometric_normal_sq}, "
            f"tangent={geometric_tangent_sq}"
        )

    omega = abs(rotation)
    return {
        "mesh_h": mesh_size,
        "boundary_measure": boundary_measure,
        "rotation": rotation,
        "initial_wall_ramp_value": 1.0,
        "l2_full_velocity_trace_mismatch":
            omega * math.sqrt(max(geometric_full_sq, 0.0)),
        "l2_normal_velocity_violation":
            omega * math.sqrt(max(geometric_normal_sq, 0.0)),
        "l2_tangential_velocity_mismatch":
            omega * math.sqrt(max(geometric_tangent_sq, 0.0)),
        "rms_full_velocity_trace_mismatch":
            omega * math.sqrt(max(geometric_full_sq, 0.0) / boundary_measure),
        "rms_normal_velocity_violation":
            omega * math.sqrt(max(geometric_normal_sq, 0.0) / boundary_measure),
        "rms_tangential_velocity_mismatch":
            omega * math.sqrt(max(geometric_tangent_sq, 0.0) / boundary_measure),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--meshes", nargs="+", type=float, default=(0.2, 0.15, 0.1),
    )
    parser.add_argument("--rotation", type=float, default=0.02)
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args()
    if PETSc.COMM_WORLD.getSize() != 1:
        parser.error("Run this independent mesh diagnostic on exactly one rank")
    if any(not math.isfinite(h) or h <= 0 for h in args.meshes):
        parser.error("All mesh sizes must be positive finite numbers")
    if not math.isfinite(args.rotation):
        parser.error("Rotation must be finite")

    results = [mismatch_for_mesh(h, args.rotation) for h in args.meshes]
    report = {
        "method": "GetFEM asm_generic exact boundary quadrature; no PDE solve",
        "geometry": "SCRIMP Disk(R=1), polygonal boundary",
        "results": results,
    }
    encoded = json.dumps(report, indent=2) + "\n"
    print(encoded, end="", flush=True)
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(encoded)


if __name__ == "__main__":
    main()
