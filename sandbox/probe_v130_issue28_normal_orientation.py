# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""#28 native GetFEM boundary-normal orientation / initial p0 trace audit.

Read-only diagnostic on EXACT production Saint-Venant weak-form model.
At t0, with original p0, evaluate GetFEM's ACTUAL boundary Normal,
Tangent=Rot*Normal, and both tangent-sign mismatches via
GetFEM asm_generic. This is not a DAE projection or physical TS run.

In particular, Gmsh-file geometry calculations are NOT evidence
that GetFEM's internal Normal has the same orientation; this is.
"""

import argparse
import json
from pathlib import Path
from unittest.mock import patch

import getfem as gf
import numpy as np
from mpi4py import MPI
from petsc4py import PETSc

import scrimp as S
from viscous_swe_2d_mpi import build_system, REGION_WALL


class Completed(Exception):
    """Expected stop after native operator setup, before initializer."""


def audit(system, mpi, geometry_order):
    t0 = float(system.time_scheme["t_0"])
    model = system.gf_model
    initial = np.asarray(model.from_variables(), dtype=float).copy()
    model.set_time(t0)
    mim = system.domain._int_method[0]

    # All these scalars use native GetFEM GWFL evaluation on the
    # SAME region and geometric mapping as the physical boundary bricks.
    # omega=0.2 and t=0 (before the decreasing control ramp).
    expressions = {
        "perimeter": "1",
        "normal_radial_flux": "X.Normal",
        "velocity_squared": "(p/rho).(p/rho)",
        "velocity_normal_squared": "((p/rho).Normal)*((p/rho).Normal)",
        "aligned_velocity_tangent": "(p/rho).Tangent",
        "aligned_error_squared":
            "((p/rho)-(0.2)*Tangent).((p/rho)-(0.2)*Tangent)",
        "reversed_error_squared":
            "((p/rho)+(0.2)*Tangent).((p/rho)+(0.2)*Tangent)",
    }
    measured = {}
    for name, expression in expressions.items():
        raw = np.asarray(
            gf.asm_generic(mim, 0, expression, REGION_WALL, model)
        )
        if raw.size != 1:
            raise ValueError(
                f"Unexpected scalar GetFEM boundary assembly {name}: {raw.shape}"
            )
        val = float(np.real(raw.reshape(-1)[0]))
        if not np.isfinite(val):
            raise FloatingPointError(f"Nonfinite GetFEM integral {name}")
        measured[name] = val

    baseline = measured["velocity_squared"]
    if baseline <= 0 or measured["perimeter"] <= 0:
        raise ValueError("Non-positive original rotation energy/perimeter")
    for key in ("aligned_error_squared", "reversed_error_squared",
                "velocity_normal_squared"):
        if measured[key] < -1e-10 * baseline:
            raise ValueError(f"Negative boundary norm squared {key}")
    scale = max(1e-25, baseline)
    aligned = np.sqrt(max(0., measured["aligned_error_squared"]) / scale)
    reversed_ = np.sqrt(max(0., measured["reversed_error_squared"]) / scale)
    normal = np.sqrt(max(0., measured["velocity_normal_squared"]) / scale)
    orientation = (
        "ALIGNED" if aligned < 0.5 * reversed_
        else "REVERSED" if reversed_ < 0.5 * aligned
        else "AMBIGUOUS"
    )
    after = np.asarray(model.from_variables(), dtype=float)
    unchanged = bool(np.array_equal(initial, after))
    all_unchanged = bool(mpi.allreduce(unchanged, op=MPI.LAND))
    if not all_unchanged:
        raise RuntimeError("Native boundary observation changed the state")

    # The MPI GetFEM scalar assembly may differ in global-reduction
    # convention; ratios remain meaningful for each rank.
    all_ranks = mpi.allgather({
        "rank": int(mpi.rank),
        "relative_aligned_trace_error": float(aligned),
        "relative_reversed_trace_error": float(reversed_),
        "relative_normal_leakage": float(normal),
        "normal_orientation_result": orientation,
        "normal_radial_flux": measured["normal_radial_flux"],
        "perimeter": measured["perimeter"],
    })
    if any(x["normal_orientation_result"] != orientation for x in all_ranks):
        raise RuntimeError("GetFEM Normal orientation differs between MPI ranks")
    report = {
        "schema": "scrimp-v130-issue28-getfem-normal-sign-v1",
        "source": "production Saint-Venant GetFEM boundary integrals",
        "geometry_order": geometry_order,
        "rotation": 0.2,
        "mpi_ranks": int(mpi.size),
        "native_getfem_integrals_rank0": measured if mpi.rank == 0 else None,
        "relative_aligned_trace_error": float(aligned),
        "relative_reversed_trace_error": float(reversed_),
        "relative_normal_leakage": float(normal),
        "normal_orientation_result": orientation,
        "rank_local_measurements": all_ranks,
        "initial_state_bitwise_preserved": all_unchanged,
        "nonlinear_initializer_executed": False,
        "physical_time_step_executed": False,
        "conclusion": (
            "Sign validation uses actual GetFEM Normal; "
            "this does not establish complete DAE compatibility."
        ),
    }
    if mpi.rank == 0:
        print(
            f"ISSUE28 NATIVE NORMAL order={geometry_order} MPI={mpi.size}: "
            f"{orientation}; aligned={aligned:.9g}, reversed={reversed_:.9g}, "
            f"normal_leakage={normal:.9g}; "
            f"integral X.Normal={measured['normal_radial_flux']:.9g}",
            flush=True,
        )
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--geometry-order", type=int, choices=(1, 2), required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    mpi = PETSc.COMM_WORLD.tompi4py()

    system = build_system(
        gamma=0.0, mesh=0.2, geometry_order=args.geometry_order,
        final_time=.02, save_dt=.002, max_step=.002,
        ramp_start=.003, ramp_end=.009, init_step=True,
        ts_type="bdf", adapt_type="none",
        init_dt=1e-9, rotation=0.2,
        initial_profile="tangent", wall_velocity_profile="facet_tangent",
        algebraic_initial_guess="zero", init_method="bdf",
        height_order=2,
    )
    output = {}
    def intercept(self):
        output["result"] = audit(self, mpi, args.geometry_order)
        raise Completed
    with patch.object(S.DPHS, "init_step", intercept):
        try:
            system.solve()
        except Completed:
            pass
        else:
            raise RuntimeError("An initializer or time step was unexpectedly run")
    if "result" not in output:
        raise RuntimeError("Native normal audit did not execute")
    if mpi.rank == 0:
        args.out.mkdir(parents=True, exist_ok=True)
        (args.out / "normal_orientation.json").write_text(
            json.dumps(output["result"], indent=2) + "\n", encoding="utf-8"
        )
    mpi.Barrier()


if __name__ == "__main__":
    main()
