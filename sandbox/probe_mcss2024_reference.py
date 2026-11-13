# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Short, faithful numerical smoke of published MCSS 2024 Saint-Venant.

DO NOT modify the published reference source or copy its forms here.
Build the actual validated experiment 3 with the archived Python code;
temporarily intercept its solve call *before its 500-second simulation*
and postprocessing, and then run the identical assembled DPHS through
the CURRENT branch's numerical backend for a short test.

This is NOT a rerun of the historical 2024 environment: PETSc,
GetFEM and SCRIMP runtime are those installed today. A failure
therefore does not disprove the published calculations.

Example (published BDF2/SuperLU setup, first 3 physical steps):
  SCRIMP_GETFEM_RHS_MODE=divide_by_mpi_size mpiexec -n 1 \\
    python -u sandbox/probe_mcss2024_reference.py \\
    --final-time 0.0003 --save-dt 0.0001 --solver superlu
"""

import argparse
import importlib.util
import json
import os
from pathlib import Path
from unittest.mock import patch

import numpy as np
from petsc4py import PETSc
import scrimp as S


class _CapturedBeforePublishedSolve(Exception):
    def __init__(self, system):
        self.system = system


def build_published_model():
    """Run EXACT published model construction, stopping before solve()."""
    source = (
        Path(__file__).resolve().parents[1]
        / "examples" / "published" / "MCSS_2024" / "shallow_water.py"
    )
    spec = importlib.util.spec_from_file_location(
        "_mcss2024_published_shallow_water", source
    )
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot import published source: {source}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    def _trap_solve(system):
        raise _CapturedBeforePublishedSolve(system)

    with patch.object(S.DPHS, "solve", _trap_solve):
        try:
            module.shallow_water(experiment=3, formulation="grad")
        except _CapturedBeforePublishedSolve as captured:
            return captured.system, str(source)
    raise RuntimeError("Published code did not call DPHS.solve()")


def describe_fields(system):
    times = np.asarray(system.solution["t"], dtype=float)
    values = np.asarray(system.solution["z"])
    if not times.size or not values.size:
        return {}
    last = np.asarray(values[-1]).ravel()
    first = np.asarray(values[0]).ravel()
    result = {}
    for name in ("h", "p", "e_h", "e_p", "U_n", "Y_n", "U", "Y"):
        try:
            start, count = system.gf_model.interval_of_variable(name)
        except (KeyError, ValueError, RuntimeError):
            continue
        start, end = int(start), int(start + count)
        old, new = np.real(first[start:end]), np.real(last[start:end])
        if new.size:
            result[name] = {
                "first_norm": float(np.linalg.norm(old)),
                "final_norm": float(np.linalg.norm(new)),
                "final_maxabs": float(np.max(np.abs(new))),
                "final_finite": bool(np.all(np.isfinite(new))),
            }
            if name == "h":
                result[name]["final_min"] = float(np.min(new))
                result[name]["final_max"] = float(np.max(new))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--final-time", type=float, default=0.0003)
    parser.add_argument("--save-dt", type=float, default=0.0001)
    parser.add_argument(
        "--solver", choices=("superlu", "mumps"), default="superlu",
        help="Published SuperLU (default) or MUMPS on same reference model",
    )
    parser.add_argument(
        "--no-init-step", action="store_true",
        help="DIAGNOSTIC ONLY; bypass published pseudo initialization",
    )
    parser.add_argument(
        "--newton-block-trace", action="store_true",
        help="Optional per-block trace of initial nonlinear Newton steps",
    )
    args = parser.parse_args()
    if (not np.isfinite(args.final_time) or args.final_time <= 0
            or not np.isfinite(args.save_dt) or args.save_dt <= 0):
        parser.error("final-time and save-dt must be finite and positive")

    if args.newton_block_trace:
        os.environ["SCRIMP_NEWTON_BLOCK_TRACE"] = "1"

    swe, path = build_published_model()
    # Keep the published finite-element forms, parameters, initial
    # fields, control ports, BDF2 and pseudo-startup EXACTLY as they
    # were constructed. Override only output horizon and factor solver
    # (when explicitly requested). Published dt=1e-4, init_dt=1e-8.
    swe.set_time_scheme(
        t_f=args.final_time,
        dt_save=args.save_dt,
        pc_factor_mat_solver_type=args.solver,
        init_step=not args.no_init_step,
    )
    rank = PETSc.COMM_WORLD.getRank()
    options_to_show = {
        key: swe.time_scheme.getString(key)
        for key in (
            "ts_type", "ts_bdf_order", "ksp_type", "pc_type",
            "pc_factor_mat_solver_type", "dt", "dt_save",
            "init_step", "init_step_ts_type",
            "init_step_nb_iter", "init_step_dt", "t_f",
        )
    }
    if rank == 0:
        print(json.dumps({
            "source": path,
            "experiment": 3,
            "formulation": "grad",
            "note": "Published model, CURRENT SCRIMP backend; "
                    "no physics or control-brick edits",
            "petsc_version": PETSc.Sys.getVersion(),
            "petsc_scalar_complex": bool(
                np.issubdtype(PETSc.ScalarType, np.complexfloating)
            ),
            "mpi_size": PETSc.COMM_WORLD.getSize(),
            "time_scheme": options_to_show,
        }, indent=2), flush=True)

    swe.solve()
    times = np.asarray(swe.solution["t"], dtype=float)
    # The monitor is NOT an authoritative record of the actual TS
    # endpoint: it samples on dt_save, potentially skipping the last
    # PETSc step. A documented run advanced at t=0.0001 and 0.0003,
    # while the last saved monitor sample was only t=0.0001.
    saved_time = float(times[-1]) if times.size else float("-inf")
    reached = float(getattr(swe, "last_ts_time", float("-inf")))
    converged_reason = int(
        getattr(swe, "last_ts_converged_reason", 0)
    )
    if rank == 0:
        print(json.dumps({
            "actual_final_time": reached,
            "last_saved_time": saved_time,
            "ts_converged_reason": converged_reason,
            "ts_accepted_steps": getattr(swe, "last_ts_steps", None),
            "ts_rejected_steps": getattr(swe, "last_ts_rejections", None),
            "ts_snes_failures": getattr(swe, "last_ts_snes_failures", None),
            "requested_final_time": args.final_time,
            "saved_steps": int(times.size),
            "state_fields": describe_fields(swe),
        }, indent=2), flush=True)
    if (converged_reason <= 0 or not np.isfinite(reached) or (
        reached < args.final_time - 1e-10 * max(1.0, args.final_time)
    )):
        raise RuntimeError(
            "Published MCSS2024 model (on current SCRIMP) did not "
            f"reach requested t_f={args.final_time}: PETSc t={reached}, "
            f"reason={converged_reason}, last saved t={saved_time}"
        )


if __name__ == "__main__":
    main()
