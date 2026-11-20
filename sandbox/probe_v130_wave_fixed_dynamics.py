# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Test discrete algebraic feasibility while preserving the wave q,p initial data.

This is a DIAGNOSTIC. It never changes DPHS production code or the
accepted pseudo-TS state used by the following original CN solve.
At t0 it evaluates actual IFunction algebraic rows with zero z_dot,
extracts the assembled PETSc stiffness subblock K_aa, and asks whether
a correction *only to algebraic unknowns* can remove those residuals.

A failed LSMR solve with large condition estimate is INCONCLUSIVE,
not proof that the wave's prescribed momentum is incompatible. Actual
IFunction is re-evaluated at the proposed state, verifying the
Jacobian extraction and whether the physical constraints are satisfied.
"""
import argparse
import json
import os
from pathlib import Path
import types
from unittest.mock import patch

from mpi4py import MPI
import numpy as np
from scipy.sparse import csr_matrix, vstack
import scrimp as S

from audit_v130_wave_dae_initial import _build_wave, _intervals
from solve_v130_dae_algebraic_block import solve_algebraic_block


# Diagnostic comparison: the published left-boundary velocity is zero
# at t=0, whereas the initial Gaussian momentum has nonzero trace.
# This control variant preserves EXACTLY the same q0/p0 and PDE but
# matches U_L(0,y)=p0(0,y)/rho(0,y). Never modify examples/wave.py.
def _build_scenario(mesh, scenario):
    if scenario == "published":
        return _build_wave(mesh)
    if scenario != "matched_left_control":
        raise ValueError(f"Unknown physical diagnostic scenario {scenario}")
    set_control = S.DPHS.set_control
    observed = []

    def replace_left_control(system, name, expression):
        if name == "Boundary control (left)":
            observed.append(expression)
            # rho=3-x, so rho(0,y)=3. The Gaussian's trace is
            # 3**(-20*(0.25+(y-0.5)**2)) at the left boundary.
            expression = (
                "3**(-20*(0.25+(y-0.5)*(y-0.5)))/3."
            )
        return set_control(system, name, expression)

    with patch.object(S.DPHS, "set_control", replace_left_control):
        system = _build_wave(mesh)
    if len(observed) != 1:
        raise RuntimeError("Expected one and only one left velocity control")
    return system


def _algebraic_indices(algebraic):
    result = np.concatenate([
        np.arange(start, start + count, dtype=np.int64)
        for start, count in algebraic.values()
    ])
    return np.sort(result)


def _owned_residual_and_block(system, state, global_indices, mpi,
                              *, extract_matrix=False):
    """Collect actual algebraic residual and optional global K_aa on rank 0."""
    t0 = float(system.time_scheme["t_0"])
    vec = system.stiffness.createVecRight()
    lo, hi = map(int, vec.getOwnershipRange())
    if hi > lo:
        vec.setValues(
            range(lo, hi),
            np.asarray(state[lo:hi], dtype=np.float64),
        )
    vec.assemble()
    zd = vec.duplicate()
    zd.set(0)
    residual = vec.duplicate()
    residual.set(0)
    system.gf_model.set_time(t0)
    system.gf_model.to_variables(
        np.array(state, dtype=np.float64, copy=True)
    )

    try:
        system.IFunction(None, t0, vec, zd, residual)
        owned = global_indices[
            (global_indices >= lo) & (global_indices < hi)
        ]
        vals = np.asarray(
            residual.getArray(readonly=True)[owned - lo]
        ).copy()
        extracted = None
        if extract_matrix:
            ptr, cols, coefficients = system.stiffness.getValuesCSR()
            K = csr_matrix(
                (coefficients, cols, ptr),
                shape=(hi - lo, int(vec.getSize())),
            )
            extracted = K[owned - lo, :][:, global_indices].tocsr()
        pieces = mpi.gather((owned, vals, extracted), root=0)
    finally:
        vec.destroy()
        zd.destroy()
        residual.destroy()

    if mpi.Get_rank() != 0:
        return None, None
    rows = np.concatenate([p[0] for p in pieces])
    values = np.concatenate([p[1] for p in pieces])
    sort = np.argsort(rows)
    if not np.array_equal(rows[sort], global_indices):
        raise RuntimeError("Incomplete/overlapping PETSc owned algebraic rows")
    r = values[sort]
    block = None
    if extract_matrix:
        blocks = [p[2] for p in pieces]
        block = vstack(blocks, format="csr")[sort, :].tocsr()
        if block.shape != (len(global_indices), len(global_indices)):
            raise RuntimeError("Incorrect extracted algebraic block")
    return r, block


def _collect_probe(system, before, accepted, dynamic, algebraic, mpi):
    aidx = _algebraic_indices(algebraic)
    r_initial, Kaa = _owned_residual_and_block(
        system, before, aidx, mpi, extract_matrix=True
    )
    if mpi.Get_rank() == 0:
        try:
            correction, solve = solve_algebraic_block(
                Kaa, r_initial
            )
            imag = float(
                np.max(np.abs(np.imag(correction)))
            )
            if imag > 1e-10:
                raise RuntimeError(
                    "Complex algebraic correction for real physical model: "
                    f"{imag:.6g}"
                )
            correction = np.asarray(
                np.real(correction), dtype=np.float64
            )
            result = {
                "error": None,
                "correction": correction,
                "solve": solve,
                "linear_prediction": r_initial + Kaa @ correction,
                "baseline_residual": r_initial,
            }
        except Exception as exc:
            result = {"error": f"{type(exc).__name__}: {exc}"}
    else:
        result = None
    result = mpi.bcast(result, root=0)
    if result["error"]:
        raise RuntimeError(result["error"])
    candidate = np.asarray(before, dtype=np.float64).copy()
    candidate[aidx] += result["correction"]

    r_candidate, _ = _owned_residual_and_block(
        system, candidate, aidx, mpi
    )
    r_accepted, _ = _owned_residual_and_block(
        system, accepted, aidx, mpi
    )
    report = None
    if mpi.Get_rank() == 0:
        linear_prediction = result["linear_prediction"]
        native_norm = float(np.linalg.norm(r_candidate))
        prediction_mismatch = float(np.linalg.norm(
            r_candidate - linear_prediction
        ))
        scale = max(
            1.0, float(np.linalg.norm(result["baseline_residual"]))
        )
        prediction_ok = prediction_mismatch <= 1e-8 * scale
        actual_compatible = bool(
            native_norm <= result["solve"]["requested_tolerance"]
            and prediction_ok
        )
        delta = {
            name: float(np.linalg.norm(
                candidate[start:start + count]
                - before[start:start + count]
            ))
            for name, (start, count) in dynamic.items()
        }
        if any(value != 0.0 for value in delta.values()):
            raise RuntimeError(
                "Algebraic-only projection modified prescribed dynamic data"
            )
        groups = {}
        initial_residual = result["baseline_residual"]
        for name, (start, count) in algebraic.items():
            mask = (aidx >= start) & (aidx < start + count)
            if int(np.count_nonzero(mask)) != count:
                raise RuntimeError(f"Invalid algebraic residual block {name}")
            values = {
                "before": initial_residual[mask],
                "projected": r_candidate[mask],
                "pseudo": r_accepted[mask],
            }
            groups[name] = {
                "rows": count,
                "residual_l2": {
                    stage: float(np.linalg.norm(array))
                    for stage, array in values.items()
                },
                "residual_max_abs": {
                    stage: float(np.max(np.abs(array)))
                    for stage, array in values.items()
                },
            }
        report = {
            "schema": "scrimp-v130-fixed-dynamics-algebraic-probe-v1",
            "image_source_revision": os.environ.get(
                "SCRIMP_IMAGE_REVISION", ""
            ),
            "mpi_ranks": mpi.Get_size(),
            "global_dofs": len(before),
            "algebraic_dofs": int(len(aidx)),
            "algebraic_residual_by_unknown": groups,
            "algebraic_intervals": {
                k: [int(v[0]), int(v[0] + v[1])]
                for k, v in algebraic.items()
            },
            "dynamic_data_change_l2": delta,
            "linear_least_squares": result["solve"],
            "independent_native_residual_after_project_l2": native_norm,
            "native_residual_after_pseudo_l2": float(
                np.linalg.norm(r_accepted)
            ),
            "jacobian_prediction_mismatch_l2": prediction_mismatch,
            "linearization_validated": prediction_ok,
            "fixed_dynamics_actual_algebraic_feasibility": (
                "PASS" if actual_compatible else "UNRESOLVED"
            ),
            "caution": (
                "UNRESOLVED can reflect inconsistent data, a rank-"
                "deficient or ill-conditioned algebraic block, or "
                "insufficient iterative accuracy. It is not a "
                "mathematical incompatibility certificate."
            ),
        }
    # Production must see the exact accepted pseudo-TS model variables.
    system.gf_model.set_time(float(system.time_scheme["t_0"]))
    system.gf_model.to_variables(
        np.array(accepted, dtype=np.float64, copy=True)
    )
    return report


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--mesh", type=float, default=0.2)
    p.add_argument(
        "--scenario", choices=("published", "matched_left_control"),
        default="published",
        help="Second scenario changes only U_L(t) to match p0/rho at t=0; "
             "all PDE bricks and initial q/p are unchanged",
    )
    p.add_argument("--output-dir", type=Path, required=True)
    args = p.parse_args()
    if not (0 < args.mesh <= 1.0):
        p.error("mesh must lie in (0,1]")
    mpi = MPI.COMM_WORLD
    os.environ.setdefault("SCRIMP_INIT_FACTOR_SOLVER_TYPE", "mumps")

    system = _build_scenario(args.mesh, args.scenario)
    system.set_time_scheme(
        ts_type="cn", ts_adapt_type="none",
        ts_exact_final_time="matchstep",
        dt=0.002, dt_save=0.002, t_f=0.01,
        ksp_type="preonly", pc_type="lu",
        pc_factor_mat_solver_type="mumps",
    )

    original_init = system.init_step
    capture = {}

    def diagnostic_initial_step(self):
        _, dynamic, algebraic = _intervals(self)
        before = np.array(
            self.gf_model.from_variables(),
            dtype=np.float64, copy=True,
        )
        original_init()
        accepted = np.array(
            self.gf_model.from_variables(),
            dtype=np.float64, copy=True,
        )
        result = _collect_probe(
            self, before, accepted, dynamic, algebraic, mpi
        )
        if mpi.Get_rank() == 0:
            capture["result"] = result

    system.init_step = types.MethodType(
        diagnostic_initial_step, system
    )
    system.solve()
    if (int(system.last_ts_converged_reason) <= 0
            or int(system.last_ts_steps) != 5
            or system.last_ts_rejections
            or system.last_ts_snes_failures):
        raise RuntimeError("Ordinary main CN integration failed")
    if mpi.Get_rank() == 0:
        if "result" not in capture:
            raise RuntimeError("Missing algebraic projection diagnostic")
        report = capture["result"]
        report.update({
            "scenario": args.scenario,
            "boundary_velocity_at_initial_time": (
                "published U_L(0,y)=0" if args.scenario == "published"
                else "matched U_L(0,y)=p0(0,y)/rho(0,y)"
            ),
            "mesh_h": args.mesh, "cn_steps": 5,
            "time_step": 0.002,
            "main_ts_converged_reason": int(
                system.last_ts_converged_reason
            ),
        })
        args.output_dir.mkdir(parents=True, exist_ok=True)
        suffix = (
            "" if args.scenario == "published"
            else f"_{args.scenario}"
        )
        destination = args.output_dir / (
            f"wave_algebraic_only_{mpi.Get_size()}r{suffix}.json"
        )
        destination.write_text(
            json.dumps(report, indent=2) + "\n",
            encoding="utf-8",
        )
        print(
            f"DAE ALGEBRAIC-ONLY {mpi.Get_size()}r "
            f"scenario={args.scenario}: "
            f"initial_r={report['linear_least_squares']['raw_constraint_residual_l2']:.8g}, "
            f"projected_native_r={report['independent_native_residual_after_project_l2']:.8g}, "
            f"accepted_pseudo_r={report['native_residual_after_pseudo_l2']:.8g}, "
            f"linearization_error={report['jacobian_prediction_mismatch_l2']:.8g}, "
            f"fixed_dynamic_feasibility="
            f"{report['fixed_dynamics_actual_algebraic_feasibility']}; "
            f"dual_lower_bound_R="
            f"{report['linear_least_squares']['dual_lower_bound_for_corrections_within_radius_l2']:.9g}",
            flush=True,
        )


if __name__ == "__main__":
    main()
