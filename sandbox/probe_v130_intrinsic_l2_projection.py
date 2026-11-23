# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
"""Native diagnosis of intrinsic-vector FE initial data.

MeshFem.eval is documented for Lagrange nodal interpolation; RT/BDM/
Nedelec basis coefficients are not vector samples at those locations.
Here we independently test a global GetFEM L2 mass projection of the
exact constant vector [1,2] on each true intrinsic FE family.

THIS IS NOT a production initializer implementation or a PDE/DAE pass.
It tests only the native GetFEM mass/projection algebra on a small mesh.
"""

import argparse
import json
import traceback
from pathlib import Path

import getfem as gf
import numpy as np
from petsc4py import PETSc

from scrimp.dphs import DPHS


CASES = {
    "rt0": "FEM_RTK(2,0)",
    "rt1": "FEM_RTK(2,1)",
    "rt2": "FEM_RTK(2,2)",
    "bdm1": "FEM_BDMK(2,1)",
    "bdm2": "FEM_BDMK(2,2)",
    "nedelec_lowest": "FEM_NEDELEC(2)",
}


class Port:
    def __init__(self, fem, region):
        self.fem, self.region = fem, region

    def get_flow(self):
        return "q"

    def get_effort(self):
        return "q"

    def get_region(self):
        return self.region

    def get_fem(self):
        return self.fem


class ProductionAdapter:
    """Execute the exact SCRIMP production set_initial_value method."""

    def __init__(self, model, fem, region):
        self.gf_model = model
        self.ports = {"q": Port(fem, region)}
        self.initial_value_set = {"q": False}
        self.lumped_variables = set()

    def set_from_vector(self, name, data):
        DPHS.set_from_vector(self, name, data)


def run_case(mesh, mim, region, name, descriptor):
    row = {"descriptor": descriptor}
    print(f"PROJECTION START {name} {descriptor}", flush=True)
    fem = gf.MeshFem(mesh, 2)
    try:
        fem.set_fem(gf.Fem(descriptor))
    except (AssertionError, ValueError, RuntimeError) as exc:
        row.update(status="UNAVAILABLE_IN_QUALIFIED_GETFEM",
                   reason=f"{type(exc).__name__}: {exc}")
        print(f"PROJECTION UNAVAILABLE {name}: {exc}", flush=True)
        return row

    dofs = np.asarray(fem.basic_dof_on_region(region)).reshape(-1)
    row.update(
        meshfem_nbdof=int(fem.nbdof()),
        meshfem_qdim=int(fem.qdim()),
        filtered_region_dofs=int(len(dofs)),
        filtered_dofs_mod_qdim=int(len(dofs) % fem.qdim()),
    )
    if fem.qdim() != 2 or not len(dofs):
        raise ValueError("Unexpected FEM qdim or no filtered DOFs")
    try:
        row["getfem_is_lagrangian"] = bool(fem.is_lagrangian())
    except Exception as exc:
        row["lagrangian_query_error"] = repr(exc)

    # Capture the genuine raw result; do NOT interpret it as coefficients.
    try:
        raw = np.asarray(fem.eval("[1.,2.]", globals(), locals()))
        row["meshfem_eval_shape"] = list(raw.shape)
        row["meshfem_eval_size"] = int(raw.size)
        row["meshfem_eval_matches_nbdof"] = bool(raw.size == fem.nbdof())
    except Exception as exc:
        row["meshfem_eval_error"] = f"{type(exc).__name__}: {exc}"

    # Isolate actual SCRIMP's current behavior without bypassing its guard.
    production_model = gf.Model("real")
    production_model.add_filtered_fem_variable("q", fem, region)
    adapter = ProductionAdapter(production_model, fem, region)
    try:
        DPHS.set_initial_value(adapter, "q", "[1.,2.]")
        row["production_set_initial_value"] = "RETURNED"
        row["production_coefficient_count"] = int(
            np.asarray(production_model.variable("q")).size)
    except Exception as exc:
        row["production_set_initial_value"] = "REJECTED"
        row["production_error"] = f"{type(exc).__name__}: {exc}"

    # Independent GetFEM GLOBAL L2 projection on the true intrinsic space:
    # find coefficients c with integral (q_h - [1;2]).v_h = 0 for all v_h.
    # This does not assume nodal DOFs or derive coefficients from eval().
    model = gf.Model("real")
    model.add_filtered_fem_variable("q", fem, region)
    model.add_linear_term(mim, "q.Test_q", region)
    model.add_source_term(mim, "[1;2].Test_q", region)
    model.assembly()
    A = np.asarray(model.tangent_matrix().full(), dtype=float)
    b = np.asarray(model.rhs(), dtype=float).reshape(-1)
    n = int(model.nbdof())
    if A.shape != (n, n) or b.shape != (n,) or n != len(dofs):
        raise RuntimeError(
            f"GetFEM assembled dimensions: matrix {A.shape}, rhs {b.shape}, "
            f"filtered DOFs {len(dofs)}")
    if not (np.isfinite(A).all() and np.isfinite(b).all()):
        raise FloatingPointError("Non-finite original GetFEM L2 system")
    symmetry = float(np.linalg.norm(A - A.T) /
                     max(np.linalg.norm(A), 1e-30))
    c = np.linalg.solve(A, b)
    relative_residual = float(np.linalg.norm(A @ c - b) /
                              max(np.linalg.norm(b), 1e-30))
    if not np.isfinite(relative_residual):
        raise FloatingPointError("Non-finite mass-projection residual")
    model.set_variable("q", c)
    stored = np.asarray(model.variable("q")).reshape(-1)
    storage_error = float(np.max(np.abs(stored - c)))
    row.update(
        projected_coefficient_count=int(c.size),
        getfem_mass_matrix_symmetry_rel=symmetry,
        getfem_mass_projection_residual_rel=relative_residual,
        projected_coefficient_norm_l2=float(np.linalg.norm(c)),
        projected_stored_coefficient_error_max=storage_error,
    )
    if (symmetry >= 1e-11 or relative_residual >= 1e-10 or
            storage_error >= 1e-12 or not np.linalg.norm(c) > 0):
        raise RuntimeError("Intrinsic L2 projection failed algebraic checks")
    row["status"] = "PASS_NATIVE_GETFEM_MASS_PROJECTION_ALGEBRA_ONLY"
    print(
        f"PROJECTION PASS {name} dofs={n} "
        f"moment_residual={relative_residual:.3g} "
        f"production={row['production_set_initial_value']}",
        flush=True)
    return row


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if PETSc.COMM_WORLD.getSize() != 1:
        raise RuntimeError("Small dense GetFEM projection diagnostic is MPI1 ONLY")
    if not args.out.parent.is_dir() or args.out.exists():
        raise RuntimeError("Expected a new file in a preexisting output mount")

    mesh = gf.Mesh("regular simplices",
                   np.linspace(0.0, 1.0, 5),
                   np.linspace(0.0, 1.0, 4))
    region = 52
    mesh.set_region(region, np.asarray(mesh.cvid()).reshape(-1)[::2])
    mim = gf.MeshIm(mesh, gf.Integ("IM_TRIANGLE(6)"))
    rows = {}
    for name, descriptor in CASES.items():
        try:
            rows[name] = run_case(mesh, mim, region, name, descriptor)
        except Exception as exc:
            rows[name] = {
                "status": "FAIL_NATIVE_GETFEM_PROJECTION",
                "descriptor": descriptor,
                "error": f"{type(exc).__name__}: {exc}",
                "traceback": traceback.format_exc(),
            }
            print(f"PROJECTION FAIL {name}: {exc}", flush=True)

    passes = [k for k, v in rows.items() if
              v["status"] == "PASS_NATIVE_GETFEM_MASS_PROJECTION_ALGEBRA_ONLY"]
    failures = [k for k, v in rows.items() if
                v["status"] == "FAIL_NATIVE_GETFEM_PROJECTION"]
    baseline_pass = "rt0" in passes
    output = {
        "schema": "scrimp-v130-native-intrinsic-l2-projection-diagnostic-v1",
        "status": ("PASS_NATIVE_GETFEM_MASS_PROJECTION_DIAGNOSTIC_ONLY"
                   if baseline_pass and not failures else "FAILED_DIAGNOSTIC"),
        "cases": rows,
        "available_cases_passed": passes,
        "failed_cases": failures,
        "production_initialization_qualified": False,
        "physical_field_reconstruction_qualified": False,
        "DAE_and_MPI_qualified": False,
        "source_expression_scope": "GWFL constant vector [1;2] only",
        "dense_mass_solve_not_a_production_scalability_strategy": True,
        "original_published_sources_unchanged": True,
    }
    args.out.write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps(output, indent=2), flush=True)
    if not baseline_pass or failures:
        raise RuntimeError("L2 projection diagnostic did not pass all available cases")
    print("PASS GETFEM INTRINSIC L2 DIAGNOSTIC ONLY", flush=True)


if __name__ == "__main__":
    main()
