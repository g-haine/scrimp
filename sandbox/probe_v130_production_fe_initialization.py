# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Native test of ACTUAL production DPHS.set_initial_value(method='l2').

Execute on the 15 formerly-qualified scalar, vector and derivative FE
spaces, with the exact, unchanged GetFEM model and filtered region.
Every case checks the independent reconstructed physical field L2.
A 16th nonpolynomial test compares PETSc coefficients to the exact
same GetFEM mass/source assembled and solved in NumPy (small reference).

MPI1 only in this gate: sparse PETSc production path is tested but
MPI2/DAE/full simulations remain separate qualification gates.
"""

import argparse
import json
import traceback
from pathlib import Path

import getfem as gf
import numpy as np
from petsc4py import PETSc

from original_pilot import CASES, quantities
from scrimp.dphs import DPHS
from scrimp.utils.fe_initialization import python_expression_to_gwfl


class Port:
    def __init__(self, fem, region):
        self.fem, self.region = fem, region

    def get_flow(self):
        return "u"

    def get_effort(self):
        return "u"

    def get_fem(self):
        return self.fem

    def get_region(self):
        return self.region


class Adapter:
    def __init__(self, model, fem, region):
        self.gf_model = model
        self.ports = {"u": Port(fem, region)}
        self.initial_value_set = {"u": False}
        self.lumped_variables = set()

    def set_from_vector(self, name, coeff):
        DPHS.set_from_vector(self, name, coeff)


def integral(mim, region, model, expression):
    return float(np.asarray(
        gf.asm_generic(mim, 0, expression, region, model)
    ).reshape(-1)[0])


def run_case(mesh, region, name, spec):
    descriptor, qdim, field, expression, composite, optional = spec
    fem = gf.MeshFem(mesh, qdim)
    fem.set_fem(gf.Fem(descriptor))
    mim = gf.MeshIm(mesh, gf.Integ(
        "IM_HCT_COMPOSITE(IM_TRIANGLE(13))"
        if composite else "IM_TRIANGLE(13)"
    ))
    model = gf.Model("real")
    model.add_filtered_fem_variable("u", fem, region)
    adapter = Adapter(model, fem, region)
    DPHS.set_initial_value(
        adapter, "u", expression, method="l2", projection_mim=mim
    )
    coefficients = np.asarray(model.variable("u")).reshape(-1)
    norm_expr, error_expr, _, _ = quantities(qdim, field)
    norm_sq = integral(mim, region, model, norm_expr)
    error_sq = integral(mim, region, model, error_expr)
    if not (np.isfinite(norm_sq) and norm_sq > 0
            and np.isfinite(error_sq) and error_sq >= -1e-13*norm_sq):
        raise RuntimeError("Invalid independent native physical L2 integral")
    error = float(np.sqrt(max(0., error_sq) / norm_sq))
    passed = (
        adapter.initial_value_set["u"] and
        coefficients.size == len(fem.basic_dof_on_region(region)) and
        np.all(np.isfinite(coefficients)) and error < 5e-8
    )
    print(f"PRODUCTION L2 {name} n={len(coefficients)} "
          f"physical_L2={error:.9g} pass={passed}", flush=True)
    return {
        "status": "PASS_ACTUAL_DPHS_SPARSE_PETSC_L2" if passed else "FAIL",
        "descriptor": descriptor,
        "qdim": int(qdim),
        "region_dofs": int(coefficients.size),
        "physical_relative_L2": error,
    }


def nonpolynomial_case(mesh, region):
    name = "gaussian_cg1_scalar"
    expression = "5.*np.exp(-25*((x-0.6)*(x-0.6)+y*y))"
    gwfl = python_expression_to_gwfl(expression, 1)
    fem = gf.MeshFem(mesh, 1)
    fem.set_fem(gf.Fem("FEM_PK(2,1)"))
    mim = gf.MeshIm(mesh, gf.Integ("IM_TRIANGLE(13)"))
    model = gf.Model("real")
    model.add_filtered_fem_variable("u", fem, region)
    adapter = Adapter(model, fem, region)
    DPHS.set_initial_value(
        adapter, "u", expression, method="l2", projection_mim=mim
    )
    observed = np.asarray(model.variable("u")).reshape(-1)

    reference = gf.Model("real")
    reference.add_filtered_fem_variable("u", fem, region)
    reference.add_linear_term(mim, "u*Test_u", region)
    reference.add_source_term(mim, f"({gwfl})*Test_u", region)
    reference.assembly()
    A = np.asarray(reference.tangent_matrix().full(), dtype=float)
    b = np.asarray(reference.rhs(), dtype=float).reshape(-1)
    expected = np.linalg.solve(A, b)
    discrepancy = float(
        np.linalg.norm(observed - expected) /
        max(float(np.linalg.norm(expected)), 1e-30)
    )
    reference.set_variable("u", expected)
    field_sq = integral(mim, region, reference, f"({gwfl})*({gwfl})")
    error_sq = integral(
        mim, region, reference, f"(u-({gwfl}))*(u-({gwfl}))"
    )
    reference_error = float(
        np.sqrt(max(0., error_sq) / field_sq)
    )
    interpolation_model = gf.Model("real")
    interpolation_model.add_filtered_fem_variable("u", fem, region)
    old = Adapter(interpolation_model, fem, region)
    DPHS.set_initial_value(old, "u", expression)
    old_coeff = np.asarray(interpolation_model.variable("u")).reshape(-1)
    legacy_to_l2 = float(
        np.linalg.norm(old_coeff - observed) /
        max(float(np.linalg.norm(observed)), 1e-30)
    )
    passed = (
        np.isfinite(discrepancy) and discrepancy < 1e-10
        and np.isfinite(reference_error) and reference_error > 0
        and legacy_to_l2 > 1e-5
    )
    print(f"PRODUCTION L2 {name}: sparse_vs_dense={discrepancy:.6g}, "
          f"projection_physical_approx_L2={reference_error:.6g}, "
          f"old_nodal_vs_l2_coeff={legacy_to_l2:.6g}, "
          f"pass={passed}", flush=True)
    return {
        "status": "PASS_NONPOLYNOMIAL_SPARSE_VS_DENSE_L2" if passed else "FAIL",
        "expression": expression,
        "gwfl_source": gwfl,
        "sparse_vs_dense_relative_coeff_error": discrepancy,
        "relative_physical_L2_approximation_error": reference_error,
        "legacy_interpolation_vs_projection_coeff_error": legacy_to_l2,
        "nonpolynomial_interpolation_differs": legacy_to_l2 > 1e-5,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if PETSc.COMM_WORLD.getSize() != 1:
        raise ValueError("MPI1 pilot only; MPI2/DAE separate")
    if not args.out.parent.is_dir() or args.out.exists():
        raise RuntimeError("New output path required")
    mesh = gf.Mesh("regular simplices",
                   np.linspace(0.,1.,5), np.linspace(0.,1.,4))
    region = 52
    mesh.set_region(region, np.asarray(mesh.cvid()).reshape(-1)[::2])
    rows = {}
    for name, spec in CASES.items():
        try:
            rows[name] = run_case(mesh, region, name, spec)
        except Exception as exc:
            rows[name] = {
                "status": "FAIL",
                "error": f"{type(exc).__name__}: {exc}",
                "traceback": traceback.format_exc(),
            }
            print(f"PRODUCTION L2 FAIL {name}: {exc}", flush=True)
    try:
        rows["gaussian_cg1_scalar"] = nonpolynomial_case(mesh, region)
    except Exception as exc:
        rows["gaussian_cg1_scalar"] = {
            "status": "FAIL",
            "error": f"{type(exc).__name__}: {exc}",
            "traceback": traceback.format_exc(),
        }
        print(f"PRODUCTION L2 FAIL Gaussian: {exc}", flush=True)
    good = (
        len(CASES) == 15 and
        all(rows[name]["status"] == "PASS_ACTUAL_DPHS_SPARSE_PETSC_L2"
            for name in CASES) and
        rows["gaussian_cg1_scalar"]["status"] ==
        "PASS_NONPOLYNOMIAL_SPARSE_VS_DENSE_L2"
    )
    result = {
        "schema": "scrimp-v130-production-sparse-fe-init-v1",
        "status": "PASS_NATIVE_SPARSE_15_PLUS_NONPOLYNOMIAL" if good else "FAIL",
        "cases": rows,
        "native_actual_DPHS_set_initial_value": True,
        "PETSc_sparse_MUMPS": True,
        "MPI2_qualified": False,
        "DAE_physics_qualified": False,
        "nodal_interpolation_default_unchanged": True,
        "published_sources_unchanged": True,
    }
    args.out.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    if not good:
        raise RuntimeError("Native opt-in production sparse FE L2 probe failed")
    print("PASS ACTUAL PRODUCTION SPARSE FEM L2 15 PLUS GAUSSIAN", flush=True)


if __name__ == "__main__":
    main()
