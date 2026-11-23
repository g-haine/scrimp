# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Independent exact nonconstant intrinsic FE L2-projection pilot.

RT/BDM: affine radial vector [X(1);X(2)].
Lowest Nedelec: affine rotation [-X(2);X(1)].
Both are exactly representable in their chosen spaces. This tests
GetFEM assembly and reconstructed-field quadrature on MPI1, not
production set_initial_value, general GWFL, or time integration.
"""

import argparse
import json
import traceback
from pathlib import Path

import getfem as gf
import numpy as np
from petsc4py import PETSc


CASES = {
    "rt0": ("FEM_RTK(2,0)", "[X(1);X(2)]"),
    "rt1": ("FEM_RTK(2,1)", "[X(1);X(2)]"),
    "rt2": ("FEM_RTK(2,2)", "[X(1);X(2)]"),
    "bdm1": ("FEM_BDMK(2,1)", "[X(1);X(2)]"),
    "bdm2": ("FEM_BDMK(2,2)", "[X(1);X(2)]"),
    "nedelec_lowest": ("FEM_NEDELEC(2)", "[-X(2);X(1)]"),
}


def test_case(mesh, mim, region, name, descriptor, expression):
    fem = gf.MeshFem(mesh, 2)
    fem.set_fem(gf.Fem(descriptor))
    if fem.qdim() != 2 or fem.is_lagrangian():
        raise RuntimeError("Expected intrinsic, non-Lagrangian vector FE")
    model = gf.Model("real")
    model.add_filtered_fem_variable("q", fem, region)
    model.add_linear_term(mim, "q.Test_q", region)
    model.add_source_term(mim, f"({expression}).Test_q", region)
    model.assembly()
    mass = np.asarray(model.tangent_matrix().full(), dtype=float)
    rhs = np.asarray(model.rhs(), dtype=float).reshape(-1)
    n = int(model.nbdof())
    if mass.shape != (n,n) or rhs.shape != (n,) or n <= 0:
        raise RuntimeError("GetFEM intrinsic DOF/mass/RHS mismatch")
    c = np.linalg.solve(mass, rhs)
    residual = float(np.linalg.norm(mass @ c - rhs) /
                     max(np.linalg.norm(rhs), 1e-30))
    model.set_variable("q", c)
    target_sq = float(np.asarray(gf.asm_generic(
        mim, 0, f"({expression}).({expression})", region, model
    )).reshape(-1)[0])
    err_sq = float(np.asarray(gf.asm_generic(
        mim, 0, f"(q-({expression})).(q-({expression}))", region, model
    )).reshape(-1)[0])
    if not (np.isfinite(target_sq) and target_sq > 0
            and np.isfinite(err_sq) and err_sq >= -1e-13*target_sq):
        raise RuntimeError("Invalid original GetFEM L2 field integral")
    error = float(np.sqrt(max(0.,err_sq)/target_sq))
    result = {
        "status": "PASS_NATIVE_AFFINE_FIELD_L2" if (
            residual < 1e-10 and error < 5e-8
        ) else "FAIL_NATIVE_AFFINE_FIELD_L2",
        "fem_descriptor": descriptor,
        "expression": expression,
        "filtered_DOFs": n,
        "mass_relative_residual": residual,
        "field_relative_L2_error": error,
        "target_L2_squared": target_sq,
        "field_L2_error_squared": err_sq,
    }
    print(f"AFFINE {name}: DOFs={n}, mass={residual:.4g}, "
          f"physical_L2={error:.4g}, status={result['status']}", flush=True)
    return result


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()
    if PETSc.COMM_WORLD.getSize() != 1:
        raise RuntimeError("This small dense affine pilot is MPI1 ONLY")
    if not a.out.parent.is_dir() or a.out.exists():
        raise RuntimeError("Refuse missing output directory or overwrite")
    mesh = gf.Mesh("regular simplices", np.linspace(0.,1.,5),
                   np.linspace(0.,1.,4))
    region = 52
    mesh.set_region(region, np.asarray(mesh.cvid()).reshape(-1)[::2])
    mim = gf.MeshIm(mesh, gf.Integ("IM_TRIANGLE(6)"))
    results = {}
    for name, (descriptor, expr) in CASES.items():
        try:
            results[name] = test_case(mesh, mim, region, name,
                                      descriptor, expr)
        except Exception as exc:
            results[name] = {
                "status":"FAIL_NATIVE_AFFINE_FIELD_L2",
                "error":f"{type(exc).__name__}: {exc}",
                "traceback":traceback.format_exc(),
            }
            print(f"AFFINE {name} FAIL: {exc}", flush=True)
    passed = all(case["status"]=="PASS_NATIVE_AFFINE_FIELD_L2"
                 for case in results.values())
    output = {
        "schema":"scrimp-v130-intrinsic-affine-physical-l2-v1",
        "status":"PASS_NATIVE_AFFINE_SIX_FAMILY" if passed else "FAIL",
        "cases":results,
        "MPI_qualified":False,
        "general_expression_qualified":False,
        "production_initialization_changed":False,
        "dense_diagnostic_not_production_implementation":True,
    }
    a.out.write_text(json.dumps(output,indent=2)+"\n",encoding="utf-8")
    if not passed:
        raise RuntimeError("Native exact affine FE projection failed")
    print("PASS NATIVE SIX INTRINSIC AFFINE PHYSICAL L2", flush=True)


if __name__=="__main__":
    main()
