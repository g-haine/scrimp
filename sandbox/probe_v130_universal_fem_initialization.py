# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Native pilot: one mass-projection algorithm across scalar/vector FE spaces.

Compare the same GetFEM L2 projection against SCRIMP's ACTUAL historical
DPHS.set_initial_value for nodal CG/DG, Hermite, Argyris, Morley, HCT,
nonconforming P1, and intrinsic vector RT/BDM/Nedelec.

Exactly representable polynomial fields are validated by independent
GetFEM quadrature of the reconstructed physical field, not by DOF count
or by the linear-system residual alone.

MPI1, small dense diagnostic ONLY: no change to production initialization,
no claim of canonical Hermite/Argyris derivative interpolation, no
non-polynomial projection promise, and no MPI/PDE/DAE qualification.
"""

import argparse
import json
import traceback
from pathlib import Path

import getfem as gf
import numpy as np
from petsc4py import PETSc

from scrimp.dphs import DPHS


# case name: (GetFEM descriptor, qdim, GWFL expression, old Python syntax,
#             special composite integration, optional availability)
CASES = {
    "cg1_scalar": (
        "FEM_PK(2,1)", 1,
        "1+X(1)-2*X(2)", "1.+x-2.*y", False, False,
    ),
    "cg2_scalar": (
        "FEM_PK(2,2)", 1,
        "1+X(1)*X(1)+X(2)*X(2)", "1.+x*x+y*y", False, False,
    ),
    "dg0_scalar": (
        "FEM_PK_DISCONTINUOUS(2,0)", 1,
        "2", "2.", False, False,
    ),
    "dg2_scalar": (
        "FEM_PK_DISCONTINUOUS(2,2)", 1,
        "1+X(1)*X(1)+X(2)*X(2)", "1.+x*x+y*y", False, False,
    ),
    "cg1_vector": (
        "FEM_PK(2,1)", 2,
        "[X(1)+X(2);2*X(2)-X(1)]", "[x+y,2.*y-x]", False, False,
    ),
    "dg1_vector": (
        "FEM_PK_DISCONTINUOUS(2,1)", 2,
        "[X(1)+X(2);2*X(2)-X(1)]", "[x+y,2.*y-x]", False, False,
    ),
    "hermite2_scalar": (
        "FEM_HERMITE(2)", 1,
        "1+X(1)+X(1)*X(1)+X(2)*X(2)*X(2)",
        "1.+x+x*x+y*y*y", False, True,
    ),
    "argyris2_scalar": (
        "FEM_ARGYRIS", 1,
        "1+X(1)*X(1)*X(2)+X(1)*X(1)*X(1)*X(1)*X(1)",
        "1.+x*x*y+x**5", False, True,
    ),
    "morley2_scalar": (
        "FEM_MORLEY(2)", 1,
        "1+X(1)*X(1)+X(2)*X(2)", "1.+x*x+y*y", False, True,
    ),
    "hct2_scalar": (
        "FEM_HCT_TRIANGLE", 1,
        "1+X(1)+X(2)*X(2)*X(2)",
        "1.+x+y*y*y", True, True,
    ),
    "reduced_hct2_scalar": (
        "FEM_REDUCED_HCT_TRIANGLE", 1,
        "1+X(1)-2*X(2)", "1.+x-2.*y", True, True,
    ),
    "nonconforming_p1_scalar": (
        "FEM_P1_NONCONFORMING", 1,
        "1+X(1)-2*X(2)", "1.+x-2.*y", False, True,
    ),
    "rt0_vector": (
        "FEM_RTK(2,0)", 2,
        "[X(1);X(2)]", "[x,y]", False, False,
    ),
    "bdm1_vector": (
        "FEM_BDMK(2,1)", 2,
        "[X(1);X(2)]", "[x,y]", False, False,
    ),
    "nedelec_vector": (
        "FEM_NEDELEC(2)", 2,
        "[-X(2);X(1)]", "[-y,x]", False, False,
    ),
}


class PortAdapter:
    def __init__(self, fem, region):
        self.fem = fem
        self.region = region

    def get_flow(self):
        return "u"

    def get_effort(self):
        return "u"

    def get_fem(self):
        return self.fem

    def get_region(self):
        return self.region


class ModelAdapter:
    def __init__(self, model, fem, region):
        self.gf_model = model
        self.ports = {"initial": PortAdapter(fem, region)}
        self.initial_value_set = {"u": False}
        self.lumped_variables = set()

    def set_from_vector(self, name, coefficients):
        DPHS.set_from_vector(self, name, coefficients)


def physical_integral(mim, region, model, expression):
    return float(np.asarray(
        gf.asm_generic(mim, 0, expression, region, model)
    ).reshape(-1)[0])


def quantities(qdim, field):
    if qdim == 1:
        return (f"({field})*({field})",
                f"(u-({field}))*(u-({field}))",
                "u*Test_u", f"({field})*Test_u")
    if qdim == 2:
        return (f"({field}).({field})",
                f"(u-({field})).(u-({field}))",
                "u.Test_u", f"({field}).Test_u")
    raise ValueError("Unsupported pilot qdim")


def run_case(mesh, region, name, spec):
    descriptor, qdim, field, python_expr, composite, optional = spec
    print(f"UNIVERSAL START {name} {descriptor} qdim={qdim}",
          flush=True)
    fem = gf.MeshFem(mesh, qdim)
    try:
        fem.set_fem(gf.Fem(descriptor))
    except (AssertionError, ValueError, RuntimeError) as exc:
        if not optional:
            raise RuntimeError(
                f"Mandatory GetFEM descriptor unavailable: {descriptor}"
            ) from exc
        print(f"UNIVERSAL UNAVAILABLE {name}: {exc}", flush=True)
        return {
            "status": "UNAVAILABLE_IN_PINNED_GETFEM",
            "descriptor": descriptor,
            "error": str(exc)[:500],
        }

    mim = gf.MeshIm(mesh, gf.Integ(
        "IM_HCT_COMPOSITE(IM_TRIANGLE(13))"
        if composite else "IM_TRIANGLE(13)"
    ))
    full = int(fem.nbdof())
    region_indices = np.asarray(
        fem.basic_dof_on_region(region)
    ).reshape(-1)
    norm_expr, error_expr, mass_expr, rhs_expr = quantities(qdim,field)
    model = gf.Model("real")
    model.add_filtered_fem_variable("u", fem, region)
    model.add_linear_term(mim, mass_expr, region)
    model.add_source_term(mim, rhs_expr, region)
    model.assembly()
    mass = np.asarray(model.tangent_matrix().full(),dtype=float)
    rhs = np.asarray(model.rhs(),dtype=float).reshape(-1)
    dofs = int(model.nbdof())
    if (dofs != len(region_indices) or mass.shape != (dofs,dofs)
            or rhs.shape != (dofs,) or dofs == 0):
        raise ValueError("True filtered GetFEM coefficient dimensions disagree")
    c = np.linalg.solve(mass,rhs)
    relative_residual = float(np.linalg.norm(mass @ c-rhs) /
                              max(np.linalg.norm(rhs),1e-30))
    model.set_variable("u",c)
    target_squared = physical_integral(mim,region,model,norm_expr)
    error_squared = physical_integral(mim,region,model,error_expr)
    if (not np.isfinite(target_squared) or target_squared <= 0
            or not np.isfinite(error_squared)
            or error_squared < -1e-13*target_squared):
        raise ValueError("Invalid physical quadrature from GetFEM")
    l2_error = float(np.sqrt(max(0.,error_squared)/target_squared))
    row = {
        "fem":descriptor,
        "qdim":qdim,
        "is_lagrangian":bool(fem.is_lagrangian()),
        "expression_gwfl":field,
        "expression_python":python_expr,
        "integration":("IM_HCT_COMPOSITE(IM_TRIANGLE(13))"
                       if composite else "IM_TRIANGLE(13)"),
        "meshfem_total_dofs":full,
        "region_filtered_dofs":dofs,
        "mass_relative_residual":relative_residual,
        "physical_L2_error_relative":l2_error,
        "physical_L2_target_squared":target_squared,
    }

    # The same production function is called in all 15 cases; this
    # does NOT replace the old method or interpret mere return as PASS.
    old_model = gf.Model("real")
    old_model.add_filtered_fem_variable("u",fem,region)
    old = ModelAdapter(old_model,fem,region)
    try:
        DPHS.set_initial_value(old,"u",python_expr)
        row["old_production_returned"] = True
        old_squared = physical_integral(mim,region,old_model,error_expr)
        if not np.isfinite(old_squared) or old_squared < -1e-13*target_squared:
            raise ValueError("Invalid old production physical integral")
        row["old_production_L2_error_relative"] = float(
            np.sqrt(max(0.,old_squared)/target_squared)
        )
    except Exception as exc:
        row["old_production_returned"] = False
        row["old_production_L2_error_relative"] = None
        row["old_production_exception"] = (
            f"{type(exc).__name__}: {exc}"
        )[:500]

    # Exactly representable polynomial: reconstructed physical field
    # must agree, not just number of DOFs or normal-equation residual.
    row["status"] = (
        "PASS_UNIVERSAL_MASS_PROJECTION_AND_PHYSICAL_L2"
        if relative_residual < 1e-10 and l2_error < 5e-8
        else "FAIL_PHYSICAL_PROJECTION"
    )
    # The old production implementation is a SEPARATE diagnostic,
    # never a precondition for the universal L2 projection to pass.
    # A matching coefficient count does not certify a physical field:
    # real CG1/DG1 vector meshes silently returned wrong fields.
    # Preserve a scalar-only baseline to catch unintended regressions
    # of existing scalar nodal evaluation.
    old_error = row["old_production_L2_error_relative"]
    row["old_production_physical_field_valid"] = (
        bool(old_error < 5e-8) if old_error is not None else False
    )
    row["old_production_silent_physical_failure"] = (
        bool(row["old_production_returned"] and old_error is not None
             and old_error >= 5e-8)
    )
    if name in {"cg1_scalar","cg2_scalar","dg0_scalar","dg2_scalar"}:
        if not row["old_production_physical_field_valid"]:
            row["status"]="FAIL_LEGACY_SCALAR_LAGRANGIAN_BASELINE"
    print(f"UNIVERSAL {name}: qdim={qdim} n={dofs} "
          f"mass={relative_residual:.3g} physical={l2_error:.3g} "
          f"old={row['old_production_L2_error_relative']} "
          f"status={row['status']}",flush=True)
    return row


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out",required=True,type=Path)
    args = parser.parse_args()
    if PETSc.COMM_WORLD.getSize() != 1:
        raise RuntimeError("Dense exploratory native pilot is MPI1 only")
    if not args.out.parent.is_dir() or args.out.exists():
        raise RuntimeError("Output must be new file in mounted host folder")
    mesh=gf.Mesh("regular simplices",
                 np.linspace(0.,1.,5),np.linspace(0.,1.,4))
    region=52
    mesh.set_region(region,np.asarray(mesh.cvid()).reshape(-1)[::2])
    results={}
    for name,spec in CASES.items():
        try:
            results[name]=run_case(mesh,region,name,spec)
        except Exception as exc:
            results[name]={
                "status":"FAIL_NATIVE_GETFEM_CASE",
                "error":f"{type(exc).__name__}: {exc}",
                "traceback":traceback.format_exc(),
            }
            print(f"UNIVERSAL FAIL {name}: {exc}",flush=True)
    mandatory=[n for n,desc in CASES.items() if not desc[-1]]
    success=all(results[n]["status"]==
                "PASS_UNIVERSAL_MASS_PROJECTION_AND_PHYSICAL_L2"
                for n in mandatory)
    advanced=[n for n,desc in CASES.items() if desc[-1]]
    success=success and all(
        results[n]["status"] in (
            "PASS_UNIVERSAL_MASS_PROJECTION_AND_PHYSICAL_L2",
            "UNAVAILABLE_IN_PINNED_GETFEM",
        ) for n in advanced
    )
    silent_legacy_failures = [
        n for n, r in results.items()
        if r.get("old_production_silent_physical_failure", False)
    ]
    output={
        "schema":"scrimp-v130-universal-fem-projection-v2",
        "status":"PASS_REQUIRED_AND_ALL_AVAILABLE_FEM" if success else "FAIL",
        "cases":results,
        "mandatory_cases":mandatory,
        "optional_high_regular_cases":advanced,
        "available_passed":[n for n,r in results.items() if r["status"]==
                            "PASS_UNIVERSAL_MASS_PROJECTION_AND_PHYSICAL_L2"],
        "unavailable":[n for n,r in results.items() if r["status"]==
                       "UNAVAILABLE_IN_PINNED_GETFEM"],
        "old_production_silent_physical_failures": silent_legacy_failures,
        "projection_and_legacy_acceptance_separate": True,
        "production_changed":False,
        "physical_scope":"Exactly represented polynomial on selected filtered 2D region",
        "MPI_qualified":False,
        "arbitrary_expression_qualified":False,
        "canonical_derivative_moment_interpolation_qualified":False,
        "dense_diagnostic_not_production_strategy":True,
    }
    args.out.write_text(json.dumps(output,indent=2)+"\n",encoding="utf-8")
    if not success:
        raise RuntimeError("One or more available FE spaces failed physical acceptance")
    print("PASS NATIVE UNIVERSAL SCALAR VECTOR HIGH REGULARITY PROJECTION",flush=True)


if __name__=="__main__":
    main()
