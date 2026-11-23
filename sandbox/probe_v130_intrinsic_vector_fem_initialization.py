# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Native GetFEM intrinsic vector FE interpolation + SCRIMP initialization.

This is a SMALL FEM-COEFFICIENT COMPATIBILITY TEST, not a physics,
curl-operator, H(curl)/H(div) trace or Maxwell solver qualification.
The only built-in Nedelec element offered by GetFEM's documented
classic descriptor is FEM_NEDELEC(P), the lowest-order edge element.
There is NO claim that a higher-order or second-kind Nedelec family
has been implemented by GetFEM or exposed by SCRIMP.FEM.
"""

import argparse
import json
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


class OriginalPort:
    def __init__(self, mf, region):
        self.mf = mf
        self.region = region

    def get_flow(self):
        return "q"

    def get_effort(self):
        return "q"

    def get_region(self):
        return self.region

    def get_fem(self):
        return self.mf


class OriginalModelAdapter:
    """Provide only fields used by the actual DPHS initializer method."""

    def __init__(self, model, fem, region):
        self.gf_model = model
        self.ports = {"intrinsic": OriginalPort(fem, region)}
        self.initial_value_set = {"q": False}
        self.lumped_variables = set()

    def set_from_vector(self, name, array):
        # This is the actual production implementation.
        DPHS.set_from_vector(self, name, array)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if PETSc.COMM_WORLD.getSize()!=1:
        raise ValueError("Source-level FE family pilot deliberately MPI1")
    if not args.out.parent.is_dir():
        raise FileNotFoundError("Expected mounted host-owned output folder")
    if args.out.exists():
        raise FileExistsError("Never overwrite native element evidence")

    mesh = gf.Mesh(
        "regular simplices",
        np.linspace(0.,1.,5),
        np.linspace(0.,1.,4),
    )
    # Select a genuine, nontrivial subset of convexes, just as the
    # original annular wave field is restricted to its mesh region.
    region=52
    mesh.set_region(region, np.asarray(mesh.cvid()).reshape(-1)[::2])
    results={}
    for name, descriptor in CASES.items():
        print(f"FEM START {name}: {descriptor}",flush=True)
        mf=gf.MeshFem(mesh,2)
        try:
            mf.set_fem(gf.Fem(descriptor))
        except (AssertionError,ValueError,RuntimeError) as exc:
            results[name]={
                "status":"UNAVAILABLE_IN_THIS_GETFEM_IMAGE",
                "descriptor":descriptor,
                "error":type(exc).__name__+": "+str(exc)[:300],
            }
            print("FEM UNAVAILABLE",name,str(exc)[:180],flush=True)
            continue

        if mf.qdim()!=2:
            raise RuntimeError(f"{name}: intrinsic vector qdim != 2")
        model=gf.Model("real")
        model.add_filtered_fem_variable("q",mf,region)
        dofs=np.asarray(mf.basic_dof_on_region(region)).reshape(-1)
        full=np.asarray(mf.eval("[1.,2.]",globals(),locals())).reshape(-1)
        if full.size!=mf.nbdof():
            raise RuntimeError(f"{name}: bad true GetFEM evaluation size")
        expected=full[dofs]
        if expected.size==0 or np.linalg.norm(expected)<=0:
            raise RuntimeError(f"{name}: physically nontrivial field vanished")
        original=OriginalModelAdapter(model,mf,region)
        # THIS is the actual production DPHS.set_initial_value() method.
        DPHS.set_initial_value(original,"q","[1.,2.]")
        observed=np.asarray(model.variable("q")).reshape(-1)
        if not original.initial_value_set["q"]:
            raise RuntimeError(f"{name}: initialization flag not set")
        if observed.shape!=expected.shape:
            raise RuntimeError(
                f"{name}: actual filtered GetFEM model DOFs {observed.size} "
                f"!= MeshFem DOFs {expected.size}"
            )
        max_error=float(np.max(abs(observed-expected)))
        if not np.isfinite(max_error) or max_error>1e-12:
            raise RuntimeError(f"{name}: coefficient orientation/index error "
                               f"{max_error:.6g}")
        results[name]={
            "status":"PASS_ACTUAL_DPHS_SET_INITIAL_VALUE_ON_GETFEM",
            "descriptor":descriptor,
            "intrinsic_qdim":int(mf.qdim()),
            "mesh_fem_full_dofs":int(mf.nbdof()),
            "filtered_region_dofs":int(dofs.size),
            "filtered_dofs_mod_qdim":int(dofs.size%mf.qdim()),
            "actual_model_variable_dofs":int(observed.size),
            "max_abs_coefficient_error":max_error,
        }
        print(
            f"FEM PASS {name} full={mf.nbdof()} region={len(dofs)} "
            f"qdim={mf.qdim()} exact_coeff_error={max_error:.6g}",flush=True
        )
    if results["rt0"]["status"]!="PASS_ACTUAL_DPHS_SET_INITIAL_VALUE_ON_GETFEM":
        raise RuntimeError("Baseline RT0 unexpectedly unavailable")
    output={
        "schema":"scrimp-v130-intrinsic-vector-family-native-pilot-v1",
        "native_status":"PASS_AVAILABLE_GETFEM_INTRINSIC_FE_COEFFICIENTS",
        "cases":results,
        "Nedelec_second_kind_supported_or_tested":False,
        "higher_order_Nedelec_supported_or_tested":False,
        "H_div_H_curl_PDE_and_interface_traces_qualified":False,
        "MPI_qualified":False,
        "published_snapshots_unchanged":True,
    }
    args.out.write_text(json.dumps(output,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(output,indent=2),flush=True)


if __name__=="__main__":
    main()
