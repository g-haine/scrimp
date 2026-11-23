# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Native MPI1/2 probe of the default sparse GetFEM/PETSc FE projection.

Uses the actual DPHS.set_initial_value() with no method keyword.
Runs selected scalar, vector, intrinsic, derivative-DOF and
nonpolynomial FE fields on the SAME original 2-D filtered mesh.

Measures: rank agreement, reconstructed physical relative L2,
MPI2-vs-MPI1 *coefficient* agreement, actual GetFEM RHS
calibration and preservation of numeric lumped variables.

This is not a DAE/Heat-Wave physics or general MPI scalability claim.
"""

import argparse
import json
import traceback
from pathlib import Path

import getfem as gf
import numpy as np
from mpi4py import MPI
from petsc4py import PETSc

from original_pilot import CASES, quantities
from scrimp.dphs import DPHS, _getfem_mpi_rhs_scale
from scrimp.utils.fe_initialization import python_expression_to_gwfl


SELECTED = (
    "cg2_scalar", "cg1_vector", "rt0_vector", "bdm1_vector",
    "nedelec_vector", "argyris2_scalar",
)
GAUSSIAN = "5.*np.exp(-25*((x-0.6)*(x-0.6)+y*y))"


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


class LumpedModel:
    def __init__(self):
        self.value = None

    def set_variable(self, name, value):
        if name != "lumped":
            raise ValueError("Wrong lumped name")
        self.value = np.asarray(value, dtype=np.float64).copy()


class LumpedAdapter:
    def __init__(self):
        self.gf_model = LumpedModel()
        self.initial_value_set = {"lumped": False}
        self.lumped_variables = {"lumped"}


def integral(mim, region, model, expression):
    return float(np.asarray(
        gf.asm_generic(mim, 0, expression, region, model)
    ).reshape(-1)[0])


def fe_case(mesh, region, key, spec, mpi_comm):
    descriptor, qdim, field, expression, composite, _ = spec
    fem = gf.MeshFem(mesh, qdim)
    fem.set_fem(gf.Fem(descriptor))
    mim = gf.MeshIm(
        mesh, gf.Integ(
            "IM_HCT_COMPOSITE(IM_TRIANGLE(13))"
            if composite else "IM_TRIANGLE(13)"
        ),
    )
    model = gf.Model("real")
    model.add_filtered_fem_variable("u", fem, region)
    obj = Adapter(model, fem, region)
    DPHS.set_initial_value(
        obj, "u", expression, projection_mim=mim
    )
    c = np.asarray(model.variable("u"), dtype=np.float64).reshape(-1)
    if not obj.initial_value_set["u"] or not np.all(np.isfinite(c)):
        raise RuntimeError("Production default failed to initialize")
    ref, diff, _, _ = quantities(qdim, field)
    ref_sq = integral(mim, region, model, ref)
    diff_sq = integral(mim, region, model, diff)
    if not np.isfinite(ref_sq) or ref_sq <= 0 or (
        not np.isfinite(diff_sq) or diff_sq < -1e-13*ref_sq
    ):
        raise RuntimeError("Invalid physical field GetFEM quadrature")
    physical = float(np.sqrt(max(0., diff_sq)/ref_sq))
    gathered = mpi_comm.allgather(c)
    if any(x.shape != c.shape for x in gathered):
        raise RuntimeError("MPI ranks have different global FE basis size")
    rank_error = max(
        float(np.linalg.norm(other-c)/max(np.linalg.norm(c),1e-30))
        for other in gathered
    )
    if physical >= 5e-8 or rank_error >= 1e-10:
        raise RuntimeError(
            f"Physical or MPI rank error {physical}, {rank_error}"
        )
    print(f"MPI{mpi_comm.size} R{mpi_comm.rank} {key} "
          f"DOFs={c.size} physical_L2={physical:.4g} "
          f"rank_error={rank_error:.4g}",flush=True)
    return {
        "status": "PASS_NATIVE_MPI_FE_FIELD",
        "descriptor": descriptor,
        "expression": expression,
        "global_coefficients": c.tolist(),
        "physical_relative_L2_error": physical,
        "getfem_raw_target_L2_squared": ref_sq,
        "rank_coefficient_rel_error": rank_error,
    }


def gaussian_case(mesh, region, mpi_comm):
    fem = gf.MeshFem(mesh, 1)
    fem.set_fem(gf.Fem("FEM_PK(2,1)"))
    mim = gf.MeshIm(mesh, gf.Integ("IM_TRIANGLE(13)"))
    model = gf.Model("real")
    model.add_filtered_fem_variable("u", fem, region)
    obj = Adapter(model, fem, region)
    DPHS.set_initial_value(obj, "u", GAUSSIAN, projection_mim=mim)
    c = np.asarray(model.variable("u"), dtype=np.float64).reshape(-1)
    source = python_expression_to_gwfl(GAUSSIAN, 1)
    norm = integral(mim, region, model, f"({source})*({source})")
    diff = integral(
        mim, region, model, f"(u-({source}))*(u-({source}))"
    )
    if not np.isfinite(norm) or norm <= 0 or not np.isfinite(diff):
        raise RuntimeError("Invalid Gaussian GetFEM physical quadrature")
    rel_physical = float(np.sqrt(max(0., diff)/norm))
    others = mpi_comm.allgather(c)
    rel_rank = max(
        float(np.linalg.norm(other-c)/max(np.linalg.norm(c),1e-30))
        for other in others
    )
    if not obj.initial_value_set["u"] or rel_rank >= 1e-10:
        raise RuntimeError("MPI rank mismatch for Gaussian")
    print(f"MPI{mpi_comm.size} R{mpi_comm.rank} Gaussian "
          f"DOFs={c.size} physical_L2={rel_physical:.4g} "
          f"rank_error={rel_rank:.4g}",flush=True)
    return {
        "status": "PASS_NATIVE_MPI_GAUSSIAN",
        "expression": GAUSSIAN,
        "global_coefficients": c.tolist(),
        "physical_relative_L2_error": rel_physical,
        "getfem_raw_target_L2_squared": norm,
        "rank_coefficient_rel_error": rel_rank,
    }


def lumped_case():
    adapter = LumpedAdapter()
    DPHS.set_initial_value(adapter, "lumped", "2.5")
    if (not adapter.initial_value_set["lumped"]
            or adapter.gf_model.value.shape != (1,)
            or float(adapter.gf_model.value[0]) != 2.5):
        raise RuntimeError("Numeric lumped initialization regression")
    return {
        "status": "PASS_DEFAULT_LUMPED_NUMERIC",
        "value": 2.5,
    }


def main():
    arg = argparse.ArgumentParser()
    arg.add_argument("--outdir", type=Path, required=True)
    a = arg.parse_args()
    comm = MPI.COMM_WORLD
    if (comm.size not in (1, 2) or
            comm.size != PETSc.COMM_WORLD.getSize()):
        raise RuntimeError("Explicit MPI1/2 native gate only")
    if int(gf.util_mpi_parallelism_level()) != 2:
        raise RuntimeError("GetFEM must be built with MPI level 2")
    if not a.outdir.is_dir():
        raise RuntimeError("Expected writable output mount")

    factor = float(_getfem_mpi_rhs_scale())
    if not np.isfinite(factor) or factor <= 0:
        raise RuntimeError("Bad GetFEM MPI RHS calibration")
    mesh = gf.Mesh("regular simplices",
                   np.linspace(0.,1.,5),np.linspace(0.,1.,4))
    region = 52
    mesh.set_region(region,np.asarray(mesh.cvid()).reshape(-1)[::2])
    rows = {}
    for name in SELECTED:
        try:
            rows[name] = fe_case(mesh,region,name,CASES[name],comm)
        except Exception as exc:
            rows[name] = {
                "status": "FAIL",
                "error": f"{type(exc).__name__}: {exc}",
                "traceback": traceback.format_exc(),
            }
            print(f"MPI{comm.size} R{comm.rank} FAIL {name}: {exc}",
                  flush=True)
            print(traceback.format_exc(),flush=True)
    try:
        rows["gaussian_cg1_scalar"] = gaussian_case(mesh,region,comm)
    except Exception as exc:
        rows["gaussian_cg1_scalar"] = {
            "status": "FAIL",
            "error": f"{type(exc).__name__}: {exc}",
            "traceback":traceback.format_exc(),
        }
        print(traceback.format_exc(),flush=True)
    try:
        rows["lumped"] = lumped_case()
    except Exception as exc:
        rows["lumped"] = {
            "status":"FAIL", "error":str(exc),
            "traceback":traceback.format_exc(),
        }
    pass_all = (
        len(rows) == len(SELECTED)+2 and
        all(r["status"].startswith("PASS_") for r in rows.values())
    )
    result = {
        "schema":"scrimp-v130-default-fe-l2-mpi-native-v1",
        "status":"PASS_RANK" if pass_all else "FAIL_RANK",
        "mpi_ranks":comm.size,
        "rank":comm.rank,
        "getfem_rhs_multiplicative_scale":factor,
        "getfem_matrix_mode":"replicated",
        "cases":rows,
        "MPI2_qualified":False,
        "DAE_and_full_physics_qualified":False,
        "published_sources_unchanged":True,
    }
    path=a.outdir/f"rank-{comm.rank}.json"
    if path.exists():
        raise RuntimeError("Refuse to overwrite native rank evidence")
    path.write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    if not pass_all:
        raise RuntimeError(f"MPI{comm.size} rank {comm.rank} FE gate failed")
    print(f"PASS DEFAULT SPARSE FE INITIALIZATION MPI{comm.size} "
          f"RANK {comm.rank}",flush=True)


if __name__=="__main__":
    main()
