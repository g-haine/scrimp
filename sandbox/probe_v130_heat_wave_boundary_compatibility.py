# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Native Heat-Wave boundary-compatible initial-data experiment, MPI1.

The original published Heat-Wave model is imported UNCHANGED.
Intercept its first DPHS.init_step, replace ONLY the projected wave
velocity p coefficients to enforce exact zero trace on outer region 20,
then invoke the ACTUAL strict algebraic initializer with unchanged
tolerances. Do not execute TS, export graphics, or relax differential
initial data inside the initializer.

Two explicitly different candidate corrections:
  clamp: set outer-boundary CG1 DOFs to 0, keep all interior values;
  constrained_l2: closest zero-boundary CG1 field in the actual
                  GetFEM physical L2 mass norm.

The constrained minimizer solves M_FF delta_F = M_FB p_B via the
existing sparse GetFEM -> PETSc COO bridge, where p_B are boundary
DOFs of the *unconstrained* L2-projected field. No dense matrix solve.
This test is NOT a new automatic global DAE policy.
"""

import argparse
import json
import math
import traceback
from pathlib import Path

import getfem as gf
from mpi4py import MPI
import numpy as np
from petsc4py import PETSc

from examples.heat_wave import heat_wave_eq
from scrimp.dphs import DPHS
from scrimp.utils.dae_initialization import InconsistentDAEInitialData
from scrimp.utils.fe_initialization import real_getfem_projection_coefficients
from scrimp.utils.linalg import convert_gmm_to_petsc


class DAEProbeComplete(Exception):
    """Stop original top-level example AFTER real initialization, before TS."""


def wave_boundary_partition(fem, wave_region=2, outer_region=20, *, field_size):
    """GetFEM filtered coefficient ordering, never geometric-node guessing."""
    domain=np.asarray(fem.basic_dof_on_region(wave_region),
                      dtype=np.int64).reshape(-1)
    boundary=np.asarray(fem.basic_dof_on_region(outer_region),
                        dtype=np.int64).reshape(-1)
    if domain.size!=field_size or np.unique(domain).size!=domain.size:
        raise RuntimeError("Wave filtered GetFEM DOFs differ from model p coefficients")
    if not boundary.size or np.unique(boundary).size!=boundary.size:
        raise RuntimeError("Nonempty unique boundary DOFs required")
    if np.setdiff1d(boundary,domain).size:
        raise RuntimeError("Outer boundary contains DOFs outside wave region")
    # Model.add_filtered_fem_variable stores the basis in GetFEM region order.
    # Never assume a vector field has qdim interleaved nodal components.
    constrained=np.flatnonzero(np.isin(domain,boundary))
    free=np.flatnonzero(~np.isin(domain,boundary))
    if not len(constrained) or not len(free):
        raise RuntimeError("Expected nontrivial free and boundary wave DOFs")
    return constrained,free,domain,boundary


def constrained_l2_field(system, wave, boundary, free):
    """Best mass-norm correction subject to boundary coefficients = 0.

    If c0 is the unconstrained mass projection, stationarity yields
    delta_B=-c0_B, M_FF delta_F=M_FB c0_B. Exactly CG1 essential
    boundary conditions, no informal coefficient clipping.
    """
    fem=system.ports["p"].get_fem()
    mim=system.domain._int_method[0]
    model=gf.Model("real")
    model.add_filtered_fem_variable("corrected_wave",fem,2)
    model.add_linear_term(mim,"corrected_wave*Test_corrected_wave",2)
    model.assembly()
    n=wave.size
    if int(model.nbdof())!=n:
        raise RuntimeError("Wrong physical CG1 L2 mass dimension")
    M=PETSc.Mat().createAIJ(size=(n,n),comm=PETSc.COMM_WORLD)
    x_b=Mb=None
    subset=Mff=rhs=delta=ksp=None
    x0=x1=grad=None
    try:
        # The exact same production distributed COO matrix bridge.
        convert_gmm_to_petsc(model.tangent_matrix(),M,comm=PETSc.COMM_WORLD,
                             distributed_coo_eligible=True)
        x_b=M.createVecRight()
        x_b.set(0.)
        x_b.setValues(boundary.astype(PETSc.IntType),
                      wave[boundary],addv=PETSc.InsertMode.INSERT_VALUES)
        x_b.assemble()
        Mb=M.createVecLeft()
        M.mult(x_b,Mb)
        subset=PETSc.IS().createGeneral(free.astype(PETSc.IntType),
                                        comm=PETSc.COMM_WORLD)
        Mff=M.createSubMatrix(subset,subset)
        rhs=Mff.createVecLeft()
        rhs.setValues(np.arange(len(free),dtype=PETSc.IntType),
                      Mb.getValues(free.astype(PETSc.IntType)),
                      addv=PETSc.InsertMode.INSERT_VALUES)
        rhs.assemble()
        delta=Mff.createVecRight()
        ksp=PETSc.KSP().create(comm=PETSc.COMM_WORLD)
        ksp.setOptionsPrefix("scrimp_hw_boundary_l2_")
        ksp.setOperators(Mff)
        ksp.setType(PETSc.KSP.Type.PREONLY)
        ksp.getPC().setType(PETSc.PC.Type.LU)
        ksp.getPC().setFactorSolverType("mumps")
        ksp.setFromOptions()
        ksp.solve(rhs,delta)
        if int(ksp.getConvergedReason())<=0:
            raise RuntimeError("Constrained L2 interior sparse solve failed")
        free_delta=real_getfem_projection_coefficients(
            np.asarray(delta.getArray(readonly=True))
        )
        corrected=wave.copy()
        corrected[boundary]=0.
        corrected[free]+=free_delta
        if not np.all(np.isfinite(corrected)) or np.any(corrected[boundary]):
            raise RuntimeError("Constrained correction violates zero boundary trace")
        x0=M.createVecRight()
        x0.setValues(np.arange(n,dtype=PETSc.IntType),
                     corrected-wave,addv=PETSc.InsertMode.INSERT_VALUES)
        x0.assemble()
        x1=M.createVecLeft()
        M.mult(x0,x1)
        gradient=real_getfem_projection_coefficients(
            np.asarray(x1.getArray(readonly=True))
        )
        stationary=float(np.linalg.norm(gradient[free])/
                         max(float(np.linalg.norm(gradient)),1e-30))
        physical_sq=float(np.vdot(x0.getArray(readonly=True),
                                 x1.getArray(readonly=True)).real)
        if not np.isfinite(physical_sq) or physical_sq < -1e-12:
            raise RuntimeError("Nonphysical mass norm of constrained L2 correction")
        if stationary>1e-8:
            raise RuntimeError(
                f"Constrained L2 free-DOF stationarity failed: {stationary}"
            )
        return corrected,{
            "mass_norm_correction":math.sqrt(max(0.,physical_sq)),
            "relative_stationarity_on_free_dofs":stationary,
            "interior_mass_solve_KSP_reason":int(ksp.getConvergedReason()),
            "interior_mass_solve_KSP_iterations":int(ksp.getIterationNumber()),
        }
    finally:
        for obj in (x1,x0,ksp,delta,rhs,Mff,subset,Mb,x_b,M):
            if obj is not None:
                obj.destroy()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate",
                        choices=("clamp","constrained_l2"),required=True)
    parser.add_argument("--out",type=Path,required=True)
    args=parser.parse_args()
    if MPI.COMM_WORLD.size!=1 or PETSc.COMM_WORLD.getSize()!=1:
        raise RuntimeError("First candidate physical DAE experiment: MPI1 only")
    if args.out.exists() or not args.out.parent.is_dir():
        raise RuntimeError("New writable output evidence file required")

    recorded={}
    original=DPHS.init_step

    def attempt_compatible_state(system):
        original_T=np.asarray(system.gf_model.variable("T")).reshape(-1).copy()
        original_q=np.asarray(system.gf_model.variable("q")).reshape(-1).copy()
        original_p=np.asarray(system.gf_model.variable("p"),
                              dtype=np.float64).reshape(-1).copy()
        fem=system.ports["p"].get_fem()
        B,F,region_dofs,boundary_dofs=wave_boundary_partition(
            fem,field_size=original_p.size
        )
        boundary_before=float(np.linalg.norm(original_p[B]))
        if not boundary_before>1e-6:
            raise RuntimeError("Wave CG1 exterior Gaussian was unexpectedly zero")
        metrics={}
        if args.candidate=="clamp":
            corrected=original_p.copy()
            corrected[B]=0.
        else:
            corrected,metrics=constrained_l2_field(
                system,original_p,B,F
            )
        if (np.any(corrected[B]) or
                not np.array_equal(corrected[~np.isin(
                    np.arange(original_p.size),B
                )],original_p[F]) and args.candidate=="clamp"):
            raise RuntimeError("Boundary clamp incorrectly changed interior DOFs")
        system.gf_model.set_variable("p",np.ascontiguousarray(
            corrected,dtype=np.float64
        ))
        if not np.array_equal(
            np.asarray(system.gf_model.variable("T")).reshape(-1),original_T
        ) or not np.array_equal(
            np.asarray(system.gf_model.variable("q")).reshape(-1),original_q
        ):
            raise RuntimeError("Boundary correction touched heat or stress states")
        recorded["candidate"]=args.candidate
        recorded["wave_region_DOFs"]=int(original_p.size)
        recorded["boundary_global_meshfem_DOFs"]=boundary_dofs.tolist()
        recorded["boundary_filtered_p_DOFs"]=B.tolist()
        recorded["free_filtered_p_DOFs"]=int(len(F))
        recorded["before_boundary_Euclidean_coeff_norm"]=boundary_before
        recorded["after_boundary_Euclidean_coeff_norm"]=float(
            np.linalg.norm(corrected[B])
        )
        recorded["wave_coefficient_relative_change"]=float(
            np.linalg.norm(corrected-original_p)/
            max(float(np.linalg.norm(original_p)),1e-30)
        )
        recorded["T_q_unchanged"]=True
        recorded["field_adjustment_metrics"]=metrics
        try:
            original(system)
        except InconsistentDAEInitialData as exc:
            recorded["status"]="UNRESOLVED_COMPATIBILITY"
            recorded["infeasibility_error"]=str(exc)
            recorded["native_initializer_report"]=system.initialization_report
        else:
            report=system.initialization_report
            if report["status"]!="PASS":
                raise RuntimeError("Native initializer returned without PASS")
            if not report["differential_initial_exact"] or report["pseudo_ts_executed"]:
                raise RuntimeError("Native initializer changed dynamic initial data")
            if not np.array_equal(
                np.asarray(system.gf_model.variable("p")).reshape(-1),corrected
            ):
                raise RuntimeError("Algebraic initializer changed corrected p")
            if not np.array_equal(
                np.asarray(system.gf_model.variable("T")).reshape(-1),original_T
            ) or not np.array_equal(
                np.asarray(system.gf_model.variable("q")).reshape(-1),original_q
            ):
                raise RuntimeError("Native DAE initializer modified T or q")
            recorded["status"]="PASS_DAE_COMPATIBLE_BOUNDARY_CORRECTION"
            recorded["native_initializer_report"]=report
        # This is a DAE-initial-state proof; intentionally do not run
        # CN TS or public postprocessing with altered physical initial data.
        raise DAEProbeComplete

    DPHS.init_step=attempt_compatible_state
    try:
        try:
            heat_wave_eq(
                heat_region=1,wave_region=2,
                final_time=.004,time_step=.001,save_step=.001,
                factor_solver="mumps",
                initial_wave_boundary="unconstrained",
            )
        except DAEProbeComplete:
            pass
        else:
            raise RuntimeError("Native compatibility probe did not intercept DAE")
    finally:
        DPHS.init_step=original

    if "status" not in recorded:
        raise RuntimeError("No true native DAE compatibility result")
    recorded.update({
        "schema":"scrimp-v130-original-heat-wave-boundary-compatibility-v1",
        "published_model_source_unchanged":True,
        "original_analytic_initial_gaussian_unchanged":True,
        "default_L2_projection_unchanged":True,
        "actual_PETSc_DAE_initializer_executed":True,
        "no_pseudo_time_or_tolerance_relaxation":True,
        "full_CN_TS_executed":False,
        "energy_power_qualified":False,
        "MPI2_qualified":False,
    })
    args.out.write_text(json.dumps(recorded,indent=2)+"\n",
                        encoding="utf-8")
    print(
        f"BOUNDARY CANDIDATE {args.candidate}: {recorded['status']}, "
        f"p_B={recorded['before_boundary_Euclidean_coeff_norm']:.6g}->"
        f"{recorded['after_boundary_Euclidean_coeff_norm']:.6g}, "
        f"relative_p_change={recorded['wave_coefficient_relative_change']:.6g}",
        flush=True
    )
    report=recorded["native_initializer_report"]
    print(f"  original IFunction algebraic norm="
          f"{report['native_final_residual_l2']:.9g}, "
          f"tolerance={report['feasibility_tolerance']:.9g}",flush=True)
    if recorded["status"]=="UNRESOLVED_COMPATIBILITY":
        for block in report["dominant_algebraic_variables"]:
            print(f"  unresolved {block['variable']}: "
                  f"{block['residual_l2']:.9g}",flush=True)


if __name__=="__main__":
    main()
