# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Native STRICT DAE initializer for a corrected two-sided Heat-Wave interface.

First run original published Heat-Wave with existing doubled region 10
and boundary-clamped p0 as a CONTROL. Compare with two experimental
one-sided facet runs, clamp and mass-minimal constrained_l2. All three
invoke the real, unrelaxed PETSc algebraic LSQR initializer on freshly
assembled production GetFEM/PETSc matrices.

The source example is imported unchanged. A reversible opt-in runtime
constructor wrapper touches only the original six interface bricks and
the two interface control ports; it constructs 30/31 using the original
P1 parent triangles. It stops IMMEDIATELY after DAE initializer, even
if algebraic compatibility is rejected: no TS, no CN, no Hamiltonian
and no physics qualification. Rejecting incompatible differential
states is a scientifically useful result, NOT a reason to relax DAE.
"""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.sparse import csr_matrix
from mpi4py import MPI
from petsc4py import PETSc

import scrimp as S
from scrimp.dphs import DPHS
from scrimp.utils.dae_initialization import InconsistentDAEInitialData
from examples.heat_wave import heat_wave_eq

from v130_hw_interface_geometry import new_one_sided_regions
from probe_v130_heat_wave_boundary_compatibility import (
    wave_boundary_partition, constrained_l2_field,
)


HEAT_BRICKS = frozenset(("M_Y_T", "B_T", "B_T^T"))
WAVE_BRICKS = frozenset(("M_Y_w", "B_w", "B_w^T"))
INTERFACE_CONTROLS={"Interface Heat":30,"Interface Wave":31}
BLOCKS={
    "heat_normal_flux_trace":("Y_T","e_Q"),
    "heat_control_to_normal_flux":("e_Q","U_T"),
    "wave_velocity_trace":("Y_w","p"),
    "wave_control_into_p":("p","U_w"),
    "wave_output_to_heat_control":("U_T","Y_w"),
    "heat_output_to_wave_control":("U_w","Y_T"),
}


class DAEOnlyComplete(Exception):
    """Stop before original CN TS regardless of DAE feasibility."""


def inspect_native_stiffness(system, variant):
    if system.gf_model.nbdof()!=7689:
        raise RuntimeError("Original 7689 Heat-Wave DOFs changed")
    ia,ja,val=system.stiffness.getValuesCSR()
    K=csr_matrix(
        (np.asarray(val).copy(),np.asarray(ja).copy(),
         np.asarray(ia).copy()),shape=(7689,7689)
    )
    spans={}
    for name in ("T","p","q","e_Q","U_T","Y_T",
                 "U_w","Y_w","U_bnd","Y_bnd"):
        start,count=map(int,system.gf_model.interval_of_variable(name))
        spans[name]=(start,count)
    result={}
    for label,(row,col) in BLOCKS.items():
        ri,rn=spans[row]
        ci,cn=spans[col]
        blk=K[ri:ri+rn,ci:ci+cn]
        result[label]={
            "Frobenius_norm":float(np.linalg.norm(blk.data)),
            "numeric_nnz":int(np.count_nonzero(blk.data)),
        }
    heat_norms=(result["heat_normal_flux_trace"]["Frobenius_norm"],
                result["heat_control_to_normal_flux"]["Frobenius_norm"])
    if variant=="baseline":
        if any(x>1e-9 for x in heat_norms):
            raise RuntimeError(
                "Original baseline no longer reproduces doubled-normal cancellation"
            )
    else:
        if any(x<=1e-4 for x in heat_norms):
            raise RuntimeError("Dual-sided heat flux coupling missing")
        if any(result[name]["Frobenius_norm"]<=1e-4
               for name in BLOCKS):
            raise RuntimeError("Dual-sided interface feedback disappeared")
    return {
        "DOFs":int(K.shape[0]),
        "stored_nnz":int(K.nnz),
        "numeric_nnz":int(np.count_nonzero(K.data)),
        "blocks":result,
        "assembly_gate":"PASS_EXPECTED_BASELINE_OR_CORRECTED_BLOCKS",
    }


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--variant",choices=("baseline","dual_sided"),
                   required=True)
    p.add_argument("--candidate",choices=("clamp","constrained_l2"),
                   required=True)
    p.add_argument("--out",type=Path,required=True)
    args=p.parse_args()

    if MPI.COMM_WORLD.size!=1 or PETSc.COMM_WORLD.getSize()!=1:
        raise RuntimeError("Strict interface DAE diagnostic MPI1 only")
    if args.out.exists() or not args.out.parent.is_dir():
        raise RuntimeError("Native DAE report path must not exist")

    original_domain=S.Domain
    original_control_port=S.Control_Port
    original_brick=S.Brick
    original_init=DPHS.init_step
    callback_count=0
    topology=[]
    modified={"ports":{},"bricks":{}}
    result={}

    def construct_domain(*pos,**named):
        domain=original_domain(*pos,**named)
        mesh=domain.get_mesh()[0]
        geo=new_one_sided_regions(mesh)
        interface=np.asarray(mesh.region(10),dtype=np.int64)
        heat=np.asarray(mesh.region(30),dtype=np.int64)
        wave=np.asarray(mesh.region(31),dtype=np.int64)
        as_set=lambda arr:set(map(tuple,arr.T.tolist()))
        if (interface.shape!=(2,156)
                or heat.shape!=(2,78) or wave.shape!=(2,78)
                or as_set(interface)!=(as_set(heat)|as_set(wave))
                or (as_set(heat)&as_set(wave))):
            raise RuntimeError("Actual original 78+78 facet topology changed")
        topology.append(geo)
        return domain

    def construct_control(*pos,**named):
        if pos and pos[0] in INTERFACE_CONTROLS and args.variant=="dual_sided":
            name=pos[0]
            if named.get("region")!=10:
                raise RuntimeError("Unexpected original interface port region")
            named["region"]=INTERFACE_CONTROLS[name]
            modified["ports"][name]=named["region"]
        return original_control_port(*pos,**named)

    def construct_brick(*pos,**named):
        if pos and pos[0] in HEAT_BRICKS|WAVE_BRICKS and args.variant=="dual_sided":
            name=pos[0]
            side=30 if name in HEAT_BRICKS else 31
            if len(pos)<3 or list(pos[2])!=[10]:
                raise RuntimeError("Unexpected original brick region: "+name)
            pos=list(pos)
            pos[2]=[side]
            pos=tuple(pos)
            modified["bricks"][name]=side
        return original_brick(*pos,**named)

    def native_DAE_gate(system):
        nonlocal callback_count
        callback_count+=1
        if callback_count!=1 or len(topology)!=1:
            raise RuntimeError("Original GetFEM assembly interception failed")
        if args.variant=="dual_sided" and (
            set(modified["ports"])!=set(INTERFACE_CONTROLS)
            or set(modified["bricks"])!=(HEAT_BRICKS|WAVE_BRICKS)):
            raise RuntimeError("Original interface constructor interception incomplete")
        assembly=inspect_native_stiffness(system,args.variant)
        old_state={
            name:np.asarray(system.gf_model.variable(name),
                            dtype=np.float64).reshape(-1).copy()
            for name in ("T","p","q")
        }
        p_orig=old_state["p"]
        b,free,_,_=wave_boundary_partition(
            system.ports["p"].get_fem(),field_size=len(p_orig)
        )
        if len(b)!=126:
            raise RuntimeError("Original 126 exterior p DOFs changed")
        metrics={}
        if args.candidate=="clamp":
            p_corrected=p_orig.copy()
            p_corrected[b]=0.
            if not np.array_equal(p_corrected[free],p_orig[free]):
                raise RuntimeError("Boundary clamp changed interior wave p")
        else:
            p_corrected,metrics=constrained_l2_field(
                system,p_orig,b,free
            )
        if (not np.all(np.isfinite(p_corrected))
                or np.any(p_corrected[b])):
            raise RuntimeError("Corrected p is not exterior-boundary compatible")
        system.gf_model.set_variable(
            "p",np.ascontiguousarray(p_corrected,dtype=np.float64)
        )
        result.update({
            "schema":"scrimp-v130-HW-dual-interface-strict-DAE-v1",
            "status":"DAE_NOT_YET_EXECUTED",
            "variant":args.variant,
            "candidate":args.candidate,
            "assembly":assembly,
            "topology":topology[0],
            "remapped_original_constructors":modified,
            "unmodified_initial_Gaussian_and_L2_projection":True,
            "original_p_boundary_DOFs":int(len(b)),
            "p_boundary_before_l2":float(np.linalg.norm(p_orig[b])),
            "p_boundary_after_l2":float(np.linalg.norm(p_corrected[b])),
            "relative_p_correction":float(
                np.linalg.norm(p_corrected-p_orig)/
                max(float(np.linalg.norm(p_orig)),1.e-30)
            ),
            "constrained_l2_metrics":metrics,
            "strict_initializer_called":True,
            "no_pseudo_TS_or_tolerance_relaxation":True,
            "original_source_and_weak_form_text_unchanged":True,
            "MPI":1,
            "main_CN_TS_executed":False,
            "physical_power_balance_qualified":False,
        })
        incompatible=None
        try:
            original_init(system)
        except InconsistentDAEInitialData as exc:
            incompatible=str(exc)
        native=system.initialization_report
        if native is None:
            raise RuntimeError(
                "Real native PETSc DAE report unavailable after strict initialization"
            )
        if (not native["differential_initial_exact"]
                or native["pseudo_ts_executed"]):
            raise RuntimeError("Strict algebraic initialization changed dynamic data")
        for name,initial in (
            ("T",old_state["T"]),("q",old_state["q"]),
            ("p",p_corrected)
        ):
            actual=np.asarray(system.gf_model.variable(name)).reshape(-1)
            if not np.array_equal(actual,initial):
                raise RuntimeError(
                    "Actual DAE initializer changed differential "+name
                )
        result["native_initialization"]=native
        result["incompatibility_exception"]=incompatible
        norm=float(native["native_final_residual_l2"])
        threshold=float(native["feasibility_tolerance"])
        if incompatible is not None:
            if native["status"]!="UNRESOLVED" or norm<=threshold:
                raise RuntimeError(
                    "Inconsistent DAE exception but recorded algebraic gate PASS"
                )
            result["status"]="INCOMPATIBLE_INITIAL_DAE_UNALTERED_PHYSICS"
        else:
            if (native["status"]!="PASS" or norm>threshold
                    or not native["accuracy_target_met"]):
                raise RuntimeError(
                    "Native DAE returned without successful strict feasibility"
                )
            result["status"]="PASS_NATIVE_STRICT_DAE_FEASIBILITY_ONLY"
        args.out.write_text(
            json.dumps(result,indent=2)+"\n",encoding="utf-8"
        )
        print(
            f"NATIVE STRICT DAE {args.variant}/{args.candidate}: "
            f"{result['status']} residual={norm:.10g} "
            f"tol={threshold:.10g} "
            f"p_boundary={result['p_boundary_before_l2']:.6g}->0 "
            f"prescribed_T_p_q_exact=True",flush=True
        )
        if incompatible is not None:
            for x in native.get("dominant_algebraic_variables",[]):
                print(
                    f"  INCOMPATIBLE {x['variable']} "
                    f"residual={x['residual_l2']:.10g}",flush=True
                )
        raise DAEOnlyComplete

    S.Domain=construct_domain
    S.Control_Port=construct_control
    S.Brick=construct_brick
    DPHS.init_step=native_DAE_gate
    try:
        try:
            heat_wave_eq(
                heat_region=1,wave_region=2,
                final_time=.004,time_step=.001,save_step=.001,
                factor_solver="mumps",
            )
        except DAEOnlyComplete:
            pass
        else:
            raise RuntimeError(
                "Original Heat-Wave passed DAE hook and unexpectedly ran TS"
            )
    finally:
        S.Domain=original_domain
        S.Control_Port=original_control_port
        S.Brick=original_brick
        DPHS.init_step=original_init
    if callback_count!=1 or not args.out.is_file():
        raise RuntimeError("Original assembled strict DAE not audited exactly once")


if __name__=="__main__":
    main()
