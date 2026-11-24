# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Sandbox-only original Heat-Wave double-sided interface assembly contrast.

Compare EXACTLY the unchanged original example assembled stiffness on
(1) its existing two-sided region 10, with (2) disk-side region 30
for HEAT forms/ports and annulus-side region 31 for WAVE forms/ports.

All changes are instance-local to a reversible Python runtime wrapper
on scrimp's exported Domain/Control_Port/Brick constructor functions.
The saved original source files are NEVER edited. There is NO main TS,
NO DAE initialization, and NO physics PASS in this first stage.

This specifically qualifies whether the two normal-flux stiffness
blocks reappear AND whether the feedback and wave blocks are intact,
before attempting constrained initial data or time integration.
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from scipy.sparse import csr_matrix
from petsc4py import PETSc

import scrimp as S
from scrimp.dphs import DPHS
from examples.heat_wave import heat_wave_eq
from v130_hw_interface_geometry import new_one_sided_regions


HEAT_BRICKS = frozenset(("M_Y_T", "B_T", "B_T^T"))
WAVE_BRICKS = frozenset(("M_Y_w", "B_w", "B_w^T"))
PORT_REGIONS = {"Interface Heat": 30, "Interface Wave": 31}
BLOCKS = {
    "heat_flux_into_interface_output": ("Y_T","e_Q"),
    "heat_boundary_input_into_flux": ("e_Q","U_T"),
    "wave_velocity_into_interface_output": ("Y_w","p"),
    "wave_boundary_input_into_velocity": ("p","U_w"),
    "feedback_wave_output_into_heat_input": ("U_T","Y_w"),
    "feedback_heat_output_into_wave_input": ("U_w","Y_T"),
    "heat_flux_into_heat_bulk": ("T","e_Q"),
    "wave_velocity_into_wave_stress": ("q","p"),
}
NAMES=("T","p","q","e_Q","U_T","Y_T","U_w","Y_w","U_bnd","Y_bnd")


class AssemblyOnlyComplete(Exception):
    """Expected end-of-experiment, before DAE and PETSc TS."""


def compare_block_norms(system, variant, topology, modifications, out):
    n=int(system.gf_model.nbdof())
    stiffness=system.stiffness
    if stiffness.getSize()!=(n,n):
        raise RuntimeError("Missing actual original GetFEM/PETSc COO stiffness")
    p,j,v=stiffness.getValuesCSR()
    K=csr_matrix(
        (np.asarray(v).copy(),np.asarray(j).copy(),np.asarray(p).copy()),
        shape=(n,n)
    )
    spans={}
    covered=[]
    for name in NAMES:
        start,size=map(int,system.gf_model.interval_of_variable(name))
        if start<0 or size<1 or start+size>n:
            raise RuntimeError("Invalid native interval for "+name)
        spans[name]={"start":start,"count":size}
        covered.extend(range(start,start+size))
    if sorted(covered)!=list(range(n)):
        raise RuntimeError("Native intervals omit or duplicate variables")
    blocks={}
    for label,(row,col) in BLOCKS.items():
        ri,rn=spans[row]["start"],spans[row]["count"]
        ci,cn=spans[col]["start"],spans[col]["count"]
        b=K[ri:ri+rn,ci:ci+cn]
        vals=np.asarray(b.data)
        blocks[label]={
            "row":row,"col":col,"shape":list(b.shape),
            "stored_entries":int(b.nnz),
            "numeric_nonzeros":int(np.count_nonzero(vals)),
            "Frobenius_norm":float(np.linalg.norm(vals)),
            "largest_absolute_entry":float(np.max(np.abs(vals)))
                if len(vals) else 0.,
        }
    for label in ("heat_flux_into_heat_bulk",
                  "wave_velocity_into_wave_stress"):
        if blocks[label]["Frobenius_norm"]<1e-3:
            raise RuntimeError("Native volume physics incorrectly assembled")
    result={
        "schema":"scrimp-v130-HW-original-two-sided-vs-one-sided-assembly-v1",
        "variant":variant,
        "status":"DIAGNOSTIC_ASSEMBLY_ONLY_NO_DAE_NO_TS",
        "stiffness_matrix_full_DOFS":n,
        "stiffness_stored_nnz":int(K.nnz),
        "stiffness_numeric_nnz":int(np.count_nonzero(K.data)),
        "variables":spans,
        "blocks":blocks,
        "mesh_interface_topology":topology,
        "actual_exported_scrimp_constructor_calls":modifications,
        "original_example_source_modified":False,
        "original_weak_forms_text_unchanged":True,
        "HEAT_on_disk_side_region30":variant=="dual_sided",
        "WAVE_on_annulus_side_region31":variant=="dual_sided",
        "mass_matrix_assembled":True,
        "DAE_initialized":False,
        "time_integration_executed":False,
        "physics_energy_power_qualified":False,
    }
    out.write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    print("ACTUAL ORIGINAL GETFEM+PETSC STIFFNESS ASSEMBLY "
          f"variant={variant} dofs={n} nnz={K.nnz}",flush=True)
    for label,info in blocks.items():
        print(
            f"  BLOCK {label} nnz={info['numeric_nonzeros']} "
            f"frob={info['Frobenius_norm']:.12g}",
            flush=True
        )
    print("DIAGNOSTIC_ASSEMBLY_ONLY_NO_DAE_NO_TS",flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--variant",choices=("baseline","dual_sided"),required=True)
    p.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    if PETSc.COMM_WORLD.getSize()!=1:
        raise RuntimeError("Heat-Wave corrected interface assembly MPI1 only")
    if args.out.exists():
        raise RuntimeError("Refuse to overwrite original native evidence")
    if not args.out.parent.is_dir():
        raise RuntimeError("Native archive directory does not exist")

    original_Domain=S.Domain
    original_Control_Port=S.Control_Port
    original_Brick=S.Brick
    original_init=DPHS.init_step
    topology_capture=[]
    modifications={
        "new_regions_created":[],
        "control_ports_remapped":{},
        "bricks_remapped":{},
    }
    callback_count=0

    def domain_wrapper(*a,**kw):
        domain=original_Domain(*a,**kw)
        mesh=domain.get_mesh()[0]
        geometry=new_one_sided_regions(mesh)
        disk=np.asarray(mesh.region(30),dtype=np.int64)
        wave=np.asarray(mesh.region(31),dtype=np.int64)
        old=np.asarray(mesh.region(10),dtype=np.int64)
        def cells(pairs):
            return set(map(tuple,pairs.T.tolist()))
        if (disk.shape!=(2,78) or wave.shape!=(2,78)
                or old.shape!=(2,156)
                or cells(old)!=(cells(disk)|cells(wave))
                or cells(disk)&cells(wave)):
            raise RuntimeError("Original 156-face interface != disk/wave disjoint halves")
        topology_capture.append(geometry)
        modifications["new_regions_created"]=[30,31]
        print("ORIGINAL REGION 10 = 78 DISK + 78 ANNULUS PARENT FACES",
              flush=True)
        return domain

    def control_wrapper(*a,**kw):
        if a and a[0] in PORT_REGIONS and args.variant=="dual_sided":
            name=a[0]
            if kw.get("region")!=10:
                raise RuntimeError("Original interface control region changed")
            kw["region"]=PORT_REGIONS[name]
            modifications["control_ports_remapped"][name]=kw["region"]
        return original_Control_Port(*a,**kw)

    def brick_wrapper(*a,**kw):
        if a and args.variant=="dual_sided":
            name=a[0]
            if name in HEAT_BRICKS or name in WAVE_BRICKS:
                region=30 if name in HEAT_BRICKS else 31
                if len(a)<3 or list(a[2])!=[10]:
                    raise RuntimeError(f"Unexpected original brick region {name}")
                a=list(a)
                a[2]=[region]
                a=tuple(a)
                modifications["bricks_remapped"][name]=region
        return original_Brick(*a,**kw)

    def before_native_DAE(system):
        nonlocal callback_count
        callback_count+=1
        if callback_count!=1 or len(topology_capture)!=1:
            raise RuntimeError("Wrong original Heat-Wave assembly callback count")
        if args.variant=="dual_sided" and (
            set(modifications["bricks_remapped"])!=
                HEAT_BRICKS|WAVE_BRICKS
            or set(modifications["control_ports_remapped"])!=
                set(PORT_REGIONS)
        ):
            raise RuntimeError("Did not intercept all exact original interface forms")
        compare_block_norms(
            system,args.variant,topology_capture[0],modifications,args.out
        )
        # No alteration to production initializer: stop BEFORE invoking it.
        raise AssemblyOnlyComplete

    S.Domain=domain_wrapper
    S.Control_Port=control_wrapper
    S.Brick=brick_wrapper
    DPHS.init_step=before_native_DAE
    try:
        try:
            heat_wave_eq(
                heat_region=1,wave_region=2,
                factor_solver="mumps",
                final_time=.004,time_step=.001,save_step=.001,
            )
        except AssemblyOnlyComplete:
            pass
        else:
            raise RuntimeError("Expected assembly-only stop before DAE/TS")
    finally:
        S.Domain=original_Domain
        S.Control_Port=original_Control_Port
        S.Brick=original_Brick
        DPHS.init_step=original_init
    if callback_count!=1 or not args.out.is_file():
        raise RuntimeError("Actual GetFEM/PETSc assembly was not reached")


if __name__=="__main__":
    main()
