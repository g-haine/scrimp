# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Diagnose GetFEM's 1D region-10 Normal against true one-sided 2D faces.

Only one real *unchanged* MUMPS original Heat-Wave trajectory (clamp,
strict DAE, default L2, COO, CN four steps). The original physics
collector is called first without modification. Then separate new
mesh regions 30 (faces of original disk 2D triangles) and 31
(inner faces of original annular 2D triangles) are added ONLY after
the solve, for independent FE postprocessing.

GetFEM's documented 'gmsh_with_lower_dim_elt' can retain lower-
dimensional convexes. The original region 10 may not encode a one-sided
normal of either adjoining 2D subdomain. Check actual region topology,
normals and integral fluxes BEFORE claiming a physical mesh correction.

Any failure of optional cross-trace expressions is preserved verbatim,
never silently interpreted as zero.
"""

import argparse
import json
from pathlib import Path
import sys

import numpy as np

import probe_v130_main_2d_examples as original
import probe_v130_heat_wave_compatible_cn as cn


FACE_HEAT = 30
FACE_WAVE = 31
REGION_ORIGINAL = 10
QUANTITIES = {
    "heat_normal_flux_sq":"(e_Q.Normal)*(e_Q.Normal)",
    "wave_normal_stress_sq":"(q.Normal)*(q.Normal)",
    "heat_flux_normal":"e_Q.Normal",
    "wave_stress_normal":"q.Normal",
    "heat_interface_flux_power":"U_T*(e_Q.Normal)",
    "wave_interface_stress_power":"Y_w*(q.Normal)",
    "heat_trace_flux_observation":"Y_T*(e_Q.Normal)",
    "heat_trace_temperature":"T*T",
    "wave_trace_velocity":"p*p",
    "effort_heat_port":"U_T*U_T",
    "effort_wave_port":"U_w*U_w",
    "output_heat_port":"Y_T*Y_T",
    "output_wave_port":"Y_w*Y_w",
}


def region_description(mesh, region):
    ids=np.asarray(mesh.region(region),dtype=np.int64)
    if ids.ndim!=2 or ids.shape[0]!=2:
        raise RuntimeError(f"Unexpected GetFEM region format {region}")
    return {
        "region":region,
        "convex_face_pairs":int(ids.shape[1]),
        "whole_convexes":int(np.count_nonzero(ids[1,:]<0)),
        "parent_convex_faces":int(np.count_nonzero(ids[1,:]>=0)),
        "first_convex_faces":ids[:,:5].tolist(),
    }


def new_one_sided_regions(mesh):
    bulk={}
    boundaries={}
    for number in (1,2):
        ids=np.asarray(mesh.region(number),dtype=np.int64)
        if ids.ndim!=2 or ids.shape[0]!=2:
            raise RuntimeError("GetFEM volume region does not use (cv,face) IDs")
        volume=np.asarray(ids[0,:],dtype=np.int64)
        if not len(volume) or np.any(ids[1,:]>=0):
            raise RuntimeError("Bulk volume must contain complete convexes")
        bulk[number]=volume
        faces=np.asarray(mesh.outer_faces(2,volume),dtype=np.int64)
        if faces.ndim!=2 or faces.shape[0]!=2:
            raise RuntimeError("Invalid GetFEM outer-face array")
        boundaries[number]=faces

    # Disk has only one boundary at the physical heat-wave interface.
    disk=boundaries[1]
    wave=boundaries[2]
    if disk.shape[1]<10 or wave.shape[1]<disk.shape[1]:
        raise RuntimeError("Unexpected Concentric bulk facet topology")
    if not np.all(np.isin(disk[0,:],bulk[1])):
        raise RuntimeError("A disk boundary face has the wrong parent")
    if not np.all(np.isin(wave[0,:],bulk[2])):
        raise RuntimeError("An annulus boundary face has the wrong parent")
    # Partition the annulus boundary faces by radial coordinate. 0.8 is
    # safely between the concentric radii 0.6 and 1.0, not a tolerance.
    annulus_inner=[]
    annulus_outer=[]
    for j in range(wave.shape[1]):
        face=wave[:,j:j+1]
        pids=np.asarray(mesh.pid_in_faces(face),dtype=np.int64)
        xyz=np.asarray(mesh.pts(pids),dtype=np.float64)
        if xyz.shape[0]<2 or xyz.shape[1]<2:
            raise RuntimeError("Expected 2D annulus face endpoints")
        r=float(np.linalg.norm(np.mean(xyz[:2,:],axis=1)))
        if r<0.8:annulus_inner.append(j)
        else:annulus_outer.append(j)
    if len(annulus_inner)!=disk.shape[1] or not annulus_outer:
        raise RuntimeError(
            "Concentric interface facet count differs across 2D subdomains: "
            f"disk={disk.shape[1]} annulus={len(annulus_inner)}, "
            f"exterior={len(annulus_outer)}"
        )
    inner=wave[:,annulus_inner]
    if any(int(r) in np.asarray(mesh.regions()).reshape(-1)
           for r in (FACE_HEAT,FACE_WAVE)):
        raise RuntimeError("Audit-only region 30 or 31 unexpectedly exists")
    mesh.set_region(FACE_HEAT, disk)
    mesh.set_region(FACE_WAVE, inner)

    normals={}
    for name,face in (("disk",disk),("annulus",inner)):
        vectors=np.asarray(mesh.normal_of_faces(face),dtype=float)
        if vectors.shape!=(2,face.shape[1]):
            raise RuntimeError(
                f"Unexpected {name} outward 2D face normal shape: {vectors.shape}"
            )
        normal_lens=np.linalg.norm(vectors,axis=0)
        if np.any(~np.isfinite(vectors)) or not np.allclose(
            normal_lens,1.,rtol=1e-3,atol=1e-3
        ):
            raise RuntimeError("2D parent mesh faces lack valid normals")
        normals[name]={
            "count":int(face.shape[1]),
            "unit_norm_min":float(np.min(normal_lens)),
            "unit_norm_max":float(np.max(normal_lens)),
            "mean_components":np.mean(vectors,axis=1).tolist(),
        }
    return {
        "original_1d_region":region_description(mesh,REGION_ORIGINAL),
        "disk_parent_region":region_description(mesh,1),
        "annulus_parent_region":region_description(mesh,2),
        "disk_interface_2d_face_region":region_description(mesh,FACE_HEAT),
        "annulus_interface_2d_face_region":region_description(mesh,FACE_WAVE),
        "outer_annulus_face_count":len(annulus_outer),
        "parent_face_normals":normals,
        "no_volume_region_modified":True,
    }


def measure_regions(system,times):
    mesh=system.domain.get_mesh()[0]
    topology=new_one_sided_regions(mesh)
    integrals={}
    failures={}
    for label,region in (
        ("original_region10",REGION_ORIGINAL),
        ("heat_faces_region30",FACE_HEAT),
        ("wave_faces_region31",FACE_WAVE),
    ):
        integrals[label]={}
        for key,expression in QUANTITIES.items():
            try:
                a=np.asarray(system.get_quantity(
                    expression,region=region,CN=True
                ),dtype=float)
                if a.shape!=np.asarray(times).shape or not np.all(np.isfinite(a)):
                    raise ValueError("Invalid GetFEM CN FE quantity")
                integrals[label][key]=a.tolist()
            except Exception as exc:
                failures[f"{label}:{key}"] = (
                    f"{type(exc).__name__}: {str(exc)[:500]}"
                )
                integrals[label][key]=None

    return topology,integrals,failures


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    if args.out.exists() and any(args.out.iterdir()):
        raise RuntimeError("Refuse to reuse native evidence")
    original_collector=original.collect_heat_wave_physics
    saved_argv=sys.argv
    called=0
    def observer(system,times,H,*a,**kw):
        nonlocal called
        called+=1
        traces,physics=original_collector(system,times,H,*a,**kw)
        topology, integrals,errors=measure_regions(system,times)
        report={
            "schema":"scrimp-v130-HW-one-sided-2D-facet-normal-audit-v1",
            "status":"DIAGNOSTIC_ONLY_NATIVE_FACET_FLUX_OBSERVATION",
            "initial_condition":"clamp",
            "original_source_and_weak_forms_unmodified":True,
            "original_physical_energy_collector_called_first":True,
            "true_CN_times":np.asarray(times).tolist(),
            "topology":topology,
            "integrals_CN":integrals,
            "unavailable_cross_FEM_traces":errors,
            "diagnostic_regions_added_after_actual_TS_solve":True,
            "original_interface_coupling_matrix_not_reassembled":True,
            "physical_interface_coupling_qualified":False,
            "main_real_CN_steps":4,
            "MPI":1,
        }
        out=args.out/"one_sided_facet_normal_diagnostic.json"
        if out.exists():
            raise RuntimeError("Refuse to overwrite native normal evidence")
        out.write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
        print("ORIGINAL 1D REGION "+json.dumps(topology["original_1d_region"]),
              flush=True)
        print("PARENT 2D FACE REGIONS: "+json.dumps({
          "disk":topology["disk_interface_2d_face_region"],
          "annulus":topology["annulus_interface_2d_face_region"],
        }),flush=True)
        for name,obs in integrals.items():
            print(f"FACET TRACE {name}: "+json.dumps({
              k:obs[k] for k in (
                "heat_normal_flux_sq","wave_normal_stress_sq",
                "heat_interface_flux_power",
                "heat_trace_flux_observation",
                "effort_heat_port","output_heat_port"
              )
            }),flush=True)
        if errors:
            print("OPTIONAL FACET TRACE ERRORS: "+
                  json.dumps(errors),flush=True)
        return traces,physics

    original.collect_heat_wave_physics=observer
    try:
        sys.argv=[
            str(Path(__file__).name),
            "--candidate","clamp","--backend","mumps",
            "--out",str(args.out),
        ]
        cn.main()
    finally:
        original.collect_heat_wave_physics=original_collector
        sys.argv=saved_argv
    if called!=1:
        raise RuntimeError("Real physical CN collector not reached")
    print("PASS_CAPTURED_TRUE_HEAT_WAVE_FACET_NORMAL_DIAGNOSTIC",
          flush=True)


if __name__=="__main__":
    main()
