# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Dependency-free geometric interface selection for original Concentric P1.

This module only requires NumPy, not the original Heat-Wave CN pilot,
physical collectors or GetFEM/PETSc imports. The interface is taken
from *oriented parent triangle faces* on opposite sides. The functions
are kept identical to those native-validated in the one-sided
diagnostic 28495b42; no mesh or production source is modified here.
"""

import numpy as np


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


def bulk_2d_convexes(mesh, region):
    """Select actual 2-D parents, not GetFEM region entries' face flags.

    Imported Gmsh region entries can be (convex_id, local_face_id)
    or (convex_id,-1). The second row is NOT a geometric-dimension
    indicator and must not be required to equal -1. The original
    Concentric Gmsh mesh is P1: disk/annulus are triangles (3 nodes)
    whereas independently imported boundary arcs are 2-node lines.
    This is a diagnostic for that 2-D mesh, not a general mixed-
    dimensional GetFEM mesh selector for arbitrary high-order cells.
    """
    entries=np.asarray(mesh.region(region),dtype=np.int64)
    if entries.ndim!=2 or entries.shape[0]!=2 or not entries.shape[1]:
        raise RuntimeError(f"Invalid or empty GetFEM region {region}")
    candidates=np.unique(entries[0,:])
    pids, offsets=mesh.pid_from_cvid(candidates)
    pids=np.asarray(pids,dtype=np.int64).reshape(-1)
    offsets=np.asarray(offsets,dtype=np.int64).reshape(-1)
    if (len(offsets)!=len(candidates)+1 or offsets[0]!=0
            or offsets[-1]!=len(pids)
            or np.any(np.diff(offsets)<2)):
        raise RuntimeError(
            f"Invalid GetFEM convex point connectivity for region {region}"
        )
    node_counts=np.diff(offsets)
    # This original P1 Concentric geometry has triangles and lines.
    if np.any(~np.isin(node_counts,(2,3))):
        raise RuntimeError(
            f"Unexpected non-P1 Concentric convex connectivity {region}: "
            +str(np.unique(node_counts).tolist())
        )
    volume=np.ascontiguousarray(candidates[node_counts==3],dtype=np.int64)
    if not len(volume):
        raise RuntimeError(
            f"No geometric triangles in native imported region {region}"
        )
    details={
        "region":int(region),
        "raw_getfem_region_entry_count":int(entries.shape[1]),
        "raw_whole_convex_entries":int(np.count_nonzero(entries[1,:]<0)),
        "raw_face_entries":int(np.count_nonzero(entries[1,:]>=0)),
        "unique_region_convex_ids":int(len(candidates)),
        "geometric_2D_triangle_parents":int(len(volume)),
        "excluded_1D_segment_convexes":int(np.count_nonzero(node_counts==2)),
        "selection_based_on_geometric_connectivity_not_face_id":True,
    }
    return volume,details


def new_one_sided_regions(mesh):
    bulk={}
    boundaries={}
    region_classification={}
    for number in (1,2):
        volume,details=bulk_2d_convexes(mesh,number)
        bulk[number]=volume
        region_classification[str(number)]=details
        faces=np.asarray(mesh.outer_faces(2,volume),dtype=np.int64)
        if faces.ndim!=2 or faces.shape[0]!=2 or not faces.shape[1]:
            raise RuntimeError(
                f"No one-sided 2-D outer faces for Concentric region {number}"
            )
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
        "bulk_region_geometric_classification":region_classification,
        "parent_face_normals":normals,
        "no_volume_region_modified":True,
    }


