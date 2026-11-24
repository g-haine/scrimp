# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Native two-adjacent-rectangle GetFEM interface topology regression.

Conforming Cartesian Q1 mesh: left [0,1]x[0,1] and right [1,2]x[0,1].
The common interface occupies only a SMALL PART of each rectangle's
boundary. A missing labeled shared facet must fail closed.
No simulation, temporal initializer, or PETSc matrix assembly.
"""

import json
import numpy as np
import getfem as gf

from scrimp.interface import (
    InterfaceTopologyError, resolve_conforming_interface,
)


def _region_of_cells(mesh, ids):
    ids = np.asarray(ids, dtype=np.int64)
    return np.vstack((ids, -np.ones(ids.size, dtype=np.int64)))


def _key(mesh, face):
    return tuple(sorted(int(x) for x in
                        np.asarray(mesh.pid_in_faces(face)).reshape(-1)))


def _face_map(mesh, convexes):
    facets = np.asarray(mesh.outer_faces(2, convexes), dtype=np.int64)
    if facets.shape[0] != 2:
        raise AssertionError("GetFEM outer_faces must be 2-by-N")
    return {
        _key(mesh, facets[:, i:i+1]): facets[:, i].copy()
        for i in range(facets.shape[1])
    }


def main():
    # Grid has two quadrangles in x on EACH side of the interface
    # and two in y, so the interface has precisely two shared facets.
    mesh = gf.Mesh(
        "cartesian", np.linspace(0., 2., 5), np.linspace(0., 1., 3),
    )
    ids = np.asarray(mesh.cvid(), dtype=np.int64).reshape(-1)
    left, right = [], []
    for cv in ids:
        points, offsets = mesh.pid_from_cvid(np.asarray([cv], dtype=np.int64))
        pids = np.asarray(points, dtype=np.int64).reshape(-1)
        coords = np.asarray(mesh.pts(pids), dtype=float)
        x = float(np.mean(coords[0, :]))
        if x < 1.:
            left.append(int(cv))
        elif x > 1.:
            right.append(int(cv))
        else:
            raise AssertionError("Cell centroid exactly on interface")
    if len(left) != 4 or len(right) != 4:
        raise AssertionError(f"Unexpected Cartesian cell split: {len(left)}, {len(right)}")
    mesh.set_region(1, _region_of_cells(mesh, left))
    mesh.set_region(2, _region_of_cells(mesh, right))

    left_boundaries = _face_map(mesh, np.asarray(left, dtype=np.int64))
    right_boundaries = _face_map(mesh, np.asarray(right, dtype=np.int64))
    shared = set(left_boundaries).intersection(right_boundaries)
    if len(shared) != 2:
        raise AssertionError(f"Expected two matching interface facets, got {len(shared)}")
    if not (len(left_boundaries) > len(shared)
            and len(right_boundaries) > len(shared)):
        raise AssertionError("Interface must be PROPER subset of both boundaries")

    # A logical interface may be represented by the two sets of parent
    # faces. This mimics existing Gmsh regions holding both oriented sides.
    logical = np.column_stack(
        [mapping[key] for key in sorted(shared)
         for mapping in (left_boundaries, right_boundaries)]
    )
    mesh.set_region(10, np.asarray(logical, dtype=np.int64))
    initial = np.asarray(mesh.region(10), dtype=np.int64).copy()
    initial_regions = np.asarray(mesh.regions()).reshape(-1).copy()
    sides = resolve_conforming_interface(mesh, 2, 10, (1, 2))
    swapped = resolve_conforming_interface(mesh, 2, 10, (2, 1))
    for subdomain in (1, 2):
        if not np.array_equal(sides.for_subdomain(subdomain),
                              swapped.for_subdomain(subdomain)):
            raise AssertionError("Subdomain ordering changed interface faces")
    if not np.array_equal(mesh.region(10), initial):
        raise AssertionError("Original interface region was modified")
    if not np.array_equal(np.asarray(mesh.regions()).reshape(-1), initial_regions):
        raise AssertionError("Resolver allocated extra GetFEM regions")

    # Missing a real shared facet from region 10 is a configuration error,
    # not permission to integrate on a partial or doubled interface.
    first = sorted(shared)[0]
    mesh.set_region(11, np.column_stack(
        (left_boundaries[first], right_boundaries[first])
    ).reshape(2, 2))
    try:
        resolve_conforming_interface(mesh, 2, 11, (1, 2))
    except InterfaceTopologyError as exc:
        if "other shared facets" not in str(exc):
            raise
    else:
        raise AssertionError("Unlabeled shared face was silently accepted")

    print(json.dumps({
        "schema": "scrimp-v130-native-two-rectangles-interface-v1",
        "status": "PASS_PARTIAL_BOUNDARY_INTERFACE_TOPOLOGY_ONLY",
        "dimension": 2,
        "element_family": "Cartesian quadrilateral",
        "left_cells": len(left),
        "right_cells": len(right),
        "shared_facets": len(shared),
        "left_boundary_facets": len(left_boundaries),
        "right_boundary_facets": len(right_boundaries),
        "proper_subset_of_both_boundaries": True,
        "unlabeled_shared_facets_rejected": True,
        "permutation_invariant": True,
        "region_10_unchanged": True,
        "no_new_regions_during_resolution": True,
        "weak_forms_assembled": False,
        "PETSc_TS_executed": False,
    }, indent=2))


if __name__ == "__main__":
    main()
