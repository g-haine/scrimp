# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp


"""Native bridge from generic Domain.resolve_interface to reference facets.

Strictly assembly PREPARATION: the generic Domain resolver constructs the
same oriented parent-convex face sets as the previously qualified, legacy
Concentric-only sandbox helper. No bulk FEM or DPHS weak form is changed.
This probe is not a proof of automatic interface assembly.
"""

import json

import numpy as np

import scrimp as S
from v130_hw_interface_geometry import new_one_sided_regions


def _keys(face_array):
    return set(map(tuple, np.asarray(face_array, dtype=np.int64).T.tolist()))


def main():
    domain = S.Domain("Concentric", {"R": 1.0, "r": 0.6, "h": 0.05})
    mesh = domain.get_mesh()[0]
    original_10 = np.asarray(mesh.region(10), dtype=np.int64).copy()

    sides = domain.resolve_interface(10, (1, 2))
    left = sides.for_subdomain(1)
    right = sides.for_subdomain(2)
    swapped = domain.resolve_interface(10, (2, 1))
    if not (np.array_equal(left, swapped.for_subdomain(1))
            and np.array_equal(right, swapped.for_subdomain(2))):
        raise AssertionError("Resolver is not invariant under subdomain swap")

    if _keys(left) & _keys(right):
        raise AssertionError("A parent face appears on both oriented sides")
    if left.shape != (2, 78) or right.shape != (2, 78):
        raise AssertionError("Unexpected original Heat-Wave face count")

    # Compare the FULL convex-face IDs, not just node-based signatures
    # or counts. The reference uses fixed radius logic; new production
    # resolver must not depend on it.
    previous = new_one_sided_regions(mesh)
    old_left = np.asarray(mesh.region(30), dtype=np.int64)
    old_right = np.asarray(mesh.region(31), dtype=np.int64)
    if _keys(left) != _keys(old_left) or _keys(right) != _keys(old_right):
        raise AssertionError("Generic and historical face maps differ")
    if not np.array_equal(np.asarray(mesh.region(10)), original_10):
        raise AssertionError("Original interface region was modified")

    print(json.dumps({
        "schema": "scrimp-v130-generic-to-reference-facets-v1",
        "status": "PASS_GENERIC_TO_ORIGINAL_ORIENTED_FACET_PARITY_ONLY",
        "logical_interface": 10,
        "subdomain_ids": [1, 2],
        "facets_per_side": [left.shape[1], right.shape[1]],
        "exact_convex_face_parity": True,
        "permutation_invariant": True,
        "legacy_radius_based_helper_used_only_as_reference": True,
        "original_logical_region_unchanged": True,
        "legacy_reference_regions_created": [30, 31],
        "previous_reference_geometry": {
            "disk_facets": previous["disk_interface_2d_face_region"]["convex_face_pairs"],
            "annulus_facets": previous["annulus_interface_2d_face_region"]["convex_face_pairs"],
        },
        "getfem_bricks_assembled": False,
        "dae_or_cn_executed": False,
        "energy_physics_qualified": False,
        "MPI2_qualified": False,
    }, indent=2))


if __name__ == "__main__":
    main()
