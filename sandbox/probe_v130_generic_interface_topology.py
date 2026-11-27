# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp


"""Native GetFEM smoke test for the production Domain interface resolver.

Run in the qualified SCRIMP environment / v1.3.0 Docker image:
    python sandbox/probe_v130_generic_interface_topology.py

The probe does not modify region 10, create extra regions, register FEM
unknowns, assemble weak forms, or claim a physical Heat-Wave PASS.
"""

import json

import numpy as np
import scrimp as S


def _compare_faces(a, b):
    return (
        np.array_equal(a.for_subdomain(1), b.for_subdomain(1))
        and np.array_equal(a.for_subdomain(2), b.for_subdomain(2))
    )


def main():
    domain = S.Domain("Concentric", {"R": 1.0, "r": 0.6, "h": 0.05})
    mesh = domain.get_mesh()[0]
    before = np.asarray(mesh.region(10), dtype=np.int64).copy()
    original_regions = tuple(int(x) for x in np.asarray(mesh.regions()).reshape(-1))

    sides = domain.resolve_interface(10, (1, 2))
    inverted = domain.resolve_interface(10, (2, 1))
    if not _compare_faces(sides, inverted):
        raise AssertionError("Subdomain permutation changed parent face selection")

    disk = sides.for_subdomain(1)
    annulus = sides.for_subdomain(2)
    if disk.shape != annulus.shape or disk.ndim != 2 or disk.shape[0] != 2:
        raise AssertionError("Inconsistent oriented side shapes")
    if disk.shape[1] == 0:
        raise AssertionError("Empty oriented interface")

    normal_disk = np.asarray(mesh.normal_of_faces(disk), dtype=float)
    normal_annulus = np.asarray(mesh.normal_of_faces(annulus), dtype=float)
    n_disk = normal_disk / np.linalg.norm(normal_disk, axis=0)
    n_annulus = normal_annulus / np.linalg.norm(normal_annulus, axis=0)
    max_defect = float(np.max(np.linalg.norm(n_disk + n_annulus, axis=0)))
    if max_defect > 1e-6:
        raise AssertionError(f"Opposite-normal error: {max_defect}")

    if not np.array_equal(mesh.region(10), before):
        raise AssertionError("Logical interface region 10 was modified")
    if tuple(int(x) for x in np.asarray(mesh.regions()).reshape(-1)) != original_regions:
        raise AssertionError("Resolver unexpectedly created GetFEM regions")

    print(json.dumps({
        "schema": "scrimp-v130-generic-interface-topology-native-v1",
        "status": "PASS_ORIENTED_TOPOLOGY_ONLY",
        "logical_interface_region": 10,
        "subdomains": [1, 2],
        "facets_per_side": int(disk.shape[1]),
        "logical_region_unchanged": True,
        "no_extra_getfem_regions": True,
        "permutation_invariant": True,
        "max_opposite_normal_vector_defect": max_defect,
        "weak_forms_assembled": False,
        "DAE_initialized": False,
        "energy_physics_qualified": False,
        "MPI2_qualified": False,
    }, indent=2))


if __name__ == "__main__":
    main()
