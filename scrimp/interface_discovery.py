# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Discover conforming physical interface labels without geometric heuristics.

Read-only discovery. Numerical region IDs are never inferred from names,
radius, element families, or ordering. The complete shared boundary of
two adjacent bulk subdomains must have exactly one logical label in this
initial implementation. Multiple separately labeled patches along the
same pair require a later generalized partition policy and fail closed.
"""

from itertools import combinations

from scrimp.interface import (
    InterfaceTopologyError,
    _boundary_map,
    _bulk_convexes,
    _region_keys,
    resolve_conforming_interface,
)


def discover_conforming_interfaces(domain):
    """Return {(mesh_id, logical_region): InterfaceSides}, no mutations.

    Candidate internal regions are identified by positive overlap
    between the region's facet signatures and the common boundary of
    two distinct bulk subdomains. All candidates then undergo the
    strict conforming, opposite-normal verification.
    """
    detected = {}
    for mesh_id, mesh in enumerate(domain.get_mesh()):
        dimension = int(domain.get_dim()[mesh_id])
        bulks = sorted({int(r)
                        for r in domain.get_subdomains()[mesh_id].values()})
        labels = sorted({int(r)
                         for r in domain.get_boundaries()[mesh_id].values()})
        if len(bulks) < 2 or not labels:
            continue
        boundary_keys = {
            region: set(_boundary_map(
                mesh, dimension, _bulk_convexes(mesh, dimension, region),
                region,
            )) for region in bulks
        }
        label_keys = {
            region: _region_keys(mesh, dimension, region)
            for region in labels
        }
        for a, b in combinations(bulks, 2):
            shared = boundary_keys[a].intersection(boundary_keys[b])
            if not shared:
                continue
            matches = [
                label for label, keys in label_keys.items()
                if shared.intersection(keys)
            ]
            if not matches:
                continue
            if len(matches) != 1:
                raise InterfaceTopologyError(
                    f"Subdomains {a},{b} have multiple logical interface "
                    f"labels {matches}; shared-interface partitioning "
                    "is not yet supported"
                )
            label = matches[0]
            sides = resolve_conforming_interface(
                mesh, dimension, label, (a, b)
            )
            key = (mesh_id, label)
            if key in detected:
                raise InterfaceTopologyError(
                    f"Logical interface {label} has more than two "
                    "adjacent subdomains"
                )
            detected[key] = sides
    return detected
