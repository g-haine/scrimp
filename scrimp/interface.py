# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Conforming, oriented two-sided interface topology on one GetFEM mesh.

A public interface region remains a single logical region. This module
resolves matching parent-convex faces on each adjacent bulk region, without
changing the original mesh region, assuming a named geometry, or inserting
problem-dependent physical signs. Nonconforming and ambiguous topology is
rejected; it must never silently lead to a doubled boundary integral.
"""

from dataclasses import dataclass

import numpy as np


class InterfaceTopologyError(ValueError):
    """A logical interface cannot be assigned two conforming outward sides."""


@dataclass(frozen=True)
class InterfaceSides:
    """Resolved faces in the same geometric-facet order on both sides.

    Each faces entry is a read-only 2-by-N array of GetFEM (convex, face)
    pairs. Normal orientation follows the parent convex, not the order
    of the logical interface or of the input subdomain tuple.
    """

    region: int
    subdomains: tuple
    faces: tuple
    facet_keys: tuple

    def for_subdomain(self, region):
        try:
            i = self.subdomains.index(int(region))
        except ValueError as exc:
            raise InterfaceTopologyError(
                f"Region {region} is not adjacent to interface {self.region}"
            ) from exc
        return self.faces[i].copy()


def _pairs(value, context):
    arr = np.asarray(value, dtype=np.int64)
    if arr.ndim != 2 or arr.shape[0] != 2:
        raise InterfaceTopologyError(
            f"{context}: expected a 2-by-N array of GetFEM region entries"
        )
    if np.any(arr[1] < -1):
        raise InterfaceTopologyError(f"{context}: invalid face indices")
    return arr


def _bulk_convexes(mesh, dimension, region):
    raw = _pairs(mesh.region(region), f"bulk {region}")
    convex_ids = np.unique(raw[0])
    if not convex_ids.size:
        raise InterfaceTopologyError(f"Empty bulk region {region}")
    # A GetFEM Gmsh-with-lower-dim-elt import can include embedded line
    # elements. Node counts cannot determine dimensions for Qk, prisms,
    # or curved geometric elements. Inspect each geometric transformation.
    transformations = mesh.geotrans(convex_ids)
    if len(transformations) != len(convex_ids):
        raise InterfaceTopologyError("GetFEM geotrans lookup size mismatch")
    result = np.asarray(
        [cv for cv, transform in zip(convex_ids, transformations)
         if int(transform.dim()) == dimension],
        dtype=np.int64,
    )
    if not result.size:
        raise InterfaceTopologyError(
            f"No dimension-{dimension} convexes in subdomain {region}"
        )
    return result


def _facet_key(mesh, face):
    points = np.asarray(mesh.pid_in_faces(face), dtype=np.int64).reshape(-1)
    if not points.size or np.unique(points).size != points.size:
        raise InterfaceTopologyError("Empty or repeated facet point indices")
    return tuple(sorted(int(p) for p in points))


def _region_keys(mesh, dimension, region):
    entries = _pairs(mesh.region(region), f"logical interface {region}")
    if not entries.shape[1]:
        raise InterfaceTopologyError(f"Empty logical interface {region}")
    ids = np.unique(entries[0])
    transforms = dict(zip((int(x) for x in ids), mesh.geotrans(ids)))
    keys = set()
    for i in range(entries.shape[1]):
        cv, local_face = (int(v) for v in entries[:, i])
        transform = transforms[cv]
        if local_face >= 0:
            if int(transform.dim()) != dimension:
                raise InterfaceTopologyError(
                    "Logical interface includes a face of wrong dimension"
                )
            key = _facet_key(mesh, entries[:, i:i + 1])
        else:
            # Physical groups can also be imported as full (d-1)-cells.
            # GetFEM returns (points, offsets) for convex connectivity.
            if int(transform.dim()) != dimension - 1:
                raise InterfaceTopologyError(
                    "Logical interface contains a non-facet convex"
                )
            pids, offsets = mesh.pid_from_cvid(np.asarray([cv], dtype=np.int64))
            pids = np.asarray(pids, dtype=np.int64).reshape(-1)
            offsets = np.asarray(offsets, dtype=np.int64).reshape(-1)
            if offsets.size != 2 or offsets[0] != 0 or offsets[-1] != pids.size:
                raise InterfaceTopologyError("Invalid lower-dimensional connectivity")
            key = tuple(sorted(int(p) for p in pids))
            if not key:
                raise InterfaceTopologyError("Empty lower-dimensional facet")
        keys.add(key)
    return keys


def _boundary_map(mesh, dimension, convex_ids, region):
    faces = _pairs(
        mesh.outer_faces(dimension, convex_ids),
        f"subdomain {region} boundary",
    )
    result = {}
    for i in range(faces.shape[1]):
        face = faces[:, i:i + 1]
        key = _facet_key(mesh, face)
        if key in result:
            raise InterfaceTopologyError(
                f"Subdomain {region} has duplicate boundary facet {key}"
            )
        result[key] = faces[:, i].copy()
    return result


def resolve_conforming_interface(mesh, dimension, interface, subdomains):
    """Return outward-facing parent facets for a two-domain interface.

    This is only a *geometric* operation. It does not infer the physics,
    automatically rewrite weak forms, or provide nonconforming/mortar
    coupling. Both adjacent subdomains must be on the SAME conforming
    GetFEM mesh, with shared point IDs along their interface.
    """
    dimension = int(dimension)
    if dimension not in (1, 2, 3):
        raise InterfaceTopologyError("Bulk dimension must be 1, 2, or 3")
    subdomains = tuple(int(x) for x in subdomains)
    if len(subdomains) != 2 or len(set(subdomains)) != 2:
        raise InterfaceTopologyError("Exactly two distinct subdomains are required")
    convexes = tuple(
        _bulk_convexes(mesh, dimension, region) for region in subdomains
    )
    if set(convexes[0]).intersection(convexes[1]):
        raise InterfaceTopologyError("Subdomains overlap in bulk convexes")

    expected = _region_keys(mesh, dimension, int(interface))
    boundaries = tuple(
        _boundary_map(mesh, dimension, cvs, region)
        for cvs, region in zip(convexes, subdomains)
    )
    shared = set(boundaries[0]).intersection(boundaries[1])
    if shared != expected:
        raise InterfaceTopologyError(
            f"Interface {interface}: unmatched logical facets="
            f"{len(expected - shared)}, other shared facets="
            f"{len(shared - expected)}; no nonconforming fallback"
        )

    ordered = tuple(sorted(expected))
    result_faces = []
    normals = []
    for boundary in boundaries:
        arr = np.ascontiguousarray(
            np.column_stack([boundary[key] for key in ordered]),
            dtype=np.int64,
        )
        values = np.asarray(mesh.normal_of_faces(arr), dtype=float)
        if values.shape != (dimension, len(ordered)) or not np.all(
            np.isfinite(values)
        ):
            raise InterfaceTopologyError("Invalid GetFEM face normals")
        lengths = np.linalg.norm(values, axis=0)
        if np.any(lengths <= 0):
            raise InterfaceTopologyError("Zero GetFEM face normal")
        normals.append(values / lengths[None, :])
        arr.setflags(write=False)
        result_faces.append(arr)
    if not np.allclose(
        np.sum(normals[0] * normals[1], axis=0),
        -1.0, atol=1e-6, rtol=0,
    ):
        raise InterfaceTopologyError("Interface outward normals are not opposite")
    return InterfaceSides(
        int(interface), subdomains, tuple(result_faces), ordered,
    )
