# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Dependency-light topology regression tests; no GetFEM/PETSc import."""

import importlib.util
from pathlib import Path
import unittest

import numpy as np


MODULE = Path(__file__).resolve().parents[1] / "scrimp" / "interface.py"
SPEC = importlib.util.spec_from_file_location("scrimp_interface_topology_test", MODULE)
interface_module = importlib.util.module_from_spec(SPEC)
import sys
sys.modules[SPEC.name] = interface_module
SPEC.loader.exec_module(interface_module)

resolve = interface_module.resolve_conforming_interface
TopologyError = interface_module.InterfaceTopologyError


class Transform:
    def __init__(self, dimension):
        self.dimension = dimension

    def dim(self):
        return self.dimension


class TwoCellMesh:
    """Two same-dimensional convexes sharing exactly one facet.

    Dimension 1/2/3 uses respectively 1/2/3 shared mesh point IDs.
    A detached (d-1)-cell may represent the imported physical interface.
    """

    def __init__(self, dimension=2, lower_dim_interface=False):
        self.dimension = dimension
        self.lower_dim_interface = lower_dim_interface
        self.shared = tuple(range(10, 10 + dimension))
        self.other_a = tuple(range(20, 20 + dimension))
        self.other_b = tuple(range(30, 30 + dimension))
        self.original = np.asarray([[1, 2], [0, 0]], dtype=np.int64)
        if lower_dim_interface:
            self.original = np.asarray([[90], [-1]], dtype=np.int64)
        self._normal_flip = False
        self._bad_point = False

    def region(self, number):
        if number == 1:
            return np.asarray([[1, 90], [-1, -1]], dtype=np.int64)
        if number == 2:
            return np.asarray([[2], [-1]], dtype=np.int64)
        if number == 10:
            return self.original
        raise KeyError(number)

    def geotrans(self, cvids=None):
        return [
            Transform(self.dimension - 1 if int(cv) == 90 else self.dimension)
            for cv in np.asarray(cvids).reshape(-1)
        ]

    def pid_from_cvid(self, cvids):
        if np.asarray(cvids).tolist() != [90]:
            raise AssertionError("Only lower-dimensional cell is supported")
        return np.asarray(self.shared), np.asarray([0, len(self.shared)])

    def outer_faces(self, dimension, cvids):
        assert dimension == self.dimension
        cv = int(np.asarray(cvids).reshape(-1)[0])
        assert cv in (1, 2)
        return np.asarray([[cv, cv], [0, 1]], dtype=np.int64)

    def pid_in_faces(self, faces):
        cv, face = [int(x) for x in np.asarray(faces)[:, 0]]
        if face == 0:
            if cv == 2 and self._bad_point:
                return np.asarray((*self.shared[:-1], 999), dtype=np.int64)
            return np.asarray(self.shared[::-1] if cv == 2 else self.shared)
        if (cv, face) == (1, 1):
            return np.asarray(self.other_a)
        if (cv, face) == (2, 1):
            return np.asarray(self.other_b)
        raise AssertionError((cv, face))

    def normal_of_faces(self, faces):
        cvs = np.asarray(faces)[0]
        result = np.zeros((self.dimension, len(cvs)), dtype=float)
        result[0, :] = [1. if int(cv) == 1 or self._normal_flip else -1.
                        for cv in cvs]
        return result


class InterfaceTopologyTests(unittest.TestCase):
    def test_dimension_independent_oriented_sides(self):
        for dim in (1, 2, 3):
            with self.subTest(dim=dim):
                mesh = TwoCellMesh(dim)
                result = resolve(mesh, dim, 10, (1, 2))
                self.assertEqual(result.region, 10)
                self.assertEqual(result.subdomains, (1, 2))
                self.assertEqual(result.facet_keys, (mesh.shared,))
                self.assertEqual(result.for_subdomain(1).tolist(), [[1], [0]])
                self.assertEqual(result.for_subdomain(2).tolist(), [[2], [0]])
                np.testing.assert_array_equal(mesh.region(10), mesh.original)

    def test_permutation_does_not_change_outward_parent_faces(self):
        mesh = TwoCellMesh()
        result = resolve(mesh, 2, 10, (2, 1))
        self.assertEqual(result.for_subdomain(1).tolist(), [[1], [0]])
        self.assertEqual(result.for_subdomain(2).tolist(), [[2], [0]])

    def test_mixed_dimensional_imported_interface(self):
        mesh = TwoCellMesh(lower_dim_interface=True)
        result = resolve(mesh, 2, 10, (1, 2))
        self.assertEqual(result.facet_keys, (mesh.shared,))
        self.assertEqual(mesh.region(10).tolist(), [[90], [-1]])

    def test_non_matching_interface_is_rejected(self):
        mesh = TwoCellMesh()
        mesh._bad_point = True
        with self.assertRaisesRegex(TopologyError, "unmatched logical facets"):
            resolve(mesh, 2, 10, (1, 2))

    def test_incoherent_normals_are_rejected(self):
        mesh = TwoCellMesh()
        mesh._normal_flip = True
        with self.assertRaisesRegex(TopologyError, "normals are not opposite"):
            resolve(mesh, 2, 10, (1, 2))

    def test_invalid_adjacency_is_rejected(self):
        with self.assertRaisesRegex(TopologyError, "distinct subdomains"):
            resolve(TwoCellMesh(), 2, 10, (1, 1))

    def test_returned_face_copy_cannot_mutate_cache(self):
        mesh = TwoCellMesh()
        resolved = resolve(mesh, 2, 10, (1, 2))
        returned = resolved.for_subdomain(1)
        returned[0, 0] = 99
        self.assertEqual(resolved.for_subdomain(1)[0, 0], 1)


if __name__ == "__main__":
    unittest.main()
