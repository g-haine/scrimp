# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Dependency-light discovery/deferred-compiler contract tests.

No GetFEM/PETSc dependency and no actual weak-form assembly: the real
numerical operator must still be qualified in native Docker.
"""

import sys
from pathlib import Path
import types
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
_fake_scrimp = types.ModuleType("scrimp")
_fake_scrimp.__path__ = [str(ROOT / "scrimp")]
_previous = sys.modules.get("scrimp")
sys.modules["scrimp"] = _fake_scrimp
try:
    from scrimp.interface_discovery import discover_conforming_interfaces
    from scrimp.interface_deferred import DeferredInterfaceCompiler
    from scrimp.interface import InterfaceTopologyError
finally:
    if _previous is None:
        del sys.modules["scrimp"]
    else:
        sys.modules["scrimp"] = _previous

from test_v130_interface_topology import TwoCellMesh


class Mesh(TwoCellMesh):
    def __init__(self):
        super().__init__(dimension=2)
        self.regions_map = {}
        self._labels = [1, 2, 10, 20]

    def region(self, number):
        if number in self.regions_map:
            return self.regions_map[number].copy()
        if number == 20:
            return np.asarray([[1], [1]], dtype=np.int64)
        return super().region(number)

    def regions(self):
        return np.asarray(self._labels, dtype=np.int64)

    def set_region(self, region, faces):
        number = int(region)
        if number in self._labels:
            raise AssertionError("Region ID must never be reused")
        self._labels.append(number)
        self.regions_map[number] = np.asarray(faces).copy()


class Domain:
    def __init__(self):
        self.mesh = Mesh()

    def get_mesh(self):
        return [self.mesh]

    def get_dim(self):
        return [2]

    def get_subdomains(self):
        return [{"A": 1, "B": 2}]

    def get_boundaries(self):
        return [{"Interface": 10, "External": 20}]


class Port:
    def __init__(self, name, region):
        self.name = name
        self.region = region
        self.effective = None

    def get_mesh_id(self):
        return 0

    def get_region(self):
        return self.region

    def set_assembly_region(self, number):
        if self.effective is not None:
            raise AssertionError("Twice bound")
        self.effective = number

    def get_isSet(self):
        return True

    def get_parameters(self):
        return []


class Brick:
    def __init__(self, name, region):
        self.name = name
        self.regions = [region]
        self.effective = None

    def get_mesh_id(self):
        return 0

    def get_regions(self):
        return self.regions

    def set_assembly_regions(self, regions):
        self.effective = list(regions)


class System:
    def __init__(self, domain):
        self.domain = domain
        self.ports = {"port-A": Port("port-A", 10),
                      "port-B": Port("port-B", 10)}
        self.bricks = {"A": Brick("A", 10), "B": Brick("B", 10)}
        self.events = []

    def audit_interface_brick_plan(self, logical, bulk, mesh_id=0):
        assert logical == 10 and bulk == (1, 2)
        return {
            "port_owners": {"port-A": 1, "port-B": 2},
            "brick_plan": [
                {"brick": "A", "owner": 1},
                {"brick": "B", "owner": 2},
            ],
        }

    def add_FEM(self, fem):
        self.events.append(("FEM", fem))

    def add_brick(self, brick):
        self.events.append(("Brick", brick.name))

    def init_parameter(self, name, port):
        raise AssertionError("No fake port has parameters")


class CompilerTests(unittest.TestCase):
    def test_discovers_shared_but_not_exterior_boundary(self):
        domain = Domain()
        found = discover_conforming_interfaces(domain)
        self.assertEqual(set(found), {(0, 10)})
        self.assertEqual(found[(0, 10)].subdomains, (1, 2))
        self.assertEqual(len(found[(0, 10)].facet_keys), 1)
        self.assertNotIn((0, 20), found)

    def test_deferred_replay_keeps_public_regions_and_order(self):
        domain = Domain()
        before = np.asarray(domain.mesh.region(10)).copy()
        system = System(domain)
        compiler = DeferredInterfaceCompiler(domain)
        compiler.record_fem("port-A")
        compiler.record_fem("port-B")
        compiler.record_brick(system.bricks["A"])
        compiler.record_brick(system.bricks["B"])
        report = compiler.compile(system)
        self.assertTrue(compiler.completed)
        self.assertEqual(report["fem_count"], 2)
        self.assertEqual(report["brick_count"], 2)
        self.assertEqual(system.events, [
            ("FEM", "port-A"), ("FEM", "port-B"),
            ("Brick", "A"), ("Brick", "B"),
        ])
        self.assertEqual(system.ports["port-A"].get_region(), 10)
        self.assertEqual(system.ports["port-B"].get_region(), 10)
        self.assertEqual(system.bricks["A"].get_regions(), [10])
        self.assertEqual(system.bricks["B"].get_regions(), [10])
        region_a = system.ports["port-A"].effective
        region_b = system.ports["port-B"].effective
        self.assertNotEqual(region_a, region_b)
        self.assertEqual(system.bricks["A"].effective, [region_a])
        self.assertEqual(system.bricks["B"].effective, [region_b])
        self.assertTrue(region_a > 20 and region_b > 20)
        np.testing.assert_array_equal(domain.mesh.region(10), before)
        self.assertEqual(len(domain.mesh.regions()), 6)
        self.assertIs(compiler.compile(system), report)

    def test_transparent_material_interface_replays_global_forms_unchanged(self):
        # 2C coenergy: q and p use the same global FEM on regions 1+2.
        # The concentric geometry has a conforming physical interface
        # (label 10), but NO interface-coupling ports or forms.
        domain = Domain()
        before = np.asarray(domain.mesh.region(10)).copy()
        system = System(domain)
        system.ports = {
            "bulk": Port("bulk", 1),
            "exterior": Port("exterior", 20),
        }
        system.bricks = {
            "bulk-1": Brick("bulk-1", 1),
            "bulk-2": Brick("bulk-2", 2),
        }
        def must_not_audit(*args, **kwargs):
            raise AssertionError(
                "No orientation audit allowed for transparent material interface"
            )
        system.audit_interface_brick_plan = must_not_audit
        compiler = DeferredInterfaceCompiler(domain)
        compiler.record_fem("bulk")
        compiler.record_fem("exterior")
        compiler.record_brick(system.bricks["bulk-1"])
        compiler.record_brick(system.bricks["bulk-2"])
        result = compiler.compile(system)
        self.assertEqual(
            result["status"], "SKIPPED_TRANSPARENT_MATERIAL_INTERFACE"
        )
        self.assertEqual(result["interfaces"], [])
        self.assertEqual(len(result["transparent_interfaces"]), 1)
        self.assertEqual(result["transparent_interfaces"][0]["logical_region"], 10)
        self.assertEqual(system.events, [
            ("FEM", "bulk"), ("FEM", "exterior"),
            ("Brick", "bulk-1"), ("Brick", "bulk-2"),
        ])
        self.assertEqual(domain.mesh.regions().tolist(), [1, 2, 10, 20])
        np.testing.assert_array_equal(domain.mesh.region(10), before)
        self.assertTrue(compiler.completed)
        self.assertIs(compiler.compile(system), result)

    def test_interface_declarations_never_silently_skipped(self):
        # Presence of EITHER an interface port OR an interface brick
        # must still go through fail-closed ownership and assembly checks.
        for declared in ("port", "brick"):
            with self.subTest(declared=declared):
                domain = Domain()
                system = System(domain)
                if declared == "port":
                    system.bricks = {"bulk": Brick("bulk", 1)}
                else:
                    system.ports = {"bulk": Port("bulk", 1)}
                def invalid(*args, **kwargs):
                    raise ValueError("unresolved declared interface ownership")
                system.audit_interface_brick_plan = invalid
                compiler = DeferredInterfaceCompiler(domain)
                compiler.record_fem("some-fem")
                with self.assertRaisesRegex(
                    ValueError, "unresolved declared interface ownership"
                ):
                    compiler.compile(system)
                self.assertFalse(compiler.completed)
                self.assertEqual(domain.mesh.regions().tolist(), [1, 2, 10, 20])
                self.assertEqual(system.events, [])

    def test_no_internal_label_does_not_enable_automatic_mode(self):
        domain = Domain()
        domain.get_boundaries = lambda: [{"External": 20}]
        self.assertEqual(discover_conforming_interfaces(domain), {})

    def test_invalid_logical_interface_rejected_before_mesh_change(self):
        domain = Domain()
        # Replace physical interface label with an external boundary.
        domain.mesh.original = np.asarray([[1], [1]], dtype=np.int64)
        self.assertEqual(discover_conforming_interfaces(domain), {})

    def test_preflight_rejects_ambiguous_owners_without_new_regions(self):
        domain = Domain()
        system = System(domain)
        compiler = DeferredInterfaceCompiler(domain)
        def invalid(*args, **kwargs):
            raise ValueError("interface weak forms have conflicting owners")
        system.audit_interface_brick_plan = invalid
        with self.assertRaisesRegex(ValueError, "conflicting owners"):
            compiler.compile(system)
        self.assertEqual(domain.mesh.regions().tolist(), [1, 2, 10, 20])


if __name__ == "__main__":
    unittest.main()
