# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Dependency-light face-region regression for the actual CN collector.

This verifies ONLY region argument routing and metadata with synthetic
five-sample traces. It does not certify numerical physics, DAE, or PETSc.
The physical gate remains the pinned-container original Heat-Wave run.
"""

import importlib.util
from pathlib import Path
import unittest

import numpy as np


SOURCE = (Path(__file__).resolve().parents[1] / "sandbox" /
          "measure_v130_hw_dual_facet_physics.py")
SPEC = importlib.util.spec_from_file_location("hw_cn_facet_collector", SOURCE)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class FakeMesh:
    def __init__(self, regions):
        self.region_ids = np.asarray(regions, dtype=np.int64)

    def regions(self):
        return self.region_ids


class FakeDomain:
    def __init__(self, regions):
        self.mesh = FakeMesh(regions)

    def get_mesh(self):
        return [self.mesh]


class FakeSystem:
    def __init__(self, regions):
        self.domain = FakeDomain(regions)
        self.calls = []

    def get_quantity(self, expression, *, region, CN=False):
        region = int(region)
        self.calls.append((expression, region, bool(CN)))
        if region == 10:
            raise AssertionError("Power must never use doubled logical region")
        if expression == ".5*T*T":
            value = 1.0
        elif expression in (".5*p*p", ".5*q.q"):
            value = 0.5
        elif expression == "e_Q.e_Q":
            value = 0.1
        elif expression in ("U_T*Y_T", "U_w*Y_w"):
            value = 0.2
        elif expression == "U_bnd*Y_bnd":
            value = 0.0
        else:
            value = 0.0
        return np.full(5, value, dtype=float)


class CollectorRegionTests(unittest.TestCase):
    def run_collector(self, regions, **kwargs):
        system = FakeSystem(regions)
        times = np.arange(5, dtype=float) * 0.001
        energy = np.full(5, 2.0, dtype=float)
        traces, report = MODULE.collect_heat_wave_dual_facet_physics(
            system, times, energy, **kwargs,
        )
        self.assertEqual(traces["H_original"].shape, (5,))
        self.assertLess(report["Hamiltonian_terms_max_difference"], 1e-14)
        self.assertTrue(report["old_doubled_interface_region_excluded_from_physical_power"])
        return system, report

    def test_default_reference_regions_preserved(self):
        system, report = self.run_collector([1, 2, 10, 20, 30, 31])
        self.assertEqual(report["interface_heat_region"], 30)
        self.assertEqual(report["interface_wave_region"], 31)
        self.assertIn(("U_T*Y_T", 30, True), system.calls)
        self.assertIn(("U_w*Y_w", 31, True), system.calls)

    def test_dynamic_compiled_regions_and_correct_metadata(self):
        system, report = self.run_collector(
            [1, 2, 10, 20, 21, 22], heat_faces=21, wave_faces=22,
        )
        self.assertEqual(report["interface_heat_region"], 21)
        self.assertEqual(report["interface_wave_region"], 22)
        self.assertIn(("U_T*Y_T", 21, True), system.calls)
        self.assertIn(("U_w*Y_w", 22, True), system.calls)
        self.assertNotIn(30, [r for _, r, _ in system.calls])
        self.assertNotIn(31, [r for _, r, _ in system.calls])
        self.assertEqual(
            report["CN_interface_power_balance_weak_form_expected_sign"],
            "-heat_region21 + wave_region22",
        )

    def test_dynamic_call_cannot_change_reference_defaults(self):
        self.run_collector(
            [1, 2, 10, 20, 21, 22], heat_faces=21, wave_faces=22,
        )
        _, report = self.run_collector([1, 2, 10, 20, 30, 31])
        self.assertEqual((report["interface_heat_region"],
                          report["interface_wave_region"]), (30, 31))

    def test_reject_invalid_or_missing_facets(self):
        times = np.arange(5, dtype=float) * 0.001
        energy = np.full(5, 2.0)
        for kwargs, pattern in (
            ({"heat_faces": 21, "wave_faces": 21}, "distinct"),
            ({"heat_faces": 10, "wave_faces": 22}, "distinct"),
            ({"heat_faces": 21, "wave_faces": 31}, "missing"),
        ):
            with self.subTest(kwargs=kwargs):
                with self.assertRaisesRegex(RuntimeError, pattern):
                    MODULE.collect_heat_wave_dual_facet_physics(
                        FakeSystem([1, 2, 10, 20, 21, 22]),
                        times, energy, **kwargs,
                    )


if __name__ == "__main__":
    unittest.main()
