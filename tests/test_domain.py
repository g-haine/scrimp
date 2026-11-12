# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Regression tests for existing built-in SCRIMP domains."""

import unittest

from scrimp import Domain


class TestDomain(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.domain = Domain("Interval", {"L": 2.0, "h": 0.5})

    def test_interval_domain_is_initialized(self):
        self.assertEqual(self.domain.get_name(), "Interval")
        self.assertTrue(self.domain.get_isSet())
        self.assertEqual(self.domain.get_dim(), [1])
        self.assertEqual(len(self.domain.get_mesh()), 1)

    def test_interval_has_expected_physical_regions(self):
        self.assertEqual(self.domain.get_subdomains()[0]["Omega"], 1)
        boundaries = self.domain.get_boundaries()[0]
        self.assertEqual(boundaries["Gamma_Left"], 10)
        self.assertEqual(boundaries["Gamma_Right"], 11)


if __name__ == "__main__":
    unittest.main()
