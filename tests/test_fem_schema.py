# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""FEM regression tests for the currently supported imperative API.

The obsolete FEMFieldSchema / infer_dim interface is not available in
SCRIMP 1.2.0; the concrete FEM implementation is tested instead.
"""

import unittest

from scrimp import Domain
from scrimp.fem import FEM


class TestFEM(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.domain = Domain("Interval", {"L": 1.0, "h": 0.5})

    def test_default_continuous_lagrange_metadata(self):
        fem = FEM("velocity", 1)
        self.assertEqual(fem.get_name(), "velocity")
        self.assertEqual(fem.get_order(), 1)
        self.assertEqual(fem.get_type(), "CG")
        self.assertIsNone(fem.get_mesh())
        self.assertFalse(fem.get_isSet())

    def test_cg_finite_element_registration(self):
        fem = FEM("velocity", 1, FEM="CG")
        fem.set_mesh(self.domain.get_mesh()[0])
        fem.set_dim(1)
        fem.set_fem()
        self.assertTrue(fem.get_isSet())
        self.assertEqual(fem.get_dim(), 1)
        self.assertIsNotNone(fem.get_fem())

    def test_dg_finite_element_registration(self):
        fem = FEM("pressure", 0, FEM="DG")
        fem.set_mesh(self.domain.get_mesh()[0])
        fem.set_dim(1)
        fem.set_fem()
        self.assertTrue(fem.get_isSet())
        self.assertEqual(fem.get_type(), "DG")
        self.assertIsNotNone(fem.get_fem())


if __name__ == "__main__":
    unittest.main()
