# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Regression tests for the current Port/Parameter API."""

import unittest

from scrimp.port import Parameter, Port


class TestPort(unittest.TestCase):
    def setUp(self):
        self.port = Port(
            "resistor", "flow", "effort", "scalar-field",
            mesh_id=0, algebraic=True, substituted=False,
            dissipative=True, region=1,
        )

    def test_port_metadata(self):
        self.assertEqual(self.port.get_name(), "resistor")
        self.assertEqual(self.port.get_flow(), "flow")
        self.assertEqual(self.port.get_effort(), "effort")
        self.assertEqual(self.port.get_kind(), "scalar-field")
        self.assertEqual(self.port.get_mesh_id(), 0)
        self.assertEqual(self.port.get_region(), 1)
        self.assertTrue(self.port.get_algebraic())
        self.assertFalse(self.port.get_substituted())
        self.assertTrue(self.port.get_dissipative())

    def test_no_fem_and_not_registered_initially(self):
        self.assertIsNone(self.port.get_fem())
        self.assertFalse(self.port.get_isSet())
        self.port.set_isSet()
        self.assertTrue(self.port.get_isSet())

    def test_parameters(self):
        parameter = Parameter(
            "resistance", "Electrical resistance", "scalar-field",
            "10.0", self.port.get_name(),
        )
        self.assertTrue(self.port.add_parameter(parameter))
        self.assertIs(self.port.get_parameter("resistance"), parameter)
        self.assertEqual(self.port.get_parameters(), [parameter])
        # A caller must not be able to mutate the internal list.
        copy = self.port.get_parameters()
        copy.clear()
        self.assertEqual(self.port.get_parameters(), [parameter])

    def test_missing_parameter_raises(self):
        with self.assertRaises(AssertionError):
            self.port.get_parameter("missing")

    def test_wrong_parameter_type_raises(self):
        with self.assertRaises(TypeError):
            self.port.add_parameter(None)

    def test_string_representation_of_registered_port(self):
        self.port.set_isSet()
        description = str(self.port)
        self.assertIn("resistor, flow, effort, scalar-field, 0", description)
        self.assertIn("True, False, True", description)

    def test_power_computation_flag(self):
        self.assertFalse(self.port.get_is_computed())
        self.port.set_is_computed()
        self.assertTrue(self.port.get_is_computed())


if __name__ == "__main__":
    unittest.main()
