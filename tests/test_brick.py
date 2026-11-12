# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Regression tests for the supported Brick API (without obsolete structures)."""

import unittest

from scrimp.brick import Brick


class TestBrick(unittest.TestCase):
    def setUp(self):
        self.brick = Brick(
            "constitutive_law", "e*u*Test_u", [1, 2],
            linear=False, dt=False, position="constitutive",
            explicit=True, mesh_id=3,
        )

    def test_identity_and_weak_form(self):
        self.assertEqual(self.brick.get_name(), "constitutive_law")
        self.assertEqual(self.brick.get_form(), "e*u*Test_u")
        self.assertEqual(self.brick.get_regions(), [1, 2])
        self.assertEqual(self.brick.get_mesh_id(), 3)

    def test_brick_flags(self):
        self.assertFalse(self.brick.get_linear())
        self.assertFalse(self.brick.get_dt())
        self.assertTrue(self.brick.get_explicit())
        self.assertEqual(self.brick.get_position(), "constitutive")

    def test_flow_brick_has_time_derivative(self):
        brick = Brick("M", "p*Test_p", [1], dt=True, position="flow")
        self.assertTrue(brick.get_linear())
        self.assertTrue(brick.get_dt())
        self.assertFalse(brick.get_explicit())
        self.assertEqual(brick.get_position(), "flow")
        self.assertEqual(brick.get_mesh_id(), 0)

    def test_brick_ids_are_registered_and_returned_as_copy(self):
        self.assertEqual(self.brick.get_id_bricks(), [])
        self.brick.add_id_brick_to_list(7)
        self.brick.add_id_brick_to_list(11)
        ids = self.brick.get_id_bricks()
        self.assertEqual(ids, [7, 11])
        ids.append(99)
        self.assertEqual(self.brick.get_id_bricks(), [7, 11])


if __name__ == "__main__":
    unittest.main()
