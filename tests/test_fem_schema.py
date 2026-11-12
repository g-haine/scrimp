# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

import unittest

from scrimp.fem import FEM
from scrimp.io.schema_loader import FEMFieldSchema


class TestFEMSchema(unittest.TestCase):
    def test_vector_schema_infer_dim(self):
        schema = FEMFieldSchema(
            name="velocity",
            order=1,
            family="CG",
            value_type="vector",
            components=2,
        )
        fem = FEM(schema)
        self.assertEqual(fem.infer_dim(mesh_dim=3, kind="vector-field"), 2)

    def test_tensor_schema_shape(self):
        schema = FEMFieldSchema(
            name="stress",
            order=1,
            family="CG",
            value_type="tensor",
            shape=(3, 3),
        )
        fem = FEM(schema)
        self.assertEqual(fem.infer_dim(mesh_dim=3, kind="tensor-field"), 9)

    def test_tensor_default_from_mesh(self):
        fem = FEM("sigma", 1, value_type="tensor")
        self.assertEqual(fem.infer_dim(mesh_dim=2, kind="tensor-field"), 4)


if __name__ == "__main__":
    unittest.main()
