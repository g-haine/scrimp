# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

# SCRIMP v1.3.0 — source-only binary GetFEM VTK decoding regression.
# Native 2D FE data is not synthesized by this test; actual qualification
# remains in the MPI pilots. This synthetic fixture guards the byte format.
import base64
import importlib.util
from pathlib import Path
import struct
import unittest
from xml.etree import ElementTree as ET

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT/"scripts"/"render-v130-wave2d-native.py"
spec = importlib.util.spec_from_file_location("scrimp_wave2d_renderer", SCRIPT)
renderer = importlib.util.module_from_spec(spec)
spec.loader.exec_module(renderer)


def data_array(values, *, vtk_type="Float32", components=1, field=False,
               native_components=None):
    dtype={"Float32":"<f4","Int32":"<i4"}[vtk_type]
    payload=np.asarray(values,dtype=dtype).tobytes()
    if field:
        # GetFEM encodes the number of *original* components in Float32:
        # 2D vectors are padded to 3D in the emitted VTK point array.
        size=np.float32(native_components if native_components is not None
                        else len(values)).tobytes()
    else:
        size=struct.pack("<I",len(payload))
    node=ET.Element("DataArray",type=vtk_type,format="binary",
                    NumberOfComponents=str(components))
    node.text=base64.b64encode(size+payload).decode("ascii")
    return node


class VtkBinaryGetfemTests(unittest.TestCase):
    def test_mesh_header_is_byte_count(self):
        pts=np.array([[0.,0.,0.],[1.,0.,0.],[0.,1.,0.]],dtype=np.float32)
        node=data_array(pts.flatten(),components=3)
        np.testing.assert_array_equal(renderer.vtk_binary(node),pts)
        connect=data_array([0,1,2],vtk_type="Int32")
        np.testing.assert_array_equal(renderer.vtk_binary(connect),[0,1,2])

    def test_scalar_point_header_is_float_count_not_payload_size(self):
        values=np.array([.1,-.3,.7],dtype=np.float32)
        field=data_array(values,field=True)
        got=renderer.vtk_binary(field,expected_count=3)
        np.testing.assert_array_equal(got,values)
        self.assertEqual(got.shape,(3,))

    def test_two_component_vector_padded_three_in_vtk(self):
        xyz=np.array([[1.,2.,0.],[3.,4.,0.],[5.,6.,0.]],dtype=np.float32)
        node=data_array(xyz.flatten(),components=3,field=True,
                        native_components=6)
        np.testing.assert_array_equal(
            renderer.vtk_binary(node,expected_count=9),xyz)
        with self.assertRaises(ValueError):
            renderer.vtk_binary(node,expected_count=12)

    def test_conservative_and_dissipative_balance_not_recomputed(self):
        t=np.array([0.,.1,.2])
        w=np.array([0.,2.,3.])
        loss=np.array([0.,.5,.7])
        H=np.array([1.,2.5,3.3])
        traces=renderer.normalized_balance({
            "t":t,"H":H,"signed_work":w,
            "cumulative_loss":loss,
            "power_left":np.array([0.,20.,10.]),
            "dissipated_power":np.array([0.,5.,2.]),
        })
        np.testing.assert_allclose(traces["residual"],0.,atol=1e-14)

    def test_substituted_coenergy_internal_and_feedback_losses(self):
        t=np.array([0.,.1,.2])
        Li=np.array([0.,.1,.25])
        Lb=np.array([0.,.05,.1])
        H=np.array([2.,1.85,1.65])
        traces=renderer.normalized_balance({
            "t":t,"H":H,
            "cumulative_internal_loss":Li,
            "cumulative_feedback_loss":Lb,
            "power_boundary":np.array([0.,-.5,-.5]),
            "dissipation_internal":np.array([0.,1.,1.5]),
        })
        np.testing.assert_allclose(traces["residual"],0.,atol=1e-14)
        self.assertTrue(np.all(traces["W"]==0))


if __name__=="__main__":
    unittest.main()
