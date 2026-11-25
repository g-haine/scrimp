# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
"""Isolated source/protocol tests: real long-time Heat–Wave visuals.

These lightweight tests never substitute native GetFEM/PETSc physical
qualification. Native acceptance is the dedicated archived Docker pilot.
"""
import ast
import base64
import struct
import tempfile
import unittest
from pathlib import Path
from xml.etree import ElementTree as ET

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
EXAMPLE=ROOT/"examples/heat_wave.py"
RUNNER=ROOT/"scripts/run-v130-heat-wave-long-showcase.sh"
PHYSICAL=ROOT/"scripts/run-v130-heat-wave-physical-pilot.sh"


class HeatWaveLongShowcaseSource(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.example=EXAMPLE.read_text(encoding="utf-8")
        cls.runner=RUNNER.read_text(encoding="utf-8")
        cls.tree=ast.parse(cls.example,filename=str(EXAMPLE))

    def test_vtk_array_reader_against_getfem_style_binary_and_ascii(self):
        funcs=[
            n for n in self.tree.body
            if isinstance(n,ast.FunctionDef)
            and n.name=="_read_native_vtk_data_array"
        ]
        self.assertEqual(len(funcs),1)
        m=ast.Module(body=funcs,type_ignores=[])
        ast.fix_missing_locations(m)
        env={"np":np}
        exec(compile(m,str(EXAMPLE),"exec"),env)
        reader=env["_read_native_vtk_data_array"]
        samples=np.array((.5,-1.,0.,3.),dtype="<f4")
        raw=samples.tobytes()
        for header in (len(raw),1):  # GetFEM scalar header variants
            elem=ET.Element("DataArray",type="Float32",format="binary")
            elem.text=base64.b64encode(
                struct.pack("<I",header)+raw
            ).decode("ascii")
            np.testing.assert_array_equal(reader(elem),samples)
        elem=ET.Element("DataArray",type="Float32",format="ascii")
        elem.text="0.5 -1. 0. 3."
        np.testing.assert_array_equal(reader(elem),samples)
        elem.text="nan"
        with self.assertRaisesRegex(ValueError,"Nonfinite"):
            reader(elem)

    def test_native_triangle_and_scalar_velocity_no_resampling(self):
        funcs=[
            n for n in self.tree.body
            if isinstance(n,ast.FunctionDef)
            and n.name in (
                "_read_native_vtk_data_array","_native_wave_vtk_frame"
            )
        ]
        self.assertEqual(len(funcs),2)
        m=ast.Module(body=funcs,type_ignores=[])
        ast.fix_missing_locations(m)
        env={"np":np,"Path":Path}
        exec(compile(m,str(EXAMPLE),"exec"),env)
        reader=env["_native_wave_vtk_frame"]

        def a(values,typ):
            values=np.asarray(values,dtype=typ)
            binary=values.tobytes()
            return base64.b64encode(
                struct.pack("<I",len(binary))+binary
            ).decode("ascii")

        xy=np.array([[.6,0,0],[1.,0,0],[.6,.4,0]],dtype="<f4")
        scalar=np.array((0.,1.,-.5),dtype="<f4")
        def vtu(cell_type):
            return f"""<VTKFile><UnstructuredGrid>
<Piece NumberOfPoints="3" NumberOfCells="1">
<Points><DataArray type="Float32" format="binary">{a(xy,"<f4")}</DataArray></Points>
<Cells>
<DataArray Name="connectivity" type="Int32" format="binary">{a([0,1,2],"<i4")}</DataArray>
<DataArray Name="offsets" type="Int32" format="binary">{a([3],"<i4")}</DataArray>
<DataArray Name="types" type="Int32" format="binary">{a([cell_type],"<i4")}</DataArray>
</Cells>
<PointData><DataArray Name="p" NumberOfComponents="1"
type="Float32" format="binary">{a(scalar,"<f4")}</DataArray></PointData>
</Piece></UnstructuredGrid></VTKFile>"""
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/"native.vtu"
            path.write_text(vtu(5))
            vertices,tri,p=reader(Path(tmp),path.name)
            np.testing.assert_array_equal(vertices,xy[:,:2])
            np.testing.assert_array_equal(tri,[[0,1,2]])
            np.testing.assert_array_equal(p,scalar)
            path.write_text(vtu(10))
            with self.assertRaisesRegex(ValueError,"triangular"):
                reader(Path(tmp),path.name)

    def test_gallery_is_opt_in_and_native_truth_only(self):
        for needle in (
            'initial_wave_profile="gaussian"',
            'initial_wave_boundary="clamp"',
            'if initial_wave_profile == "gallery_packet":',
            '"np.cos(14.*np.arctan2(y,x))"',
            '"0."',
            'long_showcase=False',
            'render_heat_wave_long_showcase(',
            'hw.export_to_pv("q")',
            'hw.export_to_pv("p")',
            'hw.export_to_pv("T")',
            '"heat_wave_FE_H_long.npz"',
            '"whispering_gallery_eigenmode_certified":False',
            '"full_field_energy_localization_certified":False',
            '"no_time_interpolation":True',
            'gallery/"wave_annulus_long.gif"',
            'gallery/"wave_outer_shell_kinetic_proxy.png"',
            'ax.loglog(times[positive],H[positive]',
            'parser.add_argument("--long-showcase",action="store_true")',
        ):
            self.assertIn(needle,self.example)

    def test_long_native_verifier_python_syntax_and_fail_closed_scope(self):
        marker="python3 - \"$out\" <<'PY'\n"
        self.assertEqual(self.runner.count(marker),1)
        verifier=self.runner.split(marker,1)[1].split("\nPY\n",1)[0]
        ast.parse(verifier,filename=str(RUNNER))
        for needle in (
            "FAIL_OR_UNQUALIFIED_LONG_HEAT_WAVE_SHOWCASE",
            "PASS_LONG_NATIVE_VISUALS_NOT_WG_MODE_PROOF",
            'raise SystemExit(0 if status["status"].startswith("PASS_") else 1)',
            'whispering_gallery_mode_qualified',
            'true_FE_long_physics_gates_qualified',
            'heat_wave_FE_H_long.npz',
            'wave_annulus_long.gif',
            'png',
            'np.all(np.diff(t)>0.)',
            'native.host.log',
        ):
            self.assertIn(needle,self.runner)
        self.assertTrue(PHYSICAL.exists())
        self.assertIn('scripts/check-v130-published-immutable.py',self.runner)
        self.assertIn('python -u -m examples.heat_wave',self.runner)
        self.assertIn('initial-wave-boundary clamp',self.runner)


if __name__=="__main__":
    unittest.main()
