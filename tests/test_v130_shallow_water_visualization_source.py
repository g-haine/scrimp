# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Source and GetFEM-like VTU parser fixtures for native SWE visuals.

Artificial XML arrays in isolated parser unit tests are NOT physical
proof. Only the real Docker visual driver can qualify GetFEM h/p
archives, and the separate native original BDF2 physics proof stands.
"""
import ast
import base64
import importlib.util
from pathlib import Path
import struct
import tempfile
import unittest
from xml.etree import ElementTree as ET

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
CANONICAL=ROOT/"examples/shallow_water.py"
HELPER=CANONICAL
RUNNER=ROOT/"scripts/run-v130-shallow-water-showcase.sh"


def array_bytes(a):
    arr=np.asarray(a)
    payload=arr.tobytes()
    return base64.b64encode(
        struct.pack("<I",len(payload))+payload
    ).decode("ascii")


def fixture_vtu(field,vector=False,badtype=False):
    xyz=np.array([[0,0,0],[1,0,0],[0,.5,0]],dtype="<f4")
    cell=np.array([0,1,2],dtype="<i4")
    off=np.array([3],dtype="<i4")
    typ=np.array([10 if badtype else 5],dtype="<i4")
    if vector:
        values=np.array([[0,0,0],[10,20,0],[-4,6,0]],dtype="<f4")
        components=' NumberOfComponents="3"'
    else:
        values=np.array([3,2.8,2.6],dtype="<f4")
        components=""
    return (
        '<VTKFile type="UnstructuredGrid"><UnstructuredGrid>'
        '<Piece NumberOfPoints="3" NumberOfCells="1">'
        '<Points><DataArray type="Float32" format="binary">'
        +array_bytes(xyz)+
        '</DataArray></Points><Cells>'
        '<DataArray Name="connectivity" type="Int32" format="binary">'
        +array_bytes(cell)+
        '</DataArray>'
        '<DataArray Name="offsets" type="Int32" format="binary">'
        +array_bytes(off)+
        '</DataArray>'
        '<DataArray Name="types" type="Int32" format="binary">'
        +array_bytes(typ)+
        '</DataArray></Cells>'
        '<PointData><DataArray Name="'+field+'" type="Float32" '
        'format="binary"'+components+'>'
        +array_bytes(values)+
        '</DataArray></PointData></Piece></UnstructuredGrid></VTKFile>'
    )


def fixture_native_quadratic_vtu(*, corrupt_mids=False,
                                 corrupt_field_header=False):
    """GetFEM 5.5 genuine shape: VTK22 CG2, six original nodes.

    Field header is float32(point count), unlike the standard UInt32
    byte header used by coordinates and connectivity in this exporter.
    """
    xyz=np.array([
        [0.,0.,0.],[1.,0.,0.],[0.,.5,0.],
        [.5,.0,0.],[.5,.25,0.],[0.,.25,0.],
    ],dtype="<f4")
    if corrupt_mids:
        xyz[3,0]=.87
    values=np.array([3.,2.8,2.6,2.9,2.7,2.8],dtype="<f4")
    header=float(len(values))+(1 if corrupt_field_header else 0)
    field=base64.b64encode(
        struct.pack("<f",header)+values.tobytes()
    ).decode("ascii")
    return (
        '<VTKFile type="UnstructuredGrid"><UnstructuredGrid>'
        '<Piece NumberOfPoints="6" NumberOfCells="1">'
        '<Points><DataArray type="Float32" format="binary">'
        +array_bytes(xyz)+
        '</DataArray></Points><Cells>'
        '<DataArray Name="connectivity" type="Int32" format="binary">'
        +array_bytes(np.arange(6,dtype="<i4"))+
        '</DataArray>'
        '<DataArray Name="offsets" type="Int32" format="binary">'
        +array_bytes(np.array([6],dtype="<i4"))+
        '</DataArray>'
        '<DataArray Name="types" type="Int32" format="binary">'
        +array_bytes(np.array([22],dtype="<i4"))+
        '</DataArray></Cells>'
        '<PointData><DataArray Name="h" type="Float32" format="binary">'
        +field+
        '</DataArray></PointData></Piece></UnstructuredGrid></VTKFile>'
    )


def fixture_native_DG1_vector_vtu():
    xyz=np.array([[0.,0.,0.],[1.,0.,0.],[0.,.5,0.]],dtype="<f4")
    val=np.array([[0.,0.,0.],[10.,20.,0.],[-4.,6.,0.]],dtype="<f4")
    field=base64.b64encode(
        struct.pack("<f",float(3*2))+val.tobytes()
    ).decode("ascii")
    return (
        '<VTKFile type="UnstructuredGrid"><UnstructuredGrid>'
        '<Piece NumberOfPoints="3" NumberOfCells="1">'
        '<Points><DataArray type="Float32" format="binary">'
        +array_bytes(xyz)+
        '</DataArray></Points><Cells>'
        '<DataArray Name="connectivity" type="Int32" format="binary">'
        +array_bytes(np.array([0,1,2],dtype="<i4"))+
        '</DataArray>'
        '<DataArray Name="offsets" type="Int32" format="binary">'
        +array_bytes(np.array([3],dtype="<i4"))+
        '</DataArray>'
        '<DataArray Name="types" type="Int32" format="binary">'
        +array_bytes(np.array([5],dtype="<i4"))+
        '</DataArray></Cells>'
        '<PointData><DataArray Name="p" NumberOfComponents="3"'
        ' type="Float32" format="binary">'
        +field+
        '</DataArray></PointData></Piece></UnstructuredGrid></VTKFile>'
    )


class ShallowWaterNativeVisualizationSource(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model=CANONICAL.read_text(encoding="utf-8")
        cls.helper=HELPER.read_text(encoding="utf-8")
        cls.runner=RUNNER.read_text(encoding="utf-8")
        for path in (CANONICAL,HELPER):
            ast.parse(path.read_text(),filename=str(path))
        spec=importlib.util.spec_from_file_location("original_swe_vtk",HELPER)
        cls.module=importlib.util.module_from_spec(spec)
        spec.loader.exec_module(cls.module)

    def test_single_file_example_contract_and_independent_offline_tools(self):
        """The user-facing tutorial must NEVER split its code into examples/."""
        self.assertEqual(HELPER,CANONICAL)
        self.assertFalse(
            (ROOT/"examples/shallow_water_visualization.py").exists()
        )
        module=ast.parse(self.model)
        expected={
            "shallow_water_eq",
            "vtk_array",
            "native_triangular_field",
            "true_pvd_files",
            "render_original_shallow_water_spatial_VTU",
            "render_native_shallow_water",
        }
        functions={
            n.name for n in module.body if isinstance(n,ast.FunctionDef)
        }
        self.assertTrue(expected.issubset(functions))
        self.assertNotIn("shallow_water_visualization",self.model)
        for utility in (
            "scripts/recover-v130-shallow-water-showcase.py",
            "scripts/verify-v130-shallow-water-showcase.py",
        ):
            source=(ROOT/utility).read_text(encoding="utf-8")
            self.assertIn("from examples.shallow_water import (",source)
            self.assertNotIn("shallow_water_visualization",source)
        self.assertNotIn(
            "shallow_water_visualization.py",self.runner
        )
        # No GetFEM/PETSc import is forced by standalone VTK helpers.
        # True native simulation imports SCRIMP only when invoked.
        self.assertNotIn("import scrimp as S",
                         self.model.split("def shallow_water_eq(",1)[0])
        self.assertIn("    import scrimp as S",self.model)

    def test_running_example_without_options_produces_all_visuals(self):
        """Execute the REAL CLI AST with a stub solver, not a string test.

        Never run GetFEM in a lightweight CI unit test. The actual
        rendering pipeline is exercised separately with VTK fixtures,
        and the native driver validates the real PDE and GetFEM fields.
        """
        import sys
        from unittest.mock import patch
        parsed=ast.parse(self.model,filename=str(CANONICAL))
        main=next(
            node for node in parsed.body
            if isinstance(node,ast.If)
            and isinstance(node.test,ast.Compare)
            and isinstance(node.test.left,ast.Name)
            and node.test.left.id=="__name__"
        )
        body=ast.Module(body=[main],type_ignores=[])
        ast.fix_missing_locations(body)
        actual=[]
        with patch.object(sys,"argv",["shallow_water.py"]):
            exec(compile(body,str(CANONICAL),"exec"),{
                "__name__":"__main__",
                "shallow_water_eq":lambda **kwargs:actual.append(kwargs),
            })
        self.assertEqual(len(actual),1)
        self.assertEqual(actual[0],{
            "factor_solver":None,
            "final_time":0.12,
            "time_step":0.0001,
            "save_step":0.002,
            "showcase":True,
            "showcase_output":None,
        })
        # The switch is explicit, but never required for normal users.
        self.assertIn("parser.set_defaults(showcase=True)",self.model)
        self.assertIn('swe.shallow_water_visuals=render_native_shallow_water(',self.model)
        self.assertIn("saint_venant_surface_vitesse.gif",self.model)
        self.assertIn("saint_venant_H_masse_bilan.png",self.model)
        self.assertIn('Path.cwd()/"outputs"/"showcase"/"shallow_water"',self.model)
        self.assertNotIn("shallow_water_visualization",self.model)

    def test_opt_in_legacy_and_custom_visual_cli_are_reproducible(self):
        import contextlib
        import io
        import sys
        from unittest.mock import patch
        parsed=ast.parse(self.model,filename=str(CANONICAL))
        main=next(
            n for n in parsed.body
            if isinstance(n,ast.If)
            and isinstance(n.test,ast.Compare)
            and isinstance(n.test.left,ast.Name)
            and n.test.left.id=="__name__"
        )
        body=ast.Module(body=[main],type_ignores=[])
        ast.fix_missing_locations(body)
        def dispatch(args):
            calls=[]
            with patch.object(sys,"argv",["shallow_water.py",*args]):
                exec(compile(body,str(CANONICAL),"exec"),{
                    "__name__":"__main__",
                    "shallow_water_eq":lambda **kwargs:calls.append(kwargs),
                })
            self.assertEqual(len(calls),1)
            return calls[0]
        legacy=dispatch(["--no-showcase"])
        self.assertFalse(legacy["showcase"])
        self.assertEqual(legacy["final_time"],.5)
        self.assertEqual(legacy["save_step"],.01)
        custom=dispatch(["--showcase","--final-time","0.2",
                         "--save-step","0.004",
                         "--showcase-output","/tmp/swe-vis"])
        self.assertTrue(custom["showcase"])
        self.assertEqual(custom["final_time"],.2)
        self.assertEqual(custom["save_step"],.004)
        self.assertEqual(custom["showcase_output"],"/tmp/swe-vis")
        with patch.object(sys,"argv",[
            "shallow_water.py","--no-showcase",
            "--showcase-output","/tmp/unused",
        ]):
            with contextlib.redirect_stderr(io.StringIO()):
                with self.assertRaises(SystemExit) as err:
                    exec(compile(body,str(CANONICAL),"exec"),{
                        "__name__":"__main__","shallow_water_eq":lambda **kw:None,
                    })
            self.assertEqual(err.exception.code,2)

    def test_original_pseudo_default_and_optional_visualization(self):
        for needle in (
            'showcase=False, showcase_output=None',
            'init_step_ts_type="pseudo"',
            'if showcase:',
            'init_step_ts_type="nonlinear_affine"',
            '"nonlinear_affine"',
            'swe.export_to_pv("h",path=str(output))',
            'swe.export_to_pv("p",path=str(output))',
            'swe.shallow_water_visuals=render_native_shallow_water(',
            'swe.export_to_pv("h")',
            'ts_bdf_order=2',
            'S.Brick("M_p", "h * p . Test_p"',
            'explicit=True, position="effort"',
            'modes.add_argument(\n        "--showcase", dest="showcase", action="store_true"',
        ):
            self.assertIn(needle,self.model,msg=needle)
        self.assertNotIn('viscous_swe_2d_mpi.py',self.model)

    def test_native_GetFEM_binary_fields_and_distinct_h_p_spaces(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp)
            for field,vector in (("h",False),("p",True)):
                path=root/(field+".vtu")
                path.write_text(fixture_vtu(field,vector=vector))
                xy,tri,data=self.module.native_triangular_field(path,field)
                self.assertEqual(xy.shape,(3,2))
                np.testing.assert_array_equal(tri,[[0,1,2]])
                if vector:
                    self.assertEqual(data.shape,(3,2))
                    np.testing.assert_allclose(data[1],[10,20])
                else:
                    self.assertEqual(data.shape,(3,))
                    np.testing.assert_allclose(data,[3,2.8,2.6])
            path=root/"h.vtu"
            path.write_text(fixture_vtu("h",badtype=True))
            with self.assertRaisesRegex(ValueError,"triangles"):
                self.module.native_triangular_field(path,"h")
            path.write_text(fixture_vtu("p",vector=True))
            with self.assertRaisesRegex(ValueError,"scalar"):
                self.module.native_triangular_field(path,"h")

    def test_GENUINE_GetFEM_CG2_quadratic_and_DG1_field_headers(self):
        # Source-14181a94 original VTU evidence: 248 quadratic
        # VTK22 cells with 6 nodes for h, and 248 VTK5 cells
        # with 3 nodes for the independent discontinuous p space.
        with tempfile.TemporaryDirectory() as tmp:
            out=Path(tmp)/"native.vtu"
            out.write_text(fixture_native_quadratic_vtu())
            xyz,tri,h=self.module.native_triangular_field(out,"h")
            self.assertEqual(xyz.shape,(6,2))
            self.assertEqual(h.shape,(6,))
            np.testing.assert_allclose(h,[3,2.8,2.6,2.9,2.7,2.8])
            self.assertEqual(tri.shape,(4,3))
            np.testing.assert_array_equal(tri,[
                [0,3,5],[3,1,4],[5,4,2],[3,4,5]
            ])
            # Deliberately corrupted geometry and binary count MUST
            # fail, not silently show a physically different surface.
            out.write_text(fixture_native_quadratic_vtu(corrupt_mids=True))
            with self.assertRaisesRegex(ValueError,"midside"):
                self.module.native_triangular_field(out,"h")
            out.write_text(
                fixture_native_quadratic_vtu(corrupt_field_header=True)
            )
            with self.assertRaisesRegex(ValueError,"field header"):
                self.module.native_triangular_field(out,"h")
            out.write_text(fixture_native_DG1_vector_vtu())
            _,p_tri,p=self.module.native_triangular_field(out,"p")
            np.testing.assert_array_equal(p_tri,[[0,1,2]])
            np.testing.assert_array_equal(p[1],[10.,20.])

    def test_original_spatial_zip_recovery_is_separate_from_physical_gate(self):
        path=ROOT/"scripts/recover-v130-shallow-water-showcase.py"
        source=path.read_text(encoding="utf-8")
        ast.parse(source,filename=str(path))
        for text in (
            "PASS_RESCUED_NATIVE_SWE_SPATIAL_VTK_ONLY",
            'old["status"]=="FAIL_OR_UNQUALIFIED_NATIVE_SWE_VISUALS"',
            "original_FE_Hamiltonian_quadratures_recovered",
            "original_FE_mass_quadratures_recovered",
            "render_original_shallow_water_spatial_VTU(",
            "Steps: 1200 (0 rejected, 0 Nonlinear solver failures)",
            "mesh_h_original_quadratic_cells",
            "no_time_interpolation",
        ):
            self.assertTrue(text in source or text in self.helper,msg=text)
        self.assertIn("def render_original_shallow_water_spatial_VTU(",self.helper)
        self.assertNotIn("interpolate",self.helper.split(
            "def render_original_shallow_water_spatial_VTU(",1
        )[1].split("def render_native_shallow_water(",1)[0].lower())

    def test_pvd_strict_FE_clock_xml_attributes(self):
        with tempfile.TemporaryDirectory() as td:
            base=Path(td)
            pvd=base/"pv"/"h"
            pvd.mkdir(parents=True)
            (pvd/"h_00000.vtu").write_text(fixture_vtu("h"))
            (pvd/"h_00001.vtu").write_text(fixture_vtu("h"))
            (pvd/"h.pvd").write_text(
                '<VTKFile><Collection>'
                '<DataSet timestep="0.0" file="h_00000.vtu"/>'
                '<DataSet timestep="0.001" file="h_00001.vtu"/>'
                '</Collection></VTKFile>'
            )
            true=[0.,.001]
            files=self.module.true_pvd_files(base,"h",true)
            self.assertEqual(len(files),2)
            with self.assertRaisesRegex(ValueError,"time differs"):
                self.module.true_pvd_files(base,"h",[0.,.005])
            with self.assertRaisesRegex(ValueError,"lengths disagree"):
                self.module.true_pvd_files(base,"h",[0.])
            (pvd/"h.pvd").write_text(
                '<VTKFile><Collection>'
                '<DataSet timestep="0" file="../escape.vtu"/>'
                '<DataSet timestep=".001" file="h_00001.vtu"/>'
                '</Collection></VTKFile>'
            )
            with self.assertRaisesRegex(ValueError,"Unsafe"):
                self.module.true_pvd_files(base,"h",true)

    def test_renders_only_genuine_VTU_v_with_correct_units_and_fixed_scales(self):
        for token in (
            'get_quantity("h",region=1)',
            'get_quantity("0.5*h*p.p/rho",region=1)',
            'get_quantity("0.5*rho*g*h*h",region=1)',
            'true_pvd_files(folder,"h",t)',
            'true_pvd_files(folder,"p",t)',
            'native_triangular_field(h_files[int(idx)],"h")',
            'native_triangular_field(p_files[int(idx)],"p")',
            'vel=momentum[ptri].mean(axis=1)/rho',
            'fixed_height_color_range',
            'fixed_velocity_color_range',
            '"no_time_interpolation":True',
            '"viscous_rotating_research_included":False',
            '"H_exact_conservation_claimed":False',
            'saint_venant_surface_vitesse.gif',
            'saint_venant_H_masse_bilan.png',
            'physical_FE_traces.npz',
        ):
            self.assertIn(token,self.helper,msg=token)

    def test_isolated_renderer_smoke_uses_only_fixture_VTU(self):
        # This is a graphical API smoke test with explicitly synthetic
        # 5-frame XML fixtures. It is NOT native GetFEM scientific proof.
        # Catch backend/3D surface/quiver errors before the expensive
        # native BDF2 example is executed on the pinned Docker image.
        import matplotlib
        matplotlib.use("Agg")
        from PIL import Image
        with tempfile.TemporaryDirectory() as temp:
            base=Path(temp)
            times=np.arange(5,dtype=float)*.001
            for name,vector in (("h",False),("p",True)):
                directory=base/"pv"/name
                directory.mkdir(parents=True)
                records=[]
                for i,t in enumerate(times):
                    filename=f"{name}_{i:05d}.vtu"
                    (directory/filename).write_text(
                        fixture_vtu(name,vector=vector)
                    )
                    records.append(
                        f'<DataSet timestep="{t}" file="{filename}"/>'
                    )
                (directory/(name+".pvd")).write_text(
                    "<VTKFile><Collection>"+ "".join(records)
                    +"</Collection></VTKFile>"
                )
            class FakePurePlottingData:
                solution={"t":list(times)}
                def get_Hamiltonian(self):
                    return 5.-.01*times
                def get_quantity(self,expression,region):
                    if expression=="h":
                        return 2.+np.zeros_like(times)
                    if expression=="0.5*h*p.p/rho":
                        return .2+.01*times
                    if expression=="0.5*rho*g*h*h":
                        return self.get_Hamiltonian()-(.2+.01*times)
                    if expression.startswith("U_") and "*Y_" in expression:
                        return np.zeros_like(times)
                    raise ValueError("Unexpected physical quantity fixture")
            output=self.module.render_native_shallow_water(
                FakePurePlottingData(),base,rho=1000.,max_gif_frames=5
            )
            self.assertEqual(output["gif_original_FE_frames"],5)
            self.assertEqual(output["true_FE_states"],5)
            self.assertFalse(output["H_exact_conservation_claimed"])
            images=sorted((base/"figures").glob("saint_venant_champs_t*.png"))
            self.assertEqual(len(images),4)
            self.assertTrue((base/"figures"/"saint_venant_H_masse_bilan.png").is_file())
            self.assertTrue((base/"physical_FE_traces.npz").is_file())
            with Image.open(base/"figures"/"saint_venant_surface_vitesse.gif") as gif:
                self.assertEqual(gif.n_frames,5)

    def test_launcher_is_native_original_and_fail_closed(self):
        for token in (
            'scrimp-mpi:v1.3.0-rc3',
            '1d6d0c90115213298047b28237a7f38087563a35',
            'scripts/check-v130-published-immutable.py',
            'python -u -m examples.shallow_water',
            '--showcase-output /artifacts',
            'examples/shallow_water.py',
            'scripts/verify-v130-shallow-water-showcase.py',
            'with zipfile.ZipFile(archive,"x"',
        ):
            self.assertIn(token,self.runner,msg=token)
        checker=(ROOT/"scripts/verify-v130-shallow-water-showcase.py").read_text()
        self.assertIn('"FAIL_OR_UNQUALIFIED_NATIVE_SWE_VISUALS"',checker)
        self.assertIn('"PASS_NATIVE_INVISCID_SWE_VISUALS"',checker)
        ast.parse(checker)
        # CI / native driver must exercise the default direct example:
        # --showcase is not a prerequisite for figures or GIF.
        self.assertNotIn(
            '--showcase --showcase-output /artifacts',self.runner
        )
        self.assertIn('parser.set_defaults(showcase=True)',self.model)
        self.assertNotIn("viscous_swe_2d_mpi.py",self.runner)


if __name__=="__main__":
    unittest.main()
