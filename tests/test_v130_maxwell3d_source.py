# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
"""Source-only contract for new 3D Nedelec/RT Maxwell tutorial.

Does not imply native GetFEM/PETSc acceptance: only the pinned
MUMPS 1/2 runner can establish the physical energy and VTK gates.
"""
import ast
import base64
import struct
import tempfile
from pathlib import Path
from xml.etree import ElementTree as ET
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
EXAMPLE = ROOT / "examples" / "maxwell_3d.py"
RUNNER = ROOT / "scripts" / "run-v130-maxwell3d-pilot.sh"
TUTORIAL = ROOT / "docs" / "source" / "examples" / "maxwell_3d.rst"


class Maxwell3DSourceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.source = EXAMPLE.read_text(encoding="utf-8")
        cls.runner = RUNNER.read_text(encoding="utf-8")
        cls.doc = TUTORIAL.read_text(encoding="utf-8")
        cls.tree = ast.parse(cls.source)

    def test_actual_3d_geometry_and_both_patched_intrinsic_families(self):
        for needle in (
            'S.Domain("Ball", {"R": 1., "h": float(mesh_size)})',
            'fem_E = S.FEM("E", -1, "FEM_NEDELEC(3)")',
            'fem_B = S.FEM("B", 0, "RT")',
            'GetFEM FEM_RTK(3,0)',
            '"FEM_NEDELEC(3)"',
            'actual_type = fem.get_type()',
            'actual_order = fem.get_order()',
            'mesh_fem = fem.get_fem()',
            'if qdim != 3 or ndof <= 0:',
        ):
            self.assertIn(needle, self.source)
        self.assertNotIn('S.FEM("E", 1, "CG")', self.source)
        self.assertNotIn('S.FEM("B", 1, "CG")', self.source)

    def test_intrinsic_fem_guard_checks_registered_scrimp_type_and_order(self):
        """Test the exact FEM.set_fem() CG fallback policy, not wrapper arrays."""
        node = next(
            n for n in self.tree.body
            if isinstance(n, ast.FunctionDef)
            and n.name == "assert_patched_intrinsic_fem"
        )
        module = ast.Module(body=[node], type_ignores=[])
        ast.fix_missing_locations(module)
        env = {}
        exec(compile(module, str(EXAMPLE), "exec"), env)
        guard = env["assert_patched_intrinsic_fem"]

        class FakeMeshFem:
            def __init__(self, qdim=3, ndof=222):
                self.dimension = qdim
                self.ndof = ndof

            def qdim(self):
                return self.dimension

            def nbdof(self):
                return self.ndof

        class FakeScrimpFEM:
            def __init__(self, family="FEM_NEDELEC(3)", order=-1,
                         native=None):
                self.family = family
                self.order = order
                self.native = FakeMeshFem() if native is None else native

            def get_type(self):
                return self.family

            def get_order(self):
                return self.order

            def get_fem(self):
                return self.native

        self.assertEqual(guard(FakeScrimpFEM(), "FEM_NEDELEC(3)", -1), 222)
        self.assertEqual(
            guard(FakeScrimpFEM("RT", 0, FakeMeshFem(ndof=275)), "RT", 0),
            275,
        )
        for family, order in (
            ("CG", 1),  # SCRIMP FEM.set_fem's exact fallback mutation
            ("FEM_NEDELEC(3)", 1),
            ("RT", 0),
            ("FEM_BDMK(3,1)", -1),
        ):
            with self.subTest(family=family, order=order):
                with self.assertRaises(RuntimeError):
                    guard(FakeScrimpFEM(family, order), "FEM_NEDELEC(3)", -1)
        for qdim, ndof in ((1,222), (3,0), (0,0)):
            with self.subTest(qdim=qdim, ndof=ndof):
                with self.assertRaises(RuntimeError):
                    guard(FakeScrimpFEM(native=FakeMeshFem(qdim, ndof)),
                          "FEM_NEDELEC(3)", -1)

    def test_native_guard_uses_wrapper_metadata_not_getfem_array_layout(self):
        guard = next(
            n for n in self.tree.body
            if isinstance(n, ast.FunctionDef)
            and n.name == "assert_patched_intrinsic_fem"
        )
        # AST: do not mistake documentation/comments for executable code.
        self.assertFalse(any(
            isinstance(node, ast.Call)
            and isinstance(node.func, ast.Attribute)
            and node.func.attr == "fem"
            for node in ast.walk(guard)
        ))
        self.assertIn('model.add_FEM(fem)', self.source)
        self.assertIn('assert_patched_intrinsic_fem(fem, expected_type, expected_order)',
                      self.source)

    def test_embedded_native_verifier_is_valid_python(self):
        """Shell syntax checks miss invalid Python embedded in a heredoc."""
        sentinel = 'python3 - "$out" <<\'PY\'\n'
        self.assertEqual(self.runner.count(sentinel), 1)
        verifier = self.runner.split(sentinel, 1)[1].split("\nPY\n", 1)[0]
        ast.parse(verifier, filename=str(RUNNER))
        self.assertIn("from pathlib import Path", verifier)
        self.assertIn("with zipfile.ZipFile", verifier)
        self.assertIn("np.isin(cell_types, (10,24))", verifier)



    def test_native_vtk_parser_and_3d_time_series_without_vtk_dependency(self):
        """Reader in the executable accepts real GetFEM-style inline arrays.

        This is a *protocol fixture*, not fake Maxwell physical evidence.
        """
        names = {"_vtk_native_array", "_native_vector_history"}
        nodes = [
            n for n in self.tree.body
            if isinstance(n, ast.FunctionDef) and n.name in names
        ]
        self.assertEqual({n.name for n in nodes}, names)
        module = ast.Module(body=nodes, type_ignores=[])
        ast.fix_missing_locations(module)
        env = {"np": np, "Path": Path}
        exec(compile(module, str(EXAMPLE), "exec"), env)
        read_series = env["_native_vector_history"]

        points = np.array(((0,0,0),(1,0,0),(0,1,0),(0,0,1)),
                          dtype="<f4")
        values = np.zeros((4,3), dtype="<f4")
        fields = []
        def packed(data, *, field_header=False):
            data = np.asarray(data)
            binary = data.tobytes()
            header = 4 if field_header else len(binary)
            return base64.b64encode(
                struct.pack("<I", header) + binary
            ).decode("ascii")

        def vtu(vectors, cell_type=10):
            return f"""<VTKFile type="UnstructuredGrid">
<UnstructuredGrid><Piece NumberOfPoints="4" NumberOfCells="1">
<Points><DataArray type="Float32" format="binary"
NumberOfComponents="3">{packed(points)}</DataArray></Points>
<Cells>
<DataArray Name="connectivity" type="Int32" format="binary">{packed(np.array([0,1,2,3],dtype="<i4"))}</DataArray>
<DataArray Name="offsets" type="Int32" format="binary">{packed(np.array([4],dtype="<i4"))}</DataArray>
<DataArray Name="types" type="Int32" format="binary">{packed(np.array([cell_type],dtype="<i4"))}</DataArray>
</Cells>
<PointData><DataArray Name="E" NumberOfComponents="3" type="Float32"
format="binary">{packed(vectors,field_header=True)}</DataArray></PointData>
</Piece></UnstructuredGrid></VTKFile>"""

        with tempfile.TemporaryDirectory() as tmp:
            subdir = Path(tmp)/"pv"/"E"
            subdir.mkdir(parents=True)
            for index in (0,1):
                vectors = values + index * np.array((1,2,3),dtype="<f4")
                (subdir/f"E_{index:05d}.vtu").write_text(
                    vtu(vectors), encoding="utf-8"
                )
            (subdir/"E.pvd").write_text(
                """<VTKFile type="Collection"><Collection>
<DataSet timestep="0.0" file="E_00000.vtu"/>
<DataSet timestep="0.02" file="E_00001.vtu"/>
</Collection></VTKFile>""", encoding="utf-8"
            )
            coords, tetra, E = read_series(
                Path(tmp), "E", np.array([0., 0.02])
            )
            np.testing.assert_array_equal(coords, points)
            np.testing.assert_array_equal(tetra, [[0,1,2,3]])
            np.testing.assert_array_equal(E[0], values)
            np.testing.assert_array_equal(
                E[1], values+np.array((1,2,3),dtype="<f4")
            )
            with self.assertRaisesRegex(ValueError,"time mismatch"):
                read_series(Path(tmp), "E", np.array([0.,0.04]))
            (subdir/"E_00001.vtu").write_text(
                vtu(values,cell_type=5),encoding="utf-8"
            )
            with self.assertRaisesRegex(ValueError,"non-tetrahedral"):
                read_series(Path(tmp), "E", np.array([0.,0.02]))

    def test_export_automatically_generates_pictures_and_animation_in_rank0(self):
        for marker in (
            "def render_maxwell_visuals(output, observed, state_sha256):",
            'folder = output / "visuals"',
            'maxwell3d_echange_energie.png',
            'maxwell3d_echange_energie.svg',
            'maxwell3d_vecteurs_E_tfinal.png',
            'maxwell3d_vecteurs_B_tfinal.png',
            'maxwell3d_B_3D_animation.gif',
            'original_FE_state_SHA256',
            '"no_synthetic_field_or_time_interpolation": True',
            'if comm.rank == 0:',
            'metadata["visuals"] = render_maxwell_visuals(',
            'system.export_to_pv(name, path=str(args.output), t="All")',
        ):
            self.assertIn(marker, self.source)
        # Rendering must be opt-in on --export, AFTER the GetFEM
        # collective VTU writer. Do not create a second MPI exporter.
        main_src = self.source.split('def main():',1)[1]
        self.assertLess(
            main_src.index('system.export_to_pv(name,'),
            main_src.index('metadata["visuals"] = render_maxwell_visuals('),
        )
        for needle in (
            'maxwell3d_B_3D_animation.gif',
            'original_FE_state_SHA256',
            'gif.n_frames!=11',
            'prov.get("MPI_writer_rank")!=0',
            'rank0_meta.get("visuals")!=prov',
        ):
            self.assertIn(needle, self.runner)

    def test_true_3d_curl_orientation(self):
        for part in (
            "Grad(v)(3,2)-Grad(v)(2,3)",
            "Grad(v)(1,3)-Grad(v)(3,1)",
            "Grad(v)(2,1)-Grad(v)(1,2)",
            '"B.Curl3(Test_E)"',
            '"-Curl3(E).Test_B"',
        ):
            self.assertIn(part, self.source)

    def test_lossless_skew_maxwell_without_fake_boundary_input(self):
        for needle in (
            'S.Brick("M_E", "E.Test_E", [1], dt=True, position="flow")',
            'S.Brick("M_B", "B.Test_B", [1], dt=True, position="flow")',
            'S.Brick("C^T", "B.Curl3(Test_E)", [1], position="effort")',
            'S.Brick("-C", "-Curl3(E).Test_B", [1], position="effort")',
            'S.Term("Electric", "0.5*E.E", [1])',
            'S.Term("Magnetic", "0.5*B.B", [1])',
            'model.set_initial_value("E", "[-y, x, 0.]")',
            'model.set_initial_value("B", "[0., 0., 0.]")',
            'init_step=False',
        ):
            self.assertIn(needle, self.source)
        # Inspect the actual call AST, not comments such as "no pseudo-TS".
        scheme_calls = [
            node for node in ast.walk(self.tree)
            if isinstance(node, ast.Call)
            and isinstance(node.func, ast.Attribute)
            and isinstance(node.func.value, ast.Name)
            and node.func.value.id == "model"
            and node.func.attr == "set_time_scheme"
        ]
        self.assertEqual(len(scheme_calls), 1)
        keywords = {key.arg: key.value for key in scheme_calls[0].keywords}
        self.assertIn("init_step", keywords)
        self.assertIsInstance(keywords["init_step"], ast.Constant)
        self.assertIs(keywords["init_step"].value, False)
        self.assertNotIn("init_step_ts_type", keywords)
        for name, expected in (("ts_type", "cn"),
                               ("ts_adapt_type", "none")):
            self.assertIn(name, keywords)
            self.assertIsInstance(keywords[name], ast.Constant)
            self.assertEqual(keywords[name].value, expected)
        self.assertNotIn("set_control(", self.source)
        self.assertNotIn("add_control_port(", self.source)

    def test_physical_energy_and_magnetic_excitation_gates(self):
        for gate in (
            '"positive_initial_em_energy"',
            '"electric_and_magnetic_exchange"',
            '"independent_fe_energy_matches"',
            '"every_CN_step_conserves_energy"',
            '"cumulative_em_energy_conserved"',
            '"pure_differential_Maxwell_no_artificial_DAE"',
            '"magnetic_gauss_law"',
        ):
            self.assertIn(gate, self.source)
        self.assertIn('Trace(Grad(B))*Trace(Grad(B))', self.source)

    def test_native_pilot_exports_original_3d_vtk_and_checks_mpi(self):
        for needle in (
            'for ranks in 1 2;',
            'mumps_mpi1',
            'mumps_mpi2',
            'actual_type = fem.get_type()',
            'python3 scripts/check-v130-published-immutable.py',
            'SCRIMP_GETFEM_MATRIX_MODE=distributed',
            'zipfile.ZipFile',
            'f\"rank{rank}.npz\"',
            'summary.json',
            'pvd=folder/name/"pv"/field/f"{field}.pvd"',
            'PASS_MAXWELL3D_MUMPS_MPI1_MPI2',
        ):
            self.assertIn(needle, self.source+self.runner)

    def test_tutorial_reuses_same_canonical_source(self):
        self.assertIn("Nédélec", self.doc)
        self.assertIn("Raviart", self.doc)
        # The tutorial uses Sphinx/LaTeX H(\mathrm{curl};\Omega) notation,
        # not the informal plaintext shorthands H(curl) and H(div).
        self.assertIn(r"H(\mathrm{curl};\Omega)", self.doc)
        self.assertIn(r"H(\mathrm{div};\Omega)", self.doc)
        self.assertIn("Maxwell", self.doc)
        self.assertIn("H(t)", self.doc)
        self.assertEqual(
            self.doc.count(".. literalinclude:: ../../../examples/maxwell_3d.py"),
            7,
        )
        for part in ("M3D-01", "M3D-02", "M3D-03",
                     "M3D-04", "M3D-05", "M3D-06", "M3D-07"):
            self.assertIn(part, self.source)
            self.assertIn(f":start-after: # TUTORIAL {part}:", self.doc)


if __name__ == "__main__":
    unittest.main()
