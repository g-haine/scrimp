# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
"""Source-only contract for new 3D Nedelec/RT Maxwell tutorial.

Does not imply native GetFEM/PETSc acceptance: only the pinned
MUMPS 1/2 runner can establish the physical energy and VTK gates.
"""
import ast
from pathlib import Path
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
            'S.FEM("E", -1, "FEM_NEDELEC(3)")',
            'S.FEM("B", 0, "RT")',
            '"FEM_RTK(3,0)"',
            '"FEM_NEDELEC(3)"',
            'assigned = np.asarray(mesh_fem.fem(), dtype=object).ravel()',
            'getfem_module.Fem(requested).char()',
            'if int(mesh_fem.qdim()) != 3 or int(mesh_fem.nbdof()) <= 0:',
        ):
            self.assertIn(needle, self.source)
        self.assertNotIn('S.FEM("E", 1, "CG")', self.source)
        self.assertNotIn('S.FEM("B", 1, "CG")', self.source)

    def test_intrinsic_fem_guard_flattens_real_getfem_layout_and_rejects_CG(self):
        """Test the isolated production guard with the array returned by MeshFem.fem."""
        node = next(
            n for n in self.tree.body
            if isinstance(n, ast.FunctionDef)
            and n.name == "assert_patched_intrinsic_fem"
        )
        module = ast.Module(body=[node], type_ignores=[])
        ast.fix_missing_locations(module)
        env = {"np": np}
        exec(compile(module, str(EXAMPLE), "exec"), env)
        guard = env["assert_patched_intrinsic_fem"]

        class FakeFem:
            def __init__(self, description):
                self.description = description

            def char(self):
                return self.description

        class FakeGetFEM:
            Fem = FakeFem

        class FakeMeshFem:
            def __init__(self, values, *, qdim=3, ndof=222):
                self.values = values
                self.dimension = qdim
                self.ndof = ndof

            def fem(self):
                return self.values

            def qdim(self):
                return self.dimension

            def nbdof(self):
                return self.ndof

        expected = "FEM_NEDELEC(3)"
        nested = np.array([[FakeFem(expected), FakeFem(expected)]],
                          dtype=object)
        self.assertEqual(guard(FakeMeshFem(nested), expected, FakeGetFEM),
                         {expected})
        for values in (
            np.array([[FakeFem("FEM_PK(3,1)")]], dtype=object),
            np.array([[FakeFem(expected), FakeFem("FEM_PK(3,1)")]],
                     dtype=object),
            np.array([], dtype=object),
            np.array([[object()]], dtype=object),
        ):
            with self.subTest(values=values.shape), self.assertRaises(RuntimeError):
                guard(FakeMeshFem(values), expected, FakeGetFEM)
        for kwargs in ({"qdim": 1}, {"ndof": 0}):
            with self.subTest(kwargs=kwargs), self.assertRaises(RuntimeError):
                guard(FakeMeshFem(nested, **kwargs), expected, FakeGetFEM)

    def test_embedded_native_verifier_is_valid_python(self):
        """Shell syntax checks miss invalid Python embedded in a heredoc."""
        sentinel = 'python3 - "$out" <<\'PY\'\n'
        self.assertEqual(self.runner.count(sentinel), 1)
        verifier = self.runner.split(sentinel, 1)[1].split("\nPY\n", 1)[0]
        ast.parse(verifier, filename=str(RUNNER))
        self.assertIn("from pathlib import Path", verifier)
        self.assertIn("with zipfile.ZipFile", verifier)
        self.assertIn("np.isin(cell_types, (10,24))", verifier)


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
            'getfem_module.Fem(requested)',
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
        self.assertIn("H(curl)", self.doc)
        self.assertIn("H(div)", self.doc)
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
