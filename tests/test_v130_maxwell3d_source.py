# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
"""Source-only contract for new 3D Nedelec/RT Maxwell tutorial.

Does not imply native GetFEM/PETSc acceptance: only the pinned
MUMPS 1/2 runner can establish the physical energy and VTK gates.
"""
import ast
from pathlib import Path
import unittest

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
            'f.char() for f in model.ports[name].get_fem().fem()',
            'gf.Fem(requested).char()',
            'if fem.qdim() != 3 or fem.nbdof() == 0:',
        ):
            self.assertIn(needle, self.source)
        self.assertNotIn('S.FEM("E", 1, "CG")', self.source)
        self.assertNotIn('S.FEM("B", 1, "CG")', self.source)

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
        self.assertNotIn("pseudo", self.source.split('model.set_time_scheme(', 1)[1].split('model.solve()',1)[0])
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
        ):
            self.assertIn(gate, self.source)
        self.assertIn('Trace(Grad(B))*Trace(Grad(B))', self.source)

    def test_native_pilot_exports_original_3d_vtk_and_checks_mpi(self):
        for needle in (
            'for ranks in 1 2;',
            'mumps_mpi1',
            'mumps_mpi2',
            'gf.Fem(requested)',
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
        for part in ("M3D-01", "M3D-02", "M3D-03",
                     "M3D-04", "M3D-05", "M3D-06", "M3D-07"):
            self.assertIn(part, self.source)


if __name__ == "__main__":
    unittest.main()
