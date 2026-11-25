# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Source-only and hostile-data checks for the original inviscid SWE pilot.

These tests are NOT a substitute for native MPI GetFEM/PETSc execution.
A genuine scientific PASS requires the four original solver outputs.
"""
import ast
import importlib.util
from pathlib import Path
import tempfile
import unittest

import numpy as np


ROOT=Path(__file__).resolve().parents[1]
MODEL=ROOT/"examples/shallow_water.py"
PROBE=ROOT/"sandbox/probe_v130_issue28_shallow_short_physical.py"
CHECK=ROOT/"scripts/verify-v130-shallow-water-extended.py"
RUN=ROOT/"scripts/run-v130-shallow-water-extended-physical.sh"


class OriginalShallowExtendedProtocol(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model=MODEL.read_text(encoding="utf-8")
        cls.probe=PROBE.read_text(encoding="utf-8")
        cls.runner=RUN.read_text(encoding="utf-8")
        cls.checker=CHECK.read_text(encoding="utf-8")
        spec=importlib.util.spec_from_file_location("sw_extended_verify",CHECK)
        cls.v=importlib.util.module_from_spec(spec)
        spec.loader.exec_module(cls.v)

    def test_real_original_model_and_physical_boundary_controls(self):
        for marker in (
            'S.Brick("M_p", "h * p . Test_p"',
            'explicit=True, position="effort"',
            'swe.set_initial_value("h", "3. - (np.sign(x-0.5)+1)/3.")',
            'swe.set_initial_value("p", f"[ 0., 0.]")',
            'init_step_ts_type="pseudo"',
            'ts_type="bdf"',
            'ts_bdf_order=2',
        ):
            self.assertIn(marker,self.model)
        for i in range(4):
            self.assertIn(f'swe.set_control("Boundary control {i}", "0.")',
                          self.model)

    def test_native_probe_preserves_original_FEM_and_extends_evidence_only(self):
        for marker in (
            "def build_original_without_published_long_run():",
            'patch.object(S.DPHS, "solve", skip_solve)',
            "module.shallow_water_eq()",
            'init_step_ts_type="nonlinear_affine"',
            'ts_type="bdf", ts_bdf_order=2',
            'ts_adapt_type="none"',
            'model.solve()',
            'parser.add_argument(',
            '"--capture-full-fe"',
            'np.asarray(model.solution["z"], dtype=np.float64)',
            'f"shallow_fe_rank{mpi.rank}.npz"',
            "np.savez_compressed(",
            '"BDF2_exact_power_conservation_claimed": False',
        ):
            self.assertIn(marker,self.probe)
        self.assertNotIn("patch.object(S.DPHS, \"init_step\"",self.probe)
        self.assertNotIn("SCRIMP_NONLINEAR_MASS_JACOBIAN=legacy",self.probe)
        self.assertNotIn('init_step_ts_type="pseudo",\n        ts_type="bdf"',self.probe)

    def test_bounded_native_four_runs_and_immutable_published_guard(self):
        for marker in (
            'for sample in "20:.0001" "40:.00005"',
            "for ranks in 1 2",
            '--final-time 0.002 --capture-full-fe',
            "--volume \"$root/scrimp:/workspace/scrimp/scrimp:ro\"",
            'scrimp-mpi:v1.3.0-rc3',
            '1d6d0c90115213298047b28237a7f38087563a35',
            'scripts/check-v130-published-immutable.py',
            'with zipfile.ZipFile(archive,"x"',
            '"exact_BDF2_energy_conservation_claimed":false',
        ):
            self.assertIn(marker,self.runner)
        self.assertNotIn("examples/published/",self.runner)
        self.assertIn('"status":"FAIL_OR_UNQUALIFIED_SWE_EXTENDED_BDF2"',
                      self.checker)
        for path in (MODEL,PROBE,CHECK):
            ast.parse(path.read_text(encoding="utf-8"),filename=str(path))

    def test_full_coefficients_reject_numerical_MPI_corruption(self):
        original={
            "t":np.array([0.,1.e-4]),
            "z":np.zeros((2,257)),
            "H":np.array([1.,.99]),
            "mass":np.ones(2),
            "boundary_power":np.zeros(2),
        }
        other={k:v.copy() for k,v in original.items()}
        checks=self.v.compare_states(original,other,tag="negative-fixture",
                                     atol=1e-10,rtol=1e-10)
        self.assertEqual(checks["z"],0.)
        other["z"][1,250]=.1
        with self.assertRaisesRegex(RuntimeError,"MPI parity FAIL"):
            self.v.compare_states(original,other,tag="negative-fixture",
                                  atol=1e-10,rtol=1e-10)

    def test_no_native_artifacts_is_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            result=self.v.audit(Path(temporary))
            self.assertTrue(result["status"].startswith("FAIL_"))
            self.assertTrue(result["errors"])
            self.assertFalse(result["BDF2_exact_power_conservation_claimed"])
            self.assertFalse(result["rotating_viscous_SWE_qualified"])

    def test_comparison_reports_observed_not_presumed_order(self):
        self.assertIn('"observed_temporal_order_qualified":False',self.checker)
        self.assertIn('"energy_defect_decreases":bool(abs(fine)<abs(coarse))',
                      self.checker)
        self.assertIn("H-H[0]-work",self.checker)
        self.assertIn("and relative_H<=5e-3",self.checker)


if __name__=="__main__":
    unittest.main()
