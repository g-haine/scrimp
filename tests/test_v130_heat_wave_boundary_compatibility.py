# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Source-only native candidate contracts (not a native DAE PASS)."""

import ast
from pathlib import Path
import unittest

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
PROBE=ROOT/"sandbox/probe_v130_heat_wave_boundary_compatibility.py"
RUNNER=ROOT/"scripts/run-v130-heat-wave-boundary-compatibility.sh"


def actual_boundary_partition():
    source=PROBE.read_text(encoding="utf-8")
    node=next(n for n in ast.parse(source).body
              if isinstance(n,ast.FunctionDef) and
              n.name=="wave_boundary_partition")
    env={"np":np}
    exec(compile(ast.Module(body=[node],type_ignores=[]),
                 str(PROBE),"exec"),env)
    return env["wave_boundary_partition"]


class CompatibleBoundarySourceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.partition=staticmethod(actual_boundary_partition())

    def test_getfem_filtered_dof_order_and_boundary_selection(self):
        class FakeFem:
            def basic_dof_on_region(self,region):
                return {2:[0,1,3,5,6],20:[3,6]}[region]
        b,f,global_dofs,original_boundary=self.partition(
            FakeFem(),field_size=5
        )
        np.testing.assert_array_equal(b,[2,4])
        np.testing.assert_array_equal(f,[0,1,3])
        np.testing.assert_array_equal(global_dofs,[0,1,3,5,6])
        np.testing.assert_array_equal(original_boundary,[3,6])

    def test_invalid_boundary_mapping_rejected(self):
        class FakeFem:
            def basic_dof_on_region(self,region):
                return {2:[0,1,3],20:[1,4]}[region]
        with self.assertRaisesRegex(RuntimeError,"outside wave region"):
            self.partition(FakeFem(),field_size=3)

    def test_real_mass_optimality_equation_and_original_fixed_DAE(self):
        source=PROBE.read_text(encoding="utf-8")
        self.assertIn('model.add_linear_term(mim,"corrected_wave*Test_corrected_wave",2)',source)
        self.assertIn("Mff=M.createSubMatrix(subset,subset)",source)
        self.assertIn("corrected[free]+=free_delta",source)
        self.assertIn("stationary>1e-8",source)
        self.assertIn("original(system)",source)
        self.assertIn("except InconsistentDAEInitialData",source)
        self.assertIn("raise DAEProbeComplete",source)
        self.assertNotIn('method="interpolation"',source)

    def test_original_published_example_is_not_rewritten(self):
        source=(ROOT/"examples/heat_wave.py").read_text(encoding="utf-8")
        self.assertIn(
         'hw.set_initial_value("p", "5.*np.exp(-25*((x-0.6)*(x-0.6)+y*y))")',
         source)
        self.assertIn('hw.set_control("Boundary", "0.")',source)
        self.assertIn('init_step_ts_type="algebraic"',source)

    def test_runner_is_original_image_mpi1_and_no_full_TS(self):
        sh=RUNNER.read_text(encoding="utf-8")
        for literal in ("--candidate \"$candidate\"",
                        "constrained_l2","SCRIMP_GETFEM_MATRIX_MODE=distributed",
                        "SCRIMP_GETFEM_BRIDGE_DISTRIBUTED=coo",
                        "mpiexec --oversubscribe --bind-to none -n 1"):
            self.assertIn(literal,sh)
        self.assertIn("full_CN_TS_executed",PROBE.read_text(encoding="utf-8"))


if __name__=="__main__":
    unittest.main()
