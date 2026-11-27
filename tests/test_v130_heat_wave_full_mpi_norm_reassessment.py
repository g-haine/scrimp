# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Negative controls for read-only full 15-second Heat-Wave MPI parity.

These SMALL SYNTHETIC arrays validate the chosen coefficient-vector
norm comparison and original scalar physical checks. They are not
native GetFEM evidence and must never substitute for original 15s data.
"""
import ast
import importlib.util
from pathlib import Path
import unittest
from unittest.mock import patch

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
REASSESS=ROOT/"scripts/reassess-v130-heat-wave-full-mpi2-norm-parity.py"
OLD=ROOT/"scripts/verify-v130-heat-wave-full-mpi2.py"
PROBE=ROOT/"sandbox/probe_v130_heat_wave_mpi_extended.py"


class NormwiseHeatWaveMPI2SourceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.reassessment_source=REASSESS.read_text(encoding="utf-8")
        cls.old_source=OLD.read_text(encoding="utf-8")
        cls.probe_source=PROBE.read_text(encoding="utf-8")
        ast.parse(cls.reassessment_source,filename=str(REASSESS))
        module=importlib.util.spec_from_file_location(
            "scrimp_HW_15s_normwise_reassessment",REASSESS
        )
        cls.reassess=importlib.util.module_from_spec(module)
        module.loader.exec_module(cls.reassess)

    def test_same_absolute_relative_tolerances_without_native_model_change(self):
        self.assertEqual(self.reassess.ATOL,5e-10)
        self.assertEqual(self.reassess.RTOL,2e-8)
        for required in (
            'ORIGINAL=Path(__file__).with_name("verify-v130-heat-wave-full-mpi2.py")',
            'original.load(',
            'compare_physics_unchanged(',
            'np.allclose(a,b,atol=ATOL,rtol=RTOL)',
            '"previous_original_strict_status_preserved":True',
            'prior.get("status")!="FAIL_OR_UNQUALIFIED_FULL_HEAT_WAVE_MPI2"',
            '"full_15s_normwise_MPI1_MPI2_parity_qualified":False',
            'original_MPI2_rank0_rank1_bitwise_identical',
            'np.array_equal(mpi20[k],mpi21[k])',
            'np.array_equal(a["z"][0],a["z0"])',
            'np.array_equal(a["z"][-1],a["zf"])',
            'normwise_all_saved_states(mpi1["z"],mpi20["z"])',
        ):
            with self.subTest(required=required):
                self.assertIn(required,self.reassessment_source)
        self.assertIn("atol=5e-10,rtol=2e-8",self.old_source)
        self.assertNotIn("import examples.heat_wave",self.reassessment_source)
        self.assertNotIn("heat_wave_eq(",self.reassessment_source)
        self.assertNotIn("set_time_scheme(",self.reassessment_source)
        self.assertNotIn("SCRIMP_GETFEM_MATRIX_MODE",self.reassessment_source)

    def test_small_local_FE_difference_passes_vector_norm_but_not_old_pointwise(self):
        a=np.tile([1.,0.,0.],(4,1))
        b=a.copy()
        b[-1,1]=2.758281700265819e-9
        with patch.object(self.reassess,"STEPS",3),patch.object(
            self.reassess,"DOFS",3
        ):
            result=self.reassess.normwise_all_saved_states(a,b)
        self.assertTrue(result["qualified"])
        self.assertEqual(result["failed_infinity_norm_steps"],0)
        self.assertEqual(result["failed_euclidean_norm_steps"],0)
        self.assertEqual(
            result["original_pointwise_failures_preserved_as_diagnostic"],1
        )
        self.assertLess(result["max_relative_coefficient_l2_per_time"],2e-8)
        self.assertTrue(result["original_initial_state_passes_old_pointwise"])

    def test_large_local_discrepancy_fails_both_real_vector_norms(self):
        a=np.tile([1.,0.,0.],(4,1))
        b=a.copy()
        b[-1,1]=1e-3
        with patch.object(self.reassess,"STEPS",3),patch.object(
            self.reassess,"DOFS",3
        ):
            result=self.reassess.normwise_all_saved_states(a,b)
        self.assertFalse(result["qualified"])
        self.assertEqual(result["first_failed_infinity_step"],3)
        self.assertEqual(result["first_failed_euclidean_step"],3)
        self.assertGreater(result["maximum_normwise_infinity_tolerance_ratio"],1)

    def test_same_original_strict_physics_threshold_rejects_corruption(self):
        baseline={"H":np.array([1.,.8,.5]),
                  "phys_CN_D_heat_flux":np.array([1.,.5,.3])}
        same={key:value.copy() for key,value in baseline.items()}
        diffs=self.reassess.compare_physics_unchanged(baseline,same)
        self.assertEqual(diffs["H"],0.)
        same["phys_CN_D_heat_flux"][1]+=.001
        with self.assertRaises(RuntimeError):
            self.reassess.compare_physics_unchanged(baseline,same)

    def test_absent_original_native_directory_fails_closed(self):
        from tempfile import TemporaryDirectory
        with TemporaryDirectory() as t:
            result=self.reassess.audit(Path(t))
        self.assertTrue(result["status"].startswith("FAIL_"))
        self.assertFalse(result["full_15s_normwise_MPI1_MPI2_parity_qualified"])
        self.assertTrue(result["errors"])


if __name__=="__main__":
    unittest.main()
