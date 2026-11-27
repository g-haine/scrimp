"""Negative controls for archived original Heat-Wave MPI parity post-mortem.

These small arrays only test the host diagnostic; not GetFEM evidence.
"""
import ast
import importlib.util
from pathlib import Path
from unittest import TestCase, main
from unittest.mock import patch

import numpy as np


ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/"sandbox/diagnose_v130_heat_wave_full_mpi_parity.py"


class HeatWaveOriginalMPIDiagnosticTests(TestCase):
    @classmethod
    def setUpClass(cls):
        code=SOURCE.read_text(encoding="utf-8")
        ast.parse(code,filename=str(SOURCE))
        spec=importlib.util.spec_from_file_location("heatwave_mpi_archive_diagnostic",SOURCE)
        cls.analysis=importlib.util.module_from_spec(spec)
        spec.loader.exec_module(cls.analysis)
        cls.source=code

    def test_pure_post_mortem_without_numerical_shortcuts(self):
        for fragment in (
            "POST_MORTEM_ONLY_NOT_A_RELEASE_ACCEPTANCE",
            "MPI2_rank0_rank1_exact_identical",
            "np.array_equal(two[key],replica[key])",
            "pointwise_failed_comparisons",
            "near_zero_violations",
            "relative_l2_over_all_saved_FE_states",
            "maximum_relative_l2_per_saved_state",
            "t_final_largest_original_DOF_differences",
            "physical_and_scalar_traces",
            "for lo in range(0,len(time),64)",
        ):
            with self.subTest(fragment=fragment):
                self.assertIn(fragment,self.source)
        self.assertNotIn("from examples.heat_wave",self.source)
        self.assertNotIn("heat_wave_eq(",self.source)
        self.assertNotIn("DPHS(",self.source)
        self.assertNotIn("PASS_HEAT_WAVE_FULL_15S_MPI1_MPI2",self.source)

    def test_detects_near_zero_pointwise_failure_without_faking_a_pass(self):
        a=np.ones((4,3),dtype=float)
        a[:,1]=0.
        b=a.copy()
        t=np.arange(4)*.005
        with patch.object(self.analysis,"STATE_SHAPE",(4,3)):
            exact=self.analysis.state_mismatch_summary(a,b,t)
            self.assertEqual(exact["pointwise_failed_comparisons"],0)
            self.assertEqual(exact["maximum_absolute_state_difference"],0.)
            b[-1,1]=2.758281700265819e-9
            difference=self.analysis.state_mismatch_summary(a,b,t)
        self.assertEqual(difference["pointwise_failed_comparisons"],1)
        self.assertEqual(
            difference["failed_comparisons_with_both_values_below_1e-3"],1
        )
        self.assertEqual(difference["worst_step"],3)
        self.assertEqual(difference["worst_DOF"],1)
        self.assertGreater(difference["maximum_pointwise_tolerance_ratio"],5.)
        self.assertLess(difference["relative_l2_over_all_saved_FE_states"],1e-8)

    def test_both_large_local_errors_and_power_corruption_are_visible(self):
        a=np.ones((4,3),dtype=float)
        b=a.copy()
        b[2,2]+=0.1
        with patch.object(self.analysis,"STATE_SHAPE",(4,3)):
            report=self.analysis.state_mismatch_summary(a,b,np.arange(4)*.005)
        self.assertGreater(report["maximum_absolute_state_difference"],.09)
        self.assertGreater(report["maximum_relative_l2_per_saved_state"],.01)
        power=self.analysis.scalar_trace_comparison(
            {"H":np.array([1.,.9,.8,.7])},
            {"H":np.array([1.,.9,.8,.71])}
        )
        self.assertFalse(power["H"]["pointwise_old_allclose"])
        self.assertGreater(power["H"]["maximum_absolute_difference"],.009)


if __name__=="__main__":
    main()
