# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Host-only negative controls for original heat-wave physical diagnostics.

These are synthetic algebraic checks, NOT a native SCRIMP PDE qualification.
"""

import importlib.util
from pathlib import Path
import unittest

import numpy as np


ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/"sandbox"/"measure_v130_heat_wave_physics.py"
spec=importlib.util.spec_from_file_location("_heat_wave_physics",SOURCE)
physical=importlib.util.module_from_spec(spec)
spec.loader.exec_module(physical)


class SyntheticGetFEM:
    def __init__(self, t, *, mismatch=False, fail_optional=False):
        self.t=np.asarray(t,dtype=float)
        self.mismatch=mismatch
        self.fail_optional=fail_optional

    def get_quantity(self, expression, region=-1, CN=False):
        if self.fail_optional and expression.startswith("("):
            raise ValueError("Synthetic trace feature unavailable")
        t=self.t
        constant=lambda value: np.ones(t.shape)*value
        values={
            (".5*T*T",1):4.-2.*t,
            (".5*p*p",2):.3+.4*t,
            (".5*q.q",2):.7+.4*t,
            ("e_Q.e_Q",1):constant(1.2),
            ("U_T*Y_T",10):constant(.8),
            ("U_w*Y_w",10):constant(.9 if self.mismatch else .8),
            ("U_bnd*Y_bnd",20):constant(0.),
        }
        if (expression,region) in values:
            return values[expression,region]
        if expression.startswith("("):
            return constant(0.)
        raise AssertionError("Unexpected mock GetFEM expression "+expression)


class HeatWaveObservableTests(unittest.TestCase):
    def test_closed_coupling_exact_synthetic_identity(self):
        t=np.array([0.,.2,.8,1.])
        model=SyntheticGetFEM(t)
        H=5.-1.2*t
        v,report=physical.collect_heat_wave_physics(model,t,H)
        np.testing.assert_allclose(v["H_original"],v["H_components"])
        np.testing.assert_allclose(v["P_interface_sum"],0.,atol=1e-14)
        np.testing.assert_allclose(
            v["total_balance_nodal_trapz_defect"],0.,atol=1e-14
        )
        np.testing.assert_allclose(
            v["heat_balance_nodal_trapz_defect"],0.,atol=1e-14
        )
        np.testing.assert_allclose(
            v["wave_balance_nodal_trapz_defect"],0.,atol=1e-14
        )
        self.assertEqual(report["dissipation_min"],1.2)
        self.assertEqual(
            report["scientific_status"],
            "DIAGNOSTIC_ONLY_MODEL_NOT_YET_QUALIFIED"
        )
        self.assertLess(report["max_CN_step_H_balance_defect"],1e-13)
        self.assertLess(report["max_CN_interface_power_defect"],1e-13)

    def test_nonconserving_gyrator_is_not_silently_accepted(self):
        t=np.array([0.,.5,1.])
        model=SyntheticGetFEM(t,mismatch=True)
        _,report=physical.collect_heat_wave_physics(model,t,5.-1.2*t)
        self.assertGreater(report["max_interface_port_power_sum"],.09)
        self.assertGreater(report["max_CN_interface_power_defect"],.09)
        self.assertEqual(
            report["scientific_status"],
            "DIAGNOSTIC_ONLY_MODEL_NOT_YET_QUALIFIED"
        )

    def test_optional_trace_unsupported_is_reported_not_set_to_zero(self):
        t=np.array([0.,.5,1.])
        model=SyntheticGetFEM(t,fail_optional=True)
        traces,report=physical.collect_heat_wave_physics(model,t,5.-1.2*t)
        self.assertNotIn("interface_T_plus_p_L2_sq",traces)
        self.assertIn(
            "interface_T_plus_p_L2_sq_error",
            report["optional_surface_traces_unavailable"]
        )

    def test_reject_negative_dissipation_and_incorrect_hamiltonian(self):
        t=np.array([0.,.5,1.])
        with self.assertRaisesRegex(ValueError,"Hamiltonian"):
            physical.collect_heat_wave_physics(
                SyntheticGetFEM(t),t,5.-1.2*t+1.
            )
        class BadD(SyntheticGetFEM):
            def get_quantity(self, expression, region=-1):
                if expression=="e_Q.e_Q":
                    return -np.ones(len(self.t))
                return super().get_quantity(expression,region)
        with self.assertRaisesRegex(ValueError,"negative"):
            physical.collect_heat_wave_physics(
                BadD(t),t,5.-1.2*t
            )

    def test_no_fake_bdf_stage_power_claim(self):
        t=np.array([0.,.1,.4])
        p=np.array([1.,2.,3.])
        expected=np.array([0.,.15,.9])
        np.testing.assert_allclose(
            physical.trapezoid(t,p),expected,atol=1e-15
        )
        report=physical.collect_heat_wave_physics(
            SyntheticGetFEM(t),t,5.-1.2*t
        )[1]
        self.assertTrue(
            report["CN_step_powers_are_assembled_at_mean_FE_state"]
        )
        self.assertTrue(report["nodal_trapezoid_is_diagnostic_only"])

    def test_refuse_nonfinite_and_wrong_shape(self):
        t=np.array([0.,.5,1.])
        with self.assertRaises(ValueError):
            physical.trapezoid(t,[2.,3.])
        with self.assertRaises(ValueError):
            physical.trapezoid([0.,0.,1.],[1.,2.,3.])
        with self.assertRaises(ValueError):
            physical.as_trace([1.,np.nan,1.],t,"bad")

    def test_physics_implementation_must_not_modify_original_example(self):
        text=SOURCE.read_text(encoding="utf-8")
        self.assertIn("system.get_quantity(",text)
        self.assertIn('"U_T*Y_T"',text)
        self.assertIn('"U_w*Y_w"',text)
        self.assertIn('"e_Q.e_Q"',text)
        self.assertNotIn("hw.set_control(",text)
        self.assertNotIn("hw.add_brick(",text)


if __name__=="__main__":
    unittest.main()
