# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Host-only regression for #30 CN stage powers vs wrong trapezoid.

The standalone helper is tested without importing GetFEM, PETSc or SCRIMP.
Only original native linear wave results establish real physical fidelity.
"""

import importlib.util
from pathlib import Path
import unittest

import numpy as np


ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/"scrimp"/"utils"/"power_integration.py"
spec=importlib.util.spec_from_file_location(
    "_v130_power_integration", SOURCE
)
module=importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
integrate_port_power=module.integrate_port_power


class PowerIntegrationTests(unittest.TestCase):
    def test_cn_integrates_step_powers_not_trapezoidal_nodes(self):
        times=np.array([0.,.002,.004,.006])
        # P0 can have a completely different instantaneous value.
        raw=np.array([123.456,-.021,-.038,-.007])
        expected=np.r_[0.,np.cumsum(np.diff(times)*raw[1:])]
        got=integrate_port_power(times,raw,convention="cn_step")
        np.testing.assert_array_equal(got,expected)
        wrong=integrate_port_power(
            times,raw,convention="legacy_trapezoid"
        )
        self.assertGreater(abs(got[1]-wrong[1]),.01)

    def test_cn_preserves_native_linear_wave_signed_H_plus_port_integral(self):
        t=np.array([0.,.001,.002,.003])
        net=np.array([17.,-.1,-.2,.05])
        signed=integrate_port_power(t,net,convention="cn_step")
        # Native linear wave's original signed port orientation:
        # H[k]-H0 = -dt[k]*P[k] (measured to roundoff in #29).
        H=.5-signed
        np.testing.assert_allclose(H+signed,H[0],rtol=0,atol=1e-15)
        wrong=integrate_port_power(t,net,convention="legacy_trapezoid")
        self.assertGreater(np.max(abs(H+wrong-H[0])),.001)

    def test_cn_restart_discontinuity_in_initial_P0_does_not_change_work(self):
        t=np.arange(0.,.025,.002)
        p=-np.linspace(.01,.025,len(t))
        full=integrate_port_power(t,p,convention="cn_step")
        cut=10
        seed=integrate_port_power(
            t[:cut+1],p[:cut+1],convention="cn_step"
        )
        resumed=p[cut:].copy()
        resumed[0]=999.  # New initial instantaneous port state.
        restart=integrate_port_power(
            t[cut:],resumed,convention="cn_step"
        )
        np.testing.assert_allclose(
            np.r_[seed,seed[-1]+restart[1:]],full,
            rtol=0,atol=1e-15
        )

    def test_nonuniform_fixed_cn_step_clock(self):
        t=np.array([0.,.001,.003,.008])
        p=np.array([100.,1.,2.,3.])
        got=integrate_port_power(t,p,convention="cn_step")
        np.testing.assert_allclose(
            got,[0.,.001,.005,.020],atol=1e-15
        )

    def test_noncn_legacy_output_is_identical_to_old_trapezoid(self):
        t=np.array([0.,.25,.75,1.])
        p=np.array([10.,11.,12.,14.])
        expected=np.zeros(len(t))
        for k in range(1,len(t)):
            expected[k]=expected[k-1]+.5*(t[k]-t[k-1])*(p[k-1]+p[k])
        np.testing.assert_allclose(
            integrate_port_power(t,p,convention="legacy_trapezoid"),
            expected,rtol=0,atol=1e-15
        )

    def test_refuse_invalid_shapes_time_and_nonfinite_or_unknown_scheme(self):
        for t,p,method in [
            ([0.,1.],[1.],"cn_step"),
            ([0.,0.],[1.,2.],"cn_step"),
            ([1.,0.],[1.,2.],"cn_step"),
            ([0.,1.],[1.,float("nan")],"cn_step"),
            ([0.,float("inf")],[1.,2.],"cn_step"),
            ([0.,1.],[1.,2.],"unknown"),
            ([[0.,1.]],[1.,2.],"cn_step"),
        ]:
            with self.subTest(time=t,power=p,method=method):
                with self.assertRaises(ValueError):
                    integrate_port_power(t,p,convention=method)

    def test_public_plot_uses_cn_step_or_explicit_legacy(self):
        s=(ROOT/"scrimp"/"dphs.py").read_text(encoding="utf-8")
        self.assertIn("from scrimp.utils.power_integration import integrate_port_power",s)
        self.assertIn('scheme = self.time_scheme.getString("ts_type", "bdf").lower()',s)
        self.assertIn('"cn_step" if scheme == "cn" else "legacy_trapezoid"',s)
        self.assertIn("values = integrate_port_power(",s)
        # The old direct trapezoid must not remain in plot_powers().
        body=s[s.index("    def plot_powers("):s.index("    def spy_Dirac(")]
        self.assertNotIn("power[k - 1] + power[k]",body)
        # Original GetFEM stage semantics are preserved.
        port=(ROOT/"scrimp"/"port.py").read_text(encoding="utf-8")
        self.assertIn("0.5 * (solution[\"z\"][k] + solution[\"z\"][k - 1])",port)


if __name__=="__main__":
    unittest.main()
