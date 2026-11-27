"""Source and negative-control tests for genuine Heat-Wave CN temporal order.

The regression checks only source contracts and miniature synthetic
Richardson algebra; it does not replace original native GetFEM states.
"""
import ast
import importlib.util
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
COMPARE=ROOT/"sandbox/compare_v130_heat_wave_CN_dt_refinement.py"
RUNNER=ROOT/"scripts/run-v130-heat-wave-CN-dt-refinement.sh"
PROBE=ROOT/"sandbox/probe_v130_heat_wave_mpi_extended.py"
EXAMPLE=ROOT/"examples/heat_wave.py"


class HWCNdtRefinementV130(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.s=COMPARE.read_text(encoding="utf-8")
        cls.run=RUNNER.read_text(encoding="utf-8")
        cls.probe=PROBE.read_text(encoding="utf-8")
        cls.example=EXAMPLE.read_text(encoding="utf-8")
        ast.parse(cls.s)
        ast.parse(cls.probe)
        ast.parse(cls.example)
        spec=importlib.util.spec_from_file_location(
            "v130_original_hw_CN_temporal_checker",COMPARE
        )
        cls.verify=importlib.util.module_from_spec(spec)
        spec.loader.exec_module(cls.verify)

    def test_unchanged_canonical_mesh_and_original_CrankNicolson(self):
        self.assertIn('S.Domain("Concentric", {"R": 1.0, "r": 0.6, "h": 0.05})',
                      self.example)
        self.assertIn('ts_type="cn"',self.example)
        self.assertIn('init_step_ts_type="algebraic"',self.example)
        self.assertIn("physics_only=True",self.probe)
        self.assertIn("from examples.heat_wave import heat_wave_eq",self.probe)
        self.assertIn('initial_wave_boundary="clamp"',self.probe)
        self.assertIn('if args.store_full_states:',self.probe)
        self.assertNotIn("heat_wave_eq(",self.s)
        self.assertNotIn("S.DPHS(",self.s)

    def test_native_four_nested_real_time_grids_and_complete_states(self):
        for needle in (
            '"dt_0p01:0.01:20"',
            '"dt_0p005:0.005:40"',
            '"dt_0p0025:0.0025:80"',
            '"dt_0p00125:0.00125:160"',
            '--profile gaussian --store-full-states',
            '--workdir /scrimp-work',
            'sandbox/probe_v130_heat_wave_mpi_extended.py',
            'sandbox/compare_v130_heat_wave_CN_dt_refinement.py',
            'check-v130-published-immutable.py',
            'SCRIMP_DOCKER_NO_SUDO:-0',
            'v130-heat-wave-CN-dt-refinement',
        ):
            with self.subTest(needle=needle):
                self.assertIn(needle,self.run)
        self.assertNotIn('export_to_pv(',self.run)
        self.assertNotIn('set_time_scheme(',self.run)
        self.assertNotIn('examples/published:rw',self.run)

    def test_scientific_verifier_rejects_shortcut_data(self):
        for needle in (
            'original_model_unchanged',
            '"same_original_mesh_h":.05',
            'a.get("z",np.array([])).shape!=(steps+1,FE_DOF)',
            'np.array_equal(a["z"][0],a["z0"])',
            'np.array_equal(a["z"][-1],a["zf"])',
            'meta.get("status")=="PASS_NATIVE_STRICT_DAE_CN_ORIENTED_POWER"',
            'all(meta.get("physical_gates",{}).values())',
            '"phys_CN_per_step_heat_balance_error"',
            '"phys_CN_per_step_wave_balance_error"',
            'np.linalg.norm(z0-z1)',
            'np.linalg.norm(z1-z2)',
            '3.0<=ratio<=5.0',
            'no_spatial_convergence_claim',
        ):
            with self.subTest(needle=needle):
                self.assertIn(needle,self.s)

    def tiny(self,dt):
        t=np.arange(round(.2/dt)+1,dtype=float)*dt
        base=np.column_stack((
            np.sin(3*t),np.cos(2*t),np.sin(7*t),np.cos(5*t)
        ))
        # Manufactured inputs ONLY to exercise the diagnostic.
        # Physical acceptance requires archived native GetFEM.
        error=dt**2*np.sin(t)[:,None]*np.array([1.,2.,-1.,.5])
        return {"t":t,"z":base+error}

    def test_true_second_order_algebra_ratio_four(self):
        a,b,c=(self.tiny(dt) for dt in (.01,.005,.0025))
        with patch.object(self.verify,"FE_DOF",4):
            report=self.verify.nested_triple(a,b,c,"synthetic_negative_control")
        self.assertAlmostEqual(report["Richardson_ratio"],4.,places=6)
        self.assertAlmostEqual(report["observed_order_log2_ratio"],2.,places=6)

    def test_bad_temporal_order_or_changed_initial_data_fails(self):
        a,b,c=(self.tiny(dt) for dt in (.01,.005,.0025))
        c["z"][0,0]=.001
        with patch.object(self.verify,"FE_DOF",4):
            with self.assertRaisesRegex(RuntimeError,"initial"):
                self.verify.nested_triple(a,b,c,"synthetic_negative_control")
        c=self.tiny(.0025)
        # Malformed deliberately non-convergent medium trajectory.
        b["z"][1:,0]+=.03
        with patch.object(self.verify,"FE_DOF",4):
            with self.assertRaisesRegex(RuntimeError,"second-order"):
                self.verify.nested_triple(a,b,c,"synthetic_negative_control")

    def test_absent_native_case_fails_closed(self):
        with tempfile.TemporaryDirectory() as work:
            result=self.verify.audit(Path(work))
        self.assertTrue(result["status"].startswith("FAIL_"))
        self.assertTrue(result["errors"])


if __name__=="__main__":
    unittest.main()
