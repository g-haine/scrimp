"""Source and negative-control regression for future Heat-Wave MPI full runs.

The historical per-DOF strict verifier is deliberately unchanged.
The accepted native gate uses full-vector norms on ALL original FE
coefficients while retaining unchanged strict scalar power checks.
No synthetic array here is claimed to be physical GetFEM evidence.
"""
import ast
import importlib.util
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
SCRIPT=ROOT/"scripts/verify-v130-heat-wave-full-mpi2-normwise.py"
RUNNER=ROOT/"scripts/run-v130-heat-wave-full-mpi2-physics.sh"
OLD=ROOT/"scripts/verify-v130-heat-wave-full-mpi2.py"


class HeatWaveNormwiseReusableSourceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.source=SCRIPT.read_text(encoding="utf-8")
        ast.parse(cls.source,filename=str(SCRIPT))
        cls.runner=RUNNER.read_text(encoding="utf-8")
        spec=importlib.util.spec_from_file_location(
            "scrimp_hw_native_reusable_norm_test", SCRIPT
        )
        cls.module=importlib.util.module_from_spec(spec)
        spec.loader.exec_module(cls.module)

    def test_uses_native_independent_original_physics_not_substitute(self):
        for key in (
            '"verify-v130-heat-wave-full-mpi2.py"',
            '"reassess-v130-heat-wave-full-mpi2-norm-parity.py"',
            'previous.load(',
            'norm.compare_physics_unchanged(x,y)',
            'norm.normwise_all_saved_states(x["z"],y["z"])',
            'np.array_equal(y[field],z[field])',
            'np.array_equal(data["z0"],data["z"][0])',
            'np.array_equal(data["zf"],data["z"][-1])',
            'pre["accepted_PETSc_TS_steps"]!=3000',
            'pre["TS_step_rejections"]!=0',
            'pre["TS_SNES_failures"]!=0',
            'original_scalar_and_power_pointwise_parity_qualified',
            'per_DOF_pointwise_parity_claimed',
            'mass_weighted_FE_L2_convergence_claimed',
            '"source"',
            '"image_revision"',
        ):
            with self.subTest(key=key):
                self.assertIn(key,self.source)
        self.assertNotIn("heat_wave_eq(",self.source)
        self.assertNotIn("S.DPHS(",self.source)
        self.assertNotIn("set_time_scheme(",self.source)

    def test_small_absolute_per_DOF_failure_is_distinct_from_vector_norm(self):
        a=np.tile([1.,0.,0.],(4,1))
        b=a.copy()
        b[-1,1]=2.758281700265819e-9
        with patch.object(self.module.norm,"STEPS",3),patch.object(
            self.module.norm,"DOFS",3
        ):
            result=self.module.norm.normwise_all_saved_states(a,b)
        self.assertTrue(result["qualified"])
        self.assertEqual(result["failed_infinity_norm_steps"],0)
        self.assertEqual(result["failed_euclidean_norm_steps"],0)
        # The original historical per-DOF metric really rejects this.
        with self.assertRaisesRegex(RuntimeError,"physical FE trace differs"):
            self.module.previous.compare(
                {"z":a},{"z":b}
            )

    def test_missing_real_native_provenance_stays_fail_closed(self):
        with tempfile.TemporaryDirectory() as folder:
            outcome=self.module.audit(Path(folder))
        self.assertTrue(outcome["status"].startswith("FAIL_"))
        self.assertFalse(outcome["full_15s_normwise_MPI1_MPI2_parity_qualified"])
        self.assertTrue(outcome["errors"])

    def test_new_run_uses_normwise_checker_not_historical_pointwise(self):
        self.assertIn(
            'scripts/verify-v130-heat-wave-full-mpi2-normwise.py "$out"',
            self.runner
        )
        self.assertIn('--out "$out/summary.json"',self.runner)
        self.assertIn(
            'PASS_HEAT_WAVE_FULL_15S_MPI1_MPI2_CN_NORM_PARITY',
            self.runner
        )
        self.assertIn('check-v130-published-immutable.py',self.runner)
        self.assertIn('--store-full-states',self.runner)
        self.assertNotIn(
            'python3 scripts/verify-v130-heat-wave-full-mpi2.py "$out"',
            self.runner
        )
        self.assertIn('np.allclose(x,y,atol=5e-10,rtol=2e-8)',
                      OLD.read_text(encoding="utf-8"))


if __name__=="__main__":
    unittest.main()
