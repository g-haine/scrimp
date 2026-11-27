# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Source and negative controls: original Heat-Wave full 15s MPI1/2.

Source tests do not qualify real native GetFEM/PETSc performance.
They guard immutable original physics, exact accepted FE history,
independent long-run MPI comparison and fail-closed acceptance.
"""
import ast
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np


ROOT=Path(__file__).resolve().parents[1]
PROBE=ROOT/"sandbox/probe_v130_heat_wave_mpi_extended.py"
CHECK=ROOT/"scripts/verify-v130-heat-wave-full-mpi2.py"
RUN=ROOT/"scripts/run-v130-heat-wave-full-mpi2-physics.sh"
ORIGINAL_CHECK=ROOT/"scripts/verify-v130-heat-wave-mpi-long-physics.py"


class HeatWaveFullMPI2V130Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.native_probe=PROBE.read_text(encoding="utf-8")
        cls.verify_source=CHECK.read_text(encoding="utf-8")
        cls.runner=RUN.read_text(encoding="utf-8")
        for path in (PROBE,CHECK,ORIGINAL_CHECK):
            ast.parse(path.read_text(encoding="utf-8"),filename=str(path))
        spec=importlib.util.spec_from_file_location(
            "v130_hw_full_mpi2_independent_checker",CHECK
        )
        cls.verifier=importlib.util.module_from_spec(spec)
        spec.loader.exec_module(cls.verifier)

    def test_native_original_getfem_model_unchanged(self):
        for needle in (
            "from examples.heat_wave import heat_wave_eq",
            "physics_only=True",
            'initial_wave_boundary="clamp"',
            "initial_wave_profile=args.profile",
            "collect_heat_wave_dual_facet_physics(",
            "np.savez_compressed(filename,**bundle)",
            'if args.store_full_states:',
            'bundle["z"]=states',
            '"all_MPI_ranks_agree_exactly"',
            "system.last_ts_steps!=args.steps",
            "missing_accepted_grid_indices",
        ):
            with self.subTest(needle=needle):
                self.assertIn(needle,self.native_probe)
        self.assertNotIn("DPHS.solve =",self.native_probe)
        self.assertNotIn("DPHS.init_step =",self.native_probe)

    def test_runner_is_true_full_horizon_MPI1_MPI2_accepted_state_gate(self):
        for needle in (
            'run_case mpi2_long 2',
            'run_case mpi1_long 1',
            '--steps 3000',
            '--dt 0.005',
            '--profile gaussian --store-full-states',
            'SCRIMP_DOCKER_NO_SUDO:-0',
            'scrimp-mpi:v1.3.0-rc3',
            'org.opencontainers.image.revision',
            'scripts/check-v130-published-immutable.py',
            'scripts/verify-v130-heat-wave-full-mpi2-normwise.py',
            '--out "$out/summary.json"',
            'all 3x3001x7689 original FE coefficients',
            'with zipfile.ZipFile(',
            'z.testzip()',
        ):
            with self.subTest(needle=needle):
                self.assertIn(needle,self.runner)
        self.assertNotIn("export_to_pv",self.runner)
        self.assertNotIn("plot_Hamiltonian",self.runner)
        self.assertNotIn("examples/published/:rw",self.runner)

    def test_independent_verifier_uses_all_states_not_only_endpoints(self):
        self.assertEqual(
            self.verifier.SHORT,("mpi1_long","mpi2_long")
        )
        self.assertEqual(self.verifier.LONG,"mpi1_long")
        for needle in (
            '"z"',
            "a[\"z\"].shape!=(steps+1,7689)",
            "compare(reference,a)",
            "np.array_equal(a[k],b[k])",
            "MPI2_rank0_rank1_bitwise_identical",
            "MPI1_vs_MPI2_all_3001_FE_states_and_physical_traces",
            "complete_15s_physical_power_balance_qualified",
            "MPI2_full_horizon_qualified",
            "missing_accepted_grid_indices",
            "3000",
            "7689",
            "0.7836421211876778",
            "0.008268389052222169",
        ):
            with self.subTest(needle=needle):
                self.assertIn(needle,self.verify_source)
        self.assertNotIn(
            "PASS_HEAT_WAVE_MPI2_SHORT_AND_MPI1_LONG_CN_PHYSICAL",
            self.verify_source,
        )

    def test_native_comparator_rejects_wrong_FE_state_and_shape(self):
        a={"z":np.zeros((5,4),dtype=float),
           "H":np.linspace(1.,.9,5)}
        b={k:v.copy() for k,v in a.items()}
        good=self.verifier.compare(a,b)
        self.assertEqual(good["z"],0.)
        b["z"][2,1]=.01
        with self.assertRaises(RuntimeError):
            self.verifier.compare(a,b)
        b["z"]=np.zeros((4,4))
        with self.assertRaises(RuntimeError):
            self.verifier.compare(a,b)

    def test_full_host_verifier_fails_closed_on_missing_native_run(self):
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)
            (path/"manifest.json").write_text(json.dumps({
                "source":"synthetic_negative_control_not_native",
                "original_model_unchanged":True,
                "full_FE_states_saved_per_rank":True,
                "profile":"gaussian","dt":.005,"steps":3000,
            }))
            result=self.verifier.audit(path)
            self.assertTrue(result["status"].startswith("FAIL_"))
            self.assertFalse(result["MPI2_full_horizon_qualified"])
            self.assertTrue(result["errors"])


if __name__=="__main__":
    unittest.main()
