# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Source/synthetic negative controls for actual nonlinear SWE TS gate.

These tests cannot replace the native 1/2-rank GetFEM/PETSc solve.
"""

import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT=Path(__file__).resolve().parents[1]
PROBE=ROOT/"sandbox/probe_v130_shallow_explicit_ts.py"
COMPARE=ROOT/"sandbox/compare_v130_shallow_explicit_ts.py"
RUNNER=ROOT/"scripts/run-v130-shallow-explicit-physical-ts.sh"


class TestRealShallowExplicitTSGate(unittest.TestCase):
    def setUp(self):
        tmp=tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.out=Path(tmp.name)
        self._write_valid()

    def _save(self,case,ranks,doc):
        (self.out/f"shallow_explicit_TS_{case}_{ranks}r.json").write_text(
            json.dumps(doc),encoding="utf-8"
        )

    def _write_valid(self):
        for case in ("published","manufactured"):
            for ranks in (1,2):
                ts=[0.,1e-4,2e-4,3e-4,4e-4]
                events=[
                    {
                        "time":t,"frozen_rhs_time":t,
                        "accepted_step_count":i,
                        "frozen_gyroscopic_rhs_norm":(
                            1e-3 if case=="manufactured" else 0.),
                    } for i,t in enumerate(ts[:-1])
                ]
                d={
                    "manufactured_nonzero_G_initial_data":
                        case=="manufactured",
                    "time_scheme":"bdf",
                    "time_steps":4,"dt":1e-4,"time_end":4e-4,
                    "TS_converged_reason":1,
                    "TS_accepted_steps":4,
                    "TS_rejected_steps":0,
                    "TS_SNES_failures":0,
                    "main_explicit_prestep":events,
                    "main_IFunction_checks":5,
                    "max_frozen_explicit_RHS_drift":0.,
                    "max_nonzero_explicit_G_norm":(
                        1e-3 if case=="manufactured" else 0.),
                    "times":ts,
                    "full_state_trajectory":[
                        [1.+i*1e-6,0.3] for i in range(len(ts))
                    ],
                    "hamiltonian":[5.+i*1e-5 for i in range(len(ts))],
                    "signed_H_final_minus_H_initial":4e-5,
                    "relative_abs_H_change":8e-6,
                    "exact_H_balance_claimed":False,
                    "nonlinear_DAE_initialization_qualified":False,
                }
                self._save(case,ranks,d)

    def _run(self):
        return subprocess.run(
            [sys.executable,str(COMPARE),str(self.out)],
            text=True,capture_output=True
        )

    def test_valid_short_ts_results_pass(self):
        out=self._run()
        self.assertEqual(out.returncode,0,out.stderr)
        summary=json.loads((
            self.out/"shallow_explicit_TS_comparison.json"
        ).read_text(encoding="utf-8"))
        self.assertEqual(summary["status"],"PASS")
        self.assertFalse(summary["exact_energy_preservation_claimed"])

    def test_nonzero_manufactured_gyro_required(self):
        p=self.out/"shallow_explicit_TS_manufactured_2r.json"
        d=json.loads(p.read_text(encoding="utf-8"))
        d["max_nonzero_explicit_G_norm"]=0.
        self._save("manufactured",2,d)
        out=self._run()
        self.assertNotEqual(out.returncode,0)
        self.assertIn("vacuous",out.stderr.lower())

    def test_newton_cache_must_remain_exact(self):
        p=self.out/"shallow_explicit_TS_published_1r.json"
        d=json.loads(p.read_text(encoding="utf-8"))
        d["max_frozen_explicit_RHS_drift"]=1e-9
        self._save("published",1,d)
        out=self._run()
        self.assertNotEqual(out.returncode,0)
        self.assertIn("Invalid physical TS record",out.stderr)

    def test_mpi_final_state_discrepancy_rejected(self):
        p=self.out/"shallow_explicit_TS_manufactured_2r.json"
        d=json.loads(p.read_text(encoding="utf-8"))
        d["full_state_trajectory"][-1][0]+=1e-3
        self._save("manufactured",2,d)
        out=self._run()
        self.assertNotEqual(out.returncode,0)
        self.assertIn("MPI disparity",out.stderr)

    def test_source_and_runner_keep_g_and_distinct_initial_data(self):
        probe=PROBE.read_text(encoding="utf-8")
        runner=RUNNER.read_text(encoding="utf-8")
        self.assertIn('examples/shallow_water.py',probe)
        self.assertIn('value="[y,-x]"',probe)
        self.assertIn('if manufactured and name=="p":',probe)
        self.assertIn(
            'original_prestep(self,ts)',probe
        )
        self.assertIn(
            'state["max_main_cache_drift"]',probe
        )
        self.assertIn(
            'self._explicit_rhs.copy()',probe
        )
        self.assertIn(
            'system.get_Hamiltonian()',probe
        )
        self.assertIn('for ranks in 1 2;',runner)
        self.assertIn('for case in published manufactured;',runner)
        self.assertIn('scrimp/utils/dae_initialization.py',runner)
        self.assertIn('scrimp-mpi:v1.3.0-rc3',runner)
        self.assertNotIn('ts_type="cn"',probe)


if __name__=="__main__":
    unittest.main()
