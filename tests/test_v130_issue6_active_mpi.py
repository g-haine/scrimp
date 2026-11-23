# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Host-only negative controls; NO native GetFEM/PETSc acceptance claimed."""

import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np


ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/"sandbox"/"compare_v130_issue6_active_mpi.py"
spec=importlib.util.spec_from_file_location("_v130_issue6_active_mpi",SOURCE)
# Import comparator from sandbox/ even when tests run from repository root.
import sys
sys.path.insert(0,str(ROOT/"sandbox"))
try:
    comparator=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(comparator)
finally:
    sys.path.pop(0)


class NativeIssue6ActivePowerTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.directory=Path(self.tmp.name)
        self.times=np.array([0.,.002,.004])
        self.fingerprint="genuine_original_model_fingerprint"
        self.powers={
            "Wall velocity":[0.,.01,.02],
            "Normal boundary velocity":[0.,-.001,-.002],
        }
        (self.directory/"mpi_physical_powers.png").write_bytes(b"P"*2048)
        for rank in range(4):
            payload={
                "schema":comparator.SCHEMA,
                "rank":rank,"world_size":4,
                "source_fingerprint":self.fingerprint,
                "physical_h_weighted_balance_is_separate":True,
                "normal_feedback_active":False,
                "sample_semantics":"P[0] instantaneous; P[k>=1] midpoint-stage",
                "times":self.times.tolist(),
                "powers":self.powers,
            }
            (self.directory/f"mpi_port_powers_rank{rank}.json").write_text(
                json.dumps(payload),encoding="utf-8"
            )

    def edit(self,rank,fn):
        path=self.directory/f"mpi_port_powers_rank{rank}.json"
        d=json.loads(path.read_text(encoding="utf-8"))
        fn(d)
        path.write_text(json.dumps(d),encoding="utf-8")

    def test_real_rank_local_1_2_4_port_contract(self):
        traces,report=comparator.rank_port_traces(
            self.directory,4,self.times,self.fingerprint
        )
        self.assertEqual(report["unique_rank_outputs"],4)
        self.assertEqual(report["max_intra_MPI_generic_port_difference"],0.)
        self.assertGreater(report["max_original_generic_port_activity"],.01)
        self.assertEqual(set(traces),set(self.powers))

    def test_missing_rank_fails_closed(self):
        (self.directory/"mpi_port_powers_rank3.json").unlink()
        with self.assertRaisesRegex(ValueError,"Missing or unexpected"):
            comparator.rank_port_traces(
                self.directory,4,self.times,self.fingerprint
            )

    def test_rank_specific_port_mismatch_fails(self):
        self.edit(2,lambda d:d["powers"]["Wall velocity"].__setitem__(1,.02))
        with self.assertRaisesRegex(ValueError,"generic port"):
            comparator.rank_port_traces(
                self.directory,4,self.times,self.fingerprint
            )

    def test_wrong_control_mesh_fingerprint_fails(self):
        self.edit(3,lambda d:d.update(source_fingerprint="wrong"))
        with self.assertRaisesRegex(ValueError,"port metadata"):
            comparator.rank_port_traces(
                self.directory,4,self.times,self.fingerprint
            )

    def test_wrong_stage_clock_fails(self):
        self.edit(1,lambda d:d.update(times=[0.,.003,.004]))
        with self.assertRaisesRegex(ValueError,"physical clock"):
            comparator.rank_port_traces(
                self.directory,4,self.times,self.fingerprint
            )

    def test_generic_power_is_not_h_weighted_physical_flux(self):
        self.edit(0,lambda d:d.update(
            physical_h_weighted_balance_is_separate=False
        ))
        with self.assertRaisesRegex(ValueError,"port metadata"):
            comparator.rank_port_traces(
                self.directory,4,self.times,self.fingerprint
            )

    def test_nontrivial_physical_power_is_required(self):
        for rank in range(4):
            self.edit(rank,lambda d:d.update(powers={
                k:[0.,0.,0.] for k in d["powers"]
            }))
        with self.assertRaisesRegex(ValueError,"vacuous"):
            comparator.rank_port_traces(
                self.directory,4,self.times,self.fingerprint
            )

    def test_all_rank_specific_output_writers_are_unique(self):
        s=(ROOT/"sandbox"/"probe_v130_rotating_controlled_ic.py").read_text(
            encoding="utf-8"
        )
        runner=(ROOT/"scripts"/"run-v130-issue6-rotating-mpi4.sh").read_text(
            encoding="utf-8"
        )
        self.assertIn('system.plot_Hamiltonian(with_powers=True, save_figure=False)',s)
        self.assertIn('f"mpi_port_powers_rank{mpi.rank}.json"',s)
        self.assertIn('if mpi.rank == 0:',s)
        self.assertIn('args.out / "mpi_physical_powers.png"',s)
        self.assertIn('SCRIMP_V130_ISSUE6_AUDIT_PORTS=1',runner)
        self.assertIn('for mesh in 0.2 0.1',runner)
        self.assertIn('cont2="$(run_native',runner)
        self.assertIn('cont4="$(run_native',runner)
        self.assertIn('restart14="$(run_native',runner)
        self.assertIn('restart41="$(run_native',runner)


if __name__=="__main__":
    unittest.main()
