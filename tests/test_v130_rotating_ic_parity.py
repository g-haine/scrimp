# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Synthetic negative tests ONLY for the native-data comparator.

Actual acceptances must still come from original PETSc/GetFEM runs.
"""

import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np


ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/"sandbox"/"compare_v130_rotating_ic_parity.py"
spec=importlib.util.spec_from_file_location("_rotating_native_parity",SOURCE)
parity=importlib.util.module_from_spec(spec)
spec.loader.exec_module(parity)


class NativeRotatingCheckpointComparatorTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        root=Path(self.tmp.name)
        self.seed=root/"seed"
        self.continuous=root/"continuous"
        self.resume=root/"resume"
        z=np.array([
            [5.,0.,0.],[5.1,.1,0.],[5.2,.2,0.],[5.3,.3,0.]
        ])
        H=np.array([100.,101.,102.,103.])
        W=np.array([0.,1.2,2.3,3.4])
        D=np.array([0.,.01,.02,.03])
        T=np.array([0.,.01,.02,.03])
        V=np.full(4,5.)
        port=np.array([0.,2.,3.,4.])
        traces={
            "H":H,"volume":V,
            "p_normal":np.zeros(4),
            "p_viscous":port,
            "d_grad":np.full(4,.2),
            "d_div":np.full(4,.1),
        }
        for directory, phase, mpi, indices in (
            (self.seed,"seed",1,[0,1,2]),
            (self.continuous,"seed",2,[0,1,2,3]),
            (self.resume,"resume",2,[2,3]),
        ):
            directory.mkdir()
            t=T[indices]
            zc=z[indices]
            Hc=H[indices]
            Wc=W[indices]-W[indices][0]
            Dc=D[indices]-D[indices][0]
            balance=Hc-Hc[0]+Dc-Wc
            values={
                "times":t,"states":zc,"full_state_start":zc[0],
                "full_state_final":zc[-1],
                "boundary_work":Wc,"viscous_dissipation":Dc,
                "power_balance":balance,
                "volume_balance":np.zeros(len(t)),
            }
            values.update({k:v[indices] for k,v in traces.items()})
            np.savez_compressed(directory/"initial_state_trace.npz",**values)
            np.savez_compressed(directory/"accepted_state.npz",state=zc[-1])
            record={
                "schema":"scrimp-v130-rotating-controlled-ic-native-v1",
                "status":"ACCEPTED_COMPATIBLE_ROTATING_DAE_INITIAL_STATE_NOT_STEADY",
                "phase":phase,"qualification_mode":"restart_parity",
                "physical_relative_power_gate_applied":False,
                "physical_energy_defect_retained_not_modified":True,
                "gamma":.25,"mesh":.2,"MPI":mpi,
                "physical_time_initial":float(t[0]),
                "physical_time_final":float(t[-1]),"dt":.01,
                "ndof":3,"accepted_steps":len(t)-1,"PETSc_reason":1,
                "PETSc_rejected":0,"SNES_failures":0,
                "original_Fa_max_l2":1e-14,
                "original_TSCN_max_l2":1e-12,
                "minimum_positive_height_FE_coefficient":5.,
                "original_volume_balance_defect":0.,
                "stationary_rotating_state_certified":False,
                "original_FE_fingerprint_sha256":"identical_hash",
                "original_power_defect_over_exchange":.4,
                "initialization":{
                    "no_pseudo_BDF_or_projection":True,
                    "all_GetFEM_and_GlobalVec_coefficients_restored_bitwise":True,
                    "PETSc_TS_history_restored":False,
                },
            }
            (directory/"initial_state_report.json").write_text(
                json.dumps(record),encoding="utf-8"
            )

    def test_parity_passes_even_for_large_common_energy_defect(self):
        report=parity.compare(self.seed,self.continuous,self.resume)
        self.assertEqual(
            report["status"],"PASS_RESTART_EQUIVALENT_TO_UNINTERRUPTED"
        )
        self.assertEqual(report["maximum_full_FE_trajectory_L2_difference"],0.)
        self.assertEqual(report["physical_energy_defect_observed_not_gated"],True)
        self.assertGreater(abs(report["continuous_final_energy_defect"]),.1)

    def test_rejects_nontrivial_intermediate_state_change(self):
        path=self.continuous/"initial_state_trace.npz"
        with np.load(path,allow_pickle=False) as archive:
            data={k:archive[k] for k in archive.files}
        data["states"][1,1]+=.01
        np.savez_compressed(path,**data)
        with self.assertRaisesRegex(ValueError,"Full FE trajectory diverged"):
            parity.compare(self.seed,self.continuous,self.resume)

    def test_rejects_different_physical_work_even_when_states_equal(self):
        path=self.resume/"initial_state_trace.npz"
        with np.load(path,allow_pickle=False) as archive:
            data={k:archive[k] for k in archive.files}
        data["boundary_work"][-1]+=.02
        np.savez_compressed(path,**data)
        with self.assertRaisesRegex(ValueError,"boundary_work"):
            parity.compare(self.seed,self.continuous,self.resume)

    def test_rejects_broken_bitwise_checkpoint_at_pause(self):
        path=self.seed/"accepted_state.npz"
        with np.load(path,allow_pickle=False) as archive:
            data={k:archive[k] for k in archive.files}
        data["state"][1]+=1e-9
        np.savez_compressed(path,**data)
        with self.assertRaisesRegex(ValueError,"checkpoint differs"):
            parity.compare(self.seed,self.continuous,self.resume)

    def test_rejects_wrong_provenance(self):
        path=self.resume/"initial_state_report.json"
        record=json.loads(path.read_text())
        record["original_FE_fingerprint_sha256"]="different"
        path.write_text(json.dumps(record))
        with self.assertRaisesRegex(ValueError,"Different model/mesh"):
            parity.compare(self.seed,self.continuous,self.resume)

    def test_rejects_pseudo_initialized_restart(self):
        path=self.resume/"initial_state_report.json"
        record=json.loads(path.read_text())
        record["initialization"]["no_pseudo_BDF_or_projection"]=False
        path.write_text(json.dumps(record))
        with self.assertRaisesRegex(ValueError,"Cold resume"):
            parity.compare(self.seed,self.continuous,self.resume)


if __name__=="__main__":
    unittest.main()
