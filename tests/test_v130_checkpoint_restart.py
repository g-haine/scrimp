# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See LICENSE for license information.

"""#29 synthetic continuous-vs-split checkpoint comparator negative gates.

Host fixtures are NEVER substituted for native original PETSc evidence.
"""

import json
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"sandbox"))
from compare_v130_checkpoint_restart import compare
from validate_v130_checkpoint_case import verify_one
from v130_checkpoint_state import (
    canonical, sha256, write_checkpoint,
)


class CheckpointReplayHostTest(unittest.TestCase):
    def setUp(self):
        temp=tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.root=Path(temp.name)
        (self.root/"manifest.json").write_text(json.dumps({
            "schema":"scrimp-v130-issue29-runner-v1",
            "original_physics_untouched":True,
            "cold_state_restart_only":True,
            "source_spinup_time":.02,
            "physical_final_time":.024,
        }),encoding="utf-8")
        detail={"ndof":6,"model":"synthetic original-P2-like placeholder"}
        self.model={"detail":detail,
                    "fingerprint_sha256":sha256(canonical(detail))}
        self.z0=np.array([5.,5.,5.,0.,0.,0.])
        self.zstar=np.array([5.,5.,5.,.10,.10,.10])
        self.zend=np.array([5.,5.,5.,.12,.12,.12])
        for mesh in (.2,.1):
            label=f"h{str(mesh).replace('.', 'p')}"
            for r in (1,2):
                self.fixture(mesh,"seed",r)
                self.fixture(mesh,"continuous",r)
                path=self.root/f"{label}_source_mpi{r}"/"accepted_state.npz"
                write_checkpoint(
                    path,z=self.zstar,t=.02,
                    provenance=self.model,ts_type="cn",
                    accepted=True,
                    extra={
                        "source_mpi_ranks":r,
                        "nominal_mesh_size":mesh,
                        "physical_control_clock":"absolute",
                        "physical_rotation_target":.2,
                    }
                )
            for source in (1,2):
                for target in (1,2):
                    self.fixture(mesh,"resume",target,source)

    def fixture(self,mesh,phase,rank,source=None):
        label=f"h{str(mesh).replace('.', 'p')}"
        if phase=="resume":
            folder=self.root/f"{label}_restart_from{source}_to{rank}"
            begin=.02;end=.024;z0=self.zstar;z1=self.zend
            H0,H1=4.01,4.013;W=.0035;D=.0005
        elif phase=="seed":
            folder=self.root/f"{label}_seed_mpi{rank}"
            begin=0.;end=.02;z0=self.z0;z1=self.zstar
            H0,H1=4.,4.01;W=.012;D=.002
        else:
            folder=self.root/f"{label}_continuous_mpi{rank}"
            begin=0.;end=.024;z0=self.z0;z1=self.zend
            H0,H1=4.,4.013;W=.0155;D=.0025
        folder.mkdir(parents=True,exist_ok=True)
        n=int(round((end-begin)/.002))
        ts=np.linspace(begin,end,n+1)
        Z=np.array([z0+(z1-z0)*k/n for k in range(n+1)])
        H=np.linspace(H0,H1,n+1)
        wall=np.linspace(0.,W,n+1)
        diss=np.linspace(0.,D,n+1)
        balance=H-H0+diss-wall
        arrays={
            "times":ts,"states":Z,"H":H,
            "volume":np.full(n+1,5.),
            "volume_balance":np.zeros(n+1),
            "supplied":wall,"dissipated":diss,
            "p_normal":np.zeros(n+1),
            "p_viscous":np.zeros(n+1),
            "d_grad":np.zeros(n+1),
            "d_div":np.zeros(n+1),
            "balance":balance,
        }
        normdata={
            "maximum_actual_original_TSCN_residual":1e-12,
            "maximum_actual_original_Fa_residual":1e-13,
            "full_step_TSCN_residuals":[1e-12]*n,
            "full_step_Fa_residuals":[1e-13]*n,
        }
        initial=(
            {"status":"SANDBOX_CHECKPOINT_COLD_RESTORED_NOT_PRODUCTION",
             "getfem_globalvec_published_bitwise":True,
             "pseudo_or_bdf_initialization_performed":False,
             "original_algebraic_residual_l2":1e-13,
             "checkpoint_sha256":sha256(z0.tobytes())}
            if phase=="resume"
            else {"strategy":"sandbox_exact_stationary_spinup_only"}
        )
        record={
            "schema":"scrimp-v130-issue29-native-cold-restart-pilot-v1",
            "status":"SANDBOX_NATIVE_SEGMENT_DIAGNOSTIC_NOT_PRODUCTION",
            "phase":phase,"mesh":mesh,"mpi_ranks":rank,
            "model_fingerprint_sha256":self.model["fingerprint_sha256"],
            "complete_final_state_sha256":sha256(Z[-1].tobytes()),
            "source_time":begin,"actual_final_time":end,
            "actual_ts_steps":n,"actual_ts_reason":1,
            "actual_ts_rejections":0,"actual_ts_snes_failures":0,
            "cold_restart_not_identical_TS_history":True,
            "general_high_index_initializer_qualified":False,
            "production_checkpoint_qualified":False,
            "initialization_report":initial,
            "exact_original_equations":normdata,
            "hamiltonian_change":H1-H0,
            "original_boundary_work":W,
            "original_viscous_dissipation":D,
            "original_H_plus_D_minus_W":float(balance[-1]),
        }
        (folder/"segment_report.json").write_text(
            json.dumps(record),encoding="utf-8")
        np.savez_compressed(folder/"segment_states.npz",**arrays)
        (self.root/f"{folder.name}.log").write_text(
            f"CHECKPOINT PHASE={phase} mesh={mesh:g} MPI={rank} "
            f"steps={n} synthetic host fixture only\\n",
            encoding="utf-8",
        )
        return record,arrays

    def test_synthetic_individual_case_resumption_gate(self):
        for mesh in (.2,.1):
            for ranks in (1,2):
                verify_one(self.root,mesh,"seed",ranks)
                verify_one(self.root,mesh,"continuous",ranks)
                for source in (1,2):
                    verify_one(
                        self.root,mesh,"resume",ranks,source_ranks=source
                    )

    def test_rejects_existing_case_without_expected_success_log(self):
        label=self.root/"h0p2_seed_mpi1.log"
        label.write_text(
            "CHECKPOINT PHASE seed (old broken grep test)\\n",
            encoding="utf-8"
        )
        with self.assertRaisesRegex(ValueError,"success marker absent"):
            verify_one(self.root,.2,"seed",1)

    def test_rejects_existing_seed_from_wrong_checkpoint_source(self):
        path=self.root/"h0p1_source_mpi2"/"accepted_state.npz"
        with np.load(path, allow_pickle=False) as data:
            arr=data["state"].copy()
            manifest=json.loads(data["manifest"].tobytes().decode("utf-8"))
        manifest["extra"]["source_mpi_ranks"]=1
        np.savez_compressed(
            path,state=arr,
            manifest=np.frombuffer(canonical(manifest),dtype=np.uint8)
        )
        with self.assertRaisesRegex(ValueError,"source/provenance"):
            verify_one(self.root,.1,"seed",2)

    def test_rejects_existing_restart_without_complete_source_case(self):
        path=self.root/"h0p2_source_mpi1"/"accepted_state.npz"
        path.unlink()
        with self.assertRaisesRegex(ValueError,"accepted-state file"):
            verify_one(self.root,.2,"resume",2,source_ranks=1)

    def test_synthetic_complete_cross_rank_replay(self):
        response=compare(self.root)
        self.assertEqual(
            response["status"],
            "NATIVE_COLD_CHECKPOINT_SPLIT_VS_CONTINUOUS_PILOT_PASSED"
        )
        self.assertFalse(response["production_DPHS_checkpoint_API_qualified"])
        self.assertFalse(response["history_exact_BDF_IMEX_restoration_qualified"])

    def test_rejects_corrupted_seed_checkpoint_state(self):
        file=self.root/"h0p2_source_mpi1"/"accepted_state.npz"
        with np.load(file,allow_pickle=False) as data:
            z=data["state"].copy()
            manifest=data["manifest"].copy()
        z[0]+=1.
        np.savez_compressed(file,state=z,manifest=manifest)
        with self.assertRaisesRegex(ValueError,"Damaged checkpoint"):
            compare(self.root)

    def test_rejects_cross_MPI_restart_mutating_one_algebraic_coefficient(self):
        record,arrays=self.fixture(.1,"resume",2,1)
        arrays["states"][0,4]+=.01
        folder=self.root/"h0p1_restart_from1_to2"
        np.savez_compressed(folder/"segment_states.npz",**arrays)
        with self.assertRaisesRegex(ValueError,"silently changed"):
            compare(self.root)

    def test_rejects_changed_full_original_solution(self):
        record,arrays=self.fixture(.2,"resume",1,2)
        arrays["states"][-1,5]+=.01
        folder=self.root/"h0p2_restart_from2_to1"
        np.savez_compressed(folder/"segment_states.npz",**arrays)
        with self.assertRaisesRegex(ValueError,"full FE state mismatch"):
            compare(self.root)

    def test_rejects_wrong_absolute_restarted_time(self):
        record,arrays=self.fixture(.2,"resume",1,1)
        record["source_time"]=0.
        path=self.root/"h0p2_restart_from1_to1"/"segment_report.json"
        path.write_text(json.dumps(record),encoding="utf-8")
        with self.assertRaisesRegex(ValueError,"segment not genuinely accepted"):
            compare(self.root)

    def test_rejects_unqualified_initial_condition_projection(self):
        record,arrays=self.fixture(.1,"resume",2,2)
        record["initialization_report"][
            "pseudo_or_bdf_initialization_performed"
        ]=True
        path=self.root/"h0p1_restart_from2_to2"/"segment_report.json"
        path.write_text(json.dumps(record),encoding="utf-8")
        with self.assertRaisesRegex(ValueError,"reprojected"):
            compare(self.root)

    def test_rejects_real_endpoint_algebraic_violation(self):
        record,arrays=self.fixture(.2,"resume",2,1)
        record["exact_original_equations"]["full_step_Fa_residuals"][-1]=1e-3
        path=self.root/"h0p2_restart_from1_to2"/"segment_report.json"
        path.write_text(json.dumps(record),encoding="utf-8")
        with self.assertRaisesRegex(ValueError,"original TSCN/algebraic"):
            compare(self.root)

    def test_rejects_physical_Hamiltonian_jump_at_restart(self):
        record,arrays=self.fixture(.1,"resume",1,2)
        arrays["H"]+=.01
        arrays["balance"]=arrays["H"]-arrays["H"][0]+(
            arrays["dissipated"]-arrays["supplied"]
        )
        folder=self.root/"h0p1_restart_from2_to1"
        np.savez_compressed(folder/"segment_states.npz",**arrays)
        with self.assertRaisesRegex(ValueError,"Hamiltonian jumped"):
            compare(self.root)

    def test_rejects_missing_native_MPI_cross_rank_case(self):
        path=self.root/"h0p2_restart_from2_to1"/"segment_states.npz"
        path.unlink()
        with self.assertRaisesRegex(ValueError,"Missing real native"):
            compare(self.root)


if __name__=="__main__":
    unittest.main()
