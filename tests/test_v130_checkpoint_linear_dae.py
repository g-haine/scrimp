# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See LICENSE for license information.

"""#29 host-only synthetic guards for SECOND linear 2-D wave DAE model.

Synthetic numbers test comparison/report rejection, NEVER count as
native GetFEM or PETSc evidence.
"""

import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"sandbox"))
from compare_v130_checkpoint_linear_dae import compare, segment

_spec=importlib.util.spec_from_file_location(
    "_scrimp_linear_dae_pure_checkpoint",
    ROOT/"scrimp"/"utils"/"dae_checkpoint.py",
)
ckpt=importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(ckpt)


class LinearDAECheckpointHostTests(unittest.TestCase):
    def setUp(self):
        tmp=tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.root=Path(tmp.name)
        (self.root/"manifest.json").write_text(json.dumps({
            "schema":"scrimp-v130-issue29-linear-DAE-runner-v1",
            "original_example":"wave_2d_linear_DAE",
            "cold_state_restart_only":True,
            "physical_final_time":.024,
        }),encoding="utf-8")
        detail={"ndof":8, "synthetic_linear_DAE":"not_natively_executed"}
        self.provenance={
            "detail":detail,
            "fingerprint_sha256":ckpt.sha256(ckpt.canonical(detail)),
        }
        self.z0=np.zeros(8)
        self.zstar=np.array([.01,.02,.03,.04,.05,.06,.07,.08])
        self.zend=self.zstar*1.2
        for r in (1,2):
            self.fixture("seed",r)
            self.fixture("continuous",r)
            checkpoint=self.root/f"wave_source_mpi{r}"/"accepted_state.npz"
            ckpt.write_checkpoint(
                checkpoint,z=self.zstar,t=.02,
                provenance=self.provenance,ts_type="cn",accepted=True,
                extra={
                    "source_MPI_ranks":r,
                    "source":"original_examples_wave_linear_DAE",
                    "physical_time_origin":"absolute",
                    "control_law":"0.2*sin(25*t)*sin(pi*y)",
                },
            )
        for source in (1,2):
            for target in (1,2):
                self.fixture("resume",target,source)

    def fixture(self,phase,ranks,source=None):
        if phase=="seed":
            a,b=0.,.02
            z0,z1=self.z0,self.zstar
            h0,h1,w=0.,.003,.004
        elif phase=="continuous":
            a,b=0.,.024
            z0,z1=self.z0,self.zend
            h0,h1,w=0.,.004,.0055
        else:
            a,b=.02,.024
            z0,z1=self.zstar,self.zend
            h0,h1,w=.003,.004,.0015
        name=(f"wave_{phase}_mpi{ranks}" if phase!="resume" else
              f"wave_restart_from{source}_to{ranks}")
        path=self.root/name
        path.mkdir(parents=True,exist_ok=True)
        n=int(round((b-a)/.002))
        times=np.linspace(a,b,n+1)
        states=np.asarray([
            z0+(z1-z0)*k/n for k in range(n+1)
        ],dtype=float)
        H=np.linspace(h0,h1,n+1)
        power=np.full(n+1,w/(b-a))
        work=np.linspace(0.,w,n+1)
        result={
            "schema":"scrimp-v130-issue29-linear-DAE-public-API-v1",
            "status":"NATIVE_LINEAR_DAE_SEGMENT_NOT_YET_CROSS_MPI_COMPARED",
            "phase":phase,"MPI":ranks,"mesh_h":.5,
            "t0":a,"tf":b,"dt":.002,"step_count":n,
            "PETSc_reason":1,"PETSc_rejected":0,"SNES_failures":0,
            "ndof":8,"dynamic_ndof":4,"algebraic_ndof":4,
            "fingerprint_sha256":self.provenance["fingerprint_sha256"],
            "initial_full_sha256":ckpt.sha256(states[0].tobytes()),
            "final_full_sha256":ckpt.sha256(states[-1].tobytes()),
            "not_SWE_model":True,
            "no_integrator_history_claimed":True,
            "production_API_used":phase in ("seed","resume"),
            "original_equations":{
                "endpoint_Fa_residuals":[1e-13]*(n+1),
                "full_step_TSCN_residuals":[1e-12]*n,
            },
            "Hamiltonian_initial":h0,"Hamiltonian_final":h1,
            "boundary_work":w,"dH_minus_W":h1-h0-w,
            "initialization_report":{
                "status":"COMPLETE_DAE_STATE_COLD_RESTORED_IN_DPHS",
                "all_GetFEM_and_GlobalVec_coefficients_restored_bitwise":True,
                "no_pseudo_BDF_or_projection":True,
                "PETSc_TS_history_restored":False,
                "original_algebraic_residual_l2":1e-13,
            } if phase=="resume" else {
                "status":"ORIGINAL_ANALYTICALLY_COMPATIBLE_ZERO_INITIAL_DAE"
            },
        }
        (path/"linear_report.json").write_text(
            json.dumps(result),encoding="utf-8"
        )
        np.savez_compressed(
            path/"linear_states.npz",
            times=times,states=states,H=H,work=work,port_power=power,
        )
        (self.root/f"{name}.log").write_text(
            f"LINEAR DAE PHASE={phase} MPI={ranks} steps={n} "
            "synthetic fixture only",encoding="utf-8"
        )
        return path,result,states

    def test_synthetic_complete_four_MPI_restarts(self):
        result=compare(self.root)
        self.assertEqual(
            result["status"],"NATIVE_LINEAR_PHDae_DPHS_COLD_RESTART_PASSED"
        )
        self.assertEqual(len(result["cases"]),4)

    def test_rejects_missing_cross_MPI_run(self):
        (self.root/"wave_restart_from1_to2"/"linear_states.npz").unlink()
        with self.assertRaisesRegex(ValueError,"Missing native linear DAE"):
            compare(self.root)

    def test_rejects_restarted_algebraic_coefficient_change(self):
        path,r,states=self.fixture("resume",2,1)
        states[0,-1]+=.001
        with np.load(path/"linear_states.npz",allow_pickle=False) as x:
            data={name:np.asarray(x[name]).copy() for name in x.files}
        data["states"]=states
        np.savez_compressed(path/"linear_states.npz",**data)
        with self.assertRaisesRegex(ValueError,"full FE state checksum"):
            compare(self.root)

    def test_rejects_original_getfem_Fa_violation(self):
        path,r,_=self.fixture("resume",1,2)
        r["original_equations"]["endpoint_Fa_residuals"][-1]=.02
        (path/"linear_report.json").write_text(json.dumps(r),encoding="utf-8")
        with self.assertRaisesRegex(ValueError,"algebraic constraint failure"):
            compare(self.root)

    def test_rejects_wrong_original_boundary_power(self):
        path,r,_=self.fixture("continuous",1)
        with np.load(path/"linear_states.npz",allow_pickle=False) as x:
            data={name:np.asarray(x[name]).copy() for name in x.files}
        data["port_power"][-1]+=1.
        np.savez_compressed(path/"linear_states.npz",**data)
        with self.assertRaisesRegex(ValueError,"boundary-work integration"):
            compare(self.root)

    def test_rejects_modified_physical_checkpoint_control(self):
        path=self.root/"wave_source_mpi1"/"accepted_state.npz"
        with np.load(path,allow_pickle=False) as x:
            z=x["state"].copy()
            header=json.loads(x["manifest"].tobytes().decode("utf-8"))
        header["extra"]["control_law"]="wrong"
        np.savez_compressed(
            path,state=z,
            manifest=np.frombuffer(ckpt.canonical(header),dtype=np.uint8)
        )
        with self.assertRaisesRegex(ValueError,"physical linear DAE checkpoint"):
            compare(self.root)

    def test_rejects_wrong_model_fingerprint(self):
        path,r,_=self.fixture("resume",2,1)
        # Valid SHA-256 syntax but WRONG actual GetFEM model identity.
        r["fingerprint_sha256"]="f"*64
        (path/"linear_report.json").write_text(json.dumps(r),encoding="utf-8")
        with self.assertRaisesRegex(ValueError,"physically different model"):
            compare(self.root)

    def test_rejects_wrong_seed_model_fingerprint(self):
        path,r,_=self.fixture("seed",2)
        r["fingerprint_sha256"]="f"*64
        (path/"linear_report.json").write_text(json.dumps(r),encoding="utf-8")
        with self.assertRaisesRegex(ValueError,"physically different model"):
            compare(self.root)

    def test_rejects_projection_without_model_mismatch(self):
        # Distinguish a genuine q0 projection from changed source physics.
        path,r,_=self.fixture("resume",2,1)
        r["initialization_report"]["no_pseudo_BDF_or_projection"]=False
        (path/"linear_report.json").write_text(json.dumps(r),encoding="utf-8")
        with self.assertRaisesRegex(ValueError,"projected/changed q0"):
            compare(self.root)


if __name__=="__main__":
    unittest.main()
