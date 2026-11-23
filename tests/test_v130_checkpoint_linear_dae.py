# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See LICENSE for license information.

"""#29 host-only synthetic guards for SECOND linear 2-D wave DAE model.

Synthetic numbers test comparison/report rejection, NEVER count as
native GetFEM or PETSc evidence.
"""

import ast
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
            h0,h1,w=0.,.003,.003
        elif phase=="continuous":
            a,b=0.,.024
            z0,z1=self.z0,self.zend
            h0,h1,w=0.,.004,.004
        else:
            a,b=.02,.024
            z0,z1=self.zstar,self.zend
            h0,h1,w=.003,.004,.001
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
        # Native Port.compute records edge/CN midpoint powers for
        # k>=1, but power[0] is instantaneous at the initial state.
        # Critically, the initial instantaneous restart sample can
        # be unrelated to the previous step's power. The old v1
        # report nonetheless incorrectly integrated all samples
        # with trapezoids and the wrong sign; the comparator must
        # transparently correct it WITHOUT rewriting this evidence.
        power=np.full(n+1,-w/(b-a))
        power[0]=0.173 if phase=="resume" else 0.
        work=np.r_[0.,np.cumsum(
            .5*np.diff(times)*(power[:-1]+power[1:])
        )]
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
            "boundary_work":float(work[-1]),
            "dH_minus_W":float(h1-h0-work[-1]),
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

    def test_original_IFunction_uses_configured_PETSc_index_width(self):
        # Regression for the first real native wave-DAE run:
        # on 64-bit Python, dtype=int defaults to int64, whereas
        # the pinned PETSc image uses 32-bit PetscInt. This test is
        # HOST ONLY and checks the native probe's actual Vec.setValues
        # call sites without importing PETSc/GetFEM.
        probe=ROOT/"sandbox"/"probe_v130_checkpoint_linear_dae.py"
        tree=ast.parse(probe.read_text(encoding="utf-8"))
        funcs=[
            node for node in tree.body
            if isinstance(node,ast.FunctionDef)
            and node.name=="original_IFunction"
        ]
        self.assertEqual(len(funcs),1)
        method=funcs[0]
        owned_assignments=[
            node for node in ast.walk(method)
            if isinstance(node,ast.Assign)
            and any(
                isinstance(target,ast.Name) and target.id=="owned"
                for target in node.targets
            )
        ]
        self.assertEqual(len(owned_assignments),1)
        self.assertEqual(
            ast.unparse(owned_assignments[0].value),
            "np.arange(lo, hi, dtype=PETSc.IntType)"
        )
        setters={
            ast.unparse(node.func.value)
            for node in ast.walk(method)
            if isinstance(node,ast.Call)
            and isinstance(node.func,ast.Attribute)
            and node.func.attr=="setValues"
            and node.args and isinstance(node.args[0],ast.Name)
            and node.args[0].id=="owned"
        }
        self.assertEqual(setters,{"vec","vel"})
        self.assertEqual(
            np.arange(0,3,dtype=np.int32).dtype,np.dtype("int32")
        )

    def test_synthetic_complete_four_MPI_restarts(self):
        result=compare(self.root)
        self.assertEqual(
            result["status"],"NATIVE_LINEAR_PHDae_DPHS_COLD_RESTART_PASSED"
        )
        self.assertEqual(len(result["cases"]),4)

    def test_native_v1_power_recovered_from_original_stage_samples(self):
        # A phase-0 instantaneous port power is NOT the CN edge power
        # used to advance energy at the first resumed time-step.
        seed=segment(self.root,"seed",1)
        resumed=segment(self.root,"resume",2,source=1)
        base=segment(self.root,"continuous",2)
        self.assertGreater(
            abs(float(resumed["power"][0]-seed["power"][-1])), .01
        )
        self.assertGreater(
            abs(float(
                seed["archived_work"][-1]+resumed["archived_work"][-1]
                -base["archived_work"][-1]
            )), 1e-5
        )
        self.assertAlmostEqual(
            float(seed["work"][-1]+resumed["work"][-1]),
            float(base["work"][-1]),places=13
        )
        self.assertLess(resumed["max_step_power_defect"],1e-14)
        self.assertEqual(
            resumed["archived_work_convention"],
            "legacy_trapezoidal_endpoint_v1_INCORRECT"
        )

    def test_v2_native_power_report_uses_CN_midpoint_step_work(self):
        # A future genuine native simulation writes v2 correctly.
        name="wave_restart_from1_to2"
        folder=self.root/name
        report=json.loads(
            (folder/"linear_report.json").read_text(encoding="utf-8")
        )
        with np.load(folder/"linear_states.npz",allow_pickle=False) as x:
            data={k:np.asarray(x[k]).copy() for k in x.files}
        dt=np.diff(data["times"])
        corrected=np.r_[0.,np.cumsum(-dt*data["port_power"][1:])]
        data["work"]=corrected
        report["schema"]="scrimp-v130-issue29-linear-DAE-public-API-v2"
        report["boundary_work_quadrature"]=(
            "minus_CN_step_midpoint_port_power_times_dt"
        )
        report["maximum_discrete_energy_balance_defect"]=float(
            np.max(np.abs(
                np.diff(data["H"])+dt*data["port_power"][1:]
            ))
        )
        report["boundary_work"]=float(corrected[-1])
        report["dH_minus_W"]=float(
            data["H"][-1]-data["H"][0]-corrected[-1]
        )
        np.savez_compressed(folder/"linear_states.npz",**data)
        (folder/"linear_report.json").write_text(
            json.dumps(report),encoding="utf-8"
        )
        accepted=segment(self.root,"resume",2,source=1)
        self.assertEqual(
            accepted["archived_work_convention"],
            "corrected_CN_step_port_power_v2"
        )
        self.assertLess(accepted["max_step_power_defect"],1e-14)
        self.assertEqual(len(compare(self.root)["cases"]),4)

    def test_rejects_wrong_original_CN_step_power_even_if_old_report_updated(self):
        folder=self.root/"wave_resume_dummy" if False else (
            self.root/"wave_restart_from1_to2"
        )
        report=json.loads(
            (folder/"linear_report.json").read_text(encoding="utf-8")
        )
        with np.load(folder/"linear_states.npz",allow_pickle=False) as x:
            data={k:np.asarray(x[k]).copy() for k in x.files}
        data["port_power"][1]*=1.5
        dt=np.diff(data["times"])
        old=np.r_[0.,np.cumsum(
            .5*dt*(data["port_power"][:-1]+data["port_power"][1:])
        )]
        data["work"]=old
        report["boundary_work"]=float(old[-1])
        report["dH_minus_W"]=float(
            data["H"][-1]-data["H"][0]-old[-1]
        )
        np.savez_compressed(folder/"linear_states.npz",**data)
        (folder/"linear_report.json").write_text(
            json.dumps(report),encoding="utf-8"
        )
        with self.assertRaisesRegex(ValueError,"POWER balance failed"):
            segment(self.root,"resume",2,source=1)

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
