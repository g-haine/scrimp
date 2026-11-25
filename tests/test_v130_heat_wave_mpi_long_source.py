# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Source and NEGATIVE-control fixtures for the native Heat–Wave MPI pilot.

Synthetic fixture arrays ONLY check the rejection logic of the HOST
verifier. They are not native GetFEM/PETSc proof; that requires archived
MPI1/MPI2+long executions and original scientific quadratures.
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
COLLECTOR=ROOT/"sandbox/measure_v130_hw_dual_facet_physics.py"
RUNNER=ROOT/"scripts/run-v130-heat-wave-mpi-long-physics.sh"
CHECK=ROOT/"scripts/verify-v130-heat-wave-mpi-long-physics.py"
EXAMPLE=ROOT/"examples/heat_wave.py"


class HeatWaveNativeMPIAndLongSourceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.probe=PROBE.read_text()
        cls.collector=COLLECTOR.read_text()
        cls.runner=RUNNER.read_text()
        cls.example=EXAMPLE.read_text()
        cls.check_source=CHECK.read_text()
        spec=importlib.util.spec_from_file_location("heat_wave_cn_verifier",CHECK)
        cls.verify=importlib.util.module_from_spec(spec)
        spec.loader.exec_module(cls.verify)

    def test_real_model_and_scalar_physics_only_output_switch(self):
        for needle in (
            'physics_only=False',
            'if physics_only:',
            'return hw',
            'hw.solve()',
            'hw.hamiltonian.add_term(term)',
            'enforce_zero_dirichlet_initial_values(',
        ):
            self.assertIn(needle,self.example)
        self.assertLess(self.example.index('hw.solve()'),
                        self.example.index('if physics_only:'))
        self.assertLess(self.example.index('if physics_only:'),
                        self.example.index('hw.export_to_pv("p")'))
        self.assertIn('physics_only=True',self.probe)
        self.assertNotIn('DPHS.init_step =',self.probe)
        self.assertNotIn('DPHS.solve =',self.probe)
        self.assertIn('initial_wave_boundary="clamp"',self.probe)

    def test_consecutive_accepted_CN_required_without_subsampled_midpoints(self):
        for needle in (
            "expected_dt=0.001, expected_steps=4,",
            "times.shape!=(expected_steps+1,)",
            "np.arange(expected_steps+1)*expected_dt",
            "expected_dt=args.dt,expected_steps=args.steps",
            "save_step=args.dt",
            "system.last_ts_steps!=args.steps",
            "system.last_ts_rejections!=0",
            "CN_per_step_heat_balance_error",
            "CN_per_step_wave_balance_error",
            "CN_P_interface_physical",
            "CN_cumulative_total_balance_error",
            "monotone_total_energy_at_every_accepted_step",
        ):
            self.assertTrue(needle in self.probe or needle in self.collector,
                            msg=needle)
        self.assertIn("raise RuntimeError(",self.probe)
        self.assertIn("all_MPI_ranks_agree_exactly",self.probe)

    def test_native_runner_phases_image_and_archive(self):
        for needle in (
            'run_case mpi1_short 1 .001 4 gaussian',
            'run_case mpi2_short 2 .001 4 gaussian',
            'run_case mpi1_long 1 "$long_dt" "$long_steps" "$long_profile"',
            "SCRIMP_V130_HW_LONG_STEPS:-3000",
            "SCRIMP_V130_HW_LONG_DT:-0.005",
            "SCRIMP_V130_HW_LONG_PROFILE:-gaussian",
            'if [[ "$short_rc" == 0 ]]; then',
            "scripts/verify-v130-heat-wave-mpi-long-physics.py",
            "--store-full-states",
            "--volume \"$root/scrimp:/workspace/scrimp/scrimp:ro\"",
            "scripts/check-v130-published-immutable.py",
            'with zipfile.ZipFile(archive,"x"',
        ):
            self.assertIn(needle,self.runner)
        self.assertNotIn('examples/published/:',self.runner)
        # All three Python entry points parse under a real AST:
        for path in (PROBE,CHECK,COLLECTOR,EXAMPLE):
            ast.parse(path.read_text(),filename=str(path))

    @staticmethod
    def fixture(root,*,long=False):
        (root/"manifest.json").write_text(
            json.dumps({"source":"synthetic_fixture_not_scientific_proof",
                        "long_dt":.005,"long_steps":12,
                        "long_profile":"gaussian"})
        )
        steps=12 if long else 4
        dt=.005 if long else .001
        t=np.arange(steps+1,dtype=float)*dt
        H=1.-.1*t
        z=np.broadcast_to(np.arange(7689,dtype=float)[None,:],
                          (steps+1,7689)).copy()
        z[:,0]+=t
        arrays={
            "t":t, "H":H,"z0":z[0],"zf":z[-1],
            "phys_H_heat":.5*H,
            "phys_H_wave":.5*H,
            "phys_H_components":H,
            "phys_CN_D_heat_flux":np.full(len(t),.1),
            "phys_CN_P_interface_physical":np.zeros(len(t)),
            "phys_CN_P_outer":np.zeros(len(t)),
            "phys_CN_per_step_heat_balance_error":np.zeros(len(t)),
            "phys_CN_per_step_wave_balance_error":np.zeros(len(t)),
            "phys_CN_per_step_H_balance_error":np.zeros(len(t)),
            "phys_CN_cumulative_total_balance_error":np.zeros(len(t)),
        }
        for label,ranks in (("mpi1_short",1),("mpi2_short",2),
                            *((("mpi1_long",1),) if long else ())):
            folder=root/label
            folder.mkdir()
            (root/f"{label}.exitcode").write_text("0")
            for rank in range(ranks):
                meta={
                    "schema":"scrimp-v130-heat-wave-MPI-and-long-CN-physical-v1",
                    "status":"PASS_NATIVE_STRICT_DAE_CN_ORIENTED_POWER",
                    "rank":rank,"MPI_ranks":ranks,
                    "steps":steps if label=="mpi1_long" else 4,
                    "dt":dt if label=="mpi1_long" else .001,
                    "profile":"gaussian","DOFs":7689,
                    "saved_all_accepted_steps":True,
                    "physical_gates":{"genuine":True},
                    "parallel_ParaView_export_qualified":False,
                    "whispering_gallery_mode_qualified":False,
                    "all_rank_native_assembly_and_quadratures":True,
                    "initial_wave_boundary":{
                        "method":"clamp",
                        "boundary_compatibility_enforced":True,
                        "after_boundary_l2_coefficients":0.,
                    },
                }
                (folder/f"rank{rank}.json").write_text(json.dumps(meta))
                a=arrays if label=="mpi1_long" else dict(arrays,
                    t=np.arange(5,dtype=float)*.001,
                    H=1.-.1*np.arange(5,dtype=float)*.001,
                )
                if label!="mpi1_long":
                    a={
                        **a,
                        "z":z[:5],
                        "z0":z[0],
                        "zf":z[4],
                        "phys_H_heat":.5*a["H"],
                        "phys_H_wave":.5*a["H"],
                        "phys_H_components":a["H"],
                        **{
                            key:np.asarray(val[:5])
                            for key,val in arrays.items()
                            if key.startswith("phys_CN_")
                        },
                    }
                np.savez_compressed(folder/f"rank{rank}.npz",**a)

    def test_host_verifier_rejects_corrupt_physical_quadratures_and_missing_rank(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp)
            self.fixture(root)
            ok=self.verify.audit(root,"short")
            self.assertEqual(ok["status"],"PASS_HEAT_WAVE_MPI1_MPI2_SHORT_PHYSICAL")
            self.assertTrue(ok["MPI2_rank0_rank1_bitwise_identical"])
            # Corrupt one true FE rank, not simply its metadata.
            f=root/"mpi2_short/rank1.npz"
            with np.load(f,allow_pickle=False) as z:
                vals={k:z[k].copy() for k in z.files}
            vals["z"][2,0]+=0.01
            np.savez_compressed(f,**vals)
            fail=self.verify.audit(root,"short")
            self.assertTrue(fail["status"].startswith("FAIL_"))
            # Restore and independently corrupt the cumulative balance.
            vals["z"][2,0]-=0.01
            np.savez_compressed(f,**vals)
            f=root/"mpi1_short/rank0.npz"
            with np.load(f,allow_pickle=False) as z:
                vals={k:z[k].copy() for k in z.files}
            vals["phys_CN_cumulative_total_balance_error"][-1]=.1
            np.savez_compressed(f,**vals)
            fail=self.verify.audit(root,"short")
            self.assertTrue(fail["status"].startswith("FAIL_"))

    def test_long_verifier_rejects_incomplete_accepted_step_clock(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp)
            self.fixture(root,long=True)
            f=root/"mpi1_long/rank0.npz"
            with np.load(f,allow_pickle=False) as z:
                vals={k:z[k].copy() for k in z.files}
            vals["t"][6]+=.001
            np.savez_compressed(f,**vals)
            failed=self.verify.audit(root,"complete")
            self.assertTrue(failed["status"].startswith("FAIL_"))
            self.assertTrue(any("time" in s.lower() or "clock" in s.lower()
                                for s in failed["errors"]))


if __name__=="__main__":
    unittest.main()
