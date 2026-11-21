# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See LICENSE for license information.

"""Synthetic host regression/negative tests of native TSCN temporal refinement.

No fabricated fixture here is physical GetFEM/PETSc evidence.
"""

import ast
import json
from pathlib import Path
import tempfile
import unittest
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"sandbox"))
from compare_v130_issue28_cn_refinement import (
    compare, OCI, TF, MESHES, STEP_COUNTS, DOFS, DYNAMIC_DOFS,
)


class TestFixedPhysicalTimeCNRefinement(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.root = Path(tmp.name)
        (self.root/"manifest.json").write_text(json.dumps({
            "schema": "v130-issue28-cn-refinement-v1",
            "image_revision": OCI,
            "original_production_sources_unchanged": True,
            "fixed_physical_final_time": TF,
            "energy_convergence_claimed": False,
        }),encoding="utf-8")
        for mesh in MESHES:
            for count in STEP_COUNTS:
                for ranks in (1,2):
                    doc, traces = self.fixture(mesh,count,ranks)
                    self.store(mesh,count,ranks,doc,traces)

    @staticmethod
    def fixture(mesh,count,ranks):
        n=DOFS[mesh]
        nh=DYNAMIC_DOFS[mesh]//3
        times=np.linspace(0.,TF,count+1)
        initial=np.zeros(n)
        initial[:nh]=5.
        initial[3*nh]=.5
        initial[3*nh+1]=-.5
        last=initial.copy()
        last[nh]=.001*(.03+1/count**2)
        states=np.asarray([initial + (last-initial)*(t/TF)
                           for t in times])
        H=4.+1.2e-6*(times/TF)
        power=np.r_[0.,np.full(count, 1.205e-6/TF)]
        dgrad=np.r_[0.,np.full(count, 4e-9/TF)]
        ddiv=np.r_[0.,np.full(count, 1e-9/TF)]
        work=np.r_[0.,np.cumsum(np.diff(times)*power[1:])]
        dissipated=np.r_[0.,np.cumsum(np.diff(times)*(dgrad[1:]+ddiv[1:]))]
        balance=H-H[0]+dissipated-work
        scale=max(abs(H[-1]-H[0]),abs(work[-1]),abs(dissipated[-1]),1e-14)
        step={"full_l2":1e-13,"dynamic_l2":1e-13,"algebraic_l2":1e-17}
        metrics={
            "schema":"scrimp-v130-issue28-cn-refinement-v1",
            "status":"NATIVE_TSCN_FIXED_FINAL_TIME_DIAGNOSTIC",
            "nominal_mesh":mesh,"mpi_ranks":ranks,
            "target_rotation":.2,
            "initial_height":5.,"initial_momentum_identically_zero":True,
            "geometry_order":2,"height_order":2,
            "original_physical_facet_tangent_wall":True,
            "original_viscous_and_gyro_weak_bricks":True,
            "original_incompatible_rotating_initial_q0_untouched":True,
            "physical_scheme":"petsc_tscn_endpoint_theta_half",
            "nonlinear_mass_PETSc_endpoint_assumption_violated":True,
            "nonlinear_mass_time_accuracy_qualified":False,
            "time_order_qualified":False,
            "physical_energy_balance_convergence_qualified":False,
            "multistep_rotating_trajectory_qualified":False,
            "general_high_index_initializer_qualified":False,
            "physical_power_audit_used_h_weighted_boundary_terms":True,
            "initial_analytic_publication":{
                "strategy":"sandbox_exact_stationary_spinup_only",
                "status":"SANDBOX_ONLY_PUBLISHED_NOT_GENERAL_INITIALIZER",
                "differential_q0_preserved_bitwise":True,
                "getfem_globalvec_published_bitwise":True,
                "no_pseudo_or_BDF_initialization":True,
                "initial_native_IFunction":{"full_l2":1e-17}
            },
            "requested_steps":count,
            "actual_total_dofs":n,
            "differential_dofs":DYNAMIC_DOFS[mesh],
            "dynamic_variable_intervals":{
                "h":[0,nh],"p":[nh,2*nh]
            },
            "ts_accepted_steps":count,
            "ts_rejected_steps":0,
            "ts_snes_failures":0,
            "ts_reason":1,
            "fixed_step":TF/count,
            "final_time_fixed":TF,
            "actual_TS_time":TF,
            "all_full_native_IFunction_TSCN_step_residuals":[step]*count,
            "all_endpoint_original_algebraic_residuals":[1e-17]*count,
            "all_implicit_midpoint_OUT_OF_SCHEME_residuals":[.001]*count,
            "maximum_native_TSCN_equation_residual_l2":1e-13,
            "maximum_endpoint_native_algebraic_residual_l2":1e-17,
            "final_dynamic_state_change_l2":float(np.linalg.norm(
                states[-1,nh:nh+2*nh]-states[0,nh:nh+2*nh])),
            "minimum_height_coefficient_over_trajectory":5.,
            "initial_hamiltonian":float(H[0]),
            "final_hamiltonian":float(H[-1]),
            "final_hamiltonian_change":float(H[-1]-H[0]),
            "integrated_physical_wall_work":float(work[-1]),
            "integrated_viscous_dissipation":float(dissipated[-1]),
            "integrated_viscous_grad":float(np.sum(np.diff(times)*dgrad[1:])),
            "integrated_viscous_div":float(np.sum(np.diff(times)*ddiv[1:])),
            "integrated_normal_port_work":0.,
            "integrated_tangent_port_work":float(work[-1]),
            "final_physical_power_balance_defect":float(balance[-1]),
            "maximum_physical_power_balance_defect_l2":
                float(np.max(np.abs(balance))),
            "final_defect_over_actual_energy_exchange":
                float(balance[-1]/scale),
            "final_volume_balance_defect":0.,
            "maximum_volume_balance_defect":0.,
        }
        arrays={
            "times":times,"states":states,
            "H":H, "kinetic":np.zeros_like(H), "potential":H,
            "p_normal":np.zeros_like(H),
            "p_viscous":power,
            "d_grad":dgrad,"d_div":ddiv,
            "volume":np.full_like(H,5.),
            "volume_balance":np.zeros_like(H),
            "supplied":work, "dissipated":dissipated,
            "balance":balance,
            "relative_balance_exchange":balance/scale,
        }
        return metrics,arrays

    def store(self,mesh,count,ranks,doc,traces):
        root=self.root/(
            f"result_h{str(mesh).replace('.', 'p')}_n{count}_mpi{ranks}"
        )
        root.mkdir(exist_ok=True)
        (root/"cn_refinement.json").write_text(
            json.dumps(doc),encoding="utf-8")
        np.savez_compressed(root/"cn_refinement.npz",**traces)

    def mutate(self,mesh,count,ranks,func):
        report,traces=self.fixture(mesh,count,ranks)
        func(report,traces)
        self.store(mesh,count,ranks,report,traces)

    def test_probe_python_syntax(self):
        path=ROOT/"sandbox"/"probe_v130_issue28_cn_refinement.py"
        ast.parse(path.read_text(encoding="utf-8"))

    def test_synthetic_refinement_passes_without_claiming_order(self):
        result=compare(self.root)
        self.assertEqual(
            result["status"],
            "NATIVE_TSCN_FIXED_FINAL_TIME_REFINEMENT_DIAGNOSTIC_ONLY")
        self.assertFalse(result["nonlinear_mass_TSCN_general_accuracy_qualified"])
        for mesh in ("h0.2","h0.1"):
            orders=result["cases"][mesh][
                "observed_order_estimates_not_general_certificates"]
            self.assertTrue(all(abs(x-2.)<1e-8 for x in orders.values()))

    def test_rejects_bad_native_TSCN_one_intermediate_step(self):
        def bad(report,traces):
            report["all_full_native_IFunction_TSCN_step_residuals"][1]={
                "full_l2":1e-2,"dynamic_l2":1e-2,"algebraic_l2":1e-2
            }
        self.mutate(.1,4,1,bad)
        with self.assertRaisesRegex(ValueError,"step 2"):
            compare(self.root)

    def test_rejects_high_index_algebraic_drift(self):
        def bad(report,traces):
            report["all_endpoint_original_algebraic_residuals"][2]=1e-3
        self.mutate(.2,8,2,bad)
        with self.assertRaisesRegex(ValueError,"step 3"):
            compare(self.root)

    def test_rejects_wrong_physical_final_time(self):
        def bad(report,traces):
            report["actual_TS_time"]=.0019
        self.mutate(.2,2,1,bad)
        with self.assertRaisesRegex(ValueError,"fixed t_final"):
            compare(self.root)

    def test_rejects_MPI_full_FEM_vector_mismatch(self):
        def bad(report,traces):
            # Corrupt an ALGEBRAIC coefficient, keeping the original
            # h/p dynamic change norm and scalar reports identical.
            # This must reach the full MPI final FE vector gate.
            traces["states"][-1,DYNAMIC_DOFS[.1]+20] += 5e-4
        self.mutate(.1,8,2,bad)
        with self.assertRaisesRegex(ValueError,"full final FE coefficient"):
            compare(self.root)

    def test_rejects_false_native_physical_work(self):
        def bad(report,traces):
            report["integrated_physical_wall_work"]+=1e-3
        self.mutate(.2,1,1,bad)
        with self.assertRaisesRegex(ValueError,"bookkeeping incorrect"):
            compare(self.root)

    def test_rejects_corrupt_original_power_trace(self):
        def bad(report,traces):
            traces["p_viscous"][1] += 1e-3
        self.mutate(.1,2,1,bad)
        with self.assertRaisesRegex(ValueError,"differs from independent"):
            compare(self.root)

    def test_rejects_false_nonlinear_mass_accuracy_claim(self):
        def bad(report,traces):
            report["nonlinear_mass_time_accuracy_qualified"]=True
        self.mutate(.1,4,1,bad)
        with self.assertRaisesRegex(ValueError,"Unjustified"):
            compare(self.root)

    def test_rejects_missing_intermediate_native_MPI_result(self):
        path=self.root/"result_h0p1_n4_mpi2"/"cn_refinement.npz"
        path.unlink()
        with self.assertRaisesRegex(ValueError,"Missing native TSCN"):
            compare(self.root)


if __name__=="__main__":
    unittest.main()
