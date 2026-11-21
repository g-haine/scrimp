# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Host-only negative controls for archived physical RHS rotating #28 pilot."""

import json
from pathlib import Path
import sys
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"sandbox"))

from compare_v130_issue28_physical_rhs_y0 import compare, EXPECTED_IMAGE


class OriginalPhysicalRHSAuditTest(unittest.TestCase):
    def setUp(self):
        tmp=tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.root=Path(tmp.name)
        (self.root/"manifest.json").write_text(
            json.dumps({
                "schema":"v130-issue28-physical-rhs-Y0-pilot-v2",
                "image_revision":EXPECTED_IMAGE,
                "physical_Y0_accepted":False,
                "no_physical_time_steps":True,
            }),encoding="utf-8",
        )
        self.cases={}
        for case in ("rest","rotating"):
            for ranks in (1,2):
                d={
                    "schema":"scrimp-v130-issue28-physical-RHS-Y0-pilot-v2",
                    "status":"DIAGNOSTIC_TRUE_RHS_SOLVED_NOT_INITIALIZED",
                    "case":case,"mpi_ranks":ranks,
                    "dofs":3300,"dynamic_dofs":1452,
                    "regular_algebraic_dofs":1716,"traction_Y_dofs":132,
                    "original_P2_geometry":True,
                    "original_q0_preserved_bitwise":True,
                    "bordered_KSP_reason":4,
                    "bordered_relative_true_rhs_residual":1e-12,
                    "matnest_original_row_and_column_IS_used":True,
                    "matnest_solution_reconstructed_by_block_IS":True,
                    "matnest_physical_rhs_roundtrip_l2":0.,
                    "candidate_qdot_coefficients":[0.0]*1452,
                    "candidate_regular_algebraic_derivative_coefficients":[0.0]*1716,
                    "candidate_traction_Y0_coefficients":
                        ([0.0]*132 if case=="rest" else [0.2]+[0.0]*131),
                    "candidate_qdot_l2":0.,
                    "candidate_regular_algebraic_derivative_l2":0.,
                    "original_native_IFunction_Fd_with_candidate_derivative_Y0_l2":1e-12,
                    "original_native_IFunction_Fa_with_candidate_Y0_l2":
                      (1e-17 if case=="rest" else 4.18e-8),
                    "candidate_traction_Y0_l2":
                      (0.0 if case=="rest" else 0.2),
                    "independent_rotating_algebraic_compatibility_qualified":False,
                    "physical_rotating_Y0_accepted":False,
                    "mesh_independent_physical_constraint_check_qualified":False,
                    "true_initialization_TS_executed":False,
                    "true_physical_TS_executed":False,
                    "initial_projection":{
                        "native_final_Fa_l2":1e-17 if case=="rest" else 4.18e-8,
                        "exact_discrete_rotating_compatibility_certified":False,
                    },
                }
                self.cases[case,ranks]=d
                self.store(case,ranks)

    def store(self,case,ranks):
        path=self.root/f"result_{case}_{ranks}"
        path.mkdir(exist_ok=True)
        (path/"physical_rhs_y0_pilot.json").write_text(
            json.dumps(self.cases[case,ranks]),encoding="utf-8",
        )

    def test_complete_but_not_qualified(self):
        result=compare(self.root)
        self.assertEqual(result["status"],"COMPLETE_DIAGNOSTIC_NOT_INITIALIZER")
        self.assertFalse(result["rotating_physical_initial_Y0_qualified"])

    def test_rejects_false_accepted_rotating_y0(self):
        self.cases["rotating",2]["physical_rotating_Y0_accepted"]=True
        self.store("rotating",2)
        with self.assertRaisesRegex(ValueError,"Unproven physical data"):
            compare(self.root)

    def test_rejects_missing_rank(self):
        (self.root/"result_rest_2"/"physical_rhs_y0_pilot.json").unlink()
        with self.assertRaisesRegex(ValueError,"Missing native report"):
            compare(self.root)

    def test_rejects_corrupted_parallel_traction_vector(self):
        self.cases["rotating",2]["candidate_traction_Y0_coefficients"] = (
            [0.0, 0.2] + [0.0]*130
        )
        self.store("rotating",2)
        with self.assertRaisesRegex(ValueError, "vector parity failure"):
            compare(self.root)

    def test_rejects_unmapped_matnest_rhs(self):
        self.cases["rotating",2]["matnest_physical_rhs_roundtrip_l2"] = 0.1
        self.store("rotating",2)
        with self.assertRaisesRegex(ValueError, "RHS packing failed"):
            compare(self.root)

    def test_rejects_wrong_rest_reference(self):
        self.cases["rest",1]["initial_projection"]["native_final_Fa_l2"]=1e-3
        self.store("rest",1)
        with self.assertRaisesRegex(ValueError,"Known exact stationary"):
            compare(self.root)


if __name__=="__main__":
    unittest.main()
