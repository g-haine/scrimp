# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Synthetic provenance and numerical-gate tests for native block audit."""
import importlib.util
import copy
import json
from pathlib import Path
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "_scrimp_compare_ph_audit",
    ROOT / "sandbox" / "compare_ph_block_audits.py",
)
COMPARE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(COMPARE)


class PhBlockComparisonTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.directory = Path(self.temp.name)
        self.records = {}
        roles = ["dynamic", "coenergy", "boundary", "algebraic"]
        fields = [
            {"name": "q", "group": "dynamic",
             "start": 0, "stop": 2, "ndofs": 2},
            {"name": "e_q", "group": "coenergy",
             "start": 2, "stop": 4, "ndofs": 2},
            {"name": "Y_B", "group": "boundary",
             "start": 4, "stop": 6, "ndofs": 2},
        ]

        def matrix(kind):
            block = [[0.] * 3 for _ in range(3)]
            role = [[0.] * 4 for _ in range(4)]
            if kind == "mass":
                block[0][0] = role[0][0] = 1.
            else:
                block[0][1] = role[0][1] = 2.
                block[1][0] = role[1][0] = 2.
                block[1][1] = role[1][1] = 3.
            block_nnz = [
                [1 if x else 0 for x in row] for row in block
            ]
            role_nnz = [
                [1 if x else 0 for x in row] for row in role
            ]
            return {
                "variable_order": [f["name"] for f in fields],
                "roles": roles,
                "block_effective_nnz": block_nnz,
                "role_effective_nnz": role_nnz,
                "block_frobenius_norm": block,
                "role_frobenius_norm": role,
            }
        matrices = {
            "flow_mass_F_and_algebraic_zero_rows": matrix("mass"),
            "linear_stiffness_all_bricks": matrix("stiffness"),
            "implicit_TS_Jacobian_last_shift": matrix("stiffness"),
        }
        for ranks in (1, 2):
            segments = [(0, 6)] if ranks == 1 else [(0, 3), (3, 6)]
            self.records[ranks] = {
                "schema": "scrimp-ph-block-audit-v1",
                "source_commit": "mock",
                "physical_model": "examples/wave.py:wave_eq",
                "mesh": 0.08,
                "mpi_ranks": ranks,
                "dofs": 6,
                "dt": 0.002,
                "steps": 2,
                "ts_type": "cn",
                "final_time": 0.004,
                "linear_mass": True,
                "linear_stiffness": True,
                "mumps_policy": "never",
                "fields": copy.deepcopy(fields),
                "owned_ranges_by_rank": [
                    {"rank": i, "lo": start, "hi": stop,
                     "role_owned_ranges": {}}
                    for i, (start, stop) in enumerate(segments)
                ],
                "assembled_matrices": copy.deepcopy(matrices),
            }
        self.flush()

    def flush(self):
        for p, data in self.records.items():
            (self.directory/f"blocks_{p}r.json").write_text(
                json.dumps(data), encoding="utf-8"
            )

    def test_exact_fields_and_block_norms_pass(self):
        result = COMPARE.compare(self.directory)
        self.assertEqual(result["ranks_compared"], [1, 2])
        self.assertTrue(
            (self.directory/"ph_block_parity.json").is_file()
        )

    def test_substituted_coenergy_model_without_e_unknowns_passes(self):
        # The coenergy formulation has independent q,p and port
        # variables but NO e_q/e_p GetFEM unknowns. A nonzero
        # physical feedback 0.5*Y is registered as a nonlinear
        # GetFEM tangent, so linear_stiffness=False is expected.
        for p in (1, 2):
            data = self.records[p]
            data["physical_model"] = (
                "examples/wave_coenergy.py:wave_coenergy_eq"
            )
            data["linear_stiffness"] = False
            data["formulation"] = {
                "formulation": "substituted_coenergy",
                "solver_candidate": "mass_dirac_port_block",
            }
            data["fields"] = [
                {"name": "q", "group": "dynamic",
                 "start": 0, "stop": 2, "ndofs": 2},
                {"name": "e_r", "group": "algebraic",
                 "start": 2, "stop": 4, "ndofs": 2},
                {"name": "U", "group": "boundary",
                 "start": 4, "stop": 6, "ndofs": 2},
            ]
            for matrix in data["assembled_matrices"].values():
                matrix["variable_order"] = ["q", "e_r", "U"]
            stiffness = data["assembled_matrices"][
                "linear_stiffness_all_bricks"
            ]
            stiffness["role_frobenius_norm"] = [
                [3., 0., 2., 1.],
                [0., 0., 0., 0.],
                [2., 0., 0., 0.],
                [1., 0., 0., 1.],
            ]
        self.flush()
        report = COMPARE.compare(self.directory)
        self.assertEqual(
            report["formulation"]["formulation"],
            "substituted_coenergy",
        )

    def test_substituted_coenergy_with_fictitious_e_block_fails(self):
        self.test_substituted_coenergy_model_without_e_unknowns_passes()
        for p in (1, 2):
            self.records[p]["fields"][0]["group"] = "coenergy"
        self.flush()
        with self.assertRaisesRegex(
            ValueError, "contains coenergy DoFs"
        ):
            COMPARE.compare(self.directory)

    def test_wrong_fem_interval_rejected(self):
        self.records[2]["fields"][1]["start"] = 3
        self.flush()
        with self.assertRaisesRegex(ValueError, "Mismatch in fields"):
            COMPARE.compare(self.directory)

    def test_missing_constitutive_coupling_rejected(self):
        matrix = self.records[1]["assembled_matrices"][
            "linear_stiffness_all_bricks"
        ]
        matrix["role_frobenius_norm"][1][1] = 0
        # Keep 1-/2-rank matrix-norm consistency so the test
        # specifically exercises the constitutive equation gate.
        self.records[2]["assembled_matrices"][
            "linear_stiffness_all_bricks"
        ]["role_frobenius_norm"][1][1] = 0
        self.flush()
        with self.assertRaisesRegex(ValueError, "Required pH coupling"):
            COMPARE.compare(self.directory)

    def test_mass_on_algebraic_row_rejected(self):
        for p in (1, 2):
            mat = self.records[p]["assembled_matrices"][
                "flow_mass_F_and_algebraic_zero_rows"
            ]
            mat["role_frobenius_norm"][2][2] = 0.5
        self.flush()
        with self.assertRaisesRegex(ValueError, "unexpected mass"):
            COMPARE.compare(self.directory)


if __name__ == "__main__":
    unittest.main()
