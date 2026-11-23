# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Host-only negative controls for archived #28 weak wall measurements."""

import ast
import json
import math
from pathlib import Path
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"sandbox"))
from compare_v130_issue28_weak_wall import compare, EXPECTED_IMAGE


class WeakWallNativeComparatorTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        (self.root/"manifest.json").write_text(
            json.dumps({
                "schema": "v130-issue28-weak-wall-moments-v1",
                "image_revision": EXPECTED_IMAGE,
                "original_fixed_h0_p0": True,
                "physical_initial_Y0_accepted": False,
                "no_physical_time_steps": True,
            }), encoding="utf-8"
        )
        for mesh in (.2, .1):
            for case in ("rest", "rotating"):
                for ranks in (1, 2):
                    self.store(mesh, case, ranks, self.fixture(mesh, case, ranks))

    @staticmethod
    def fixture(mesh, case, ranks):
        perimeter = 6.28
        omega = 0. if case == "rest" else .2
        weak = 0. if case == "rest" else 2e-6*(mesh/.2)**2
        strong = 0. if case == "rest" else 3e-5*(mesh/.2)**2
        vectors = {
            "full": [weak, weak*.5],
            "normal": [weak, 0.],
            "tangential": [0., weak*.5],
        }
        moments = {}
        dual = {}
        for name, vec in vectors.items():
            norm = math.hypot(*vec)
            moments[name] = {
                "gf_output_layout": "full-model-DOF-vector",
                "mpi_vector_reduction_mode":
                    "serial" if ranks == 1 else "replicated-global",
                "per_rank_native_vector": [vec[:] for _ in range(ranks)],
                "assembled_vector": vec[:],
                "coefficient_l2": norm,
            }
            rms = 0. if case == "rest" else (
                weak if name == "full"
                else weak*.8 if name == "normal" else weak*.2
            )
            dual[name] = {
                "riesz_L2_projected_mismatch_l2": rms*math.sqrt(perimeter),
                "riesz_boundary_RMS": rms,
                "boundary_mass_solver_relative_residual": 1e-15,
                "boundary_mass_is_negative_native_J_UU": True,
            }
        return {
            "schema": "scrimp-v130-issue28-weak-wall-moments-v1",
            "status": "NATIVE_WEAK_MOMENT_DIAGNOSTIC_ONLY",
            "mpi_ranks": ranks, "nominal_mesh": mesh,
            "geometry_order": 2,
            "rotation": omega, "dofs": 100,
            "boundary_U_dofs": 2,
            "raw_per_rank_perimeter": [perimeter]*ranks,
            "boundary_perimeter": perimeter,
            "perimeter_reduction_mode":
                "serial" if ranks == 1 else "replicated-global",
            "strong_wall_vector_mismatch_rms": strong,
            "raw_per_rank_strong_wall_squared_integral":
                [strong*strong*perimeter]*ranks,
            "strong_squared_integral_reduction_mode":
                "serial" if ranks == 1 else "replicated-global",
            "moment_vectors": moments,
            "normal_tangential_decomposition_l2": 0.,
            "boundary_mass_riesz": dual,
            "original_h0_p0_Y_preserved_bitwise": True,
            "entire_state_rollback_bitwise": True,
            "no_initialized_traction_Y0": True,
            "no_initialization_or_physical_TS": True,
            "algebraic_feasibility_proven": False,
        }

    def store(self, mesh, case, ranks, data):
        folder = self.root/(
            f"result_h{str(mesh).replace('.', 'p')}_{case}_{ranks}"
        )
        folder.mkdir(exist_ok=True)
        (folder/"weak_wall.json").write_text(
            json.dumps(data), encoding="utf-8"
        )

    def test_source_python_syntax(self):
        probe = ROOT/"sandbox"/"probe_v130_issue28_weak_wall.py"
        ast.parse(probe.read_text(encoding="utf-8"))

    def test_passes_but_does_not_certify_initializer(self):
        result = compare(self.root)
        self.assertEqual(
            result["status"], "COMPLETE_WEAK_WALL_DIAGNOSTIC_NOT_INITIALIZER"
        )
        self.assertFalse(result["rotating_exact_weak_DAE_compatibility_certified"])

    def test_rejects_false_accepted_Y0(self):
        r = self.fixture(.2, "rotating", 2)
        r["no_initialized_traction_Y0"] = False
        self.store(.2, "rotating", 2, r)
        with self.assertRaisesRegex(ValueError, "Unqualified physical"):
            compare(self.root)

    def test_rejects_mismatched_MPI_wall_vector(self):
        r = self.fixture(.1, "rotating", 2)
        r["moment_vectors"]["full"]["per_rank_native_vector"][1] = [0.3, 0.2]
        self.store(.1, "rotating", 2, r)
        with self.assertRaisesRegex(ValueError, "Replicated GetFEM MPI weak vector mismatch"):
            compare(self.root)

    def test_rejects_false_Riesz_norm(self):
        r = self.fixture(.2, "rotating", 2)
        r["boundary_mass_riesz"]["full"]["riesz_boundary_RMS"] = 1.
        self.store(.2, "rotating", 2, r)
        with self.assertRaisesRegex(ValueError, "Riesz mass physical norm"):
            compare(self.root)

    def test_rejects_missing_rank(self):
        p = (self.root/"result_h0p1_rest_2"/"weak_wall.json")
        p.unlink()
        with self.assertRaisesRegex(ValueError, "Missing native report"):
            compare(self.root)


if __name__ == "__main__":
    unittest.main()
