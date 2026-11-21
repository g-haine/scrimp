# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Host negative controls for #28 native FE physical compatibility evidence."""

import ast
import json
from pathlib import Path
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"sandbox"))
from compare_v130_issue28_physical_compatibility import (
    compare, EXPECTED_IMAGE, checked_integrals, mode_against_serial,
)


class OriginalPhysicalCompatibilityTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        (self.root/"manifest.json").write_text(json.dumps({
            "schema": "v130-issue28-physical-compatibility-v1",
            "image_revision": EXPECTED_IMAGE,
            "original_geometry_P2": True, "original_fixed_h0_p0": True,
            "physical_initial_data_accepted": False,
            "no_physical_time_steps": True,
        }), encoding="utf-8")
        for h in (.2, .1):
            for case in ("rest", "rotating"):
                for ranks in (1, 2):
                    self.write(h, case, ranks, self.example(h, case, ranks))

    @staticmethod
    def integral(value, ranks, replicated):
        v = ([value] if ranks == 1 else
             [value, value] if replicated else [0.6*value, 0.4*value])
        return {"per_rank": v, "sum": sum(v)}

    @classmethod
    def sample(cls, mesh, case, ranks):
        replicated = (ranks == 2 and mesh == .2)
        factor = ranks if replicated else 1
        n = 0.0 if case == "rest" else mesh/10
        t = 0.0 if case == "rest" else mesh/8
        area, length = 3.13, 6.27
        phys = {
            "area": cls.integral(area, ranks, replicated),
            "perimeter": cls.integral(length, ranks, replicated),
            "initial_p_normal_squared":
                cls.integral(length*n*n, ranks, replicated),
            "initial_p_tangent_mismatch_squared":
                cls.integral(length*t*t, ranks, replicated),
            "initial_p_wall_vector_mismatch_squared":
                cls.integral(length*(n*n+t*t), ranks, replicated),
            "initial_p_tangent_squared":
                cls.integral(length*0.04, ranks, replicated),
        }
        return {
            "raw_integrals": phys,
            "rms": {
                key.replace("_squared", "_rms"):
                    (max(0., val["sum"] / (area*factor if key.startswith("interior_")
                                          else length*factor)))**0.5
                for key, val in phys.items()
                if key not in ("area", "perimeter")
            },
        }

    @classmethod
    def example(cls, h, case, ranks):
        baseline = cls.sample(h, case, ranks)
        fa = 1e-17 if case == "rest" else 4.18e-8
        if case == "rest":
            proj = {
                "type": "exact-analytical-rest", "native_Fa_l2": fa,
                "physical_integrals": cls.sample(h, case, ranks),
            }
        elif h == .2:
            proj = {
                "type": "distributed-rectangular-LSQR-two-tolerances",
                "trials": [
                    {"kind": kind, "native_Fa_l2": fa,
                     "normal_equation_AT_Fa_l2": 1e-10,
                     "algebraic_correction_l2": 13.,
                     "relative_coefficient_change_from_standard":
                         (0. if kind == "standard" else 1e-7),
                     "iterations": 5100,
                     "physical_integrals": cls.sample(h, case, ranks)}
                    for kind in ("standard", "tight")
                ],
            }
        else:
            proj = {"type": "not-run-refinement-boundary-only"}
        return {
            "schema": "scrimp-v130-issue28-physical-compatibility-v1",
            "status": "DIAGNOSTIC_ONLY_NOT_INITIALIZED",
            "mpi_ranks": ranks, "mesh_nominal_h": h,
            "geometry_order": 2, "height_order": 2,
            "rotation": 0. if case == "rest" else .2,
            "dofs": 100, "boundary_at_prescribed_q0": baseline,
            "projection": proj, "original_q0_preserved_bitwise": True,
            "state_rollback_bitwise": True,
            "accepted_rotating_Y0": False,
            "exact_DAE_compatibility_qualified": False,
            "initialization_TS_executed": False,
            "physical_TS_executed": False,
        }

    def write(self, h, case, ranks, data):
        p = self.root/f"result_h{str(h).replace('.', 'p')}_{case}_{ranks}"
        p.mkdir(exist_ok=True)
        (p/"physical_compatibility.json").write_text(
            json.dumps(data), encoding="utf-8"
        )

    def test_new_native_probe_is_valid_python_source(self):
        # Python syntax check without requiring PETSc/GetFEM on the host.
        source = ROOT/"sandbox"/"probe_v130_issue28_physical_compatibility.py"
        ast.parse(source.read_text(encoding="utf-8"))

    def test_complete_diagnostic_not_initializer(self):
        r = compare(self.root)
        self.assertEqual(
            r["status"], "COMPLETE_PHYSICAL_TRACE_DIAGNOSTIC_NOT_INITIALIZER"
        )
        self.assertFalse(r["exact_rotating_DAE_initial_compatibility_qualified"])
        self.assertAlmostEqual(r["normal_RMS_fine_over_coarse"], .5)

    def test_rejects_false_compatibility_acceptance(self):
        r = self.example(.2, "rotating", 2)
        r["accepted_rotating_Y0"] = True
        self.write(.2, "rotating", 2, r)
        with self.assertRaisesRegex(ValueError, "falsely claims"):
            compare(self.root)

    def test_rejects_mpi_physical_norm_mismatch(self):
        r = self.example(.1, "rotating", 2)
        v = r["boundary_at_prescribed_q0"]["raw_integrals"]
        v["initial_p_normal_squared"]["per_rank"][0] *= 2
        v["initial_p_normal_squared"]["sum"] = sum(
            v["initial_p_normal_squared"]["per_rank"]
        )
        r["boundary_at_prescribed_q0"]["rms"]["initial_p_normal_rms"] = (
            v["initial_p_normal_squared"]["sum"] / v["perimeter"]["sum"]
        )**0.5
        self.write(.1, "rotating", 2, r)
        with self.assertRaisesRegex(ValueError, "decomposition|parity|mismatch"):
            compare(self.root)

    def test_rejects_missing_native_rank(self):
        (self.root/"result_h0p2_rest_2"/"physical_compatibility.json").unlink()
        with self.assertRaisesRegex(ValueError, "Missing native evidence"):
            compare(self.root)

    def test_rejects_nonphysical_negative_integral(self):
        p = self.sample(.2, "rotating", 1)
        p["raw_integrals"]["initial_p_normal_squared"]["per_rank"] = [-1.]
        p["raw_integrals"]["initial_p_normal_squared"]["sum"] = -1.
        with self.assertRaisesRegex(ValueError, "Negative squared"):
            checked_integrals(p, 1)

    def test_mpi_scalar_assembly_mode_calibrated(self):
        serial = {"area": 3.13, "perimeter": 6.27}
        self.assertEqual(
            mode_against_serial(serial, {"area": 6.26, "perimeter": 12.54}),
            "replicated-global"
        )
        self.assertEqual(
            mode_against_serial(serial, dict(serial)), "summed-local"
        )


if __name__ == "__main__":
    unittest.main()
