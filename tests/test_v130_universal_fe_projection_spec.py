# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Source-only constraints for universal FE projection native pilot.

No GetFEM/PETSc dependency in the GitHub Actions source workflow.
These tests cannot substitute for the original native physical run.
"""

import ast
from pathlib import Path
import unittest


SOURCE = (
    Path(__file__).resolve().parents[1]
    / "sandbox/probe_v130_universal_fem_initialization.py"
)


class UniversalFESpecTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.source = SOURCE.read_text(encoding="utf-8")
        cls.tree = ast.parse(cls.source)
        decl = next(
            n for n in cls.tree.body if isinstance(n, ast.Assign)
            and any(isinstance(t, ast.Name) and t.id == "CASES"
                    for t in n.targets)
        )
        cls.cases = ast.literal_eval(decl.value)

    def test_every_polynomial_is_parseable_gwfl_without_python_power(self):
        for name, (_, _, gwfl, _, _, _) in self.cases.items():
            with self.subTest(name=name):
                self.assertNotIn("^", gwfl)
                self.assertNotIn("**", gwfl)

    def test_required_scalar_and_vector_and_derivative_families(self):
        for name in (
            "cg1_scalar", "cg2_scalar", "dg0_scalar", "dg2_scalar",
            "cg1_vector", "dg1_vector", "rt0_vector",
            "bdm1_vector", "nedelec_vector", "hermite2_scalar",
            "argyris2_scalar", "morley2_scalar", "hct2_scalar",
            "reduced_hct2_scalar", "nonconforming_p1_scalar",
        ):
            self.assertIn(name, self.cases)
        self.assertEqual(
            sum(not spec[-1] for spec in self.cases.values()), 9
        )

    def test_legacy_vector_mismatch_not_rejected_as_projection_failure(self):
        run = next(
            n for n in self.tree.body if isinstance(n, ast.FunctionDef)
            and n.name == "run_case"
        )
        run_source = ast.get_source_segment(self.source, run)
        self.assertIn(
            "old_production_silent_physical_failure", run_source
        )
        self.assertIn(
            '{"cg1_scalar","cg2_scalar","dg0_scalar","dg2_scalar"}',
            run_source,
        )
        self.assertNotIn("FAIL_LEGACY_LAGRANGIAN_COMPATIBILITY", run_source)

    def test_no_production_mutation(self):
        self.assertIn("production_changed", self.source)
        self.assertIn("canonical_derivative_moment_interpolation_qualified",
                      self.source)


if __name__ == "__main__":
    unittest.main()
