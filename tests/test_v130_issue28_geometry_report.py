# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Source-only #28 report contracts: no GetFEM, PETSc or Docker required."""

import ast
from pathlib import Path
import unittest


SOURCE = Path(__file__).resolve().parents[1] / (
    "sandbox/probe_v130_issue28_feasibility.py"
)


class Issue28GeometryReportTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        module = ast.parse(SOURCE.read_text(encoding="utf-8"))
        cls.audit = next(
            f for f in module.body
            if isinstance(f, ast.FunctionDef) and f.name == "_audit"
        )
        cls.main = next(
            f for f in module.body
            if isinstance(f, ast.FunctionDef) and f.name == "main"
        )

    def test_audit_does_not_reference_cli_args_out_of_scope(self):
        # This catches the original runtime NameError, which py_compile
        # cannot detect, without attempting a native PDE simulation.
        used = {
            x.id for x in ast.walk(self.audit)
            if isinstance(x, ast.Name) and isinstance(x.ctx, ast.Load)
        }
        self.assertNotIn("args", used)

    def test_geometry_report_reads_system_configuration(self):
        report = next(
            x for x in ast.walk(self.audit)
            if isinstance(x, ast.Dict) and any(
                isinstance(k, ast.Constant)
                and k.value == "original_polygon_facet_boundary"
                for k in x.keys
            )
        )
        fields = {
            k.value: v for k, v in zip(report.keys, report.values)
            if isinstance(k, ast.Constant)
        }
        for field in (
            "original_polygon_facet_boundary",
            "geometry_order_requested",
            "getfem_import_has_nonlinear_geometric_transform",
        ):
            self.assertIn(field, fields)
        fields_source = ast.unparse(fields["original_polygon_facet_boundary"])
        self.assertIn("system._issue28_geometry_order", fields_source)

    def test_geometry_is_set_in_main_before_audit(self):
        assignments = [
            x for x in ast.walk(self.main) if isinstance(x, ast.Assign)
            and any(
                isinstance(t, ast.Attribute)
                and t.attr == "_issue28_geometry_order"
                for t in x.targets
            )
        ]
        self.assertEqual(len(assignments), 1)


if __name__ == "__main__":
    unittest.main()
