# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Guard historical published examples ahead of algebraic default migration.

Qualifying default behavior still requires native GetFEM/PETSc
after the explicit legacy-mode source migration. This module
only checks syntax and intent, not numerical behavior.
"""

import ast
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]
LEGACY_EXAMPLES = (
    "heat.py",
    "heat_wave.py",
    "mass_spring_damper.py",
    "nonlinear_pendulum.py",
    "shallow_water.py",
    "wave.py",
    "wave_coenergy.py",
    "wave_dissipative.py",
)


class TestDefaultMigration(unittest.TestCase):
    def test_all_original_examples_preserve_explicit_pseudo(self):
        for name in LEGACY_EXAMPLES:
            with self.subTest(example=name):
                src=(ROOT / "examples" / name).read_text(
                    encoding="utf-8"
                )
                tree=ast.parse(src)
                found=[]
                for node in ast.walk(tree):
                    if not (
                        isinstance(node,ast.Call)
                        and isinstance(node.func,ast.Attribute)
                        and node.func.attr=="set_time_scheme"
                    ):
                        continue
                    vals={
                        kw.arg:kw.value for kw in node.keywords
                        if kw.arg
                    }
                    if "init_step_ts_type" in vals:
                        arg=vals["init_step_ts_type"]
                        if isinstance(arg,ast.Constant):
                            found.append(arg.value)
                        else:
                            found.append("<dynamic>")
                self.assertTrue(found,"No set_time_scheme call")
                self.assertTrue(
                    all(x=="pseudo" for x in found),
                    f"Historical example changed initializer: {found}",
                )

    def test_original_wave_gaussian_kept_unchanged(self):
        src=(ROOT / "examples/wave.py").read_text(
            encoding="utf-8"
        )
        self.assertIn(
            '3**(-20*((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5)))',
            src,
        )
        self.assertIn('wave.set_initial_value("p", p_0)', src)

    def test_default_promotion_is_not_premature(self):
        src=(ROOT / "scrimp/dphs.py").read_text(
            encoding="utf-8"
        )
        # Gate intentionally pinned to current opt-in-only source;
        # when second native model and full migration pass, this
        # test MUST be updated in the same commit that changes the
        # global default and adds the new no-option native test.
        self.assertIn(
            '"init_step_ts_type", "pseudo"', src
        )
        self.assertIn(
            'if init_type == "algebraic":', src
        )

    def test_sphinx_documents_legacy_modes(self):
        src=(ROOT / "docs/source/initialization.rst").read_text(
            encoding="utf-8"
        )
        self.assertIn('init_step_ts_type="pseudo"', src)
        self.assertIn('init_step_ts_type="algebraic"', src)
        self.assertIn('init_step_bdf_order=2', src)
        self.assertIn('init_step_from_negative_time=True', src)


if __name__=="__main__":
    unittest.main()
