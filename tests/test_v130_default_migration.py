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
import json
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]
LEGACY_EXAMPLES = (
    "heat.py",
    "mass_spring_damper.py",
    "nonlinear_pendulum.py",
    "shallow_water.py",
)
# All four 2D wave/tutorial variants have since migrated to
# independently qualified strict algebraic initializations.
MODERN_ALGEBRAIC_EXAMPLES = (
    "heat_wave.py",
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

    def test_separately_qualified_wave_tutorials_use_algebraic(self):
        for name in MODERN_ALGEBRAIC_EXAMPLES:
            with self.subTest(name=name):
                source=(ROOT/"examples"/name).read_text(encoding="utf-8")
                self.assertIn('init_step_ts_type="algebraic"',source)
                self.assertNotIn('init_step_ts_type="pseudo"',source)

    def test_heat_wave_deliberate_CN_algebraic_physics_migration(self):
        source=(ROOT/"examples"/"heat_wave.py").read_text(
            encoding="utf-8"
        )
        self.assertIn('init_step_ts_type="algebraic"',source)
        self.assertIn('ts_type="cn"',source)
        self.assertIn('ts_adapt_type="none"',source)
        self.assertIn('S.FEM("Heat flux", k, "RT")',source)
        self.assertIn('S.FEM("q", k, "RT")',source)

    def test_modern_anisotropic_wave_initial_data_match_homogeneous_boundary(self):
        src=(ROOT/"examples/wave.py").read_text(encoding="utf-8")
        self.assertIn('system.set_initial_value("p", "0.")',src)
        self.assertIn('system.set_initial_value("q", "[0., 0.]")',src)
        self.assertIn('init_step_ts_type="algebraic"',src)

    def test_candidate_default_is_algebraic_and_legacy_is_explicit(self):
        src=(ROOT / "scrimp/dphs.py").read_text(
            encoding="utf-8"
        )
        # Source candidate is promoted only to prepare independent
        # *native* no-option verification. A green source check alone
        # cannot certify the final default in a release.
        self.assertGreaterEqual(
            src.count('"init_step_ts_type", "algebraic"'), 2
        )
        # PETSc.Options materializes the default. A fallback string
        # alone is insufficient: the old pseudo assignment silently
        # defeated the intended default in native runtime.
        self.assertIn(
            'self.time_scheme["init_step_ts_type"] = "algebraic"',
            src,
        )
        self.assertNotIn(
            'self.time_scheme["init_step_ts_type"] = "pseudo"',
            src,
        )
        self.assertIn(
            'explicitly_selected_initializer = "init_step_ts_type" in kwargs',
            src,
        )
        self.assertIn(
            'if init_type == "algebraic":', src
        )
        self.assertIn(
            'if init_type not in ("pseudo", "bdf"):', src
        )
        # New nonlinear-state/affine-algebraic initializer is
        # explicitly opt-in; the existing *default* remains algebraic,
        # while pseudo/BDF remain legal and unchanged.
        self.assertIn(
            'if init_strategy not in ("algebraic", "nonlinear_affine", "pseudo", "bdf"):',
            src,
        )
        self.assertIn(
            '"init_step_ts_type must be algebraic, nonlinear_affine, "',
            src,
        )
        self.assertIn('if init_type == "nonlinear_affine":', src)
        self.assertIn('if init_type not in ("pseudo", "bdf"):', src)
        # init_step=False must remain legal if no strategy was explicitly
        # selected. The explicit but contradictory opt-in is rejected.
        self.assertIn(
            'and explicitly_selected_initializer:', src
        )
        self.assertIn(
            'when no initialization is intended', src
        )

    def test_published_maxwell_keeps_legacy_pseudo(self):
        src=(ROOT / "examples/published/COMPEL_2025/Maxwell_2D.py"
             ).read_text(encoding="utf-8")
        self.assertIn('init_step_ts_type="pseudo"', src)

    def test_legacy_notebooks_keep_pseudo_explicitly(self):
        for name, expected in (
            ("Heat_2D.ipynb", 1),
            ("Heat_Wave_2D.ipynb", 1),
            ("Shallow_water_2D.ipynb", 1),
            ("Wave_2D.ipynb", 2),
            ("Wave_2D_CoEnergy.ipynb", 1),
        ):
            with self.subTest(notebook=name):
                j=json.loads((
                    ROOT / "examples/notebooks" / name
                ).read_text(encoding="utf-8"))
                cells=[
                    "".join(cell.get("source", []))
                    for cell in j.get("cells", [])
                    if cell.get("cell_type") == "code"
                ]
                source="\n".join(cells)
                self.assertEqual(
                    source.count("set_time_scheme("), expected
                )
                self.assertEqual(
                    source.count('init_step_ts_type="pseudo"'), expected
                )

    def test_legacy_benchmark_and_lumped_sandbox_explicit(self):
        for path in (
            "sandbox/bench_acoustic_2d.py",
            "sandbox/sloshing_tank.py",
        ):
            with self.subTest(path=path):
                src=(ROOT / path).read_text(encoding="utf-8")
                self.assertIn('init_step_ts_type="pseudo"', src)

    def test_gui_script_generator_does_not_duplicate_strategy(self):
        src=(ROOT / "scrimp/utils/GUI.py").read_text(
            encoding="utf-8"
        )
        self.assertIn("has_initializer_choice",src)
        self.assertIn('file.write(\'init_step_ts_type="pseudo",\')',src)
        self.assertIn('"init_step_ts_type")',src)

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
