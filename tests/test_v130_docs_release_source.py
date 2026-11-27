"""Release-documentation parity checks; runnable without GetFEM/Sphinx.

These checks verify the *source contract*, NOT an HTML or PDF build.
The actual build is executed by scripts/run-v130-docs-build.sh and
the Docs GitHub Actions workflow; only their outputs can qualify it.
"""
from pathlib import Path
import unittest


ROOT=Path(__file__).resolve().parents[1]
DOC=ROOT/"docs/source"


class DocumentationReleaseSourceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.energy=(DOC/"energy_and_power.rst").read_text(encoding="utf-8")
        cls.index=(DOC/"index.rst").read_text(encoding="utf-8")

    def test_energy_guide_in_user_navigation(self):
        self.assertIn("    energy_and_power",self.index)
        self.assertTrue((DOC/"energy_and_power.rst").is_file())
        self.assertIn("BDF, BDF2, ARKIMEX",self.energy)
        self.assertIn("not a numerical proof of passivity",self.energy)

    def test_distinct_real_cn_and_approximate_nodal_quadrature(self):
        for phrase in (
            '"cn_step"', '"instantaneous"',
            "two consecutive end states saved",
            "P[0]", "not universally exact",
            "only a post-processing approximation",
            "Butcher weights", "stored", "time-dependent",
        ):
            with self.subTest(phrase=phrase):
                self.assertIn(phrase,self.energy)
        self.assertIn("power_sampling.py",self.energy)
        self.assertIn("power_integration.py",self.energy)

    def test_no_old_generic_output_paths_in_canonical_tutorials(self):
        pages=["heat","heat_wave","shallow_water",
               "mass_spring_damper","nonlinear_pendulum"]
        forbidden=("outputs/showcase/","outputs/pv/",
                   "outputs/png/Hamiltonian","outputs/mesh/")
        for name in pages:
            source=(DOC/"examples"/f"{name}.rst").read_text(encoding="utf-8")
            with self.subTest(page=name):
                for old in forbidden:
                    self.assertNotIn(old,source,
                        f"outdated canonical output path {old} in {name}")
                self.assertIn("outputs/",source)

    def test_heat_wave_cli_is_runnable_and_published_results_real(self):
        src=(DOC/"examples/heat_wave.rst").read_text(encoding="utf-8")
        for line in src.splitlines():
            if (line.strip().startswith("python -m examples.heat_wave")
                    or line.strip().startswith("--final-time 0.2")):
                self.assertFalse(
                    line.rstrip().endswith("\\\\"),
                    "double backslash does NOT continue a POSIX shell line"
                )
        self.assertIn("single-output-tree check",src)
        self.assertIn("full nine-case",src)
        self.assertIn("mesh/",src)
        self.assertIn(":doc:`../energy_and_power`",(DOC/"examples/shallow_water.rst").read_text())

    def test_gui_scope_is_not_an_unfinished_promise(self):
        gui=(DOC/"GUI.rst").read_text(encoding="utf-8")
        self.assertNotIn("#TO COMPLETE#",gui)
        self.assertIn("experimental",gui.lower())

    def test_sphinx_build_reproducibility_entry_points_exist(self):
        runner=(ROOT/"scripts/run-v130-docs-build.sh").read_text()
        self.assertIn("python -m sphinx",runner)
        self.assertIn("-b html",runner)
        self.assertIn("-b latex",runner)
        self.assertIn("latexmk",runner)
        workflow=(ROOT/".github/workflows/documentation.yml").read_text()
        self.assertIn("sphinx-build -b html",workflow)
        self.assertIn("sphinx-build -b latex",workflow)
        self.assertIn("latex-action",workflow)


if __name__=="__main__":
    unittest.main()
