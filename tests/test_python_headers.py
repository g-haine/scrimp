# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Ensure every tracked Python source carries the standard SCRIMP header.

Run with: python -m unittest -v tests.test_python_headers
This test depends only on Python and Git (no GetFEM/PETSc needed).
"""

import subprocess
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
HEADER = "# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics\n#\n# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3\n# \n# See the LICENSE file for license information.\n#\n# github: https://github.com/g-haine/scrimp"


class PythonHeaderTest(unittest.TestCase):
    def test_all_tracked_python_files_have_header(self):
        try:
            paths = subprocess.check_output(
                ["git", "-C", str(ROOT), "ls-files", "-z", "--", "*.py"]
            ).decode("utf-8").split("\0")
        except (FileNotFoundError, subprocess.CalledProcessError) as exc:
            self.skipTest(f"Git worktree unavailable: {exc}")
        paths = [path for path in paths if path]
        self.assertTrue(paths, "No Python files found in repository")
        invalid = []
        for path in paths:
            contents = (ROOT / path).read_text(encoding="utf-8-sig")
            if not contents.startswith(HEADER + "\n"):
                invalid.append(path)
        self.assertFalse(invalid, "Invalid SCRIMP headers:\n" + "\n".join(invalid))

    def test_gui_generator_uses_current_header(self):
        template = (ROOT / "scrimp" / "utils" / "GUI.py").read_text(
            encoding="utf-8"
        )
        self.assertIn('heading = """' + HEADER + "\n", template)


if __name__ == "__main__":
    unittest.main()
