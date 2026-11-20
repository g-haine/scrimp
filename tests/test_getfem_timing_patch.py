# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Source contract for the approved native GetFEM MPI timing patch.

The former scripts/patch-getfem-mpi-timing.py helper was replaced
by the versioned 0000 patch. This test checks the actual shipped patch;
the separate GetFEM source audit checks clean application of all patches.
No PETSc/GetFEM installation is required here.
"""

from pathlib import Path
import unittest


PATCH = (
    Path(__file__).resolve().parents[1]
    / "scripts/patches/getfem-5.5/0000-mpi-assembly-timing.patch"
)
PRINT = 'cout << "Assembly time " << MPI_Wtime()-t_ref << endl;'


class GetFEMTimingPatchTest(unittest.TestCase):
    def test_approved_patch_gates_assembly_print_on_opt_in(self):
        diff = PATCH.read_text(encoding="utf-8")
        added = [
            line[1:] for line in diff.splitlines()
            if line.startswith("+") and not line.startswith("+++")
        ]
        removed = [
            line[1:] for line in diff.splitlines()
            if line.startswith("-") and not line.startswith("---")
        ]
        self.assertIn("#include <cstdlib>", added)
        self.assertIn(
            '      const char *flag = std::getenv("SCRIMP_GETFEM_ASSEMBLY_TIMING");',
            added,
        )
        self.assertIn(
            r"      if (flag && flag[0] == '1' && flag[1] == '\0')",
            added,
        )
        self.assertEqual(sum(PRINT in line for line in added), 1)
        self.assertEqual(sum(PRINT in line for line in removed), 1)
        self.assertIn("    if (MPI_IS_MASTER()) {", added)
        self.assertIn("    }", added)

    def test_patch_changes_only_getfem_model_source(self):
        diff = PATCH.read_text(encoding="utf-8")
        self.assertEqual(
            [line for line in diff.splitlines() if line.startswith("diff --git ")],
            ["diff --git a/src/getfem_models.cc b/src/getfem_models.cc"],
        )
        self.assertIn("@@ -2879,", diff)


if __name__ == "__main__":
    unittest.main()
