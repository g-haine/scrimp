# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Regression tests for the *diagnostic-only* native GetFEM MPI patch.

No GetFEM/PETSc installation is required to run these tests.
"""
import importlib.util
from pathlib import Path
import unittest


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "patch-getfem-mpi-timing.py"
spec = importlib.util.spec_from_file_location("getfem_mpi_timing_patch", SCRIPT)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)

UPSTREAM = """#include <iomanip>
namespace getfem {
void model::assembly(build_version version) {
#if GETFEM_PARA_LEVEL > 0
  if (MPI_IS_MASTER()) cout << "Assembly time " << MPI_Wtime()-t_ref << endl;
#endif
}
}
"""


class GetFEMTimingPatchTest(unittest.TestCase):
    def test_native_assembly_print_requires_explicit_opt_in(self):
        patched = module.patch_source(UPSTREAM)
        self.assertIn("#include <cstdlib>", patched)
        self.assertIn('std::getenv("SCRIMP_GETFEM_ASSEMBLY_TIMING")', patched)
        self.assertIn("flag[0] == '1' && flag[1] == '\\0'", patched)
        self.assertIn(
            'cout << "Assembly time " << MPI_Wtime()-t_ref << endl;', patched
        )
        self.assertEqual(module.patch_source(patched), patched)

    def test_refuses_unexpected_upstream_source(self):
        with self.assertRaisesRegex(ValueError, "Expected exactly one"):
            module.patch_source(UPSTREAM.replace("Assembly time ", "Time "))

    def test_refuses_multiple_timing_prints(self):
        modified = UPSTREAM.replace(
            "#endif",
            '  if (MPI_IS_MASTER()) cout << "Assembly time " '
            '<< MPI_Wtime()-t_ref << endl;\n#endif',
        )
        with self.assertRaisesRegex(ValueError, "found 2"):
            module.patch_source(modified)


if __name__ == "__main__":
    unittest.main()
