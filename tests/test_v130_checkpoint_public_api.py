# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See LICENSE for license information.

"""#29 host preflight for production checkpoint API (no PETSc/GMsh needed).

Check exact interoperability with the previously native-qualified
sandbox checkpoint bytes/schema and fail-closed model fingerprint.
Full numerical DPHS.restore must still be verified with native PETSc.
"""

import ast
import importlib.util
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"sandbox"))
from v130_checkpoint_state import (
    model_provenance as sandbox_provenance,
    read_checkpoint as sandbox_read_checkpoint,
    write_checkpoint as sandbox_write_checkpoint,
)

# The host-only fake is loaded by exact path to avoid importing any
# unrelated installed "tests" namespace or SCRIMP's native dependencies.
fixture_path = ROOT/"tests"/"test_v130_checkpoint_state.py"
fake_spec = importlib.util.spec_from_file_location(
    "_scrimp_checkpoint_native_free_fakes", fixture_path
)
fake_module = importlib.util.module_from_spec(fake_spec)
fake_spec.loader.exec_module(fake_module)
System = fake_module.System
PATH = ROOT/"scrimp"/"utils"/"dae_checkpoint.py"
spec = importlib.util.spec_from_file_location("_scrimp_checkpoint_host", PATH)
production = importlib.util.module_from_spec(spec)
spec.loader.exec_module(production)


class CheckpointPublicAPIHostTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.path = Path(self.tmp.name)/"accepted.npz"
        self.state = np.array([5.,5.,5.,.5,.5,.5])
        self.provenance = production.model_provenance(System())

    def test_production_helper_matches_native_qualified_model_signature(self):
        self.assertEqual(self.provenance, sandbox_provenance(System()))

    def test_production_checkpoint_readable_by_original_sandbox_format(self):
        production.write_checkpoint(
            self.path, z=self.state, t=.02,
            provenance=self.provenance,
            ts_type="cn", accepted=True,
        )
        t,z,head = sandbox_read_checkpoint(
            self.path, current_provenance=sandbox_provenance(System())
        )
        self.assertEqual(t, .02)
        np.testing.assert_array_equal(z,self.state)
        self.assertEqual(head["kind"],"cold_state_only_not_ts_history")

    def test_original_native_qualified_format_readable_by_production(self):
        sandbox_write_checkpoint(
            self.path, z=self.state, t=.02,
            provenance=sandbox_provenance(System()),
            ts_type="cn", accepted=True,
        )
        t,z,head = production.read_checkpoint(
            self.path, current_provenance=self.provenance
        )
        self.assertEqual(t,.02)
        np.testing.assert_array_equal(z,self.state)

    def test_production_rejects_modified_physical_boundary_law(self):
        production.write_checkpoint(
            self.path,z=self.state,t=.02,
            provenance=self.provenance,
            ts_type="cn",accepted=True
        )
        with self.assertRaisesRegex(ValueError,"mesh, FE numbering"):
            production.read_checkpoint(
                self.path,
                current_provenance=production.model_provenance(
                    System(control="different boundary")
                ),
            )

    def test_dphs_public_methods_and_real_cold_injection_exist(self):
        tree = ast.parse(
            (ROOT/"scrimp"/"dphs.py").read_text(encoding="utf-8")
        )
        classes = [x for x in tree.body if isinstance(x,ast.ClassDef)
                   and x.name=="DPHS"]
        self.assertEqual(len(classes),1)
        methodnames={x.name for x in classes[0].body
                     if isinstance(x,(ast.FunctionDef,ast.AsyncFunctionDef))}
        self.assertTrue({
            "save_checkpoint","load_checkpoint",
            "_publish_pending_checkpoint_after_assembly","_solve_impl"
        }<=methodnames)
        impl=next(x for x in classes[0].body
                  if isinstance(x,ast.FunctionDef)
                  and x.name=="_solve_impl")
        code=ast.unparse(impl)
        self.assertIn("self._publish_pending_checkpoint_after_assembly()",code)
        self.assertIn("self.init_step()",code)
        self.assertIn("elif self.time_scheme['init_step']",code)
        self.assertIn("self.gf_model.from_variables()",code)


if __name__ == "__main__":
    unittest.main()
