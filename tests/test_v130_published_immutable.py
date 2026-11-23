# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Mandatory negative tests for immutable archived publication sources."""

import importlib.util
from pathlib import Path
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
GUARD = ROOT / "scripts" / "check-v130-published-immutable.py"
spec = importlib.util.spec_from_file_location("_published_guard", GUARD)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


class ImmutablePublishedTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.path = "examples/published/ORIGINAL_2000/source.py"
        self.file = self.root / self.path
        self.file.parent.mkdir(parents=True)
        self.file.write_bytes(b"original publication snapshot\n")
        self.sha = module.git_blob_sha(self.file.read_bytes())
        self.expected = {self.path: self.sha}
        self.index = {self.path: ("100644", self.sha)}

    def check(self, index=None, untracked=(), expected=None):
        return module.validate(
            self.root,
            self.index if index is None else index,
            untracked,
            self.expected if expected is None else expected,
        )

    def test_original_bytes_and_index_pass(self):
        self.assertEqual(self.check(), [])

    def test_modified_unstaged_worktree_fails(self):
        self.file.write_bytes(b"changed physical model\n")
        self.assertTrue(any("WORKTREE BYTES" in x for x in self.check()))

    def test_staged_modified_index_fails(self):
        newindex={self.path:("100644",module.git_blob_sha(b"modified\n"))}
        self.assertTrue(any("GIT INDEX" in x for x in self.check(index=newindex)))

    def test_deletion_fails(self):
        self.assertTrue(any("DELETED" in x for x in self.check(index={})))

    def test_replacement_symlink_fails(self):
        target=self.root/"new_script"
        target.write_bytes(b"new\n")
        self.file.unlink()
        self.file.symlink_to(target)
        self.assertTrue(any("symlinked" in x for x in self.check()))

    def test_rename_or_added_snapshot_fails(self):
        extra="examples/published/ORIGINAL_2000/corrected.py"
        index={**self.index,extra:("100644",self.sha)}
        self.assertTrue(any("ADDED/RENAMED" in x for x in self.check(index=index)))

    def test_untracked_new_code_in_published_fails(self):
        self.assertTrue(any("UNTRACKED" in x for x in self.check(
            untracked=["examples/published/ORIGINAL_2000/new.py"]
        )))

    def test_file_permission_or_symlink_index_change_fails(self):
        index={self.path:("100755",self.sha)}
        self.assertTrue(any("MODE" in x for x in self.check(index=index)))

    def test_real_snapshot_baseline_is_pinned_and_no_scripts_altered(self):
        self.assertEqual(
            module.BASELINE_COMMIT,
            "b6fe4933415fb19eb02a8448a9fc587432196227",
        )
        self.assertEqual(len(module.EXPECTED_GIT_BLOBS),8)
        self.assertTrue(all(
            path.startswith("examples/published/")
            and len(sha)==40
            for path,sha in module.EXPECTED_GIT_BLOBS.items()
        ))


if __name__ == "__main__":
    unittest.main()
