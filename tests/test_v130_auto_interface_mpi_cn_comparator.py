# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Host-side regression for the real MPI CN comparison's inline Python.

Exercises the *actual* inline comparator extracted from its shell
driver, rather than a rewritten imitation, with fake completed
per-rank reports. Does not claim to run or replace PETSc/GetFEM.
This test would have caught the historical uninitialized host_error
in the SUCCESS path, and checks archiving even when a comparison
fails or a native MPI rank exited nonzero.
"""

import contextlib
import io
import json
import os
from pathlib import Path
from types import SimpleNamespace
import sys
import tempfile
import traceback
import unittest
from unittest.mock import patch

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
DRIVER = ROOT / "scripts" / "run-v130-auto-interface-mpi-cn.sh"
START = 'python3 - "$out" <<\'PY\'\n'
STOP = '\nPY\n'


def inline_comparator():
    source = DRIVER.read_text(encoding="utf-8")
    if source.count(START) != 1:
        raise AssertionError("Expected one native-driver inline Python comparator")
    tail = source.split(START, 1)[1]
    if STOP not in tail:
        raise AssertionError("Comparator heredoc terminator is missing")
    return tail.split(STOP, 1)[0] + "\n"


def prepare_fake_evidence(root):
    """Generate minimal reports and deterministic *simulated* trace files."""
    for heat in (1, 2):
        for size in (1, 2):
            label = f"heat{heat}_mpi{size}"
            folder = root / label
            folder.mkdir()
            (root / f"{label}.exitcode").write_text("0\n")
            for rank in range(size):
                record = {
                    "MPI_ranks": size,
                    "heat_region": heat,
                    "rank": rank,
                    "status": "PASS_RANK_MPI_CN_AND_PHYSICAL_GATES",
                    "rank_consensus": True,
                    "real_PETSc_CN_steps": 4,
                    "physical_gates": {"fake_physics": True},
                    "plot_and_ParaView_export_skipped_on_all_ranks": True,
                    "physical_GetFEM_assembly_executed_on_all_ranks": True,
                    "rank_fingerprints": [{"sha": "same"}] * size,
                    "automatic_interface": {"logical_region": 10},
                }
                (folder / f"rank{rank}.json").write_text(
                    json.dumps(record, indent=2) + "\n"
                )
                np.savez_compressed(
                    folder / f"rank{rank}.npz",
                    times=np.arange(5, dtype=float) / 1000,
                    states=np.arange(20, dtype=float).reshape(5, 4),
                    H=np.array([5., 4.9, 4.8, 4.7, 4.6]),
                    phys_CN_P_interface_heat=np.ones(5),
                    phys_CN_P_interface_wave=np.ones(5),
                )


def compare(root):
    """Execute the EXACT production inline comparator in this process.

    Starting another Python interpreter is unnecessary to test an inline
    script and has been observed to terminate with SIGKILL (-9) in the
    native conda environment. Running in-process avoids a second Python
    process (and an independent NumPy/BLAS import) while still exercising
    the actual comparator source, with its real sys.argv and exit code.
    All global state and captured streams are restored on exit.
    """
    stdout, stderr = io.StringIO(), io.StringIO()
    exit_code = 0
    namespace = {"__name__": "__main__", "__file__": str(DRIVER)}
    with patch.dict(os.environ, {"SCRIMP_V130_MPI_PV_DIAGNOSTIC": "0"}):
      with patch.object(sys, "argv", [str(DRIVER), str(root)]):
        with contextlib.redirect_stdout(stdout), contextlib.redirect_stderr(stderr):
            try:
                exec(compile(inline_comparator(), str(DRIVER), "exec"), namespace)
            except SystemExit as exc:
                exit_code = int(exc.code) if isinstance(exc.code, int) else 1
            except Exception:
                traceback.print_exc()
                exit_code = 1
    return SimpleNamespace(
        returncode=exit_code, stdout=stdout.getvalue(),
        stderr=stderr.getvalue(),
    )


class HostComparatorRegression(unittest.TestCase):
    def test_actual_embedded_comparator_is_valid_python(self):
        compile(inline_comparator(), str(DRIVER), "exec")

    def test_success_writes_summary_and_archive_without_nameerror(self):
        with tempfile.TemporaryDirectory(prefix="v130-mpi-cn-success-") as d:
            root = Path(d) / "run"
            root.mkdir()
            prepare_fake_evidence(root)
            outcome = compare(root)
            self.assertEqual(outcome.returncode, 0, outcome.stderr + outcome.stdout)
            summary = json.loads((root / "comparison.json").read_text())
            self.assertEqual(
                summary["status"],
                "PASS_BOTH_AUTO_INTERFACE_ORIENTATIONS_MPI2_DAE_CN_ENERGY_PARITY",
            )
            self.assertIsNone(summary["host_comparison_failure"])
            self.assertEqual(set(summary["cases"]), {
                "heat1_wave2", "heat2_wave1"
            })
            self.assertTrue(Path(str(root) + ".zip").is_file())

    def test_real_trace_mismatch_fails_but_preserves_archive(self):
        with tempfile.TemporaryDirectory(prefix="v130-mpi-cn-diff-") as d:
            root = Path(d) / "run"
            root.mkdir()
            prepare_fake_evidence(root)
            path = root / "heat1_mpi2" / "rank0.npz"
            with np.load(path) as old:
                values = {key: np.asarray(old[key]).copy() for key in old.files}
            values["states"][1, 2] += .5
            np.savez_compressed(path, **values)
            outcome = compare(root)
            self.assertNotEqual(outcome.returncode, 0)
            summary = json.loads((root / "comparison.json").read_text())
            self.assertEqual(
                summary["status"],
                "FAIL_OR_UNQUALIFIED_AUTO_INTERFACE_MPI2_DAE_CN_PHYSICS",
            )
            self.assertEqual(
                summary["cases"]["heat1_wave2"]["status"],
                "FAIL_OR_UNQUALIFIED_MPI2_DAE_CN_PARITY",
            )
            self.assertIsNone(summary["host_comparison_failure"])
            self.assertTrue(Path(str(root) + ".zip").exists())

    def test_nonzero_native_rank_exit_fails_with_archive(self):
        with tempfile.TemporaryDirectory(prefix="v130-mpi-cn-native-fail-") as d:
            root = Path(d) / "run"
            root.mkdir()
            prepare_fake_evidence(root)
            (root / "heat2_mpi2.exitcode").write_text("1\n")
            outcome = compare(root)
            self.assertNotEqual(outcome.returncode, 0)
            summary = json.loads((root / "comparison.json").read_text())
            self.assertEqual(
                summary["cases"]["heat2_wave1"]["MPI2"]["exit_code"], 1,
            )
            self.assertIsNone(summary["host_comparison_failure"])
            self.assertTrue(Path(str(root) + ".zip").is_file())

    def test_malformed_rank_report_archived_as_host_error(self):
        with tempfile.TemporaryDirectory(prefix="v130-mpi-cn-bad-meta-") as d:
            root = Path(d) / "run"
            root.mkdir()
            prepare_fake_evidence(root)
            (root / "heat1_mpi2" / "rank1.json").write_text("not JSON\n")
            outcome = compare(root)
            self.assertNotEqual(outcome.returncode, 0)
            summary = json.loads((root / "comparison.json").read_text())
            self.assertIn("JSONDecodeError", summary["host_comparison_failure"])
            self.assertTrue(Path(str(root) + ".zip").is_file())


if __name__ == "__main__":
    unittest.main()
