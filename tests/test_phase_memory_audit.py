# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Synthetic tests for phase-resolved RSS audit (no PETSc/GetFEM)."""
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, ROOT / path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


SAMPLER = load("_scrimp_phase_memory", "scrimp/utils/phase_memory_audit.py")
ANALYZER = load("_scrimp_phase_report", "sandbox/analyze_phase_memory.py")


class PhaseMemoryAuditTests(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.root = Path(tmp.name)

    def test_reads_actual_process_rss_without_mpi(self):
        audit = SAMPLER.PhaseMemoryAudit(
            self.root, rank=0, ranks=1,
            source_commit="local-test", mode="compact",
        )
        audit.sample("before")
        audit.sample("after")
        file = audit.write()
        saved = json.loads(file.read_text(encoding="utf-8"))
        self.assertEqual(saved["schema"], "scrimp-phase-memory-v1")
        self.assertEqual(
            [r["phase"] for r in saved["samples"]],
            ["before", "after"],
        )
        self.assertGreater(saved["samples"][0]["VmRSS_kib"], 0)
        self.assertGreaterEqual(
            saved["samples"][-1]["VmHWM_kib"],
            saved["samples"][-1]["VmRSS_kib"],
        )

    def _build_fake_reports(self):
        for case in ("acoustic", "swe"):
            for mode in ("legacy", "compact"):
                for rank in (0, 1):
                    directory = self.root / case / mode
                    directory.mkdir(parents=True, exist_ok=True)
                    entries = []
                    for idx, name in enumerate(ANALYZER.PHASES):
                        rss = 100000 + rank * 1000 + idx * 100
                        if mode == "compact":
                            rss -= 1500
                        entries.append({
                            "phase": name, "VmRSS_kib": rss,
                            "VmHWM_kib": rss + 100,
                        })
                    report = {
                        "schema": "scrimp-phase-memory-v1",
                        "source_commit": "same-sha",
                        "rank": rank, "ranks": 2,
                        "matrix_memory_mode": mode,
                        "samples": entries,
                    }
                    (directory / (
                        f"memory_phase_rank{rank:04d}.json"
                    )).write_text(
                        json.dumps(report), encoding="utf-8"
                    )

    def test_comparison_of_all_phases_and_both_ranks(self):
        self._build_fake_reports()
        report = ANALYZER.analyze(self.root)
        self.assertIn("acoustic", report)
        self.assertIn("swe", report)
        self.assertLess(
            report["acoustic"]["after_mass_assembly"]["delta_mib"], 0
        )
        self.assertTrue(
            (self.root / "matrix_phase_memory_summary.json").is_file()
        )
        self.assertTrue(
            (self.root / "matrix_phase_memory_timeseries.csv").is_file()
        )

    def test_missing_checkpoint_rejected(self):
        self._build_fake_reports()
        file = self.root / "swe/compact/memory_phase_rank0001.json"
        report = json.loads(file.read_text(encoding="utf-8"))
        report["samples"] = [
            x for x in report["samples"]
            if x["phase"] != "before_main_ts_solve"
        ]
        file.write_text(json.dumps(report), encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "missing checkpoints"):
            ANALYZER.analyze(self.root)

    def test_mixed_commits_rejected(self):
        self._build_fake_reports()
        file = self.root / "acoustic/compact/memory_phase_rank0000.json"
        report = json.loads(file.read_text(encoding="utf-8"))
        report["source_commit"] = "wrong-sha"
        file.write_text(json.dumps(report), encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "source commit"):
            ANALYZER.analyze(self.root)


if __name__ == "__main__":
    unittest.main()
