# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Synthetic-only current v1.3.0 3M-DOF analyzer/plot regression.

No synthetic timings are reported as real SCRIMP measurements.
"""
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "sandbox"))
SPEC = importlib.util.spec_from_file_location(
    "v130_large_wave_analysis_test",
    ROOT / "sandbox" / "analyze_v130_large_wave_scaling.py",
)
ANALYZE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(ANALYZE)

RANKS = (1, 2, 4, 6)
SHA = "synthetic-commit-for-source-checks"


class LargeWaveAnalysisTests(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.dir = Path(temp.name)
        self.repetitions = 2
        self.dofs = 24
        self.ranks = (1, 2)
        self.manifest = {
            "schema": "scrimp-v130-current-large-wave-v1",
            "current_sha": SHA,
            "mesh": .005,
            "repetitions": self.repetitions,
            "ranks": "1,2",
            "dt": .002,
            "final_time": .02,
            "steps": 10,
            "source": "examples/wave.py",
            "variant": "current",
            "mumps_policy": "on_change",
        }
        (self.dir / "manifest.json").write_text(json.dumps(self.manifest))
        self.data = {}
        for rep in range(1, self.repetitions+1):
            for rank in self.ranks:
                label = f"current_on_change_{rank}r_rep{rep}"
                self.data[(rank, rep)] = self.make_record(rank, rep)
                self.write_record(label, self.data[(rank, rep)])
                np.savez_compressed(
                    self.dir / (label + "_state.npz"),
                    initial=np.linspace(0., 1., self.dofs),
                    final=np.linspace(1., 2., self.dofs),
                )
                (self.dir / (label + ".log")).write_text(
                    "Event Stage 1: WAVE_INITIAL_CONSISTENCY\n"
                    "MatLUFactorNum 1 1.0 .02 extra\n"
                    "Event Stage 2: WAVE_MAIN_CN\n"
                    "MatLUFactorSym 1 1.0 .03 extra\n"
                    "MatLUFactorNum 2 1.0 .04 extra\n"
                    "MatSolve 10 1.0 .01 extra\n",
                    encoding="utf-8",
                )

    def make_record(self, rank, rep):
        time = 2.4 / rank + .18 + rep*.01
        rss = [30000 + rank*200 + j for j in range(rank)]
        phase_names = (
            "setup.getfem.mass", "setup.bridge.mass",
            "setup.getfem.stiffness", "setup.bridge.stiffness",
            "initial.consistency_total", "main.petsc.ts_solve",
            "main.petsc.callback.ijacobian",
            "main.petsc.callback.ifunction",
            "main.petsc.csr_pattern_check",
        )
        return {
            "variant": "current", "source_sha": SHA,
            "repetition": rep, "mpi_ranks": rank,
            "mesh_h": .005, "dofs": self.dofs,
            "scheme": "cn", "ts_adapt_type": "none",
            "steps": 10, "rejections": 0, "snes_failures": 0,
            "ts_convergence_reason": 2, "time_final_monitor": .02,
            "dt": .002, "matrix_backend": "distributed_coo",
            "matrix_memory_mode": "compact",
            "mumps_policy": "on_change",
            "main_cn_mumps_icntl28": 2 if rank > 1 else 1,
            "wall_solve_s": time,
            "wall_build_s": .08,
            "wall_total_s": time + .08,
            "first_step_wall_s": .13 / rank,
            "steady_step_mean_s": .07 / rank,
            "step_times_max_rank_s": [.13 / rank] + [.07/rank]*9,
            "profile": {
                "enabled": True,
                "events": {name: {"seconds_rank_max": .1/rank}
                           for name in phase_names},
            },
            "rank_hwm_kib": rss,
            "max_rank_hwm_kib": max(rss),
            "sum_rank_hwm_kib": sum(rss),
            "mumps_reuse_stats": {"main": {
                "ijacobian_calls": 10, "numeric_updates": 2,
                "symbolic_updates": 1, "reused_operators": 8,
            }},
        }

    def write_record(self, label, record):
        (self.dir / (label + ".json")).write_text(
            json.dumps(record), encoding="utf-8"
        )

    def test_reports_validate_and_generate_all_png_and_svg(self):
        result = ANALYZE.analyze(
            self.dir, list(self.ranks), self.repetitions,
        )
        self.assertTrue(result["all_full_state_parity"])
        self.assertTrue(result["all_exact_mumps_factorization_gates"])
        self.assertEqual(len(result["rows"]), 2)
        self.assertEqual(len(result["figures_png"]), 12)
        self.assertEqual(len(result["figures_svg"]), 12)
        for path in result["figures_png"] + result["figures_svg"]:
            self.assertGreater((self.dir / path).stat().st_size, 1000)
        self.assertTrue(
            (self.dir / "v130_large_wave_performance.md").is_file()
        )
        self.assertTrue(
            (self.dir / "v130_large_wave_scaling.csv").is_file()
        )
        self.assertTrue(
            (self.dir / "v130_large_wave_runs.csv").is_file()
        )

    def test_wrong_petsc_numeric_lu_fails(self):
        log = self.dir / "current_on_change_2r_rep1.log"
        log.write_text(log.read_text().replace(
            "MatLUFactorNum 2 1.0 .04", "MatLUFactorNum 5 1.0 .04"
        ))
        with self.assertRaisesRegex(RuntimeError, "MUMPS"):
            ANALYZE.load_and_validate(
                self.dir, self.manifest, list(self.ranks), self.repetitions
            )

    def test_incomplete_rss_fails(self):
        self.data[(2, 1)]["rank_hwm_kib"] = [100]
        self.write_record("current_on_change_2r_rep1", self.data[(2, 1)])
        with self.assertRaisesRegex(RuntimeError, "RSS"):
            ANALYZE.load_and_validate(
                self.dir, self.manifest, list(self.ranks), self.repetitions
            )

    def test_mpi_state_mismatch_fails(self):
        path = self.dir / "current_on_change_2r_rep2_state.npz"
        with np.load(path) as z:
            initial, final = z["initial"].copy(), z["final"].copy()
        final[3] += .02
        np.savez_compressed(path, initial=initial, final=final)
        with self.assertRaisesRegex(RuntimeError, "full-state parity"):
            ANALYZE.load_and_validate(
                self.dir, self.manifest, list(self.ranks), self.repetitions
            )

    def test_no_fake_release_reference(self):
        self.assertNotIn("release_sha", self.manifest)
        self.assertNotIn("release_120", self.manifest["variant"])


if __name__ == "__main__":
    unittest.main()
