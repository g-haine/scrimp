# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Synthetic plotting gates; never presented as measured timings."""
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

SOURCE = (Path(__file__).resolve().parents[1] /
          "sandbox" / "plot_release120_coo_scaling.py")
SPEC = importlib.util.spec_from_file_location("_release_plot", SOURCE)
PLOT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(PLOT)


class Release120ComparisonPlotTests(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.dir = Path(tmp.name)
        self.runs = []
        for version, pols, ranks in (
            ("release_120", ("release_native",), (1,)),
            ("current", ("never", "on_change"), (1, 2, 4, 6)),
        ):
            for p in pols:
                for rank in ranks:
                    for rep in (1, 2):
                        filename = (f"{version}_{p}_{rank}r_rep{rep}")
                        timing = (
                            6.0 if version == "release_120"
                            else (1.3 + 0.07 * (p == "on_change")) /
                            np.sqrt(rank)
                        )
                        timing += .01 * rep
                        states = np.array([0., .1, .5, 2., 3.])
                        json_obj = {
                            "benchmark": "examples/wave.py:wave_eq()",
                            "variant": version,
                            "source_sha": (
                                "releasesha" if version == "release_120"
                                else "newsha"
                            ),
                            "source_tag": (
                                "v1.2.0" if version == "release_120" else ""
                            ),
                            "repetition": rep, "mpi_ranks": rank,
                            "mesh_h": .08, "dofs": 5,
                            "scheme": "cn", "ts_adapt_type": "none",
                            "dt": .002,
                            "time_final_requested": .04,
                            "time_final_monitor": .04,
                            "save_dt": .01, "steps": 20,
                            "rejections": 0, "snes_failures": 0,
                            "ts_convergence_reason": (
                                None if version == "release_120" else 2
                            ),
                            "matrix_backend": (
                                "release_v1.2.0_native"
                                if version == "release_120" else
                                "distributed_coo"
                            ),
                            "mumps_policy": p,
                            "wall_build_s": .04,
                            "wall_solve_s": float(timing),
                            "wall_total_s": float(timing + .04),
                            "state_file": filename + "_state.npz",
                        }
                        np.savez_compressed(
                            self.dir / (filename + "_state.npz"),
                            initial=states, final=states + .5
                        )
                        (self.dir / (filename + ".json")).write_text(
                            json.dumps(json_obj), encoding="utf-8"
                        )

    def test_full_1_2_4_6_grid_produces_all_figures(self):
        result = PLOT.plot(self.dir, 2)
        self.assertEqual(result["provenance"]["baseline_tag"], "v1.2.0")
        self.assertEqual(result["provenance"]["repetitions"], 2)
        for name in (
            "01_version120_vs_coo_serial",
            "02_strong_scaling_coo",
            "03_speedup_coo",
            "04_efficiency_coo",
        ):
            self.assertGreater(
                (self.dir / (name + ".png")).stat().st_size, 3000
            )
            self.assertGreater(
                (self.dir / (name + ".svg")).stat().st_size, 3000
            )

    def test_current_only_salvages_scaling_after_v120_state_mismatch(self):
        # Reproduce the observed failure: release v1.2.0 states do NOT
        # match current even before the first CN step. The valid current
        # MPI series still merits 3 figures, with NO release speedup.
        baseline_state = (
            self.dir /
            "release_120_release_native_1r_rep1_state.npz"
        )
        with np.load(baseline_state, allow_pickle=False) as file:
            reference_initial = file["initial"]
            reference_final = file["final"]
        np.savez_compressed(
            baseline_state,
            initial=reference_initial + .4,
            final=reference_final + 2.0,
        )
        with self.assertRaisesRegex(ValueError, "PARITY FAILED"):
            PLOT.plot(self.dir, 2)
        result = PLOT.plot(self.dir, 2, current_only=True)
        self.assertFalse(result["release_comparison_accepted"])
        self.assertEqual(result["provenance"]["baseline_tag"], None)
        self.assertIsNone(result["baseline_median_s"])
        self.assertFalse(
            (self.dir / "01_version120_vs_coo_serial.png").exists()
        )
        for stem in (
            "02_strong_scaling_coo", "03_speedup_coo",
            "04_efficiency_coo",
        ):
            self.assertGreater(
                (self.dir / (stem + ".png")).stat().st_size, 3000
            )

    def test_current_only_still_rejects_rank_dependent_state_drift(self):
        file = self.dir / "current_never_4r_rep2_state.npz"
        with np.load(file, allow_pickle=False) as data:
            initial, final = data["initial"], data["final"]
        np.savez_compressed(
            file, initial=initial, final=final + 0.3
        )
        with self.assertRaisesRegex(ValueError, "PARITY FAILED"):
            PLOT.plot(self.dir, 2, current_only=True)

    def test_missing_MPI_measurements_cannot_be_fabricated(self):
        target = self.dir / "current_on_change_6r_rep2.json"
        target.unlink()
        with self.assertRaises(FileNotFoundError):
            PLOT.plot(self.dir, 2)

    def test_rejects_physical_drift_not_just_timing(self):
        target = self.dir / "current_never_4r_rep2_state.npz"
        with np.load(target, allow_pickle=False) as f:
            initial = f["initial"]
            final = f["final"]
        final[-1] += .5
        np.savez_compressed(target, initial=initial, final=final)
        with self.assertRaisesRegex(ValueError, "PARITY FAILED"):
            PLOT.plot(self.dir, 2)

    def test_wrong_source_cannot_be_called_v120(self):
        file = self.dir / "release_120_release_native_1r_rep1.json"
        data = json.loads(file.read_text(encoding="utf-8"))
        data["source_tag"] = "main"
        file.write_text(json.dumps(data), encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "NOT the original"):
            PLOT.plot(self.dir, 2)


if __name__ == "__main__":
    unittest.main()
