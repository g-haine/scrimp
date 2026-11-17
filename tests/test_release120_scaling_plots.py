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
import subprocess
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
                        # Synthetic timing records for ten accepted CN
                        # steps. Rank 1 is deliberately slower, and
                        # step 1 carries a large startup penalty.
                        nsteps = 10
                        step_base = timing / 13.
                        step_profile = np.array(
                            [3 * step_base] +
                            [step_base * (1. + 0.004 * j)
                             for j in range(nsteps - 1)]
                        )
                        rank_steps = [
                            (step_profile * (
                                1. if rank == 1 else
                                .95 + .05 * r / (rank - 1)
                            )).tolist()
                            for r in range(rank)
                        ]
                        monitor_costs = [
                            [0.0005] * (nsteps + 1)
                            for _ in range(rank)
                        ]
                        steady = step_profile[1:]
                        half = len(steady) // 2
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
                            "time_final_requested": .02,
                            "time_final_monitor": .02,
                            "save_dt": .01, "steps": nsteps,
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
                            "step_warmup_excluded": 1,
                            "step_indices": list(range(1, nsteps + 1)),
                            "step_physical_times": (
                                np.arange(1, nsteps+1) * .002
                            ).tolist(),
                            "step_times_max_rank_s":
                                step_profile.tolist(),
                            "step_times_by_rank_s": rank_steps,
                            "monitor_times_by_rank_s": monitor_costs,
                            "first_step_wall_s": float(step_profile[0]),
                            "steady_step_count": nsteps - 1,
                            "steady_step_mean_s": float(np.mean(steady)),
                            "steady_step_median_s": float(np.median(steady)),
                            "steady_step_early_mean_s":
                                float(np.mean(steady[:half])),
                            "steady_step_late_mean_s":
                                float(np.mean(steady[-half:])),
                            "steady_step_window_size": half,
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
            "05_steady_step_time_coo",
            "06_steady_step_speedup_coo",
            "07_steady_step_efficiency_coo",
            "08_first_step_overhead_coo",
        ):
            self.assertGreater(
                (self.dir / "png" / (name + ".png")).stat().st_size, 3000
            )
            self.assertGreater(
                (self.dir / "svg" / (name + ".svg")).stat().st_size, 3000
            )
            # No visual artefacts at the root of the output folder.
            self.assertFalse((self.dir / (name + ".png")).exists())
            self.assertFalse((self.dir / (name + ".svg")).exists())

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
            (self.dir / "png" / "01_version120_vs_coo_serial.png").exists()
        )
        for stem in (
            "02_strong_scaling_coo", "03_speedup_coo",
            "04_efficiency_coo",
            "05_steady_step_time_coo",
            "06_steady_step_speedup_coo",
            "07_steady_step_efficiency_coo",
            "08_first_step_overhead_coo",
        ):
            self.assertGreater(
                (self.dir / "png" / (stem + ".png")).stat().st_size, 3000
            )
            self.assertGreater(
                (self.dir / "svg" / (stem + ".svg")).stat().st_size, 3000
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

    def test_mpi_only_shell_cli_does_not_require_release(self):
        launcher = SOURCE.with_name("run_release120_coo_scaling.sh")
        result = subprocess.run(
            ["bash", str(launcher), "--mpi-only", "--help"],
            capture_output=True, text=True, timeout=10,
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("--mpi-only", result.stdout)
        self.assertIn("No v1.2.0 tag, worktree", result.stdout)
        result = subprocess.run(
            ["bash", str(launcher), "--mpi-only", "0.025", "0"],
            capture_output=True, text=True, timeout=10,
        )
        self.assertEqual(result.returncode, 2)
        # Tag/worktree setup is guarded by !mpi_only.
        content = launcher.read_text(encoding="utf-8")
        self.assertIn("if (( mpi_only )); then", content)
        self.assertIn("if (( !mpi_only )); then", content)
        self.assertIn("plot_flags+=(--current-only)", content)
        syntax = subprocess.run(
            ["bash", "-n", str(launcher)],
            capture_output=True, text=True, timeout=10,
        )
        self.assertEqual(syntax.returncode, 0, syntax.stderr)

    def test_tampered_mpi_step_maximum_is_rejected(self):
        target = self.dir / "current_never_4r_rep1.json"
        data = json.loads(target.read_text(encoding="utf-8"))
        data["step_times_max_rank_s"][3] *= 0.5
        target.write_text(json.dumps(data), encoding="utf-8")
        with self.assertRaisesRegex(
            ValueError, "incorrect MPI-max step wall times"
        ):
            PLOT.plot(self.dir, 2, current_only=True)

    def test_tampered_steady_step_mean_is_rejected(self):
        target = self.dir / "current_on_change_2r_rep2.json"
        data = json.loads(target.read_text(encoding="utf-8"))
        data["steady_step_mean_s"] *= 2.
        target.write_text(json.dumps(data), encoding="utf-8")
        with self.assertRaisesRegex(
            ValueError, "inconsistent steady_step_mean_s"
        ):
            PLOT.plot(self.dir, 2, current_only=True)

    def test_old_campaign_keeps_total_plots_without_faked_steps(self):
        # Older JSON has no clocked accepted steps. Tsolve/N must never
        # masquerade as an observed post-startup step cost.
        prefixes = (
            "step_", "steady_step_", "first_step_", "monitor_times_"
        )
        for path in self.dir.glob("*.json"):
            data = json.loads(path.read_text(encoding="utf-8"))
            data = {
                key: val for key, val in data.items()
                if not key.startswith(prefixes)
            }
            path.write_text(json.dumps(data), encoding="utf-8")
        result = PLOT.plot(self.dir, 2, current_only=True)
        self.assertFalse(result["direct_step_timing_available"])
        self.assertIsNone(result["step_timing"])
        self.assertTrue(
            (self.dir / "png" / "03_speedup_coo.png").exists()
        )
        self.assertFalse(
            (self.dir / "png" / "06_steady_step_speedup_coo.png").exists()
        )

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
