# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Synthetic gates for short linear/nonlinear compact matrix A/B analyzer."""
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "_scrimp_matrix_memory_comparison",
    ROOT / "sandbox" / "compare_compact_matrix.py",
)
ANALYZER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(ANALYZER)


SWE_FIELDS = (
    "t", "kinetic", "potential", "H", "d_grad", "d_div",
    "p_normal", "p_viscous", "supplied", "dissipated",
    "balance", "relative_balance", "relative_balance_exchange",
    "volume", "normal_volume_outflow", "cumulative_outflow",
    "volume_balance", "relative_volume_balance", "lyapunov",
)


class CompactMatrixComparisonTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        values = np.arange(20, dtype=np.float64).reshape(4, 5) * .001
        for case in ("acoustic", "swe"):
            for mode in ("legacy", "compact"):
                folder = self.root / case / mode
                folder.mkdir(parents=True)
                compressed = {"states": values}
                if case == "acoustic":
                    compressed.update({
                        "times": np.array([0., .001, .002, .003]),
                        "hamiltonian": np.linspace(0.1, .4, 4),
                        "power_0": np.linspace(.001, .004, 4),
                    })
                    report = {
                        "mode": "memory",
                        "matrix_memory_mode": mode,
                        "source_commit": "fakeSHA",
                        "mpi_size": 2, "dofs": 5,
                        "n_snapshots": 4, "steps": 3,
                        "rejected": 0, "snes_failures": 0,
                        "final_time": .003,
                    }
                    stem = "memory_physical"
                else:
                    compressed.update({
                        field: (np.array([0., .001, .002, .003])
                                if field == "t" else
                                np.linspace(.001, .004, 4))
                        for field in SWE_FIELDS
                    })
                    report = {
                        "mode": "memory",
                        "matrix_memory_mode": mode,
                        "source_commit": "fakeSHA",
                        "ranks": 2, "dofs": 5, "snapshots": 4,
                        "mesh": .35, "dt": .001, "final_time": .003,
                        "spinup_duration": .03, "ts_type": "cn",
                        "ts_adapt_type": "none",
                        "mumps_reset_policy": "on_change",
                        "convergence": {
                            "reason": 2, "steps": 3, "rejections": 0,
                            "snes_failures": 0, "t_final": .003,
                        },
                    }
                    stem = "memory_swe"
                flags = (case == "acoustic")
                report["matrix_memory_per_rank"] = [
                    {
                        "mpi_rank": rank,
                        "matrix_memory_mode": mode,
                        "linear_mass": flags,
                        "linear_stiffness": flags,
                        "tangent_mass_is_mass": mode == "compact" and flags,
                        "tangent_stiffness_is_stiffness":
                            mode == "compact" and flags,
                        "hessian_fd_matrix_allocations":
                            int(case == "swe" and mode == "compact"),
                        "hessian_fd_probe_allocations":
                            int(case == "swe" and mode == "compact"),
                        "hessian_fd_matrix_reuses":
                            3 if case == "swe" and mode == "compact" else 0,
                        "VmHWM_kib": (120000 if mode == "legacy"
                                        else 100000) + 100 * rank,
                        "VmRSS_kib": (110000 if mode == "legacy"
                                        else 90000) + 100 * rank,
                    }
                    for rank in range(2)
                ]
                np.savez_compressed(folder / (stem + ".npz"), **compressed)
                (folder / (stem + ".json")).write_text(
                    json.dumps(report), encoding="utf-8"
                )

    def test_full_parity_and_rank_metrics_pass(self):
        report = ANALYZER.compare(self.root)
        self.assertIn("swe", report["cases"])
        self.assertLess(
            report["cases"]["acoustic"][
                "peak_hwm_delta_mib_compact_minus_legacy"
            ], 0.
        )
        self.assertTrue(
            (self.root / "matrix_memory_comparison.json").is_file()
        )

    def test_changed_swe_balance_is_rejected(self):
        target = self.root / "swe/compact/memory_swe.npz"
        with np.load(target, allow_pickle=False) as file:
            arrays = {name: file[name] for name in file.files}
        arrays["balance"][1] = 0.1
        np.savez_compressed(target, **arrays)
        with self.assertRaisesRegex(ValueError, "swe/balance"):
            ANALYZER.compare(self.root)

    def test_missing_linear_matrix_alias_is_rejected(self):
        target = self.root / "acoustic/compact/memory_physical.json"
        data = json.loads(target.read_text(encoding="utf-8"))
        data["matrix_memory_per_rank"][1][
            "tangent_mass_is_mass"
        ] = False
        target.write_text(json.dumps(data), encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "mass tangent ownership"):
            ANALYZER.compare(self.root)

    def test_nonreused_nonlinear_hessian_is_rejected(self):
        target = self.root / "swe/compact/memory_swe.json"
        data = json.loads(target.read_text(encoding="utf-8"))
        data["matrix_memory_per_rank"][0][
            "hessian_fd_matrix_reuses"
        ] = 0
        target.write_text(json.dumps(data), encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "workspace was not reused"):
            ANALYZER.compare(self.root)


if __name__ == "__main__":
    unittest.main()
