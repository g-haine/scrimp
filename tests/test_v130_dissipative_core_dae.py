# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Source-only gates for second independent dissipative linear DAE.

Synthetic fixtures test that the acceptance comparator rejects
non-physical damping, energy defects, and MPI trajectory mismatch.
They do not qualify the native GetFEM/PETSc formulation.
"""

import ast
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
COMPARE = ROOT / "sandbox/compare_v130_dissipative_core_dae.py"
PROBE = ROOT / "sandbox/probe_v130_dissipative_core_dae.py"


class V130DissipativeCoreDAETest(unittest.TestCase):
    def setUp(self):
        folder = tempfile.TemporaryDirectory()
        self.addCleanup(folder.cleanup)
        self.folder = Path(folder.name)
        self._write()

    def _write(self, *, negative_dissipation=False, high_defect=False,
               trajectory_mismatch=False):
        for ranks in (1, 2):
            times = np.arange(11) * .002
            states = np.zeros((11, 5), dtype=float)
            states[:, 0] = np.arange(11)*.004
            states[:, 1] = 2.
            if trajectory_mismatch and ranks == 2:
                states[-1, 2] += .002
            positive = 1e-5
            physical = {
                "hamiltonian": [2.+k*.001 for k in range(11)],
                "balance_defects": [0.]*11,
                "step_defects": [0.]*10,
                "absolute_dissipation_work": positive,
                "integrated_dissipation": (
                    -positive if negative_dissipation else positive
                ),
                "absolute_signed_boundary_work": 4e-5,
                "integrated_signed_boundary_work": -4e-5,
                "max_absolute_balance_defect": (
                    .1 if high_defect else 1e-14
                ),
                "midpoint_port_powers": {
                    "Damping": [.01]*10,
                    "Boundary control (left)": [.02]*10,
                },
            }
            record = {
                "schema": "scrimp-v130-dissipative-second-linear-dae-v1",
                "rank_count": ranks,
                "mesh": .2, "dt": .002, "steps": 10,
                "initial_dynamic_drift": {"q": 0., "p": 0.},
                "dae_projection": {
                    "strategy": "algebraic",
                    "status": "PASS",
                    "differential_initial_exact": True,
                    "pseudo_ts_executed": False,
                    "accuracy_target_met": True,
                    "native_final_residual_l2": 1e-13,
                    "feasibility_tolerance": 1e-9,
                },
                "physical": physical,
                "power_tolerance": 1e-5,
            }
            stem = f"dissipative_core_{ranks}r"
            (self.folder / f"{stem}.json").write_text(
                json.dumps(record), encoding="utf-8"
            )
            np.savez_compressed(
                self.folder / f"{stem}.npz",
                times=times, states=states,
                initial=np.array([0., 2., 0., 0., 0.]),
            )

    def _run(self):
        return subprocess.run(
            [sys.executable, str(COMPARE), str(self.folder)],
            capture_output=True, text=True,
        )

    def test_source_includes_original_model_dissipation(self):
        src = PROBE.read_text(encoding="utf-8")
        self.assertIn(
            'examples/wave_dissipative.py', src
        )
        self.assertIn(
            'kwargs["init_step_ts_type"] = "algebraic"', src
        )
        self.assertIn("if use_default:", src)
        self.assertIn(
            '"Damping"', src
        )
        self.assertIn(
            'port.get_flow() + "*" + port.get_effort()', src
        )
        self.assertIn(
            'cumulative_dissipation', src
        )

    def test_getfem_port_region_conventions(self):
        # Execute only the isolated pure helper through its AST,
        # not the native GetFEM/PETSc import unavailable in source CI.
        # The real failure was an unlabelled (region=None) damping
        # port passed to GetFEM asm's scalar integer region argument.
        src = PROBE.read_text(encoding="utf-8")
        tree = ast.parse(src)
        funcs = [
            node for node in tree.body
            if isinstance(node, ast.FunctionDef)
            and node.name == "_getfem_region"
        ]
        self.assertEqual(len(funcs), 1)
        isolated = ast.Module(body=funcs, type_ignores=[])
        ns = {"np": np}
        exec(compile(isolated, str(PROBE), "exec"), ns)
        get_region = ns["_getfem_region"]

        class Port:
            def __init__(self, name, region):
                self.name = name
                self.region = region

            def get_region(self):
                return self.region

            def get_name(self):
                return self.name

        self.assertEqual(
            get_region(Port("Damping", None)), -1
        )
        self.assertEqual(
            get_region(Port("Boundary control (left)", 13)), 13
        )
        self.assertEqual(
            get_region(Port("Boundary control (top)", np.int64(12))), 12
        )
        with self.assertRaisesRegex(TypeError, "noninteger region"):
            get_region(Port("invalid", [1]))
        self.assertIn(
            "region = _getfem_region(port)", src
        )
        self.assertIn(
            "], 0, exp, region, model", src
        )

    def test_synthetic_positive_mpi_parity(self):
        result = self._run()
        self.assertEqual(result.returncode, 0, result.stderr)
        summary = json.loads((
            self.folder / "second_linear_dae_comparison.json"
        ).read_text(encoding="utf-8"))
        self.assertEqual(summary["scientific_status"], "PASS")

    def test_negative_dissipation_is_rejected(self):
        self._write(negative_dissipation=True)
        result = self._run()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Non-vacuous dissipative", result.stderr)

    def test_bad_energy_identity_is_rejected(self):
        self._write(high_defect=True)
        result = self._run()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("power gate FAIL", result.stderr)

    def test_later_mpi_state_discrepancy_is_rejected(self):
        self._write(trajectory_mismatch=True)
        result = self._run()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("trajectory mismatch", result.stderr)


if __name__ == "__main__":
    unittest.main()
