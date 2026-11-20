# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Synthetic negative controls for the new true-implicit-default native gate.

No source-only test proves GetFEM/PETSc physics or memory correctness.
"""

import ast
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
VERIFY = ROOT / "sandbox/verify_v130_default_promotion.py"
CORE_PROBE = ROOT / "sandbox/probe_v130_core_dae_initial.py"
DISSIPATIVE_PROBE = ROOT / "sandbox/probe_v130_dissipative_core_dae.py"


class TestV130DefaultPromotion(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.folder = Path(tmp.name)
        self._write()

    def _json(self, name, value):
        (self.folder / name).write_text(
            json.dumps(value), encoding="utf-8",
        )

    def _write(self):
        self._json("core_dae_qualification.json", {
            "qualified_opt_in_core_DAE": "PASS",
            "prototype_full_state_parity": {
                "1r": {"complete_trajectory_max_scaled_l2": 0.},
                "2r": {"complete_trajectory_max_scaled_l2": 1e-13},
            },
        })
        self._json("second_linear_dae_comparison.json", {
            "scientific_status": "PASS",
            "max_scaled_full_CN_MPI_discrepancy": 1e-13,
        })
        for ranks in (1, 2):
            for scenario in ("polynomial_compatible", "published"):
                negative = scenario == "published"
                self._json(
                    f"core_algebraic_{scenario}_{ranks}r.json",
                    {
                        "mpi_ranks": ranks,
                        "scenario": scenario,
                        "mode": "algebraic",
                        "effective_PETSc_initial_strategy": "algebraic",
                        "status": (
                            "EXPECTED_REJECTION" if negative else "PASS"
                        ),
                        "implicit_default_requested": True,
                        "explicit_init_step_ts_type_present": False,
                        "initialization_options": {"init_step": "True"},
                        "physical_cn_executed": not negative,
                        "initial_data_preserved": True,
                        "init_report": {
                            "strategy": "algebraic",
                            "pseudo_ts_executed": False,
                            "differential_initial_exact": True,
                        },
                    },
                )
            self._json(f"dissipative_core_{ranks}r.json", {
                "implicit_default_requested": True,
                "explicit_init_step_ts_type_present": False,
                "effective_PETSc_initial_strategy": "algebraic",
                "dae_projection": {
                    "strategy": "algebraic",
                    "pseudo_ts_executed": False,
                },
            })
        for mode, explicit in (("pseudo", "pseudo"), ("bdf2", "bdf")):
            self._json(
                f"core_{mode}_polynomial_compatible_1r.json",
                {
                    "mode": mode,
                    "status": "PASS",
                    "effective_PETSc_initial_strategy": explicit,
                    "implicit_default_requested": False,
                    "explicit_init_step_ts_type_present": True,
                    "initialization_options": {
                        "init_step_ts_type": explicit,
                    },
                },
            )

    def _run(self):
        return subprocess.run(
            [sys.executable, str(VERIFY), str(self.folder)],
            capture_output=True, text=True,
        )

    def test_correct_implicit_default_fixture_passes(self):
        result = self._run()
        self.assertEqual(result.returncode, 0, result.stderr)
        report = json.loads((
            self.folder / "default_promotion_qualification.json"
        ).read_text(encoding="utf-8"))
        self.assertEqual(report["status"], "PASS")
        self.assertTrue(report["actual_default_option_omitted"])
        self.assertFalse(report["new_default_is_released"])

    def test_explicit_projection_disguised_as_default_fails(self):
        path = self.folder / "core_algebraic_polynomial_compatible_1r.json"
        data=json.loads(path.read_text(encoding="utf-8"))
        data["explicit_init_step_ts_type_present"] = True
        self._json(path.name, data)
        result=self._run()
        self.assertNotEqual(result.returncode,0)
        self.assertIn("actual default", result.stderr)

    def test_false_default_pet_sc_strategy_is_rejected(self):
        path = self.folder / "core_algebraic_polynomial_compatible_1r.json"
        data = json.loads(path.read_text(encoding="utf-8"))
        data["effective_PETSc_initial_strategy"] = "pseudo"
        self._json(path.name, data)
        result = self._run()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("actual default", result.stderr)

    def test_false_dissipative_default_pet_sc_strategy_is_rejected(self):
        path = self.folder / "dissipative_core_2r.json"
        data = json.loads(path.read_text(encoding="utf-8"))
        data["effective_PETSc_initial_strategy"] = "pseudo"
        self._json(path.name, data)
        result = self._run()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Dissipative wave", result.stderr)

    def test_gaussian_must_not_be_accepted_by_default(self):
        path=self.folder / "core_algebraic_published_2r.json"
        data=json.loads(path.read_text(encoding="utf-8"))
        data["status"]="PASS"
        self._json(path.name,data)
        result=self._run()
        self.assertNotEqual(result.returncode,0)
        self.assertIn("incompatible Gaussian",result.stderr)

    def test_explicit_legacy_bdf_must_remain_selectable(self):
        path=self.folder / "core_bdf2_polynomial_compatible_1r.json"
        data=json.loads(path.read_text(encoding="utf-8"))
        data["initialization_options"]["init_step_ts_type"]="algebraic"
        self._json(path.name,data)
        result=self._run()
        self.assertNotEqual(result.returncode,0)
        self.assertIn("historical",result.stderr.lower())

    def test_builder_defaults_are_cleared_without_touching_explicit_controls(self):
        # Execute only the pure guard function from its AST: no
        # GetFEM, PETSc, mpi4py or physical solver in source-only CI.
        tree = ast.parse(CORE_PROBE.read_text(encoding="utf-8"))
        node = next(
            x for x in tree.body
            if isinstance(x, ast.FunctionDef)
            and x.name == (
                "_clear_benchmark_initializer_for_implicit_default"
            )
        )
        code = compile(
            ast.Module(body=[node], type_ignores=[]),
            str(CORE_PROBE), "exec",
        )
        scope = {}
        exec(code, scope)
        reset = scope[
            "_clear_benchmark_initializer_for_implicit_default"
        ]

        class Options:
            def __init__(self, value=None):
                self.value = value
                self.deleted = False

            def hasName(self, name):
                if name != "init_step_ts_type":
                    raise AssertionError(name)
                return self.value is not None

            def getString(self, name, default=""):
                return self.value if self.value is not None else default

            def delValue(self, name):
                if name != "init_step_ts_type":
                    raise AssertionError(name)
                self.deleted = True
                self.value = None

        for inherited in ("pseudo", "algebraic"):
            with self.subTest(inherited=inherited):
                opts = Options(inherited)
                reset(opts, enabled=True)
                self.assertIsNone(opts.value)
                self.assertTrue(opts.deleted)
        for inherited in ("pseudo", "algebraic", "bdf"):
            with self.subTest(explicit=inherited):
                opts = Options(inherited)
                reset(opts, enabled=False)
                self.assertEqual(opts.value, inherited)
                self.assertFalse(opts.deleted)
        opts = Options()
        reset(opts, enabled=True)
        self.assertFalse(opts.deleted)
        for bad in ("bdf", "unrecognized"):
            with self.subTest(unsupported_builder=bad):
                opts = Options(bad)
                with self.assertRaisesRegex(
                    RuntimeError, "Unexpected inherited benchmark"
                ):
                    reset(opts, enabled=True)
                self.assertEqual(opts.value, bad)
                self.assertFalse(opts.deleted)

    def test_probes_actually_omit_default_option_in_native_mode(self):
        wave=CORE_PROBE.read_text(encoding="utf-8")
        diss=DISSIPATIVE_PROBE.read_text(encoding="utf-8")
        self.assertIn('if not args.use_default:',wave)
        self.assertIn(
            "_clear_benchmark_initializer_for_implicit_default(",
            wave,
        )
        self.assertIn(
            'options.delValue("init_step_ts_type")',
            wave,
        )
        self.assertIn('kwargs["init_step_ts_type"] = option',wave)
        self.assertIn('if use_default:',diss)
        self.assertIn('kwargs["init_step_ts_type"] = "algebraic"',diss)
        self.assertIn(
            'kwargs.pop("init_step_ts_type", None)',diss
        )
        self.assertIn(
            "system.time_scheme.hasName(\"init_step_ts_type\")",diss
        )


if __name__ == "__main__":
    unittest.main()
