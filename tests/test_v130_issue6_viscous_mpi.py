# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Synthetic negative contracts for #6 native MPI comparator, not physics."""

import importlib.util
import json
import tempfile
import unittest
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "sandbox/compare_v130_issue6_viscous_mpi.py"
spec = importlib.util.spec_from_file_location("v130_issue6_compare", SOURCE)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


class Issue6MPIComparatorTests(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.root = Path(tmp.name)
        self.write_valid()

    def write_valid(self):
        time = np.array([0., 0.01, 0.02])
        h = np.array([10., 9.99, 9.98])
        zero = np.zeros(3)
        energy = {
            "t": time, "hamiltonian": h,
            "kinetic": np.array([1., 0.99, 0.98]),
            "potential": np.array([9., 9., 9.]),
            "diss_grad": np.array([0., 1., 1.]),
            "diss_div": zero, "power_normal": zero,
            "power_viscous": zero, "supplied": zero,
            "dissipated": np.array([0., 0.01, 0.02]),
            "balance": zero, "relative_balance_exchange": zero,
            "volume": np.array([5., 5., 5.]),
            "cumulative_outflow": zero, "volume_balance": zero,
        }
        for case in ("open", "closed"):
            for ranks in (1, 2, 4):
                p = self.root / f"{case}_{ranks}r"
                p.mkdir()
                np.savez_compressed(p / "energies.npz", **energy)
                np.savez_compressed(p / "state_endpoints.npz",
                                    initial=np.array([5., 0.]),
                                    final=np.array([5., 0.02]))
                (p / "metrics.json").write_text(json.dumps({
                    "mpi_ranks": ranks, "petsc_ts_converged_reason": 1,
                    "petsc_ts_accepted_steps": 2, "petsc_ts_rejected_steps": 0,
                    "petsc_ts_snes_failures": 0,
                    "minimum_height_coefficient": 5.,
                    "maximum_absolute_balance_error_over_actual_exchange": 0.,
                    "maximum_absolute_relative_volume_balance_error": 0.,
                }), encoding="utf-8")
                for rank in range(ranks):
                    (p / f"mpi_port_powers_rank{rank}.json").write_text(
                        json.dumps({
                            "rank": rank, "world_size": ranks,
                            "times": time.tolist(),
                            "powers": {"algebraic": [0., 0.1, 0.2]},
                        }), encoding="utf-8")

    def compare(self, expect_rest=False):
        return module.compare(self.root, ["open", "closed"], 0.02,
                              0.05, 1e-4, not expect_rest,
                              expect_rest=expect_rest)

    def write_rest(self):
        """True stationary synthetic records with deliberately meaningless
        100% exchange-scaled residuals; still enforce absolute physics."""
        for case in ("open", "closed"):
            for ranks in (1, 2, 4):
                p = self.root / f"{case}_{ranks}r"
                with np.load(p / "energies.npz") as f:
                    data = {k: f[k].copy() for k in f.files}
                for key in ("hamiltonian", "potential"):
                    data[key][:] = 10.
                for key in ("kinetic", "diss_grad", "diss_div",
                            "power_normal", "power_viscous", "supplied",
                            "dissipated", "balance", "relative_balance_exchange"):
                    data[key][:] = 0.
                np.savez_compressed(p / "energies.npz", **data)
                metrics_file = p / "metrics.json"
                m = json.loads(metrics_file.read_text(encoding="utf-8"))
                m["maximum_absolute_balance_error_over_actual_exchange"] = 1.
                m["minimum_height_coefficient"] = 5.
                m["maximum_height_coefficient"] = 5.
                metrics_file.write_text(json.dumps(m), encoding="utf-8")
                for rank in range(ranks):
                    power_file = p / f"mpi_port_powers_rank{rank}.json"
                    power_doc = json.loads(power_file.read_text(encoding="utf-8"))
                    power_doc["powers"]["algebraic"] = [0., 0., 0.]
                    power_file.write_text(json.dumps(power_doc), encoding="utf-8")

    def test_stationary_equilibrium_ignores_ill_conditioned_exchange_ratio(self):
        self.write_rest()
        result = self.compare(expect_rest=True)
        self.assertEqual(result["gate"], "rest_equilibrium")
        self.assertFalse(result["exchange_scaled_metric_applicable"])
        self.assertIsNone(
            result["comparisons"]["open"]["max_exchange_scaled_power_balance_defect"]
        )

    def test_rest_rejects_hidden_kinetic_activity(self):
        self.write_rest()
        p = self.root / "open_2r" / "energies.npz"
        with np.load(p) as f:
            data = {k: f[k].copy() for k in f.files}
        data["kinetic"][-1] = 0.01
        np.savez_compressed(p, **data)
        with self.assertRaisesRegex(ValueError, "rest equilibrium violation"):
            self.compare(expect_rest=True)

    def test_rest_rejects_energy_drift_despite_zero_exchange(self):
        self.write_rest()
        p = self.root / "open_1r" / "energies.npz"
        with np.load(p) as f:
            data = {k: f[k].copy() for k in f.files}
        data["hamiltonian"][-1] += 0.01
        data["balance"][-1] += 0.01
        np.savez_compressed(p, **data)
        with self.assertRaisesRegex(ValueError, "rest equilibrium violation"):
            self.compare(expect_rest=True)

    def test_valid_synthetic_contract(self):
        result = self.compare()
        self.assertEqual(result["status"], "PASS")
        self.assertEqual(set(result["comparisons"]), {"open", "closed"})

    def test_missing_rank_trace_fails(self):
        (self.root / "open_4r" / "mpi_port_powers_rank3.json").unlink()
        with self.assertRaises(FileNotFoundError):
            self.compare()

    def test_rank_local_power_disagreement_fails(self):
        p = self.root / "open_4r" / "mpi_port_powers_rank1.json"
        d = json.loads(p.read_text(encoding="utf-8"))
        d["powers"]["algebraic"][-1] += 0.5
        p.write_text(json.dumps(d), encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "generic port"):
            self.compare()

    def test_serial_vs_mpi_endpoint_disagreement_fails(self):
        p = self.root / "closed_2r" / "state_endpoints.npz"
        np.savez_compressed(p, initial=np.array([5., 0.]),
                            final=np.array([4.5, 0.02]))
        with self.assertRaisesRegex(ValueError, "full final"):
            self.compare()

    def test_nonpositive_height_fails(self):
        p = self.root / "open_2r" / "metrics.json"
        d = json.loads(p.read_text(encoding="utf-8"))
        d["minimum_height_coefficient"] = -0.01
        p.write_text(json.dumps(d), encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "nonpositive sampled"):
            self.compare()

    def test_energy_exchange_error_fails(self):
        p = self.root / "closed_1r" / "metrics.json"
        d = json.loads(p.read_text(encoding="utf-8"))
        d["maximum_absolute_balance_error_over_actual_exchange"] = 0.4
        p.write_text(json.dumps(d), encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "physical power defect"):
            self.compare()

    def test_inconsistent_signed_power_balance_fails(self):
        p = self.root / "open_1r" / "energies.npz"
        with np.load(p) as f:
            data = {k: f[k] for k in f.files}
        data["balance"] = np.array([0., 0., 0.2])
        np.savez_compressed(p, **data)
        with self.assertRaisesRegex(ValueError, "signed physical power identity"):
            self.compare()


if __name__ == "__main__":
    unittest.main()
