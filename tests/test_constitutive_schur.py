# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""No-PETSc tests of the uncondensed Schur identity and MPI report gates."""
import copy
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "_compare_constitutive_schur",
    ROOT / "sandbox" / "compare_constitutive_mumps_schur.py",
)
CHECK = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(CHECK)


class ConstitutiveSchurTests(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.directory = Path(temp.name)
        self.base = {
            "schema": "scrimp-constitutive-schur-probe-v1",
            "source_commit": "abcdef",
            "mpi_ranks": 1,
            "mesh_h": .08,
            "dofs": 12,
            "coenergy_dofs": 5,
            "complement_dofs": 7,
            "coenergy_fields": ["e_q", "e_p"],
            "complement_fields": [
                "Y_B", "U_B", "q", "p",
            ],
            "constitutive_Q_frobenius_norm": 5.0,
            "constitutive_Q_symmetric": True,
            "constitutive_Q_min_diagonal": 1.0,
            "constitutive_Q_manufactured_solution_relative_error": 1e-14,
            "Q_solver": "PETSc KSPPREONLY + PCLU + MUMPS",
            "Q_factor_operator_configuration_count": 1,
            "Q_reused_across_shifts_and_rhs": True,
            "reference_pde": {
                "mumps_full_system_unchanged": True,
                "converged_reason": 2,
                "accepted_steps": 2,
            },
            "shift_results": [
                {
                    "sigma": val,
                    "relative_Aee_plus_Q": 1e-15,
                    "Q_MUMPS_manufactured_relative_residual": 2e-15,
                    "exact_schur_action_relative_error": 3e-15,
                    "full_DAE_algebraic_relative_residual": 4e-15,
                    "KSP_iterations": 1,
                    "KSP_converged_reason": 4,
                    "MUMPS_Q_operator_reused": True,
                }
                for val in (1000., 2000.)
            ],
        }
        self.reports = {
            1: copy.deepcopy(self.base),
            2: {**copy.deepcopy(self.base), "mpi_ranks": 2},
        }
        self.flush()

    def flush(self):
        for p, doc in self.reports.items():
            (self.directory / f"schur_{p}r.json").write_text(
                json.dumps(doc), encoding="utf-8"
            )

    def test_docker_mounts_scrimp_writable_default_outputs(self):
        # SCRIMP's import-time config.set_paths() creates
        # outputs/{log,png,mesh,pv} below the container workdir.
        # A /artifacts-only mount is not sufficient as the source
        # package/example/sandbox mounts are intentionally readonly.
        runner = (
            ROOT / "sandbox" / "run_constitutive_mumps_schur.sh"
        ).read_text(encoding="utf-8")
        self.assertIn(
            '-v "$out:/workspace/scrimp/outputs"', runner
        )
        self.assertIn('-v "$out:/artifacts"', runner)

    def test_block_schur_sign_and_uncondensed_residual(self):
        # Two *different* physical TS shifts. Original equations
        # with constitutive row -Q e + Cx=0, not e=Qz.
        Q = np.array([[3., .2], [.2, 2.]])
        C = np.array([[1., 4., 0.], [2., 0., .5]])
        D = np.array([[.1, 1.], [3., 2.], [1., -1.]])
        dynamic_mass = np.diag([2., 4., 1.])
        K = np.array([[2., .5, 0.], [0., 1., 1.], [1., 0., 3.]])
        x = np.array([.5, -.7, 1.2])
        e = np.linalg.solve(Q, C @ x)
        for sigma in (1000., 2000.):
            G = K+sigma*dynamic_mass
            A = np.block([[-Q, C], [D, G]])
            full = A @ np.r_[e, x]
            Sx = G @ x + D @ np.linalg.solve(Q, C @ x)
            np.testing.assert_allclose(full[:2], 0, atol=1e-12)
            np.testing.assert_allclose(
                full[2:], Sx, rtol=1e-12, atol=1e-12
            )
            np.testing.assert_allclose(A[:2, :2], -Q)

    def test_native_report_mpi_parity_passes(self):
        out = CHECK.compare(self.directory)
        self.assertEqual(out["MPI_ranks"], [1, 2])
        self.assertEqual(len(out["shift_parity"]), 2)
        self.assertTrue(
            (self.directory / "constitutive_schur_parity.json").is_file()
        )

    def test_changed_shift_constitutive_block_fails(self):
        self.reports[2]["shift_results"][1][
            "relative_Aee_plus_Q"
        ] = .01
        self.flush()
        with self.assertRaisesRegex(ValueError, "relative_Aee_plus_Q"):
            CHECK.compare(self.directory)

    def test_rank_dependent_constitutive_norm_fails(self):
        self.reports[2]["constitutive_Q_frobenius_norm"] = 6.
        self.flush()
        with self.assertRaisesRegex(ValueError, "Q norms differ"):
            CHECK.compare(self.directory)

    def test_rank_dependent_field_sizes_fail(self):
        self.reports[2]["coenergy_dofs"] = 4
        self.flush()
        with self.assertRaisesRegex(ValueError, "coenergy_dofs"):
            CHECK.compare(self.directory)

    def test_mumps_operator_reset_not_accepted(self):
        self.reports[2]["Q_factor_operator_configuration_count"] = 2
        self.flush()
        with self.assertRaisesRegex(ValueError, "provenance"):
            CHECK.compare(self.directory)

    def test_missing_variable_block_fails(self):
        self.reports[1]["complement_dofs"] = 1
        self.reports[2]["complement_dofs"] = 1
        self.flush()
        with self.assertRaisesRegex(ValueError, "cover all"):
            CHECK.compare(self.directory)


if __name__ == "__main__":
    unittest.main()
