# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Pure SciPy tests for finding whether native GetFEM or COO loses mass."""
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np
import scipy.sparse as sp

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "_scrimp_bridge_localizer",
    ROOT / "sandbox" / "compare_raw_getfem_petsc.py",
)
LOCALIZER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(LOCALIZER)


class GetfemBridgeLocalizationTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.folder = Path(self.temp.name)
        self.matrix = sp.csr_matrix(
            [[5., 0., 1., 0.],
             [0., 4., 0., 0.],
             [1., 0., 3., 0.],
             [0., 0., 0., 2.]],
            dtype=np.float64,
        )
        fields = [
            {"name": "q", "start": 0, "stop": 2},
            {"name": "p", "start": 2, "stop": 4},
        ]
        for n in (1, 2):
            (self.folder / f"blocks_{n}r.json").write_text(
                json.dumps({
                    "dofs": 4, "mesh": .12, "fields": fields,
                    "physical_model": "mock-coenergy",
                    "capture_matrices": True,
                    "mesh_sha256": "same-test-mesh",
                    "getfem_mpi_level": 2,
                    "getfem_module_path": "mock-getfem",
                }),
                encoding="utf-8",
            )

    def store_raw(self, label, n, rank, matrix):
        matrix = sp.csr_matrix(matrix)
        np.savez_compressed(
            self.folder / f"raw_getfem_{label}_{n}r_rank{rank}.npz",
            shape=matrix.shape, indptr=matrix.indptr,
            indices=matrix.indices, data=matrix.data,
        )

    def store_petsc(self, label, n, matrix):
        matrix = sp.csr_matrix(matrix)
        nr = matrix.shape[0]
        ranges = [(0, nr)] if n == 1 else [(0, 2), (2, nr)]
        for rank, (lo, hi) in enumerate(ranges):
            block = matrix[lo:hi].tocsr()
            np.savez_compressed(
                self.folder / f"petsc_{label}_{n}r_rank{rank}.npz",
                lo=lo, hi=hi,
                shape=matrix.shape, indptr=block.indptr,
                indices=block.indices, data=block.data,
            )

    def make_reports(self, *, raw_2_sum_scale, petsc_2_scale):
        for label in ("mass", "stiffness"):
            self.store_raw(label, 1, 0, self.matrix)
            self.store_petsc(label, 1, self.matrix)
            for rank in (0, 1):
                # Two disjoint contributions are not required for this
                # synthetic test; each rank contributes half the
                # prescribed native global sum.
                self.store_raw(
                    label, 2, rank,
                    .5 * raw_2_sum_scale * self.matrix,
                )
            self.store_petsc(
                label, 2, petsc_2_scale * self.matrix
            )

    def test_native_getfem_wrong_but_coo_preserves_sum(self):
        self.make_reports(raw_2_sum_scale=.5, petsc_2_scale=.5)
        data = LOCALIZER.analyze(self.folder)
        self.assertTrue(data["same_mesh_sha256"])
        mass = data["matrices"]["mass"]
        self.assertAlmostEqual(
            mass["native_sum_2r_vs_native_1r"][
                "relative_fro_error"
            ], .5, places=12
        )
        self.assertAlmostEqual(
            mass["petsc_2r_vs_native_sum_2r"][
                "relative_fro_error"
            ], 0., places=12
        )

    def test_getfem_raw_correct_but_coo_drops_half(self):
        self.make_reports(raw_2_sum_scale=1., petsc_2_scale=.5)
        data = LOCALIZER.analyze(self.folder)
        mass = data["matrices"]["mass"]
        self.assertAlmostEqual(
            mass["native_sum_2r_vs_native_1r"][
                "relative_fro_error"
            ], 0., places=12
        )
        self.assertAlmostEqual(
            mass["petsc_2r_vs_native_sum_2r"][
                "relative_fro_error"
            ], .5, places=12
        )

    def test_missing_native_snapshot_rejected(self):
        self.make_reports(raw_2_sum_scale=1., petsc_2_scale=1.)
        (self.folder / "raw_getfem_mass_2r_rank1.npz").unlink()
        with self.assertRaises(FileNotFoundError):
            LOCALIZER.analyze(self.folder)

    def test_different_mesh_hash_warns_not_parity(self):
        self.make_reports(raw_2_sum_scale=1., petsc_2_scale=1.)
        path = self.folder / "blocks_2r.json"
        data = json.loads(path.read_text())
        data["mesh_sha256"] = "other-mesh"
        path.write_text(json.dumps(data))
        result = LOCALIZER.analyze(self.folder)
        self.assertFalse(result["same_mesh_sha256"])


if __name__ == "__main__":
    unittest.main()
