# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Pure-Python correctness gates for pH-DAE field/CSR classification."""
import importlib.util
from pathlib import Path
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "_scrimp_ph_block_layout",
    ROOT / "scrimp" / "utils" / "ph_block_layout.py",
)
PH = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(PH)


class _FakePort:
    def __init__(self, flow, effort, algebraic):
        self.flow = flow
        self.effort = effort
        self.algebraic = algebraic

    def get_flow(self):
        return self.flow

    def get_effort(self):
        return self.effort

    def get_algebraic(self):
        return self.algebraic


class _FakeGetFEM:
    def __init__(self):
        # Deliberately interleaved state/coenergy and boundary indices.
        self.intervals = {
            "q": (0, 1), "e_q": (1, 1),
            "Y_B": (2, 1), "p": (3, 1),
            "U_B": (4, 1), "e_p": (5, 1),
        }

    def nbdof(self):
        return 6

    def interval_of_variable(self, name):
        if name not in self.intervals:
            raise RuntimeError("Variable not registered")
        return self.intervals[name]


class PhBlockLayoutTests(unittest.TestCase):
    def setUp(self):
        self.model = _FakeGetFEM()
        self.kwargs = {
            "states": {"q": object(), "p": object()},
            "costates": {"e_q": object(), "e_p": object()},
            "ports": {
                "q": _FakePort("q", "e_q", False),
                "p": _FakePort("p", "e_p", False),
                "B": _FakePort("Y_B", "U_B", True),
            },
            "controls": {
                "B": _FakePort("Y_B", "U_B", True)
            },
        }

    def fields(self):
        return PH.discover_layout(self.model, **self.kwargs)

    def test_exact_global_coverage_and_noncontiguous_groups(self):
        fields = self.fields()
        self.assertEqual(
            [x["name"] for x in fields],
            ["q", "e_q", "Y_B", "p", "U_B", "e_p"],
        )
        self.assertEqual(
            [x["group"] for x in fields],
            ["dynamic", "coenergy", "boundary",
             "dynamic", "boundary", "coenergy"],
        )
        a = PH.owned_group_ranges(fields, 0, 3)
        b = PH.owned_group_ranges(fields, 3, 6)
        self.assertEqual(a["coenergy"], [[1, 2]])
        self.assertEqual(b["coenergy"], [[5, 6]])
        self.assertEqual(a["boundary"], [[2, 3]])
        self.assertEqual(b["boundary"], [[4, 5]])
        self.assertEqual(sum(
            y-x for groups in (a, b)
            for ranges in groups.values()
            for x, y in ranges
        ), 6)

    def test_missing_boundary_multipliers_rejected(self):
        self.kwargs["controls"] = {}
        self.kwargs["ports"].pop("B")
        with self.assertRaisesRegex(ValueError, "leave a gap"):
            self.fields()

    def test_role_conflict_rejected(self):
        self.kwargs["costates"]["q"] = object()
        with self.assertRaisesRegex(ValueError, "multiple pH roles"):
            self.fields()

    def test_owned_sparse_blocks_equal_full_dense_reference(self):
        fields = self.fields()
        a = np.array([
            [3., 2., 0., 0., 0., 1.],
            [1., 6., 0., 0., 0., 0.],
            [0., 0., 7., 2., 0., 0.],
            [0., 0., 2., 8., 0., 0.],
            [0., 0., 0., 0., 5., 1.],
            [0., 0., 0., 0., 1., 9.],
        ])
        sums = None
        for lo, hi in ((0, 3), (3, 6)):
            m = a[lo:hi]
            rows, cols = np.nonzero(m)
            # Build a normal owned-row CSR including explicit stored
            # structural zeros (PETSc can keep zero diagonal entries).
            rows_full = []
            cols_full = []
            vals_full = []
            ptr = [0]
            for r in range(m.shape[0]):
                for col in range(6):
                    if m[r, col] or col == r + lo:
                        rows_full.append(r)
                        cols_full.append(col)
                        vals_full.append(m[r, col])
                ptr.append(len(vals_full))
            s, e, f = PH.classify_owned_csr(
                fields, lo, hi,
                (ptr, cols_full, vals_full),
            )
            if sums is None:
                sums = [s, e, f]
            else:
                sums = [v+w for v, w in zip(sums, (s, e, f))]
        stored, effective, fro2 = sums
        self.assertEqual(int(effective.sum()), int(np.count_nonzero(a)))
        self.assertGreaterEqual(int(stored.sum()), int(effective.sum()))
        np.testing.assert_allclose(
            float(np.sqrt(fro2.sum())),
            float(np.linalg.norm(a)),
            atol=1e-12,
        )
        # Couplings cannot be silently reshuffled by field grouping.
        expected = np.zeros_like(fro2)
        for i in range(6):
            for j in range(6):
                expected[i, j] = a[i, j]**2
        np.testing.assert_allclose(fro2, expected)
        group = PH.aggregate_role_pairs(fields, effective)
        np.testing.assert_equal(
            int(group.sum()), int(effective.sum())
        )

    def test_global_columns_and_csr_layout_enforced(self):
        fields = self.fields()
        with self.assertRaisesRegex(ValueError, "GLOBAL PETSc"):
            PH.classify_owned_csr(
                fields, 0, 1,
                ([0, 1], [6], [1.])
            )
        with self.assertRaisesRegex(ValueError, "Invalid owned CSR"):
            PH.classify_owned_csr(
                fields, 0, 2,
                ([0, 1], [0], [1.])
            )


if __name__ == "__main__":
    unittest.main()
