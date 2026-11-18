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


class _FakeCoState:
    def __init__(self, *, substituted=False):
        self.substituted = substituted

    def get_substituted(self):
        return self.substituted


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

    def test_native_audit_uses_user_owned_docker_workspace(self):
        script = (
            ROOT / "sandbox" / "run_ph_block_audit.sh"
        ).read_text(encoding="utf-8")
        self.assertIn(
            'SCRIMP_DOCKER_WORKSPACE:-$HOME/Simulations/workspace',
            script,
        )
        # Main runs, optional raw-matrix comparison, and final
        # block-parity comparison each run in the owned workspace.
        self.assertEqual(
            script.count("--workdir /scrimp-work"), 3
        )
        self.assertEqual(
            script.count('-v "$host_workspace:/scrimp-work"'), 3
        )
        self.assertIn('-e PYTHONPATH=/workspace/scrimp', script)
        self.assertIn('model="${2:-wave}"', script)
        self.assertIn('--model "$model"', script)
        self.assertIn(
            'python3 sandbox/diagnose_ph_block_parity.py "$out"',
            script,
        )

    def test_concentric_region_exclusion_is_sandbox_only(self):
        # The problematic Gmsh interface is preserved; the opt-in
        # diagnostic subtracts its region from the two volume regions.
        # No defaults or production geometry generator are modified.
        source = (
            ROOT / "sandbox" / "audit_ph_blocks.py"
        ).read_text(encoding="utf-8")
        runner = (
            ROOT / "sandbox" / "run_ph_block_audit.sh"
        ).read_text(encoding="utf-8")
        self.assertIn(
            "native_mesh.region_subtract(1, 10)", source
        )
        self.assertIn(
            "native_mesh.region_subtract(2, 10)", source
        )
        self.assertIn(
            'args.volume_only_regions', source
        )
        self.assertIn(
            'SCRIMP_PH_VOLUME_ONLY_REGIONS:-0', runner
        )
        self.assertIn(
            'capture_args+=(--volume-only-regions)', runner
        )
        self.assertIn(
            '--volume-only-regions only supports wave_coenergy',
            source,
        )

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

    def test_implicit_wave_selects_constitutive_schur_candidate(self):
        fields = self.fields()
        plan = PH.classify_ph_formulation(
            fields, costates=self.kwargs["costates"]
        )
        self.assertEqual(plan["formulation"], "implicit_constitutive")
        self.assertEqual(
            plan["solver_candidate"], "coenergy_schur_mumps"
        )

    def test_substituted_wave_has_no_independent_coenergy_block(self):
        # As in examples/wave_coenergy.py: q,p are the
        # coenergy variables; e_q/e_p are substituted *away*
        # and appear only through their constitutive mass bricks.
        self.model.intervals = {
            "q": (0, 2),
            "p": (2, 3),
            "e_r": (5, 1),
            "U": (6, 1),
            "Y": (7, 1),
        }
        self.model.nbdof = lambda: 8
        self.kwargs["costates"] = {
            "e_q": _FakeCoState(substituted=True),
            "e_p": _FakeCoState(substituted=True),
        }
        self.kwargs["ports"] = {
            "q": _FakePort("q", "e_q", False),
            "p": _FakePort("p", "e_p", False),
            "Damping": _FakePort("e_r", "e_r", True),
            "Boundary": _FakePort("U", "Y", True),
        }
        self.kwargs["controls"] = {
            "Boundary": _FakePort("U", "Y", True),
        }
        fields = self.fields()
        self.assertFalse(
            any(f["group"] == "coenergy" for f in fields)
        )
        self.assertEqual(
            {f["name"] for f in fields},
            {"q", "p", "e_r", "U", "Y"}
        )
        self.assertEqual(
            [f["group"] for f in fields],
            ["dynamic", "dynamic", "algebraic",
             "boundary", "boundary"],
        )
        plan = PH.classify_ph_formulation(
            fields, costates=self.kwargs["costates"]
        )
        self.assertEqual(plan["formulation"], "substituted_coenergy")
        self.assertEqual(
            plan["solver_candidate"], "mass_dirac_port_block"
        )

    def test_falsely_declared_substituted_costate_is_rejected(self):
        self.kwargs["costates"] = {
            "e_q": _FakeCoState(substituted=True),
            "e_p": _FakeCoState(substituted=True),
        }
        # The skipped substituted e_q/e_p intervals leave holes in
        # the remaining GetFEM vector; full DoF coverage catches this
        # inconsistent declaration without calling GetFEM C++ on the
        # nonexistent substituted coenergy name.
        with self.assertRaisesRegex(ValueError, "leave a gap"):
            self.fields()

    def test_substituted_names_are_never_probed_in_getfem(self):
        self.model.intervals = {
            "q": (0, 2), "p": (2, 3),
            "U": (5, 1), "Y": (6, 1),
        }
        self.model.nbdof = lambda: 7
        self.kwargs["costates"] = {
            "e_q": _FakeCoState(substituted=True),
            "e_p": _FakeCoState(substituted=True),
        }
        self.kwargs["ports"] = {
            "q": _FakePort("q", "e_q", False),
            "p": _FakePort("p", "e_p", False),
            "Boundary": _FakePort("U", "Y", True),
        }
        self.kwargs["controls"] = {
            "Boundary": _FakePort("U", "Y", True)
        }
        calls = []
        original = self.model.interval_of_variable

        def recorded(name):
            calls.append(name)
            if name in ("e_q", "e_p"):
                raise AssertionError(
                    "Substituted variable must not be queried in GetFEM"
                )
            return original(name)

        self.model.interval_of_variable = recorded
        fields = self.fields()
        self.assertNotIn("e_q", calls)
        self.assertNotIn("e_p", calls)
        self.assertEqual(sum(f["ndofs"] for f in fields), 7)

    def test_non_substituted_missing_coenergy_is_rejected(self):
        del self.model.intervals["e_q"]
        with self.assertRaisesRegex(
            ValueError, "Declared non-substituted"
        ):
            self.fields()

    def test_partially_substituted_formulation_is_not_assumed_explicit(self):
        del self.model.intervals["e_p"]
        self.model.nbdof = lambda: 5
        self.kwargs["costates"] = {
            "e_q": _FakeCoState(substituted=False),
            "e_p": _FakeCoState(substituted=True),
        }
        fields = self.fields()
        plan = PH.classify_ph_formulation(
            fields, costates=self.kwargs["costates"]
        )
        self.assertEqual(
            plan["formulation"], "partially_substituted"
        )
        self.assertEqual(
            plan["solver_candidate"], "mixed_block_decomposition"
        )

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
