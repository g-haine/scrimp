# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Source-only production policy and MPI shared VTK fail-closed tests.

Do not import native GetFEM/PETSc. The real root-only replica
validator is extracted from the DPHS source with Python AST, and
executed against deterministic minimal model and MPI fakes.
Only genuine pinned native Docker runs qualify the actual export.
"""

import ast
import importlib.util
from pathlib import Path
import types
import unittest

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "scrimp" / "utils" / "production_policy.py"
SPEC = importlib.util.spec_from_file_location("v130_production_policy", SOURCE)
POLICY = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(POLICY)


def load_live_replica_validator():
    """Compile the actual method, avoiding a mock of its implementation."""
    tree = ast.parse((ROOT / "scrimp" / "dphs.py").read_text())
    cls = next(node for node in tree.body if isinstance(node, ast.ClassDef)
               and node.name == "DPHS")
    method = next(node for node in cls.body if isinstance(node, ast.FunctionDef)
                  and node.name == "_validate_root_only_pv_replication")
    module = ast.fix_missing_locations(ast.Module(body=[method], type_ignores=[]))
    return compile(module, "scrimp/dphs.py (real replica validator)", "exec")


class Model:
    def nbdof(self):
        return 4

    def mesh_fem_of_variable(self, name):
        if name != "p":
            raise KeyError(name)
        return types.SimpleNamespace(nbdof=lambda: 4)


class Mesh:
    def pts(self):
        return np.array([[0., 1., 2.], [0., 0., 0.]])

    def nbcvs(self):
        return 2

    def cvid(self):
        return np.array([0, 1], dtype=np.int64)

    def pid_from_cvid(self, ids):
        assert list(ids) == [0, 1]
        return np.array([0, 1, 1, 2]), np.array([0, 2, 4])


class Comm:
    def __init__(self, different=False):
        self.different = different
        self.seen = 0

    def tompi4py(self):
        return self

    def allgather(self, data):
        self.seen += 1
        if self.different and data[1] is None:
            broken = list(data[0])
            broken[-1] = "different rank trajectory"
            return [data, (tuple(broken), None)]
        return [data, data]


def fake_system(*, bad_size=False, nan=False):
    z = np.array([1., 2., 3., 4.])
    if bad_size:
        z = z[:3]
    if nan:
        z[1] = np.nan
    return types.SimpleNamespace(
        solution={"t": [0., .001], "z": [z.copy(), z.copy()]},
        gf_model=Model(),
        ports={"p": types.SimpleNamespace(get_mesh_id=lambda: 0)},
        domain=types.SimpleNamespace(get_mesh=lambda: [Mesh()]),
    )


class ProductionSwitchTests(unittest.TestCase):
    def test_conforming_compiler_is_enabled_by_default(self):
        self.assertTrue(POLICY.auto_interface_enabled({}))

    def test_conforming_compiler_can_be_disabled_and_reenabled(self):
        self.assertFalse(POLICY.auto_interface_enabled({
            "SCRIMP_AUTO_INTERFACE": "0",
        }))
        self.assertTrue(POLICY.auto_interface_enabled({
            "SCRIMP_AUTO_INTERFACE": "1",
        }))

    def test_old_compiler_pilot_flag_still_works(self):
        self.assertFalse(POLICY.auto_interface_enabled({
            "SCRIMP_V130_AUTO_INTERFACE_PILOT": "0",
        }))
        self.assertTrue(POLICY.auto_interface_enabled({
            "SCRIMP_V130_AUTO_INTERFACE_PILOT": "1",
        }))

    def test_production_switch_overrides_legacy_pilot(self):
        self.assertFalse(POLICY.auto_interface_enabled({
            "SCRIMP_AUTO_INTERFACE": "0",
            "SCRIMP_V130_AUTO_INTERFACE_PILOT": "1",
        }))

    def test_reject_invalid_auto_interface_flag(self):
        with self.assertRaisesRegex(ValueError, "SCRIMP_AUTO_INTERFACE"):
            POLICY.auto_interface_enabled({"SCRIMP_AUTO_INTERFACE": "maybe"})

    def test_single_rank_preserves_original_export_policy(self):
        self.assertEqual(POLICY.pv_export_mode({}, 1), "legacy")

    def test_mpi_two_defaults_to_one_export_writer(self):
        self.assertEqual(POLICY.pv_export_mode({}, 2), "root")

    def test_legacy_writer_optout(self):
        self.assertEqual(POLICY.pv_export_mode({
            "SCRIMP_PV_EXPORT_MODE": "legacy",
        }, 2), "legacy")

    def test_explicit_root_selected_in_serial_and_mpi(self):
        for n in (1, 2):
            self.assertEqual(POLICY.pv_export_mode({
                "SCRIMP_PV_EXPORT_MODE": "root",
            }, n), "root")

    def test_old_root_only_pilot_alias(self):
        self.assertEqual(POLICY.pv_export_mode({
            "SCRIMP_V130_PV_ROOT_ONLY_PILOT": "1",
        }, 2), "root")
        self.assertEqual(POLICY.pv_export_mode({
            "SCRIMP_V130_PV_ROOT_ONLY_PILOT": "0",
        }, 2), "legacy")

    def test_invalid_export_mode_rejected(self):
        with self.assertRaisesRegex(ValueError, "SCRIMP_PV_EXPORT_MODE"):
            POLICY.pv_export_mode({
                "SCRIMP_PV_EXPORT_MODE": "distributed_guess",
            }, 2)


class RootOnlyReplicatedVTKTests(unittest.TestCase):
    def setUp(self):
        self.comm = Comm()
        namespace = {"np": np, "comm": self.comm}
        exec(load_live_replica_validator(), namespace)
        self.validate = namespace["_validate_root_only_pv_replication"]

    def test_real_validator_accepts_identical_full_states_and_mesh(self):
        self.validate(fake_system(), "p", "All")
        self.assertEqual(self.comm.seen, 1)

    def test_real_validator_rejects_incomplete_global_state(self):
        with self.assertRaisesRegex(RuntimeError, "FE snapshot"):
            self.validate(fake_system(bad_size=True), "p", "All")
        self.assertEqual(self.comm.seen, 1)

    def test_real_validator_rejects_nonfinite_state(self):
        with self.assertRaisesRegex(RuntimeError, "Nonfinite"):
            self.validate(fake_system(nan=True), "p", "Final")
        self.assertEqual(self.comm.seen, 1)

    def test_real_validator_rejects_rank_disagreement(self):
        self.comm.different = True
        with self.assertRaisesRegex(RuntimeError, "identical complete mesh"):
            self.validate(fake_system(), "p", "All")
        self.assertEqual(self.comm.seen, 1)

    def test_real_validator_rejects_missing_export_field(self):
        with self.assertRaisesRegex(RuntimeError, "Unknown export variable"):
            self.validate(fake_system(), "missing", "All")

    def test_real_validator_rejects_invalid_time_selection(self):
        with self.assertRaisesRegex(RuntimeError, "Unknown ParaView time"):
            self.validate(fake_system(), "p", "Unsupported")


if __name__ == "__main__":
    unittest.main()
