# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Focused, no-linear-solve checks for MUMPS reset-policy selection/CSR.

Run in the same environment as SCRIMP:
    python -u sandbox/check_mumps_reset_policy_unit.py

This test checks the policy default, compatibility with the legacy
boolean flag, and exact stored local CSR detection. It does not prove
that PETSc/MUMPS factorization is safe; the real nonlinear MPI
validation is documented in SYMBOLIC_REUSE_VALIDATION.md.
"""

import os
import numpy as np
from scrimp.dphs import (
    DPHS,
    _all_bricks_are_linear,
    _get_mumps_reset_policy,
    _local_stored_csr_pattern,
    _local_csr_pattern_changed,
)


class FakeBrick:
    def __init__(self, linear, position="effort", dt=False):
        self.linear = linear
        self.position = position
        self.dt = dt

    def get_linear(self):
        return self.linear

    def get_position(self):
        return self.position

    def get_dt(self):
        return self.dt


def brick_linearity_checks():
    assert not _all_bricks_are_linear({})
    assert _all_bricks_are_linear({"M": FakeBrick(True)})
    assert _all_bricks_are_linear({
        "mass": FakeBrick(True),
        "source": FakeBrick(False, position="source"),
        "boundary": FakeBrick(True),
    }), "Prescribed control source is RHS-only, not a nonlinear tangent"
    assert not _all_bricks_are_linear({
        "mass": FakeBrick(True),
        "boundary": FakeBrick(False),
        "source": FakeBrick(False, position="source"),
    })
    assert not _all_bricks_are_linear({
        "mass": FakeBrick(True),
        "nonlinear_explicit": FakeBrick(
            False, position="effort",
        ),
    })

    # Isolate set_linear_flags without constructing a GetFEM domain.
    system = DPHS.__new__(DPHS)
    system.bricks = {
        "mass": FakeBrick(True, dt=True),
        "control_mass": FakeBrick(True),
        "prescribed_control": FakeBrick(False, position="source"),
    }
    system.set_linear_flags()
    assert system.linear_mass is True
    assert system.linear_stiffness is True
    system.bricks["stiffness_nl"] = FakeBrick(False)
    system.set_linear_flags()
    assert system.linear_mass is True
    assert system.linear_stiffness is False
    system.bricks["mass_nl"] = FakeBrick(False, dt=True)
    system.set_linear_flags()
    assert system.linear_mass is False
    assert system.linear_stiffness is False


class Matrix:
    def __init__(self, rowptr, columns, ownership=(0, 2), values=None):
        self.rowptr = np.array(rowptr, dtype=np.int32)
        self.columns = np.array(columns, dtype=np.int32)
        self.ownership = ownership
        self.values = np.asarray(
            np.ones(len(columns)) if values is None else values
        )

    def getValuesCSR(self):
        return self.rowptr, self.columns, self.values

    def getOwnershipRange(self):
        return self.ownership


def policy_checks():
    envvars = (
        "SCRIMP_MUMPS_RESET_POLICY",
        "SCRIMP_MUMPS_RESET_PC_ON_JACOBIAN",
    )
    saved = {key: os.environ.get(key) for key in envvars}
    try:
        for key in envvars:
            os.environ.pop(key, None)
        assert _get_mumps_reset_policy() == "on_change"
        assert _get_mumps_reset_policy(all_bricks_linear=False) == "on_change"
        assert _get_mumps_reset_policy(all_bricks_linear=True) == "never"
        os.environ["SCRIMP_MUMPS_RESET_PC_ON_JACOBIAN"] = "1"
        assert _get_mumps_reset_policy() == "always"
        assert _get_mumps_reset_policy(all_bricks_linear=True) == "always"
        os.environ["SCRIMP_MUMPS_RESET_PC_ON_JACOBIAN"] = "0"
        assert _get_mumps_reset_policy() == "never"
        assert _get_mumps_reset_policy(all_bricks_linear=False) == "never"
        os.environ["SCRIMP_MUMPS_RESET_POLICY"] = "on_change"
        assert _get_mumps_reset_policy() == "on_change"
        assert _get_mumps_reset_policy(all_bricks_linear=True) == "on_change"
        os.environ["SCRIMP_MUMPS_RESET_POLICY"] = "always"
        assert _get_mumps_reset_policy() == "always"
        os.environ["SCRIMP_MUMPS_RESET_POLICY"] = "never"
        assert _get_mumps_reset_policy() == "never"
        os.environ["SCRIMP_MUMPS_RESET_POLICY"] = "invalid"
        try:
            _get_mumps_reset_policy()
        except ValueError:
            pass
        else:
            raise AssertionError("Invalid explicit policy was accepted")
        os.environ.pop("SCRIMP_MUMPS_RESET_POLICY")
        os.environ["SCRIMP_MUMPS_RESET_PC_ON_JACOBIAN"] = "invalid"
        try:
            _get_mumps_reset_policy()
        except ValueError:
            pass
        else:
            raise AssertionError("Invalid legacy policy was accepted")
    finally:
        for key, value in saved.items():
            if value is None:
                os.environ.pop(key, None)
            else:
                os.environ[key] = value


def csr_checks():
    mat = Matrix([0, 2, 4], [0, 1, 0, 1])
    original = _local_stored_csr_pattern(mat)
    mat.values *= 10
    assert not _local_csr_pattern_changed(
        original, _local_stored_csr_pattern(mat)
    ), "A numerical-only update must not require new symbolic LU"
    # Ensure the returned snapshot is detached from PETSc CSR storage.
    mat.columns[2] = 2
    assert original[2][2] == 0
    updated = _local_stored_csr_pattern(mat)
    assert _local_csr_pattern_changed(original, updated)
    mat.columns[2] = 0
    assert not _local_csr_pattern_changed(
        original, _local_stored_csr_pattern(mat)
    )
    assert _local_csr_pattern_changed(
        original, _local_stored_csr_pattern(
            Matrix([0, 1, 4], [0, 1, 0, 1])
        )
    )
    assert _local_csr_pattern_changed(
        original, _local_stored_csr_pattern(
            Matrix([0, 2, 4], [0, 1, 0, 1], ownership=(2, 4))
        )
    )


if __name__ == "__main__":
    brick_linearity_checks()
    policy_checks()
    csr_checks()
    print(
        "PASS: MUMPS automatic never for all-linear bricks, "
        "on_change for nonlinear bricks, explicit overrides and exact CSR"
    )
