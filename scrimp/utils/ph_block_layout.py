# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Read-only GetFEM/PETSc variable-block introspection for pH-DAEs.

GetFEM's unknown vector may interleave dynamic state, co-energy and
boundary variables in arbitrary contiguous variable ranges. This
module NEVER assumes z/e are in separate leading/trailing segments,
and never changes or condenses the physical finite-element system.

No PETSc/GetFEM dependency at import time. The native MPI validation
checks the actual assembled Jacobian and rank ownership separately.
"""

import numpy as np


GROUPS = ("dynamic", "coenergy", "boundary", "algebraic")


def discover_layout(model, *, states, costates, ports, controls):
    """Return complete disjoint GetFEM variable intervals with roles.

    A full-cover check is mandatory: silently omitting one boundary
    multiplier from PCFIELDSPLIT can make the resulting Schur operator
    singular or violate the port-Hamiltonian balance.
    """
    expected = int(model.nbdof())
    if expected <= 0:
        raise ValueError("GetFEM model has no unknowns")
    roles = {}

    def register(name, role):
        if not name:
            return
        if role not in GROUPS:
            raise ValueError(f"Unknown group {role!r}")
        old = roles.get(name)
        if old is not None and old != role:
            raise ValueError(
                f"Variable {name!r} occurs in multiple pH roles: "
                f"{old}, {role}"
            )
        roles[name] = role

    # Substituted co-states remain part of the *physical model* but
    # are NOT independent GetFEM variables. In wave_coenergy.py,
    # e_q and e_p are substituted into the dynamic mass blocks.
    # Do not infer a fake -Q_c block from their names.
    substituted = {
        name for name, obj in costates.items()
        if getattr(obj, "get_substituted", lambda: False)()
    }
    for port in ports.values():
        if getattr(port, "get_substituted", lambda: False)():
            effort = port.get_effort()
            if effort != port.get_flow():
                substituted.add(effort)

    for state in states:
        register(state, "dynamic")
    for costate in costates:
        register(costate, "coenergy")
    for control in controls.values():
        register(control.get_flow(), "boundary")
        register(control.get_effort(), "boundary")
    for port in ports.values():
        role = "algebraic" if port.get_algebraic() else "dynamic"
        register(port.get_flow(), roles.get(port.get_flow(), role))
        register(
            port.get_effort(),
            roles.get(port.get_effort(), "algebraic"),
        )

    fields = []
    for name, group in roles.items():
        # The GetFEM C++ binding prints "logic_error exception caught"
        # BEFORE Python can catch a missing interval. Substituted
        # coenergies have no independent DoFs by definition: never
        # query their nonexistent GetFEM variable intervals.
        if name in substituted:
            continue
        try:
            start, count = model.interval_of_variable(name)
        except (KeyError, RuntimeError, ValueError) as exc:
            # A genuinely substituted effort has no independent
            # GetFEM variable. Missing ANY non-substituted unknown
            # is an error, even if the remaining fields cover
            # a seemingly valid contiguous subset of the vector.
            if name in substituted:
                continue
            raise ValueError(
                f"Declared non-substituted variable {name!r} is "
                "missing from the GetFEM unknown vector"
            ) from exc
        if name in substituted:
            raise ValueError(
                f"Substituted effort {name!r} unexpectedly exists as "
                "an independent GetFEM unknown"
            )
        start, count = int(start), int(count)
        if start < 0 or count <= 0:
            raise ValueError(f"Invalid GetFEM range for {name}: {start}+{count}")
        fields.append({
            "name": name, "group": group,
            "start": start, "stop": start + count,
            "ndofs": count,
        })
    fields.sort(key=lambda f: f["start"])
    cursor = 0
    for field in fields:
        if field["start"] != cursor:
            raise ValueError(
                "GetFEM block ranges overlap or leave a gap at index "
                f"{cursor}: {field['name']} begins at {field['start']}. "
                "Do not construct a preconditioner until all boundary "
                "and constitutive unknowns are identified."
            )
        cursor = field["stop"]
    if cursor != expected:
        raise ValueError(
            f"GetFEM variables cover {cursor}/{expected} DoFs. "
            "A missing flow, coenergy or port constraint must be "
            "classified before attempting block preconditioning."
        )
    if not any(f["group"] == "dynamic" for f in fields):
        raise ValueError("Missing pH dynamic variable block")
    # No coenergy DoFs is entirely valid when all co-states are
    # substituted, e.g. the co-energy wave formulation.
    if (not any(f["group"] == "coenergy" for f in fields)
            and any(name not in substituted for name in costates)):
        raise ValueError(
            "Non-substituted co-states declared, but no independent "
            "coenergy DoFs were found"
        )
    return fields


def classify_ph_formulation(fields, *, costates):
    """Determine which candidate block strategy is structurally legal.

    Only a classification: NOT an authorization to enable any new
    linear solver. In particular, a substituted co-energy model
    has no independent constitutive -Q_c block, and an explicit-
    constitutive Schur solver must NEVER be attempted on it.
    """
    independent = {
        f["name"] for f in fields if f["group"] == "coenergy"
    }
    all_costates = set(costates)
    substituted = {
        name for name, costate in costates.items()
        if getattr(costate, "get_substituted", lambda: False)()
    }
    if not all_costates:
        return {
            "formulation": "other_or_no_costates",
            "solver_candidate": "global_mumps_reference",
            "reason": "No registered co-states: inspect full DAE manually",
        }
    if independent and not substituted:
        if independent != all_costates:
            raise ValueError(
                "Explicit constitutive partition omits some co-states"
            )
        return {
            "formulation": "implicit_constitutive",
            "solver_candidate": "coenergy_schur_mumps",
            "reason": "All coenergies exist as independent GetFEM unknowns",
        }
    if not independent and substituted == all_costates:
        return {
            "formulation": "substituted_coenergy",
            "solver_candidate": "mass_dirac_port_block",
            "reason": (
                "Constitutive laws are substituted into dynamic "
                "mass/stiffness; no independent -Q_c block exists"
            ),
        }
    if independent == all_costates - substituted:
        return {
            "formulation": "partially_substituted",
            "solver_candidate": "mixed_block_decomposition",
            "reason": (
                "Some constitutive relations are explicit and "
                "some already substituted; generic block audit needed"
            ),
        }
    raise ValueError("Registered co-state GetFEM blocks are inconsistent")


def owned_group_ranges(fields, lo, hi):
    """Sorted GLOBAL index intervals of each group owned by one MPI rank.

    The intervals can be used to construct local-global PETSc IS,
    including noncontiguous coenergy indices. They are never
    materialized as global vectors in the analysis itself.
    """
    lo, hi = int(lo), int(hi)
    if not (0 <= lo <= hi <= fields[-1]["stop"]):
        raise ValueError("Invalid PETSc ownership range")
    groups = {name: [] for name in GROUPS}
    for field in fields:
        first = max(lo, field["start"])
        last = min(hi, field["stop"])
        if first < last:
            groups[field["group"]].append([first, last])
    covered = sum(
        last - first
        for intervals in groups.values()
        for first, last in intervals
    )
    if covered != hi-lo:
        raise ValueError("PETSc rank ownership not covered by GetFEM fields")
    return groups


def classify_owned_csr(fields, lo, hi, csr, *, atol=0.0):
    """Sum per-variable and per-role sparse Frobenius norms from local rows.

    The three `csr` arrays must contain GLOBAL column indices
    (PETSc MatGetValuesCSR MPIAIJ). Includes only owned rows and
    performs no global gather, dense conversion, or matrix mutation.

    Returns numpy arrays that callers can MPI.Allreduce(SUM).
    Absolute-value zero tests distinguish allocated structural
    zeros from real nonzero coupling blocks.
    """
    ptr, cols, vals = (np.asarray(a) for a in csr)
    lo, hi = int(lo), int(hi)
    nf = len(fields)
    n = fields[-1]["stop"]
    nnz = len(vals)
    if (ptr.ndim != 1 or cols.ndim != 1 or vals.ndim != 1
            or ptr.size != hi-lo+1 or cols.size != nnz
            or ptr[0] != 0 or ptr[-1] != nnz
            or np.any(np.diff(ptr) < 0)):
        raise ValueError("Invalid owned CSR row pointer and value sizes")
    if nnz and (np.min(cols) < 0 or np.max(cols) >= n):
        raise ValueError(
            "CSR columns must use GLOBAL PETSc indices. "
            "Check MPIAIJ off-process column numbering."
        )
    if not np.isfinite(atol) or atol < 0:
        raise ValueError("Absolute sparsity tolerance must be finite >= 0")

    starts = np.array([f["start"] for f in fields], dtype=np.int64)
    stops = np.array([f["stop"] for f in fields], dtype=np.int64)
    rows = np.repeat(
        np.arange(lo, hi, dtype=np.int64),
        np.diff(ptr).astype(np.int64),
    )
    row_fields = np.searchsorted(starts, rows, side="right") - 1
    col_fields = np.searchsorted(starts, cols, side="right") - 1
    if nnz and (
        np.any(row_fields < 0) or np.any(col_fields < 0)
        or np.any(rows >= stops[row_fields])
        or np.any(cols >= stops[col_fields])
    ):
        raise ValueError("CSR row/column outside classified GetFEM interval")
    ids = row_fields * nf + col_fields
    is_effective = np.abs(vals) > atol
    if nnz and not np.all(np.isfinite(vals)):
        raise FloatingPointError("Nonfinite PETSc matrix entries")
    stored = np.bincount(ids, minlength=nf*nf).reshape(nf, nf)
    effective = np.bincount(
        ids[is_effective], minlength=nf*nf,
    ).reshape(nf, nf)
    fro2 = np.bincount(
        ids, weights=np.abs(vals)**2, minlength=nf*nf,
    ).reshape(nf, nf)
    return stored.astype(np.int64), effective.astype(np.int64), fro2


def aggregate_role_pairs(fields, values):
    """Collapse per-variable matrix to (dynamic/coenergy/boundary/etc.)."""
    v = np.asarray(values)
    nf = len(fields)
    if v.shape != (nf, nf):
        raise ValueError("Expected per-variable square matrix")
    lookup = {g: i for i, g in enumerate(GROUPS)}
    ids = np.array([lookup[f["group"]] for f in fields])
    return np.bincount(
        (ids[:, None]*len(GROUPS) + ids[None, :]).reshape(-1),
        weights=v.reshape(-1),
        minlength=len(GROUPS)**2,
    ).reshape(len(GROUPS), len(GROUPS))
