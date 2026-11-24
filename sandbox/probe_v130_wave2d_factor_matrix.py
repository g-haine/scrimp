# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# See LICENSE for license information.

"""Capture the EXACT first-main-step PETSc factor matrix of wave 2D.

This is an opt-in diagnostic wrapper: it does not modify the canonical
examples/wave.py, SCRIMP's equations, the Jacobian, TS/SNES/KSP,
LU configuration, time step or DAE initialization.

The original DPHS.IJacobian callback first assembles the actual operator.
Immediately afterwards, BEFORE PC factorization, take a byte-for-byte
copy of the owned AIJ CSR buffers for the PC matrix P and its M and K
components. Record original GetFEM variable intervals for the ten zero
diagonal positions, while keeping the physical solver untouched.

This script intentionally supports only a single MPI rank. Multi-rank
matrix snapshots require an explicitly collective global-CSR format,
not a dishonest concatenation of local index spaces.
"""

import hashlib
import json
import os
from pathlib import Path
import traceback

import numpy as np
import scrimp as S
from mpi4py import MPI


def capture_first_main_jacobian():
    if MPI.COMM_WORLD.Get_size() != 1:
        raise RuntimeError("Exact global CSR capture is serial-only")
    location = os.environ.get("SCRIMP_V130_WAVE2D_CAPTURE_DIR", "")
    if not location:
        raise RuntimeError("Set SCRIMP_V130_WAVE2D_CAPTURE_DIR to /artifacts")
    target = Path(location)
    target.mkdir(parents=True, exist_ok=True)
    original = S.DPHS.IJacobian

    def wrapped(system, TS, t, z, zd, shift, A, P):
        result = original(system, TS, t, z, zd, shift, A, P)
        if (getattr(system, "_profile_phase", "") != "main"
                or getattr(system, "_v130_wave2d_captured", False)):
            return result
        if not (np.isfinite(t) and np.isfinite(shift)
                and abs(float(t)-.008) < 1.e-9
                and abs(float(shift)-250.) < 1.e-7):
            raise RuntimeError(
                "Unexpected first CN Jacobian time/shift "
                f"t={t!r}, shift={shift!r}; refusing unlabelled capture"
            )
        system._v130_wave2d_captured = True

        mats = {
            "P": P,  # EXACT PC input that PETSc will factorize
            "M": system.tangent_mass,
            "K": system.tangent_stiffness,
        }
        arrays = {}
        shapes = {}
        zero_diagnostics = {}
        for label, matrix in mats.items():
            shape = tuple(int(x) for x in matrix.getSize())
            start, stop = map(int, matrix.getOwnershipRange())
            if shape != (2230, 2230) or (start, stop) != (0, 2230):
                raise RuntimeError(
                    f"Cannot interpret non-global owned matrix {label}: "
                    f"shape={shape} ownership={(start, stop)}"
                )
            indptr, indices, values = matrix.getValuesCSR()
            indptr = np.asarray(indptr).copy()
            indices = np.asarray(indices).copy()
            values = np.asarray(values).copy()
            if (indptr.shape != (shape[0]+1,)
                    or int(indptr[0]) != 0
                    or int(indptr[-1]) != len(values)
                    or len(indices) != len(values)
                    or not np.all(np.isfinite(values))):
                raise RuntimeError(f"Malformed/nonfinite native CSR {label}")
            arrays[f"{label}_indptr"] = indptr
            arrays[f"{label}_indices"] = indices
            arrays[f"{label}_data"] = values
            shapes[label] = list(shape)
            diag = np.zeros(shape[0], dtype=float)
            structural_present = np.zeros(shape[0], dtype=bool)
            zero_rows = []
            for row in range(shape[0]):
                a, b = int(indptr[row]), int(indptr[row + 1])
                rowcol, rowval = indices[a:b], values[a:b]
                diagonal = (rowcol == row)
                if np.any(diagonal):
                    diag[row] = float(np.sum(rowval[diagonal]))
                    structural_present[row] = True
                if not np.any(rowval):
                    zero_rows.append(row)
            zero_diagnostics[label] = {
                "stored_nnz": int(values.size),
                "numerical_nnz": int(np.count_nonzero(values)),
                "numerical_zero_diagonal_ids":
                    np.flatnonzero(diag == 0).tolist(),
                "structurally_missing_diagonal_ids":
                    np.flatnonzero(~structural_present).tolist(),
                "numerical_zero_row_ids": zero_rows,
                "min_nonzero_abs_entry": float(
                    np.min(np.abs(values[values != 0]))
                ) if np.any(values != 0) else None,
            }
        variables = []
        seen = set()
        for port in system.ports.values():
            for variable in (port.get_flow(), port.get_effort()):
                if not variable or variable in seen:
                    continue
                seen.add(variable)
                try:
                    start, width = map(
                        int, system.gf_model.interval_of_variable(variable)
                    )
                except (KeyError, ValueError, RuntimeError):
                    continue
                variables.append(dict(
                    name=variable, lo=start,
                    hi=start+width, dofs=width,
                ))
        variables.sort(key=lambda d: d["lo"])
        digest = hashlib.sha256()
        for key in sorted(arrays):
            digest.update(key.encode())
            digest.update(np.ascontiguousarray(arrays[key]).tobytes())
        report = {
            "schema": "scrimp-v130-first-main-CN-exact-PC-factor-matrix-v1",
            "qualification": "DIAGNOSTIC_ONLY_NO_BACKEND_SUCCESS_CLAIM",
            "time": float(t),
            "TS_shift": float(shift),
            "source_model": "examples.wave:wave_eq",
            "selected_solver": str(
                system.time_scheme.getString(
                    "pc_factor_mat_solver_type", ""
                )
            ).lower(),
            "matrix_shapes": shapes,
            "csr_SHA256": digest.hexdigest(),
            "matrix_statistics": zero_diagnostics,
            "getfem_variable_intervals": variables,
        }
        np.savez_compressed(target/"first_main_CN_factor.npz", **arrays)
        (target/"first_main_CN_factor.json").write_text(
            json.dumps(report, indent=2)+"\n", encoding="utf-8"
        )
        print(
            "CAPTURED_FIRST_REAL_CN_PC_MATRIX "
            f"solver={report['selected_solver']} "
            f"shape=2230x2230 shift={float(shift):g} "
            "original_P_M_K_csr=yes", flush=True
        )
        return result

    S.DPHS.IJacobian = wrapped
    try:
        from examples.wave import main as original_main
        original_main()
    finally:
        S.DPHS.IJacobian = original


if __name__ == "__main__":
    try:
        capture_first_main_jacobian()
    except Exception:
        # Preserve exactly the original native solver's exit and
        # trace; the accompanying script archives this log.
        traceback.print_exc()
        raise
