# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Controlled SuperLU_DIST pivot/order sweep on ONE immutable Heat-Wave matrix.

Real, unchanged original PETSc first CN IJacobian snapshot. Each run
is a new MPI1 process. Set only documented SuperLU_DIST row/column
permutation options; NEVER replace tiny pivots, perturb the matrix or
change residual tolerances. Compare against the known exact manufactured
solution and verify the PETSc source Mat CSR hash before/after.
"""

import argparse
import json
from pathlib import Path
import time

import numpy as np
from petsc4py import PETSc

from replay_base import matrix_digest, matrix_audit


EXPECTED_CSR_SHA = (
    "c94c69c4e4307d7a5aa61d9eba186782c0495b49431cadcb7281d7a7b5bffaff"
)
ROWS = ("default", "NOROWPERM", "LargeDiag_MC64", "LargeDiag_AWPM")
COLS = ("default", "MMD_AT_PLUS_A", "MMD_ATA", "METIS_AT_PLUS_A")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--matrix", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--rowperm", choices=ROWS, required=True)
    ap.add_argument("--colperm", choices=COLS, required=True)
    args = ap.parse_args()
    if (PETSc.COMM_WORLD.getSize() != 1 or not args.matrix.is_file()
            or args.output.exists() or not args.output.parent.is_dir()):
        raise RuntimeError("Require MPI1 existing immutable matrix and NEW output")

    opts = PETSc.Options()
    option_values = {}
    if args.rowperm != "default":
        option_values["mat_superlu_dist_rowperm"] = args.rowperm
    if args.colperm != "default":
        option_values["mat_superlu_dist_colperm"] = args.colperm
    # Do not replace tiny pivots, and never modify the matrix.
    option_values["mat_superlu_dist_replacetinypivot"] = "false"
    for key, value in option_values.items():
        opts[key] = value

    viewer = PETSc.Viewer().createBinary(
        str(args.matrix), mode="r", comm=PETSc.COMM_WORLD
    )
    A = PETSc.Mat().create(comm=PETSc.COMM_WORLD)
    try:
        A.setType("aij")
        A.load(viewer)
    finally:
        viewer.destroy()

    x_true = rhs = x = check = ksp = None
    try:
        expected = matrix_digest(A)
        if expected != EXPECTED_CSR_SHA:
            raise RuntimeError(
                "Changed input CN matrix, refusing permutation comparison"
            )
        audit = matrix_audit(A)
        if audit["size"] != 7689:
            raise RuntimeError("Wrong original Heat-Wave CN dimension")

        x_true = A.createVecRight()
        rhs = A.createVecLeft()
        x = A.createVecRight()
        check = A.createVecLeft()
        lo, hi = x_true.getOwnershipRange()
        x_true.setValues(
            np.arange(lo, hi, dtype=PETSc.IntType),
            1.+0.05*np.sin(np.arange(lo, hi, dtype=np.float64)),
            addv=PETSc.InsertMode.INSERT_VALUES,
        )
        x_true.assemble()
        A.mult(x_true, rhs)
        x.set(0.)
        rhs_norm = float(rhs.norm())

        ksp = PETSc.KSP().create(comm=PETSc.COMM_WORLD)
        ksp.setOperators(A)
        ksp.setType(PETSc.KSP.Type.PREONLY)
        pc = ksp.getPC()
        pc.setType(PETSc.PC.Type.LU)
        pc.setFactorSolverType("superlu_dist")
        # Actual library option path, not an isolated Python mock.
        ksp.setFromOptions()
        report = {
            "schema": "scrimp-v130-original-CN-superlu-dist-permutations-v1",
            "status": "DIAGNOSTIC_NOT_TS_PHYSICAL_PASS",
            "rowperm_requested": args.rowperm,
            "colperm_requested": args.colperm,
            "PETSc_options_explicitly_set": option_values,
            "input_CSR_SHA256": expected,
            "matrix_size": audit["size"],
            "matrix_structural_rank_including_stored_zeros":
                audit["structural_rank"],
            "zero_diagonal_entries": audit["numeric_zero_diagonal_entries"],
            "consistent_RHS_norm": rhs_norm,
            "tiny_pivot_replacement_disabled": True,
            "original_CN_TS_qualified": False,
            "physical_Hamiltonian_or_power_qualified": False,
        }

        start = time.monotonic()
        try:
            ksp.setUp()
            ksp.solve(rhs, x)
        except PETSc.Error as exc:
            report["PETSc_error"] = str(exc)
        report["wall_seconds"] = float(time.monotonic()-start)
        report["KSP_reason"] = int(ksp.getConvergedReason())
        report["KSP_iterations"] = int(ksp.getIterationNumber())
        report["PC_failed"] = str(pc.getFailedReason())

        # Obtain factor Mat's *actual* choices, if possible.
        factor_info = args.output.with_suffix(".factor-view.txt")
        try:
            F = pc.getFactorMatrix()
            ascii_viewer = PETSc.Viewer().createASCII(
                str(factor_info), comm=PETSc.COMM_WORLD
            )
            try:
                F.view(ascii_viewer)
            finally:
                ascii_viewer.destroy()
            view = factor_info.read_text(encoding="utf-8",errors="replace")
            report["factor_PETSc_view_available"] = True
            report["selected_rowperm_verified"] = (
                args.rowperm == "default"
                or f"Row permutation {args.rowperm}" in view
            )
            report["selected_colperm_verified"] = (
                args.colperm == "default"
                or f"Column permutation {args.colperm}" in view
            )
            report["factor_options_excerpt"] = [
                line.strip() for line in view.splitlines()
                if any(key in line for key in (
                    "Row permutation", "Column permutation",
                    "tiny pivots", "Equilibration"
                ))
            ]
        except Exception as exc:
            report["factor_PETSc_view_available"] = False
            report["factor_view_error"] = repr(exc)
            report["selected_rowperm_verified"] = False
            report["selected_colperm_verified"] = False

        if report["KSP_reason"] > 0 and "PETSc_error" not in report:
            A.mult(x, check)
            check.axpy(-1.,rhs)
            diff = x.copy()
            try:
                diff.axpy(-1.,x_true)
                report["relative_linear_residual"] = float(
                    check.norm()/max(rhs_norm,1e-30)
                )
                report["relative_state_error"] = float(
                    diff.norm()/max(float(x_true.norm()),1e-30)
                )
            finally:
                diff.destroy()
        else:
            report["relative_linear_residual"] = None
            report["relative_state_error"] = None

        report["output_CSR_SHA256"] = matrix_digest(A)
        if report["output_CSR_SHA256"] != expected:
            raise RuntimeError("SuperLU_DIST modified original PETSc Mat")
        good = (
            report["KSP_reason"] > 0
            and report["relative_linear_residual"] is not None
            and report["relative_linear_residual"] < 1e-9
            and report["relative_state_error"] < 1e-6
            and report["selected_rowperm_verified"]
            and report["selected_colperm_verified"]
        )
        report["status"] = (
            "PASS_IDENTICAL_CN_MATRIX_SOLVED_WITH_THIS_ORDERING"
            if good else "FAIL_LU_OR_ACCURACY_OR_OPTION_VERIFICATION"
        )
        args.output.write_text(json.dumps(report,indent=2)+"\n",
                               encoding="utf-8")
        print(
            f"{args.rowperm}/{args.colperm} {report['status']} "
            f"PC={report['PC_failed']} KSP={report['KSP_reason']} "
            f"residual={report['relative_linear_residual']} "
            f"error={report['relative_state_error']} "
            f"actual_factor_options={report['factor_options_excerpt'] if report['factor_PETSc_view_available'] else 'unavailable'}",
            flush=True,
        )
    finally:
        if ksp is not None:
            ksp.destroy()
        for obj in (check,x,rhs,x_true,A):
            if obj is not None:
                obj.destroy()


if __name__ == "__main__":
    main()
