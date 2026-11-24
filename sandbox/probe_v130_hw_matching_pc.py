# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""EXPERIMENTAL serial PETSc Python PC: (P A)^(-1) P r.

For the unchanged original Heat-Wave DAE. The outer KSP/TS operates on
its ORIGINAL operator and original RHS. The inner PETSc KSP factors a
left-row-matched matrix computed from the actual preconditioning
operator at each PC setup, and *applies the same row map to every RHS*.

This is a sandbox-only, serial feasibility experiment. Never install
as an automatic universal solver fallback without MPI, dynamic-graph,
rank, lifetime and full trajectory qualification.
"""

import hashlib

import numpy as np
from petsc4py import PETSc

from probe_v130_hw_cn_exact_row_matching import (
    row_match, petsc_from_permuted_coo,
)
from replay_base import matrix_digest


class RowMatchedSuperLUPC:
    """PCPython protocol implementing original-order inverse exactly."""

    def __init__(self):
        self.matrix = None
        self.inner = None
        self.permuted_rhs = None
        self.check = None
        self.original = None
        self.rows = None
        self.setup_history = []
        self.apply_count = 0
        self.max_original_relative_residual = 0.0

    def reset(self, pc):
        for item in (
            self.check, self.permuted_rhs, self.inner, self.matrix
        ):
            if item is not None:
                item.destroy()
        self.check = self.permuted_rhs = self.inner = self.matrix = None
        self.original = self.rows = None

    def destroy(self, pc):
        self.reset(pc)

    def setUp(self, pc):
        if PETSc.COMM_WORLD.getSize() != 1:
            raise RuntimeError("Row-matched Python PC is MPI1-only research")
        self.reset(pc)
        original, preconditioner = pc.getOperators()
        if original.getSize() != (7689, 7689):
            raise RuntimeError("Wrong original Heat-Wave CN state dimension")
        if preconditioner.getSize() != original.getSize():
            raise RuntimeError("Incompatible PETSc KSP operator sizes")
        self.original = preconditioner
        before = matrix_digest(preconditioner)
        source, matched, rows, meta = row_match(preconditioner)
        if meta["matched_numeric_zero_diagonal"] != 0:
            raise RuntimeError("Matching failed to remove numerical diag zeros")
        self.matrix = petsc_from_permuted_coo(matched)
        self.rows = rows
        self.permuted_rhs = self.matrix.createVecLeft()
        self.check = preconditioner.createVecLeft()
        inner = PETSc.KSP().create(comm=PETSc.COMM_WORLD)
        self.inner = inner
        inner.setOptionsPrefix("scrimp_hw_rowmatched_inner_")
        inner.setOperators(self.matrix)
        inner.setType(PETSc.KSP.Type.PREONLY)
        inner_pc = inner.getPC()
        inner_pc.setType(PETSc.PC.Type.LU)
        inner_pc.setFactorSolverType("superlu_dist")
        # Same qualified settings as successful direct row-matched LU.
        opts = PETSc.Options()
        opts["mat_superlu_dist_rowperm"] = "NOROWPERM"
        opts["mat_superlu_dist_colperm"] = "MMD_AT_PLUS_A"
        opts["mat_superlu_dist_replacetinypivot"] = "false"
        inner.setFromOptions()
        inner.setUp()
        if str(inner_pc.getFailedReason()) not in ("0", "NONE"):
            raise RuntimeError(
                "Row-matched SuperLU_DIST numeric factor setup failed "
                f"{inner_pc.getFailedReason()}"
            )
        after = matrix_digest(preconditioner)
        if before != after:
            raise RuntimeError("PC setup changed the original PETSc Jacobian")
        self.setup_history.append({
            "source_matrix_CSR_SHA256": before,
            "row_matching_SHA256": meta["row_permutation_sha256"],
            "rows_moved": meta["rows_moved"],
            "source_stored_nnz": meta["original_stored_nonzeros"],
            "matched_zero_diagonal": meta["matched_numeric_zero_diagonal"],
            "min_matched_abs_diagonal":
                meta["matched_smallest_absolute_diagonal"],
            "inner_factor_backend":inner_pc.getFactorSolverType(),
            "inner_pc_failure":str(inner_pc.getFailedReason()),
        })
        print(
            "MATCHED-PC SETUP "
            f"count={len(self.setup_history)} "
            f"moved={meta['rows_moved']} "
            f"min_diag={meta['matched_smallest_absolute_diagonal']:.8g} "
            "factor=superlu_dist", flush=True
        )

    def apply(self, pc, b, x):
        if self.inner is None or self.rows is None:
            raise RuntimeError("Matched PC apply without successful setUp")
        values = np.asarray(b.getArray(readonly=True)).copy()
        if len(values) != len(self.rows):
            raise RuntimeError("PETSc PC apply ownership wrong (MPI1 only)")
        self.permuted_rhs.setValues(
            np.arange(len(values),dtype=PETSc.IntType),
            values[self.rows],
            addv=PETSc.InsertMode.INSERT_VALUES,
        )
        self.permuted_rhs.assemble()
        self.inner.solve(self.permuted_rhs,x)
        if int(self.inner.getConvergedReason()) <= 0:
            raise RuntimeError(
                "Inner matched LU solve failed: "
                f"{self.inner.getConvergedReason()}"
            )
        # Validate the returned vector in the ORIGINAL row order.
        self.original.mult(x,self.check)
        self.check.axpy(-1.,b)
        relative=float(
            self.check.norm()/max(float(b.norm()),1e-30)
        )
        if not np.isfinite(relative) or relative > 1.e-9:
            raise RuntimeError(
                "Matched PC returned a wrong ORIGINAL-equation solution "
                f"{relative:.12g}"
            )
        self.apply_count+=1
        self.max_original_relative_residual=max(
            self.max_original_relative_residual,relative
        )
        if self.apply_count<=5:
            print(
                f"MATCHED-PC APPLY {self.apply_count} "
                f"original_rel_residual={relative:.12g}",flush=True
            )

    def audit(self):
        return {
            "setup_count":len(self.setup_history),
            "apply_count":self.apply_count,
            "max_relative_residual_original_preconditioner":
                self.max_original_relative_residual,
            "setup_history":self.setup_history,
            "permutation_applied_to_matrix_and_each_RHS":True,
            "unknown_vector_unpermuted":True,
            "no_pivot_replacement_or_diagonal_regularization":True,
            "MPI2_qualified":False,
        }


if __name__=="__main__":
    raise SystemExit(
        "Use sandbox/probe_v130_hw_matched_pc_snapshot.py or "
        "sandbox/probe_v130_hw_matched_pc_cn.py"
    )
