# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Source-only guard for original Heat-Wave TS Python matched preconditioner."""

import ast
from pathlib import Path
import unittest

ROOT=Path(__file__).resolve().parents[1]
PC=ROOT/"sandbox/probe_v130_hw_matching_pc.py"
SNAP=ROOT/"sandbox/probe_v130_hw_matched_pc_snapshot.py"
CN=ROOT/"sandbox/probe_v130_hw_matched_pc_cn.py"
RUNNER=ROOT/"scripts/run-v130-hw-matched-PC-cn.sh"


class MatchedPCSourceContracts(unittest.TestCase):
    def test_real_petsc_python_pc_protocol_and_original_rhs(self):
        source=PC.read_text(encoding="utf-8")
        tree=ast.parse(source)
        cls=next(n for n in tree.body if isinstance(n,ast.ClassDef)
                 and n.name=="RowMatchedSuperLUPC")
        names={n.name for n in cls.body if isinstance(n,ast.FunctionDef)}
        self.assertGreaterEqual(len({"setUp","apply","reset","destroy"}&names),4)
        for token in (
            "original, preconditioner = pc.getOperators()",
            "row_match(preconditioner)",
            "petsc_from_permuted_coo(matched)",
            "values[self.rows]",
            "self.inner.solve(self.permuted_rhs,x)",
            "self.original.mult(x,self.check)",
            "self.check.axpy(-1.,b)",
            "relative > 1.e-9",
            'inner_pc.setFactorSolverType("superlu_dist")',
            'opts["mat_superlu_dist_replacetinypivot"] = "false"',
        ):
            self.assertIn(token,source)
        self.assertNotIn("A.shift(",source)

    def test_outer_operator_original_and_real_Ts_optin_only(self):
        s=SNAP.read_text(encoding="utf-8")
        self.assertIn("ksp.setOperators(A)",s)
        self.assertIn("outer_pc.setType(PETSc.PC.Type.PYTHON)",s)
        self.assertIn("outer_pc.setPythonContext(ctx)",s)
        t=CN.read_text(encoding="utf-8")
        self.assertIn("DPHS.IJacobian=attach_matched_PC",t)
        self.assertIn("result=original_ijacobian(",t)
        self.assertIn("pc.setPythonContext(ctx)",t)
        self.assertIn("DPHS.IJacobian=original_ijacobian",t)
        self.assertIn("original_cn.main()",t)

    def test_pinned_source_and_two_stage_native_gate(self):
        sh=RUNNER.read_text(encoding="utf-8")
        for token in (
            "check-v130-published-immutable.py",
            "first_cn_J_001.petscbin",
            "--candidate clamp",
            "PASS_SNAPSHOT_PYTHON_PC_ORIGINAL_EQUATIONS_ONLY",
            "real TS matched-PC",
            "MPLBACKEND=Agg",
        ):
            self.assertIn(token,sh)
        self.assertNotIn("LargeDiag_MC64",sh)
        self.assertNotIn("LargeDiag_AWPM",sh)


if __name__=="__main__":
    unittest.main()
