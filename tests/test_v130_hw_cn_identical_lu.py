# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Source-only guard on original identical-CN-Jacobian LU replay."""

from pathlib import Path
import unittest

ROOT=Path(__file__).resolve().parents[1]
REPLAY=ROOT/"sandbox/probe_v130_hw_cn_identical_lu.py"
RUNNER=ROOT/"scripts/run-v130-hw-cn-identical-lu.sh"

class HeatWaveIdenticalLUReplayTests(unittest.TestCase):
    def test_capture_uses_unchanged_production_first_ijacobian_hook(self):
        d=(ROOT/"scrimp/dphs.py").read_text(encoding="utf-8")
        assert 'SCRIMP_JACOBIAN_SNAPSHOT_PREFIX' in d
        assert 'P.view(viewer)' in d
        shell=RUNNER.read_text(encoding="utf-8")
        for literal in (
            "SCRIMP_JACOBIAN_SNAPSHOT_PREFIX=/artifacts/first_cn_J",
            "first_cn_J_001.petscbin",
            "--backend superlu_dist",
            "--capture-ts-failure",
            "SCRIMP_GETFEM_BRIDGE_DISTRIBUTED=coo",
            "check-v130-published-immutable.py",
        ):
            self.assertIn(literal,shell)

    def test_both_factorizers_receive_same_binary_loaded_mat(self):
        source=REPLAY.read_text(encoding="utf-8")
        for literal in (
            'for solver in ("mumps","superlu_dist")',
            'A.load(viewer)',
            'ksp.setOperators(A)',
            'pc.setFactorSolverType(solver)',
            'matrix_digest(A)',
            'structural_rank(graph)',
            'matrix_audit(A)',
            'no_diagonal_shift_or_pivot_replacement',
            'full_Crank_Nicolson_TS_completed',
        ):
            self.assertIn(literal,source)
        self.assertNotIn('A.shift(',source)
        self.assertNotIn('setPreallocationCOO',source)
        self.assertNotIn('from scrimp',source)

    def test_source_immutability(self):
        published=(ROOT/"scripts/check-v130-published-immutable.py").read_text()
        self.assertIn("EXPECTED_GIT_BLOBS",published)


if __name__=="__main__":
    unittest.main()
