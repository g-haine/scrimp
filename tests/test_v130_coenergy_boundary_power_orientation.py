# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
"""Regression of the REAL 2C Grad-Grad boundary-power orientation.

Recorded *original native FE* arrays, not manufactured model results:
  sandbox/outputs/v130-wave2d-coenergy-e1c5bc50-bRo68Qaf.zip,
  commit e1c5bc50, GetFEM/PETSc MUMPS, MPI1, dt=.008.
The original 2C pilot was correctly marked FAIL: U=-0.5Y
is ACTIVE for the Grad-Grad convention with Y=-p. The raw
GetFEM Hamiltonian and CN powers prove H-H0+Li-Lb=0,
not H-H0+Li+Lb=0, to machine precision.

This proves only the sign orientation in that historical run.
It MUST NOT count as physical acceptance of the repaired
U=+0.5Y native MPI1/2 run.
"""
import unittest
import numpy as np


# Native original data printed to 17 significant digits from rank0.npz.
H = np.array([
    0.25557110317334153, 0.25161684045952981,
    0.24777839262916021, 0.24407407658520644,
    0.24050864703669372, 0.2370750541067578,
    0.23376459829292304, 0.23057703848349603,
    0.22752342262351954, 0.22462238814172406,
    0.22189336900953158,
])
P_PORT = np.array([
    -0.0088144948076408942, -0.010403663710850459,
    -0.01444318636734819, -0.019386109771087928,
    -0.023654493182651937, -0.025904684875218208,
    -0.0260741589940265, -0.025372695141445974,
    -0.025369714915396826, -0.027175025169683183,
    -0.031049256744950689,
])
D_INNER = np.array([
    0.5095548240302572, 0.50468650293730644,
    0.49424916516360129, 0.48242561526524469,
    0.46933318674669883, 0.45510380111729781,
    0.43988113572334703, 0.42381767131985121,
    0.40707169741242311, 0.38980433539411991,
    0.37217664826897773,
])
TIME = np.arange(11)*.008


class OriginalGradGradPowerOrientationTests(unittest.TestCase):
    def test_real_fe_hamiltonian_has_negative_port_sign(self):
        dt=np.diff(TIME)
        # The ACTUAL native discrete pH balance is:
        # delta(H) = -dt*(D_inner+U*Y).
        err=np.diff(H)+dt*(D_INNER[1:]+P_PORT[1:])
        self.assertLess(np.max(np.abs(err)),1e-14)
        self.assertLess(np.max(P_PORT),0.)
        self.assertGreater(np.min(D_INNER),0.)

    def test_wrong_sign_explains_entire_old_gate_failure(self):
        dt=np.diff(TIME)
        Li=np.r_[0.,np.cumsum(dt*D_INNER[1:])]
        Lb_old=np.r_[0.,np.cumsum(-dt*P_PORT[1:])]
        expected=H-H[0]+Li-Lb_old
        wrong=H-H[0]+Li+Lb_old
        self.assertLess(np.max(np.abs(expected)),1e-14)
        self.assertTrue(np.allclose(wrong[1:],2*Lb_old[1:],
                                    atol=1e-14,rtol=0))
        self.assertGreater(np.max(np.abs(wrong)),.003)

    def test_canonical_passivity_law_is_positive_feedback_coefficient(self):
        from pathlib import Path
        root=Path(__file__).resolve().parents[1]
        code=(root/"examples"/"wave_coenergy.py").read_text(encoding="utf-8")
        self.assertIn('expressions = [f"{FEEDBACK_GAIN}*Y"]',code)
        self.assertNotIn('expressions = [f"-{FEEDBACK_GAIN}*Y"]',code)
        self.assertIn('step_boundary = dt*P[1:]',code)
        self.assertIn('"positive_nontrivial_feedback_loss"',code)
        self.assertIn('np.max(np.abs(feedback-P))',code)


if __name__=="__main__":
    unittest.main()
