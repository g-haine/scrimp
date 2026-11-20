# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Pure test for native mass-Jacobian gate; not a numerical replacement."""

import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT=Path(__file__).resolve().parents[1]
COMPARE=ROOT/"sandbox/compare_v130_shallow_water_mass_jacobian.py"
PROBE=ROOT/"sandbox/check_v130_shallow_water_mass_jacobian.py"
CORE=ROOT/"scrimp/dphs.py"


class TestMassJacobianQualification(unittest.TestCase):
    def setUp(self):
        tmp=tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.out=Path(tmp.name)
        self._write()

    def _write(self, *, wrong_default=False, missing_mass=False):
        for ranks in (1,2):
            for mode in ("default","legacy"):
                base=1e-9 if mode=="default" else 1e-2
                results=[
                    {"shift":s,"zd_norm":zd,"relative_error":base}
                    for zd,s in ((0.,0.),(0.,1.),(1.,0.),(1.,1.))
                ]
                if mode=="legacy":
                    results[0]["relative_error"]=1e-9
                    results[1]["relative_error"]=1e-9
                report={
                    "mode":mode, "ranks":ranks,
                    "nonlinear_mass":not missing_mass,
                    "main_TS_executed":False,
                    "independent_time_derivative":True,
                    "dofs":80,
                    "mass_mode_env":(
                        "legacy" if mode=="legacy" else
                        ("finite_difference" if wrong_default else None)
                    ),
                    "outcomes":results,
                    "maximum_directional_relative_error":max(
                        x["relative_error"] for x in results
                    ),
                }
                (self.out/f"shallow_{mode}_{ranks}r.json").write_text(
                    json.dumps(report),encoding="utf-8",
                )

    def _run(self):
        return subprocess.run(
            [sys.executable,str(COMPARE),str(self.out)],
            capture_output=True,text=True,
        )

    def test_qualifier_accepts_default_and_worse_legacy_fixture(self):
        x=self._run()
        self.assertEqual(x.returncode,0,x.stderr)
        data=json.loads((
            self.out/"shallow_mass_jacobian_qualification.json"
        ).read_text(encoding="utf-8"))
        self.assertEqual(data["status"],"PASS")
        self.assertFalse(data["nonlinear_DAE_initializer_qualified"])

    def test_rejects_explicit_correction_disguised_as_default(self):
        self._write(wrong_default=True)
        x=self._run()
        self.assertNotEqual(x.returncode,0)
        self.assertIn("explicitly configured",x.stderr)

    def test_rejects_missing_nonlinear_mass(self):
        self._write(missing_mass=True)
        x=self._run()
        self.assertNotEqual(x.returncode,0)
        self.assertIn("Wrong native Jacobian",x.stderr)

    def test_production_default_is_finite_difference(self):
        s=CORE.read_text(encoding="utf-8")
        self.assertIn(
            '"SCRIMP_NONLINEAR_MASS_JACOBIAN", "finite_difference"',
            s,
        )
        self.assertIn(
            "self._add_nonlinear_mass_hessian_fd(t, z, zd, target)",s
        )
        self.assertIn("target.axpy(sig, self.tangent_mass",s)
        source=PROBE.read_text(encoding="utf-8")
        self.assertIn(
            'examples/shallow_water.py',source
        )
        self.assertIn(
            "os.environ.get(\"SCRIMP_NONLINEAR_MASS_JACOBIAN\")",source
        )
        self.assertIn(
            'if args.mode == "default" and env_mode is not None:',source
        )
        self.assertIn(
            "system.IJacobian(ts, 0., z, zd, shift, m, m)",source
        )


if __name__=="__main__":
    unittest.main()
