# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
"""Check physically compatible initial-data code before native PETSc gates.

Uses GetFEM-like fake DOF ordering ONLY for isolated index mapping.
The DAE/CN/power balance still requires actual native validation.
"""
import ast
from pathlib import Path
import unittest
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
FE=ROOT/"scrimp/utils/fe_initialization.py"
EXAMPLE=ROOT/"examples/heat_wave.py"


class HeatWaveBoundaryInitialData(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.fe=FE.read_text(encoding="utf-8")
        cls.example=EXAMPLE.read_text(encoding="utf-8")

    def test_filtered_boundary_indices_are_not_guessed(self):
        node=next(
            n for n in ast.parse(self.fe).body
            if isinstance(n,ast.FunctionDef)
            and n.name=="filtered_dirichlet_dofs"
        )
        module=ast.Module(body=[node],type_ignores=[])
        ast.fix_missing_locations(module)
        env={"np":np}
        exec(compile(module,str(FE),"exec"),env)
        f=env["filtered_dirichlet_dofs"]
        class Fem:
            def basic_dof_on_region(self,region):
                # Filtered basis order deliberately differs from mesh node IDs.
                return {2: np.array([19,7,42,3]),
                        20: np.array([3,19])}[region]
        fixed,free=f(Fem(),2,20,coefficient_count=4)
        np.testing.assert_array_equal(fixed,[0,3])
        np.testing.assert_array_equal(free,[1,2])
        with self.assertRaises(ValueError):
            f(Fem(),2,20,coefficient_count=5)
        class Invalid:
            def basic_dof_on_region(self,region):
                return {2:[19,7,42,3],20:[7,21]}[region]
        with self.assertRaises(ValueError):
            f(Invalid(),2,20,coefficient_count=4)

    def test_model_uses_explicit_default_clamp_after_l2_projection(self):
        for needle in (
            'initial_wave_boundary="clamp"',
            'initial_wave_profile="gaussian"',
            'hw.set_initial_value("p", p_initial)',
            'hw.set_initial_value("q", "[0.,0.]")',
            'enforce_zero_dirichlet_initial_values(',
            'volume_region=wave_region, boundary_region=20',
            'method=initial_wave_boundary',
            'np.array_equal(',
            'init_step_ts_type="algebraic"',
            '"unconstrained_historical_diagnostic"',
            'hw.heat_wave_initialization_report = report',
        ):
            self.assertIn(needle,self.example)
        self.assertLess(
            self.example.index('hw.set_initial_value("p", p_initial)'),
            self.example.index("enforce_zero_dirichlet_initial_values("),
        )
        self.assertLess(
            self.example.index("enforce_zero_dirichlet_initial_values("),
            self.example.index('hw.set_time_scheme('),
        )

    def test_constrained_l2_is_true_sparse_mass_projection(self):
        for needle in (
            "def constrained_l2_zero_boundary_correction(",
            "M_FF (c-c0)_F = M_FB (c0)_B",
            "model.add_filtered_fem_variable(",
            "model.add_linear_term(",
            'm_free = mat.createSubMatrix(subset, subset)',
            'ksp.getPC().setFactorSolverType(factor_solver)',
            "corrected[fixed] = 0.",
            "corrected[free] += delta_free",
            "stationarity =",
            "not np.isfinite(mass_squared)",
            'if PETSc.COMM_WORLD.getSize() != 1:',
            "def enforce_zero_dirichlet_initial_values(",
            'if method not in ("clamp", "constrained_l2"):',
            'system.set_from_vector(name, corrected)',
        ):
            self.assertIn(needle,self.fe)
        self.assertNotIn('init_step_ts_type="pseudo"',self.example)

    def test_old_probe_boundary_fixes_remain_reproducible(self):
        for path in (
            "sandbox/probe_v130_heat_wave_boundary_compatibility.py",
            "sandbox/probe_v130_hw_dual_sided_dae.py",
            "sandbox/probe_v130_auto_interface_mpi_cn.py",
        ):
            source=(ROOT/path).read_text(encoding="utf-8")
            self.assertIn('initial_wave_boundary="unconstrained"',source)


if __name__=="__main__":
    unittest.main()
