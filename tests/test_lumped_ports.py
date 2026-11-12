"""Regression checks for spatially global scalar ports.

Run with: python -m unittest tests.test_lumped_ports
Requires the full SCRIMP environment (GetFEM, PETSc, Gmsh).
"""
import unittest
import numpy as np

try:
    import scrimp as S
except ModuleNotFoundError as exc:
    if exc.name in {"getfem", "petsc4py", "gmsh"}:
        raise unittest.SkipTest("SCRIMP numerical dependencies unavailable") from exc
    raise


class LumpedPortsTest(unittest.TestCase):
    def setUp(self):
        self.dphs = S.DPHS("real")
        self.dphs.set_domain(S.Domain("Interval", {"L": 2.0, "h": 0.5}))
        self.p = S.State("P", "Global total momentum", "lumped")
        self.dphs.add_state(self.p)
        self.dphs.add_costate(S.CoState("e_P", "Global velocity", self.p))

    def test_dynamic_port_has_exactly_one_dof_per_variable(self):
        self.dphs.add_lumped_port("P")
        self.assertEqual(self.dphs.gf_model.interval_of_variable("P")[1], 1)
        self.assertEqual(self.dphs.gf_model.interval_of_variable("e_P")[1], 1)
        self.dphs.set_initial_value("P", 3.0)
        self.assertTrue(self.dphs.initial_value_set["P"])
        np.testing.assert_allclose(self.dphs.gf_model.variable("P"), [3.0])
        self.assertIn("P", self.dphs.lumped_variables)

    def test_global_algebraic_port_has_one_variable(self):
        self.dphs.add_port(
            S.Port("tank_velocity", "v", "unused", "lumped",
                   algebraic=True, substituted=True)
        )
        self.dphs.add_lumped_port("tank_velocity")
        self.assertEqual(self.dphs.gf_model.interval_of_variable("v")[1], 1)
        self.assertNotIn("v", self.dphs.initial_value_set)
        self.assertIn("tank_velocity", self.dphs.ports)

    def test_invalid_port_and_duplicate_registration(self):
        with self.assertRaises(KeyError):
            self.dphs.add_lumped_port("does_not_exist")
        with self.assertRaises(ValueError):
            self.dphs.add_FEM(S.FEM("P", 1))
        self.dphs.add_lumped_port("P")
        with self.assertRaises(ValueError):
            self.dphs.add_lumped_port("P")
        self.assertIn("P", str(self.dphs.ports["P"]))

    def test_non_numeric_initial_data_is_rejected(self):
        self.dphs.add_lumped_port("P")
        with self.assertRaises(ValueError):
            self.dphs.set_initial_value("P", "sin(x)")


if __name__ == "__main__":
    unittest.main()
