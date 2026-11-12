# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Minimal regressions for DPHS construction and scalar-port bookkeeping."""

import unittest

import scrimp as S


class TestDPHS(unittest.TestCase):
    def test_new_real_system_has_no_ports(self):
        system = S.DPHS("real")
        self.assertEqual(system.ports, {})
        self.assertEqual(system.lumped_variables, set())

    def test_new_real_system_can_register_a_global_state(self):
        system = S.DPHS("real")
        system.set_domain(S.Domain("Interval", {"L": 1.0, "h": 0.5}))
        state = S.State("q", "Displacement", "lumped")
        system.add_state(state)
        system.add_costate(S.CoState("e_q", "Elastic effort", state))
        system.add_lumped_port("q")
        self.assertIn("q", system.lumped_variables)
        self.assertIn("e_q", system.lumped_variables)
        self.assertEqual(system.gf_model.interval_of_variable("q")[1], 1)


if __name__ == "__main__":
    unittest.main()
