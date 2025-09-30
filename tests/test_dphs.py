import unittest
from unittest.mock import Mock

from scrimp import CoState, DPHS, State


class TestDPHSBuilder(unittest.TestCase):
    def setUp(self):
        self.dphs = DPHS()
        self.state = State("x", "displacement", "scalar-field", region=1, mesh_id=0)
        self.costate = CoState("p", "momentum", self.state)

    def test_set_domain_through_builder(self):
        domain = Mock()
        domain.get_name.return_value = "mock-domain"

        self.dphs.set_domain(domain)

        self.assertIs(self.dphs.domain, domain)
        domain.get_name.assert_called()

    def test_register_state_and_costate(self):
        self.dphs.add_state(self.state)
        self.dphs.add_costate(self.costate)

        self.assertIn(self.state.get_name(), self.dphs.states)
        self.assertIn(self.costate.get_name(), self.dphs.costates)
        self.assertIn(self.state.get_name(), self.dphs.ports)
        self.assertIs(self.state.get_port(), self.dphs.ports[self.state.get_name()])
        self.assertIn(self.state.get_name(), self.dphs.initial_value_set)
        self.assertFalse(self.dphs.initial_value_set[self.state.get_name()])

    def test_fluent_builder_usage(self):
        domain = Mock()
        domain.get_name.return_value = "mock-domain"

        (
            self.dphs.builder.with_domain(domain)
            .add_state(self.state)
            .add_costate(self.costate)
        )

        self.assertIs(self.dphs.domain, domain)
        self.assertIn(self.state.get_name(), self.dphs.states)
        self.assertEqual(self.costate.get_state(), self.state)
        self.assertIsNotNone(self.state.get_port())
        self.assertIn(self.state.get_name(), self.dphs.ports)


if __name__ == "__main__":
    unittest.main()
