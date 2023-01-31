import unittest
from ..src.costate import CoState
from ..src.state import State
from ..src.port import Port


class TestState(unittest.TestCase):

    def test_name(self):
        state = State("name", "description", "kind", 1, 2)
        self.assertEqual(state.get_name(), "name")  # add assertion here

    def test_description(self):
        state = State("name", "description", "kind", 1, 2)
        self.assertEqual(state.get_description(), "description")  # add assertion here

    def test_kind(self):
        state = State("name", "description", "kind", 1, 2)
        self.assertEqual(state.get_kind(), "kind")  # add assertion here

    def test_region(self):
        state = State("name", "description", "kind", 1, 2)
        self.assertEqual(state.get_region(), 1)  # add assertion here

    def test_mesh_id(self):
        state = State("name", "description", "kind", 1, 2)
        self.assertEqual(state.get_mesh_id(), 2)  # add assertion here

    def test_set_costate(self):
        state = State("name", "description", "kind", 1, 2)
        # try to set a costate when the state doesn't have one yet
        costate1 = CoState("name_costate1", "descrtiption_costate", state, True)
        state.set_costate(costate1)
        self.assertEqual(costate1, state.get_costate())

        # try to set a costate when the state does have one already
        costate2 = CoState("name_costate2", "descrtiption_costate", state, True)
        state.set_costate(costate2)
        self.assertEqual(costate1, state.get_costate())

    def test_get_costate(self):
        state = State("name", "description", "kind", 1, 2)
        self.assertIsNone(state.get_costate())
        costate = CoState("name_costate", "descrtiption_costate", state, True)
        state.set_costate(costate)
        self.assertEqual(costate, state.get_costate())

    def test_set_port(self):
        state = State("name", "description", "kind", 1, 2)
        port = Port("name_port", "flow_port", "effort_port", "kind_port", 0, False, False, 0)
        state.set_port(port)
        self.assertEqual(port, state.get_port())

    def test_get_port(self):
        state = State("name", "description", "kind", 1, 2)
        self.assertIsNone(state.get_port())
        port = Port("name_port", "flow_port", "effort_port", "kind_port", 0, False, False, 0)
        state.set_port(port)
        self.assertEqual(port, state.get_port())

    def test_str(self):
        state = State("name", "description", "kind", 2, 1)
        self.assertEqual("A state variable: name, describing: description, has been initialized as a: kind, on mesh: 1, in region numbered 2",state.__str__())
        state = State("name", "description", "kind",)
        self.assertEqual("A state variable: name, describing: description, has been initialized as a: kind, on mesh: 0",state.__str__())

if __name__ == '__main__':
    unittest.main()

