import unittest
from scrimp.costate import CoState
from scrimp.state import State
from scrimp.port import Port


class TestCoState(unittest.TestCase):
    def test_str(self):
        state = State("name_state", "description_state", "kind_state", 1, 2)
        # try to set a costate when the state doesn't have one yet
        costate = CoState("name_costate", "descrtiption_costate", state, True)
        str1 = "A co-state variable: name_costate, describing: descrtiption_costate, associated to state: name_state, has been initialized as a: kind_state, on mesh: 2, in region numbered 1"
        str2 = "The constitutive relations between the state: name_state, and the co-state: name_costate, will be substituted for the resolution: variable name_costate, will not be considered as an unknown"
        self.assertEqual(str1 + "\n" + str2, costate.__str__())

    def test_get_state(self):
        state = State("name_state", "description_state", "kind_state", 1, 2)
        # try to set a costate when the state doesn't have one yet
        costate = CoState("name_costate", "descrtiption_costate", state, True)
        self.assertEqual(state, costate.get_state())

    def test_get_substituted(self):
        state = State("name_state", "description_state", "kind_state", 1, 2)
        costate = CoState("name_costate", "descrtiption_costate", state)
        self.assertEqual(False, costate.get_substituted())


if __name__ == "__main__":
    unittest.main()
