from scrimp import State


class CoState(State):
    """This class defines a Co-State."""

    def __init__(self, name: str, description: str, state: State, substituted=False):
        """This constructor creates a Co-State.

        Args:
            name (str): name of the Co-State
            description (str): description of the Co-State
            state (State): the State to which the Co-State is bounded
            substituted (bool, optional): boolean that defines whether to substitute the variable. Defaults to False.
        """
        assert isinstance(
            state, State
        ), f"State {state} must be added before its co-state"
        super().__init__(
            name, description, state.get_kind(), state.get_region(), state.get_mesh_id()
        )
        del self._costate
        self._state = state
        self._substituted = substituted

    def __str__(self):
        where = ""
        if self.get_state().get_region() is not None:
            where = ", in region numbered " + str(self.get_state().get_region())
        str1 = f"A co-state variable: {self.get_name()}, describing: {self.get_description()}, associated to state: {self.get_state().get_name()}, has been initialized as a: {self.get_kind()}, on mesh: {self.get_mesh_id()}{where}"
        str2 = f'The constitutive relations between the state: {self.get_state().get_name()}, and the co-state: {self.get_name()}, will{(not self.get_substituted()) * " not"} be substituted for the resolution: variable {self.get_name()}, will { (self.get_substituted()) * "not"} be considered as an unknown'
        return str1 + "\n" + str2

    def get_state(self) -> object:
        """This function gets the State of the Costate.

        Returns:
            object: State
        """
        return self._state

    def get_substituted(self) -> bool:
        """This function gets the boolean indicating wether to substitute the variable or not.

        Returns:
            bool: boolean indicating wether to substitute the variable
        """
        return self._substituted


if __name__ == "__main__":
    print("COSTATE module OK!")
