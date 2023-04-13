from scrimp import Port


class Control_Port(Port):
    def __init__(
        self,
        name: str,
        name_control: str,
        description_control: str,
        name_observation: str,
        description_observation: str,
        kind: str,
        region: int = None,
        position: str = "effort",
        mesh_id: int = 0,
    ):
        if position == "effort":
            flow = name_observation
            effort = name_control
        elif position == "flow":
            flow = name_control
            effort = name_observation
        else:
            raise ValueError("Position", position, "is not available for control port")

        super().__init__(
            name,
            flow,
            effort,
            kind,
            mesh_id,
            algebraic=True,
            substituted=False,
            region=region,
        )

        self.__name_control = name_control
        self.__name_observation = name_observation
        self.__description_control = description_control
        self.__description_observation = description_observation

    def get_name_control(self) -> str:
        """This function gets the name of the control.

        Returns:
            str: the name of the control
        """
        return self.__name_control

    def get_name_obervation(self) -> str:
        """This function gets the name of the obervation.

        Returns:
            str: the name of the obervation
        """
        return self.__name_obervation

    def get_description_control(self) -> str:
        """This function gets the description of the control.

        Returns:
            str: the description of the control
        """
        return self.__description_control

    def get_description_observation(self) -> str:
        """This function gets the description of the observation.

        Returns:
            str: the description of the observation
        """
        return self.__description_observation


if __name__ == "__main__":
    Control_Port(
        "Boundary control (bottom)",
        "U_B",
        "Normal force",
        "Y_B",
        "Velocity trace",
        "scalar-field",
        region=10,
    )
