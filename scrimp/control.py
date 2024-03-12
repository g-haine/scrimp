# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2024 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""
- file:             control.py
- authors:          Giuseppe Ferraro, Ghislain Haine
- date:             31 may 2023
- brief:            class for control port object
"""

from scrimp.port import Port
import logging


class Control_Port(Port):
    """This class defines a Control Port."""

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
        """_summary_

        Args:
            name (str): the name of the port, will be used mainly for plotting purpose.
            name_control (str): the name of the control, used in forms for Brick definition.
            description_control (str): the physical description of the control, used for plotting purpose.
            name_observation (str): the name of the observation, used in forms for Brick definition.
            description_observation (str): the physical description of the observation, used for plotting purpose.
            kind (str): must be `scalar-field` or `vector-field`.
            region (int, optional): the region of mesh_id where the control is applied, defaults to None (=everywhere).
            position (str, optional): says if the control is on the `flow` side or on the `effort` side of the Dirac structure. Defaults to `effort`.
            mesh_id (int, optional): the id of the mesh where the form applies. Defaults to 0.
        """

        if position == "effort":
            flow = name_observation
            effort = name_control
        elif position == "flow":
            flow = name_control
            effort = name_observation
        else:
            logging.error(f"Position {position} is not available for control port")
            raise ValueError

        super().__init__(
            name,
            flow,
            effort,
            kind,
            mesh_id,
            algebraic=True,
            substituted=False,
            dissipative=True,  #: This allows to compute correctly the balance and to be plot with an opposite sign
            region=region,
        )

        self.__name_control = name_control
        self.__name_observation = name_observation
        self.__description_control = description_control
        self.__description_observation = description_observation
        self.__position = position

    def get_name_control(self) -> str:
        """This function gets the name of the control.

        Returns:
            str: the name of the control
        """

        return self.__name_control

    def get_name_obervation(self) -> str:
        """This function gets the name of the obervation.

        Returns:
            str: the name of the observation
        """

        return self.__name_observation

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

    def get_position(self) -> str:
        """This function gets the position of the control in the Dirac structure.

        Returns:
            str: the position of the control
        """

        return self.__position


if __name__ == "__main__":
    Control_Port(
        "Boundary control (bottom)",
        "U_B",
        "Normal force",
        "Y_B",
        "Velocity trace",
        "scalar-field",
        region=10,
        position="effort",
    )
