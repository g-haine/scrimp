# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2024 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""
- file:             state.py
- authors:          Giuseppe Ferraro, Ghislain Haine
- date:             31 may 2023
- brief:            class for state object
"""

import petsc4py
import sys

petsc4py.init(sys.argv)
from petsc4py import PETSc

comm = PETSc.COMM_WORLD
rank = comm.getRank()

import logging


class State:
    """This class defines a State."""

    def __init__(
        self,
        name: str,
        description: str,
        kind: str,
        region: int = None,
        mesh_id: int = 0,
    ):
        """This constructor create a State.

        Args:
            name (str): name of the State
            description (str): description of the State
            kind (str): kind of the State
            region (int, optional): region of the State. Defaults to None.
            mesh_id (int, optional): id of the mesh. Defaults to 0.
        """

        self._name = name
        self._description = description
        self._kind = kind
        self._region = region
        self._mesh_id = mesh_id
        self._port = None
        self._costate = None
        self._mesh_id = mesh_id

    def __str__(self):
        where = ""
        if self.get_region() is not None:
            where = ", in region numbered " + str(self.get_region())
        return f"A state variable: {self.get_name()}, describing: {self.get_description()}, has been initialized as a: {self.get_kind()}, on mesh: {self.get_mesh_id()}{where}"

    def get_name(self) -> str:
        """This function gets the name of the State.

        Returns:
            str: name of the state
        """

        return self._name

    def get_description(self) -> str:
        """This function gets the description of the State.

        Returns:
            str: description of the state
        """

        return self._description

    def get_kind(self) -> str:
        """This function gets the kind of the State.

        Returns:
            str: kind of the state
        """

        return self._kind

    def get_region(self) -> int:
        """This function gets the integer number of the region of the state.

        Returns:
            int: region of the State.
        """

        return self._region

    def set_costate(self, costate):
        """This function sets a Co-state to the State.

        Args:
            costate (Costate): Co-state
        """

        if self.get_costate() is None and costate is not None:
            from scrimp.costate import CoState

            try:
                assert isinstance(costate, CoState)
            except AssertionError:
                logging.error("Bad type for costate")
                raise TypeError
            self._costate = costate
        else:
            if rank == 0:
                logging.info("A costate is already present for this state")

    def get_costate(self) -> object:
        """This function gets the Co-state of the state

        Returns:
            object: Costate
        """

        return self._costate

    def set_port(self, port):
        """This function sets a port to the state.

        Args:
            port (Port): Port
        """

        if self.get_port() is None and port is not None:
            from scrimp.port import Port

            try:
                assert isinstance(port, Port)
            except AssertionError:
                logging.error("Bad type for port")
                raise TypeError
            self._port = port
        else:
            if rank == 0:
                logging.info("A port is already present for this state")

    def get_port(self) -> object:
        """This function gets the port of the state

        Returns:
            object: Port
        """

        return self._port

    def get_mesh_id(self) -> int:
        """This funtion gets the integer number of the mesh of the state.

        Returns:
            int: id of the mesh
        """

        return self._mesh_id
