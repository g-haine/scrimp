# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2023 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""
- file:             hamiltonian.py
- authors:          Giuseppe Ferraro, Ghislain Haine
- date:             31 may 2023
- brief:            class for hamiltonian and term objects
"""

import time
import logging
import numpy as np
import matplotlib.pyplot as plt
import getfem as gf
from scrimp.domain import Domain
from petsc4py import PETSc
import petsc4py
import os
import sys

petsc4py.init(sys.argv)

comm = PETSc.COMM_WORLD
rank = comm.getRank()


module_path = os.path.join(__file__[:-22], "outputs")
print(module_path)


class Term:
    """This class defines a term for the Hamiltoninan."""

    def __init__(
        self, description: str, expression: str, regions: str, mesh_id: int = 0
    ):
        """This constructor defines the object Term for the Hamiltonian functional.

        Args:
            description (str): the name or description of the term (e.g. 'Kinetic energy')
            expression (str): the formula, using the `Model` variables, defining the term. Parameters are allowed (e.g. '0.5*q.T.q')
            regions (str): the region IDs of the mesh where the expression has to be evaluated
            mesh_id (sint): the mesh id of the mesh where the regions belong to.
        """

        self.__description = description
        self.__expression = expression
        self.__regions = regions
        self.__mesh_id = mesh_id
        self.__values = []

    def get_description(self) -> str:
        """This function gets the description of the term.

        Returns:
            str: description of the term
        """

        return self.__description

    def get_expression(self) -> str:
        """This function gets the matematical expression of the term.

        Returns:
            str: matematical expression of the term.
        """

        return self.__expression

    def get_regions(self) -> str:
        """This function gets the regions of the term.

        Returns:
            str: regions of the termthe region IDs of the mesh where the expression has to be evaluated
        """

        return self.__regions

    def get_mesh_id(self) -> int:
        """This function gets the mesh id of the mesh where the regions belong to..

        Returns:
            int: the mesh id of the mesh where the regions belong to.
        """

        return self.__mesh_id

    def get_values(self) -> list:
        """This function gets the valeus of the term.

        Returns:
            list: list of values of the term
        """

        return self.__values.copy()

    def set_value(self, value):
        """This function sets a value for the term.

        Args:
            value : a value for the term
        """

        try:
            assert isinstance(value, float)
        except AssertionError as err:
            logging.error(
                "Can't add value {value}, bad type"
            )
            raise err
        self.__values.append(value)


class Hamiltonian:
    """This class defines the Hamiltoninan."""

    def __init__(self, name: str) -> None:
        """This constructor defines the object Hamiltonian functional.

        Args:
            name (str): the name or description of the Hamiltonian
        """

        self.__terms = []
        self.__name = name
        self.__is_computed = False
        self.__n = None

    def add_term(self, term: Term):
        """This function adds a term to the term list of the Hamiltonian

        Args:
            term (Term): term for the Hamiltonian
        """

        try:
            assert isinstance(term, Term)
        except AssertionError as err:
            logging.error(
                f"Term {term} does not exist."
            )
            raise err

        self.__terms.append(term)

    def compute(self, solution: dict, gf_model: gf.Model, domain: Domain):
        """Compute each `term` constituting the Hamiltonian

        Args:
            - solutions (dict):         The solution of the dphs
            - gf_model (GetFEM Model):  The model getfem of the dphs
            - domain (Domain):          The domain of the dphs
        """

        if rank == 0:
            logging.info(
                "Start computing the Hamiltonian"
            )

        start = time.perf_counter()
        for t, _ in enumerate(solution["t"]):
            gf_model.to_variables(solution["z"][t])
            for _, term in enumerate(self):
                term_value_at_t = 0.0
                for region in term.get_regions():
                    term_value_at_t += gf.asm(
                        "generic",
                        domain._int_method[term.get_mesh_id()],
                        0,
                        term.get_expression(),
                        region,
                        gf_model,
                    )
                term.set_value(term_value_at_t)

        self.set_is_computed()

        if rank == 0:
            logging.info(
                f"Hamiltonian has been computed in {str(time.perf_counter() - start)} s"
            )

    def get_name(self) -> str:
        """This function returns the name of the Hamiltonion.

        Args:
            name (str): name of the Hamiltonian
        """

        return self.__name

    def set_name(self, name: str):
        """This function set the name for the Hamiltonion. For plotting purposes.

        Args:
            name (str): name of the Hamiltonian
        """

        self.__name = name

    def __len__(self) -> int:
        """This function returns the length of the Hamiltonian, that is the number of therms.

        Returns:
            int: _description_
        """
        return len(self.__terms)

    def __iter__(self):
        self.__n = -1
        return self

    def __next__(self):
        if self.__n < len(self) - 1:
            self.__n += 1
            return self.__terms[self.__n]
        else:
            raise StopIteration

    def get_terms(self) -> list:
        """This function returns a copy of the list of terms of the Hamiltonian.

        Returns:
            list: list of terms of the Hamiltonian.
        """

        return self.__terms.copy()

    def __getitem__(self, index) -> Term:
        """This function return the term correspondig ti the index.

        Args:
            index (_type_): _description_

        Returns:
            Term: _description_
        """

        try:
            return self.__terms[index]
        except IndexError as err:
            logging.error(
                f"index {index} out of {len(self.__terms)}"
            )
            raise err

    def set_is_computed(self):
        """This function sets the Hamiltonian as computed."""

        self.__is_computed = True

    def get_is_computed(self) -> bool:
        """This function returns True if the Hamiltonian is computed, False otherwise.

        Returns:
            bool: flag that indicates if the hamiltonian terms have been computed
        """

        return self.__is_computed
