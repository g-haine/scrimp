# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""
- file:             port.py
- authors:          Giuseppe Ferraro, Ghislain Haine
- date:             31 may 2023
- brief:            class for port and parameter objects
"""

import time
import logging
import numpy as np
import getfem as gf
from scrimp.fem import FEM
from scrimp.domain import Domain
from petsc4py import PETSc
import petsc4py
import sys

petsc4py.init(sys.argv)

comm = PETSc.COMM_WORLD
rank = comm.getRank()

class Parameter:
    """This class describes the Parameter for a Port."""

    def __init__(
        self,
        name: str,
        description: str,
        kind: str,
        expression: str,
        name_port: str,
    ):
        """This constructor defines the object Parameter for a Port.

        Args:
            name (str): name of the parameter
            description (str): description of the paramter
            kind (str): kind of the parameter
            expression (str): expression of the parameter
            name_port (str): name of the Port to apply the parameter
        """

        self.__name = name
        self.__description = description
        self.__kind = kind
        self.__expression = expression
        self.__name_port = name_port

    def get_name(self) -> str:
        """This function gets the name of the parameter.

        Returns:
            str: name of the parameter
        """

        return self.__name

    def get_description(self) -> str:
        """This function gets the description of the parameter.

        Returns:
            str: description of the parameter
        """

        return self.__description

    def get_kind(self) -> str:
        """This function gets the kind of the parameter.

        Returns:
            str: kind of the parameter
        """

        return self.__kind

    def get_expression(self) -> str:
        """This function gets the matematical expression of the parameter.

        Returns:
            str: matematical expression of the parameter.
        """

        return self.__expression

    def get_name_port(self) -> str:
        """This function gets the name of the port whose the parameter is bounded.

        Returns:
            str: name of the port whose the parameter is bounded.
        """

        return self.__name_port


class Port:
    """
    A class to handle a port of a dpHs

    It is mainly constituted of a flow variable, an effort variable, and a fem.
    """

    def __init__(
        self,
        name: str,
        flow: str,
        effort: str,
        kind: str,
        mesh_id: int = 0,
        algebraic: bool = True,
        substituted: bool = False,
        dissipative: bool = True,
        region: int = None,
    ):
        """Constructor of a `port` of a discrete port Hmiltonian system (dpHs).

        Args:
            name (str): the name of the port
            flow (str): the name of the flow variable
            effort (str): the name of the effort variable
            kind (str): the type of the variables, must be `scalar-field` or `vector-field`
            mesh_id (int): the id of the mesh where the variables belong
            algebraic (bool): if `True`, the equation associated to this port is algebraic, otherwise dynamic
            substituted (bool): if `True`, the constitutive relation is substituted and there is only a getfem variable for the effort
            dissipative (bool, optional): for post-processing purpose, indicates wether the port is dissipative-like or not. Defaults to `True`
            region (int): the integer identifying the region in the mesh where the port belong, useful for boundary ports
        """

        self.__isSet = False  #: A boolean to check if the domain has been set
        self.__name = name  #: The name of the `port`
        self.__flow = flow  #: The name of the flow variable
        self.__effort = effort  #: The name of the effort variable
        self.__kind = kind  #: The type of the variables (e.g. `scalar-field`)
        self.__mesh_id = mesh_id  #: The id of the mesh where the variables belong
        self.__algebraic = algebraic  #: If `True`, the equation associated to this port is algebraic, otherwise dynamic
        self.__substituted = substituted  #: If `True`, the getfem `Model` will only have an unknown variable for the effort: the constitutive relation is substituted into the mass matrix on the flow side
        self.__dissipative = dissipative  #: If `True`, the power associated must be taken with a negative sign (see plot_hamiltonian in DPHS class)
        self.__parameters = (
            []
        )  #: A list of parameters acting on the variables of the `port`
        self.__fem = None  #: A getfem Meshfem object to discretize the `port`
        #: If any, the int of the region of mesh_id where the flow/effort variables belong
        self.__region = region
        self.__power = list()  # : The raw signed flow/effort products
        # Sampling is part of the result, never an implicit CN convention.
        # "unknown" preserves explicit user-set_power calls without
        # misrepresenting their integration meaning.
        self.__power_sampling = "unknown"
        # : A boolean flag that indicates if the powers of the `port` has been computed
        self.__is_computed = False

    def get_isSet(self) -> bool:
        """This funcion gets the boolean value that indicates wether the port is set or not.


        Returns:
            bool: value that indicates if the port is set.
        """

        return self.__isSet

    def set_isSet(self) -> bool:
        """This funcion sets the boolean value that indicates the port is set.


        Returns:
            bool: value that indicates if the port is set.
        """

        self.__isSet = True

    def get_name(self) -> str:
        """This function gets the name of the port.

        Returns:
            str: name of the port
        """

        return self.__name

    def get_flow(self) -> str:
        """This function gets the name of the flow variable.

        Returns:
            str: name of the flow variable
        """

        return self.__flow

    def get_effort(self) -> str:
        """This function gets the name of the effort variable.

        Returns:
            str: name of the effort variable
        """

        return self.__effort

    def get_kind(self) -> str:
        """This function gets the type of the variables (e.g. `scalar-field`)

        Returns:
            str: type of the variables (e.g. `scalar-field`)
        """

        return self.__kind

    def get_mesh_id(self) -> int:
        """This function gets he id of the mesh where the variables belong

        Returns:
            int: The id of the mesh where the variables belong
        """

        return self.__mesh_id

    def get_algebraic(self) -> bool:
        """This function gets boolean value of the algebraic parameter of the port.
         If `True`, the equation associated to this port is algebraic, otherwise dynamic and the flow is derivated in time at resoltuion

        Returns:
            bool: value of the algebraic parameter
        """

        return self.__algebraic

    def get_substituted(self) -> bool:
        """This function gets boolean value of the subtituted parameter.
        If `True, the getfem `Model` will only have an unknown variable for the effort: the constitutive relation is substituted into the mass matrix on the flow side

        Returns:
            bool: value of the subtituted parameter
        """

        return self.__substituted

    def get_dissipative(self) -> bool:
        """This function gets boolean value for the dissipativeness flag of the port.
        If `True`, the power associated to the port gets a negative sign

        Returns:
            bool: value of the dissipativeness flag
        """

        return self.__dissipative

    def get_region(self) -> int:
        """This function gets the region of the mesh.
        If any, the int of the region of mesh_id where the flow/effort variables belong

        Returns:
            int: region of the mesh
        """

        return self.__region

    def get_fem(self):
        """This function returns the fetfem Meshfem object to discretize the port.

        Returns:
            _type_: the getfem Meshfem object to discretize the port
        """

        return self.__fem

    def set_fem(self, fem: FEM):
        """This function sets the Meshfem getfem object defining the finite element method to use to discretize the port.

        Args:
            fem (FEM): the FEM object to use
        """

        fem.set_fem()

        try:
            assert fem.get_isSet()
            self.__fem = fem.get_fem()
            self.__isSet = True
        except AssertionError as err:
            if rank == 0:
                logging.error(f"Error while setting the FEM of port: {self.get_name()}")
            raise err

    def get_parameter(self, name) -> Parameter:
        """This function return the parameter with a specific name.

        Args:
            name (str): the name of the parameter of interest

        Returns:
            Parameter: the desired parameter, None otherwise
        """

        for p in self.__parameters:
            if p.get_name() == name:
                return p

        logging.error(f"Parameter with name: {name} does not exit!")
        raise AssertionError

    def get_parameters(self) -> list:
        """This function returns the list of all the parameters inserted for the port.

        Returns:
            list(Parameter): list of all the parameters inserted for the port
        """

        return self.__parameters.copy()

    def add_parameter(self, parameter: Parameter) -> bool:
        """This function adds a Parameter object that is acting on the variables of the port.

        Args:
            parameter (Parameter): parameter for the port.:

        Returns:
            bool: True if the insertion has been complete correctly, False otherwise

        """

        if isinstance(parameter, Parameter):
            self.__parameters.append(parameter)
            return True
        else:
            if rank == 0:
                logging.error("Bad type of parameter")
            raise TypeError

    def init_parameter(self, name: str, expression: str):
        """This function sets the chosen parameter object for the current port by initialization in the FE basis.

        Args:
            name (str): the name of the parameter object
            expression (str):

        Returns:
            out (numpy array): the evaluation of the parameter in the fem of the port.
        """

        for p in self.__parameters:
            if p.get_name() == name:
                assert self.__isSet, (
                    "A fem must be set for port '",
                    self.__name,
                    "' before initialization",
                )
                evaluation = self.__fem.eval(expression, globals(), locals())
                return evaluation
        assert False, (f"Parameter '{name}' must be added before intialization")

    def set_power(self, power, *, sampling="unknown"):
        """Provide an externally calculated signed flow/effort series.

        The caller must explicitly supply `sampling="instantaneous"`
        or `sampling="cn_step"` before plot_powers can integrate it.
        The historical one-argument call remains valid, but the
        integration convention is deliberately UNKNOWN.
        """
        if sampling not in ("unknown", "instantaneous", "cn_step"):
            raise ValueError("Unknown port power sampling convention")
        self.__power = power
        self.__power_sampling = sampling

    def set_is_computed(self):
        """This function sets the power of the Port as computed."""

        self.__is_computed = True

    def get_is_computed(self) -> bool:
        """This function returns True if the power of the Port is computed, False otherwise.

        Returns:
            bool: flag that indicates if the power flowing through the Port has been computed
        """

        return self.__is_computed

    def compute(
        self, solution: dict, gf_model: gf.Model, domain: Domain,
        *, sampling: str = "instantaneous",
    ):
        """Compute raw signed flow/effort powers for an algebraic port.

        Sampling is explicit rather than hardcoding Crank--Nicolson for
        every PETSc time integrator:

        `instantaneous`: at every accepted (t[k], z[k]), suitable as
        nodal *samples* independent of the time scheme. Trapezoidal
        integration of these samples is approximate and not a claim
        of BDF/ARKIMEX stage quadrature or a power identity.

        `cn_step`: at (t[0], z[0]) for k=0; for k>=1, evaluate
        at (t[k], 0.5*(z[k-1]+z[k])), preserving the already-native-
        qualified SCRIMP CN convention exactly. In particular,
        time-dependent GetFEM expressions are evaluated at t[k] as
        they were historically, **not** at the physical midpoint time.

        Args:
            solution: Accepted time/FE-state history, containing
                matching `t` and `z` sequences.
            gf_model: Original GetFEM model (collective MPI assembly).
            domain: SCRIMP integration domain.
            sampling: `instantaneous` (default) or `cn_step`.

        The returned `get_power()` array preserves the original port
        algebraic orientation. `get_power_sampling()` identifies the
        numerical meaning of each stored value.
        """
        if sampling not in ("instantaneous", "cn_step"):
            raise ValueError("Unsupported port-power sampling: " + str(sampling))

        if self.get_algebraic() and not self.get_substituted():
            if rank == 0:
                logging.info(
                    f"Start computing the power flowing through {self.get_name()}"
                )

            start = time.perf_counter()
            if self.get_region() is None:
                region = -1
            else:
                region = self.get_region()

            if self.get_kind() == "scalar-field":
                multiply = "*"
            elif self.get_kind() == "vector-field":
                multiply = "."
            elif self.get_kind() == "tensor-field":
                multiply = ":"
            else:
                raise ValueError("Unknown port kind " + str(self.get_kind()))

            form = self.get_flow() + multiply + self.get_effort()
            times = np.asarray(solution["t"], dtype=float)
            states = solution["z"]
            if (times.ndim != 1 or len(times) == 0
                    or len(times) != len(states)
                    or not np.all(np.isfinite(times))
                    or not np.all(np.diff(times) > 0.)):
                raise ValueError("Invalid accepted port-power time/state history")

            # Do not leave stale data marked computed if a new assembly
            # fails. Each rank runs the SAME GetFEM assembly sequence.
            self.__is_computed = False
            self.__power_sampling = "unknown"
            powers = []
            for k, t in enumerate(times):
                # Both the FE state and model time must correspond to
                # THIS sample. Controls may explicitly depend on time.
                gf_model.set_time(float(t))
                if sampling == "cn_step" and k > 0:
                    state = 0.5 * (solution["z"][k] + solution["z"][k - 1])
                else:
                    state = solution["z"][k]
                gf_model.to_variables(state)
                value = gf.asm(
                    "generic",
                    domain._int_method[self.get_mesh_id()],
                    0, form, region, gf_model,
                )
                scalar = np.asarray(value)
                if (scalar.ndim != 0 or not np.isfinite(scalar)
                        or (np.iscomplexobj(scalar) and
                            abs(scalar.imag) > 1e-12 * max(1., abs(scalar.real)))):
                    raise FloatingPointError(
                        "Invalid original signed power in port "
                        + self.get_name()
                    )
                powers.append(value)

            self.__power = powers
            self.__power_sampling = sampling
            self.set_is_computed()
            if rank == 0:
                logging.info(
                    f"Power in {self.get_name()} computed ({sampling}) in "
                    f"{time.perf_counter() - start} s"
                )
        elif rank == 0:
            logging.info(
                f"Power cannot be computed for dynamic or substituted port {self.get_name()}"
            )

    def get_power_sampling(self) -> str:
        """Numerical meaning of the raw port-power history.

        `instantaneous` means every P[k]=P(t[k], z[k]).
        `cn_step` means P[0] is instantaneous while P[k>=1]
        uses the mean CN state for the interval ending at t[k].
        `unknown` means no integration convention was provided.
        """
        return self.__power_sampling

    def get_power(self):
        """Return signed powers with the convention in get_power_sampling()."""

        try:
            assert self.__is_computed
            return self.__power
        except AssertionError as err:
            logging.error(
                "Power for port %s has not been computed yet",
                self.__name,
            )
            raise err

    def __str__(self) -> str:
        """This function displays all the info of the port.

        Returns:
            str: detailed info of the port
        """

        try:
            assert self.__isSet
            if self.__fem is not None:
                self.__fem.display()
        except AssertionError as err:
            logging.error(
                "Port",
                self.__name,
                "has not been set yet (a fem or lumped registration is missing)",
            )
            raise err

        return f"{self.__name}, {self.__flow}, {self.__effort}, {self.__kind}, {str(self.__mesh_id)}, {str(self.__algebraic)}, {str(self.__substituted)}, {str(self.__dissipative)}, {self.__parameters}"


if __name__ == "__main__":
    p = Port("name", "flow", "effort", "kind", 0, True, True, 1)
    param = Parameter(
        "name_para", "description", "kind_param", "expression", p.get_name()
    )
    p.add_parameter(p)
    print(p)
