# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2022 Ghislain Haine
#
# See the LICENSE file in the root directory for license information.
#
# github: https://github.com/g-haine/scrimp

"""
- file:             dpHs.py
- authors:          Ghislain Haine, Florian Monteghetti, Giuseppe Ferraro
- date:             22 nov. 2022
- last modified:    13 dec. 2022
- brief:            class for distributed port-Hamiltonian system
"""

import os
import sys
import time
import getfem as gf

import petsc4py

petsc4py.init()
from petsc4py import PETSc

comm = PETSc.COMM_WORLD

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pylab as pl
import logging

from utils.linalg import extract_gmm_to_petsc, convert_PETSc_to_scipy
from src import set_default_path

from src.domain import Domain
from src.state import State
from src.costate import CoState
from src.port import Port, Parameter
from src.brick import Brick
from src.hamiltonian import Term, Hamiltonian


class DPHS:
    """
    A generic class handling distributed pHs using the GetFEM tools

    This is a wrapper in order to simplify the coding process

    Access to fine tunings is preserved as much as possible
    """

    def __init__(self, basis_field="real"):
        """
        Constructor of the `distributed port-Hamiltonian system` dpHs class

        :param basis_field: basis field for unknowns (must be `real` or `complex`)
        :type basis_field: str
        """
        # Set the log file
        logging.basicConfig(
            filename="dphs.log",
            encoding="utf-8",
            level=logging.DEBUG,
            filemode="w",
            format="",
        )
        #: The `domain` of a dpHs is an object that handle mesh(es) and dict of regions with getfem indices (for each mesh), useful to define `bricks` (i.e. forms) in the getfem syntax
        self.domain = None
        #: Clear and init `time_scheme` member, which embed PETSc TS options database
        self.get_cleared_TS_options()
        #: The dict of `states`, store many infos for display()
        self.states = dict()
        #: The dict of `costates`, store many infos for display()
        self.costates = dict()
        #: The dict of `ports`, store many infos for display()
        self.ports = dict()
        #: The dict of `bricks`, associating getfem `bricks` to petsc matrices obtained by the PFEM, store many infos for display()
        self.bricks = dict()
        #: The dict of `controls`, collecting information about control ports. Also appear in `ports` entries
        self.controls = dict()
        #: The `Hamiltonian` of a dpHs is a list of dict containing several useful information for each term
        self.hamiltonian = Hamiltonian("Hamiltonian")
        #: Store the computations of the powers f(t)*e(t) on each algebraic port
        self.powers = dict()
        #: To check if the powers have been computed
        self.powers_computed = False
        #: Tangent (non-linear + linear) mass matrix of the system in PETSc CSR format
        self.tangent_mass = PETSc.Mat().create(comm=comm)
        #: Non-linear mass matrix of the system in PETSc CSR format
        self.nl_mass = PETSc.Mat().create(comm=comm)
        #: Linear mass matrix of the system in PETSc CSR format
        self.mass = PETSc.Mat().create(comm=comm)
        #: Tangent (non-linear + linear) stiffness matrix of the system in PETSc CSR format
        self.tangent_stiffness = PETSc.Mat().create(comm=comm)
        #: Non-linear stiffness matrix of the system in PETSc CSR format
        self.nl_stiffness = PETSc.Mat().create(comm=comm)
        #: Linear stiffness matrix of the system in PETSc CSR format
        self.stiffness = PETSc.Mat().create(comm=comm)
        #: rhs of the system in PETSc Vec
        self.rhs = PETSc.Vec().create(comm=comm)
        #: To check if the initial values have been setted before time-resolution
        self.initial_value_setted = dict()
        #: To check if the PETSc TS time-integration parameters have been setted before time-integration
        self.time_scheme["isset"] = False
        #: For monitoring time in TS resolution
        self.ts_start = 0
        #: Will contain both time t and solution z
        self.solution = dict()
        #: Time t where the solution have been saved
        self.solution["t"] = list()
        #: Solution z at time t
        self.solution["z"] = list()
        #: To check if the system has been solved
        self.solve_done = False
        #: A PETSc Vec for residual computation
        self.F = PETSc.Vec().create(comm=comm)
        #: A PETSc Mat for Jacobian computation
        self.J = PETSc.Mat().create(comm=comm)
        #: A PETSc Vec buffering computation
        self.buffer = PETSc.Vec().create(comm=comm)

        #: A getfem `Model` object that is use as core for the dpHs
        self.gf_model = gf.Model(basis_field)
        #: Says to getfem that the `Model` is time-dependent
        self.gf_model.add_initialized_data("t", 0.0, sizes=1)

        # HACK: Macros x, y and z are not available for source term otherwise (why?)
        # Relative problem: source term does not accept numpy functions for the moment
        self.gf_model.add_macro("x", "X(1)")
        self.gf_model.add_macro("y", "X(2)")
        self.gf_model.add_macro("z", "X(3)")

        print("A model with", basis_field, "unknowns has been initialized")

    def set_domain(self, domain: Domain):
        """This function sets a domain for the dphp.
        TODO: If not built_in, given from a script 'name.py' or a .geo file with args in the dict 'parameters' should be able to handle several meshes e.g. for interconnections, hence the list type

        Args:
            name (str): id of the domain, either for built in, or user-defined auxiliary script
            parameters (dict): parameters for the construction, either for built in, or user-defined auxiliary script
        """
        self.domain = domain
        logging.debug(f"domain: {domain.get_name()} has been set")

    def add_state(self, state: State):
        """This functions adds a state to the dpHS.

        Args:
            state (State): the state
        """

        self.states[state.get_name()] = state
        logging.debug(f"state: {state.get_name()} has been added")

    def add_costate(self, costate: CoState):
        """This function adds a costate to the costate list of the dhps.

        Args:
            costate (Costate): the costate"""
        state = costate.get_state()
        state.set_costate(costate)
        self.costates[costate.get_name()] = costate

        port = Port(
            state.get_name(),
            state.get_name(),
            costate.get_name(),
            costate.get_kind(),
            state.get_mesh_id(),
            algebraic=False,
            substituted=costate.get_substituted(),
            region=state.get_region(),
        )
        self.add_port(port)

        state.set_port(port)
        costate.set_port(port)

        logging.debug(
            f"costate: {costate.get_name()} has been added to state: {costate.get_state().get_name()}"
        )
        logging.debug(
            f"state: {state.get_name()} has new costate: {state.get_costate().get_name()}"
        )

    def add_port(self, port: Port):
        """This function adds a `port` object to the dpHs

        Args:
        port (Port):the port"""
        self.ports[port.get_name()] = port
        logging.debug(f"port: {port.get_name()} has been added")

    def add_FEM(self, name_port, order, FEM="CG"):
        """This function defines a FEM (Finite Element Method) for the variables associated to a `port` of the dpHs

        This FEM is a member of the `port`, but it is linked to the getfem `Model` at this stage

        The `kind` member of the `port` is used to deduce the dimension of the FE (`scalar-field` = 1, `vector-field` = dimension of the `mesh_id`-th mesh in `domain`)

        TODO: handle `tensor-field`

        Args:
        name_port (str): the name of the `port` to discretize
        order (int): the order of the FE
        FEM (str): the FE to use, default `CG` for the classical Lagrange element, see port.set_FEM() for more details
        """

        assert self.domain.get_isSet(), "Domain must be setted before adding FEM"

        assert name_port in self.ports.keys(), ("Port", name_port, "does not exist")

        # Select the right dimension
        if self.ports[name_port].get_kind() == "scalar-field":
            dim = 1
        elif self.ports[name_port].get_kind() == "vector-field":
            dim = self.domain.get_dim()[self.ports[name_port].get_mesh_id()]
        else:
            raise ValueError(
                "Unknown kind of variables", self.ports[name_port].get_kind()
            )

        self.ports[name_port].set_fem(
            self.domain.get_mesh()[self.ports[name_port].get_mesh_id()], dim, order, FEM
        )

        # If region is not None, the variable is restricted to the region of mesh_id. Useful for boundary ports or interconnected dpHs on the same mesh
        if self.ports[name_port].get_region() is not None:
            self.gf_model.add_filtered_fem_variable(
                self.ports[name_port].get_flow(),
                self.ports[name_port].get_fem(),
                self.ports[name_port].get_region(),
            )
        else:
            self.gf_model.add_fem_variable(
                self.ports[name_port].get_flow(), self.ports[name_port].get_fem()
            )

        # If the port is algebraic and not substituted, we also add the effort to the list of variables
        if (
            self.ports[name_port].get_algebraic()
            and not self.ports[name_port].get_substituted()
        ):
            if self.ports[name_port].get_region() is not None:
                self.gf_model.add_filtered_fem_variable(
                    self.ports[name_port].get_effort(),
                    self.ports[name_port].get_fem(),
                    self.ports[name_port].get_region(),
                )
            else:
                self.gf_model.add_fem_variable(
                    self.ports[name_port].get_effort(), self.ports[name_port].get_fem()
                )

        # If the port is dynamic and the co-state is not substituted, the co-state must be added as an unknown variable
        if self.ports[name_port].get_effort() in self.costates.keys():
            assert not self.ports[name_port].get_algebraic(), (
                "Port",
                name_port,
                "should not be algebraic because",
                self.ports[name_port].get_effort(),
                "is a co-state",
            )
            if not self.ports[name_port].get_substituted():
                if self.ports[name_port].get_region() is not None:
                    self.gf_model.add_filtered_fem_variable(
                        self.ports[name_port].get_effort(),
                        self.ports[name_port].get_fem(),
                        self.ports[name_port].get_region(),
                    )
                else:
                    self.gf_model.add_fem_variable(
                        self.ports[name_port].get_effort(),
                        self.ports[name_port].get_fem(),
                    )

        # If the port is dynamic, the state needs initialization for time-resolution
        if not self.ports[name_port].get_algebraic():
            self.initial_value_setted[self.ports[name_port].get_flow()] = False

    def add_parameter(self, parameter: Parameter):
        """This function defines a time-independent possibly space-varying parameter (x, y and z are the space variables to use) associated to a `port` of the dpHs

        This parameter is a member of the `port`, but it is added to the dpHs at this stage

        If the FEM of the port has already been setted, the parameter is initialized in this FEM without call needed

        Args:
        parameter (Parameter): the parameter"""
        name_port = parameter.get_name_port()
        self.ports[name_port].add_parameter(parameter)
        logging.debug(
            f"parameter: {parameter.get_name()} has been added to port: {parameter.get_name_port()}"
        )

        if self.ports[name_port].get_is_set():
            self.init_parameter(parameter.get_name(), name_port)

    def init_parameter(self, name, name_port):
        """This function initializes the parameter name in the FEM of the `port` of the dpHs
        evaluate the parameter's expression in the FEM of the `port` name_port and add it to the getfem `Model`

        A parameter is a member of the `port` where it belongs, but it is initialized from
        the parent dpHs of the `port` at this stage, where it is added to the getfem `Model` object

        Args:
        name (str): the name of the parameter as defined with add_parameter()
        name_port (str): the name of the `port` where the parameter belongs
        """

        evaluation = self.ports[name_port].init_parameter(
            name, self.ports[name_port].get_parameter(name).get_expression()
        )
        logging.debug(
            f"expression: {self.ports[name_port].get_parameter(name).get_expression()} for parameter: {name} of port: {name_port}"
        )

        sizes = None
        if self.ports[name_port].get_parameter(name).get_kind() == "tensor-field":
            sizes = evaluation.shape[0]

        self.gf_model.add_initialized_fem_data(
            name, self.ports[name_port].get_fem(), evaluation, sizes=sizes
        )
        logging.debug(
            f"parameter: {name} has been initialized with the FEM of port: {name_port}"
        )

    def set_initial_value(self, name_variable, expression):
        """This function sets the initial value of the variable `name_variable` of the dpHs from an expression

        Args:
        name_variable (str): the name of the variable to set
        expression (str): the expression of the function to use
        """

        if name_variable in self.initial_value_setted.keys():
            if self.initial_value_setted[name_variable]:
                print("Warning: defining again the initial value of", name_variable)
        else:
            raise ValueError("Variable", name_variable, "does not need initialization")

        evaluation = self.gf_model.mesh_fem_of_variable(name_variable).eval(
            expression, globals(), locals()
        )

        initial_value = None
        for name_port in self.ports.keys():
            if (
                name_variable == self.ports[name_port].get_flow()
                or name_variable == self.ports[name_port].get_effort()
            ):
                if self.ports[name_port].get_region() == None:
                    initial_value = evaluation
                else:
                    nb_dofs_total = self.ports[name_port].get_fem().nbdof()
                    dofs_on_region = (
                        self.ports[name_port]
                        .get_fem()
                        .basic_dof_on_region(self.ports[name_port].get_region())
                    )
                    qdim = self.ports[name_port].get_fem().qdim()
                    size = dofs_on_region.shape[0]
                    shape = (qdim, int(size / qdim))
                    evaluation_long = np.reshape(evaluation, (nb_dofs_total,))
                    initial_value_long = evaluation_long[dofs_on_region]
                    initial_value = np.reshape(initial_value_long, shape)
                continue

        if initial_value.all() == None:
            raise ValueError(name_variable, "can not be found in ports")

        self.set_from_vector(name_variable, initial_value)
        self.initial_value_setted[name_variable] = True
        logging.debug(
            f"variable: {name_variable} has been initialized with: {expression}"
        )

    def set_from_vector(self, name_variable, x):
        """This function sets the value of the variable `name_variable` in the getfem `Model` from a numpy vector

        Args:
        name_variable (str): the name of the variable to set
        x (numpy array): the vector of values
        """

        self.gf_model.set_variable(name_variable, x)
        logging.debug(f"value: {x} has been set for variable: {name_variable}")

    def add_brick(self, brick: Brick):
        """This function adds a `brick` in the getfem `Model` thanks to a form in GWFL getfem language

        The form may be non-linear

        Args:
        brick (Brick): the brick
        """
        name_brick = brick.get_name()
        mesh_id = brick.get_mesh_id()
        position = brick.get_position()
        form = brick.get_form()

        self.bricks[name_brick] = brick

        # Flows are on the left-hand side => need a minus for fully implicit formulation in time-resolution
        if position == "flow":
            form = "-(" + form + ")"

        for region in brick.get_regions():
            if brick.get_linear():
                id_brick = self.gf_model.add_linear_term(
                    self.domain._int_method[mesh_id], form, region
                )
                s = ("Linear form '",)

            else:
                id_brick = self.gf_model.add_nonlinear_term(
                    self.domain._int_method[mesh_id], form, region
                )
                s = ("Non-linear form '",)

            logging.debug(
                f"{s},{self.bricks[name_brick].get_form()} has been added as: {position} relation on region: {region} of mesh: {mesh_id}"
            )

            self.bricks[name_brick].add_id_brick_to_list(id_brick)
            logging.debug(
                f"brick IDs: {id_brick} has been added to brick: {name_brick}"
            )

    def add_control_port(
        self,
        name,
        name_control,
        description_control,
        name_observation,
        description_observation,
        kind,
        region=None,
        position="effort",
        mesh_id=0,
    ):
        """This function adds a control `port` to the dpHs

        :param name: the name of the port
        :type name: str
        :param name_control: the name of the control variable
        :type name_control: str
        :param description_control: a physically motivated description of the control
        :type description_control: str
        :param name_observation: the name of the observation variable
        :type name_observation: str
        :param description_observation: a physically motivated description of the observation
        :type description_observation: str
        :param kind: the type of the variables, must be `scalar-field` or `vector-field`
        :type kind: str
        :param region: default=None, the id of the region where the port applies in the mesh mesh_id
        :type region: int
        :param position: default='effort', the position in the flow-effort formulation. Note that 'flow' means that the observation will be a Lagrange multiplier (see examples for more details).
        :type position: str
        :param mesh_id: the id of the mesh where the port belong
        :type mesh_id: int

        :return: appends a `port` to the dpHs and add a control `name` to the dict `controls`
        """

        if position == "effort":
            flow = name_observation
            effort = name_control
        elif position == "flow":
            flow = name_control
            effort = name_observation
        else:
            raise ValueError("Position", position, "is not available for control port")

        self.add_port(
            Port(
                name,
                flow,
                effort,
                kind,
                mesh_id,
                algebraic=True,
                substituted=False,
                region=region,
            )
        )
        self.controls[name] = {
            "name_control": name_control,
            "description_control": description_control,
            "name_observation": name_observation,
            "description_observation": description_observation,
            "kind": kind,
            "region": region,
            "position": position,
            "mesh_id": mesh_id,
            "isset": False,
        }

    def set_control(self, name, expression):
        """This function applies a source term `expression` to the control port `name`

        :param name: the name of the port
        :type name: str
        :param expression: the expression of the source term
        :type expression: str

        :return: construct the constitutive relation `M u = F` setting the control variable u using `brick` and `source` in the getfem `Model`
        """

        if self.controls[name]["kind"] == "scalar-field":
            times = "*"
        elif self.controls[name]["kind"] == "vector-field":
            times = "."
        elif self.controls[name]["kind"] == "tensor-field":
            times = ":"
        else:
            raise ValueError(
                "Unknown kind", self.controls[name]["kind"], "for control port"
            )

        u = self.controls[name]["name_control"]

        mass_form = (
            u + times + "Test_" + u
        )  # form of the mass matrix for the control variable
        self.add_brick(
            Brick(
                "M_" + u,
                mass_form,
                [self.controls[name]["region"]],
                linear=True,
                dt=False,
                position="constitutive",
                mesh_id=self.controls[name]["mesh_id"],
            )
        )

        # Construct the form
        expression_form = "(" + expression + ")" + times + "Test_" + u
        # Add the source brick
        self.add_brick(
            Brick(
                u + "_source",
                expression_form,
                [self.controls[name]["region"]],
                False,
                False,
                "source",
                self.controls[name]["mesh_id"],
            )
        )

        logging.debug(
            f"Control function has been setted to: {u} = {expression} on region: {self.controls[name]['region']} of mesh: {self.controls[name]['mesh_id']}"
        )
        self.controls[name]["isset"] = True
        logging.debug(f"control port: {self.controls[name]} has been set:")

    def assemble_mass(self):
        """This function performs the assembly of the bricks dt=True and linear=True and set the PETSc.Mat attribute `mass`"""

        for _, brick in self.bricks.items():
            # Enable the bricks that are dynamical and linear
            if brick.get_dt() and brick.get_linear():
                brick.enable_id_bricks(self.gf_model)

        self.gf_model.assembly(option="build_matrix")
        size = self.gf_model.nbdof()
        self.mass = extract_gmm_to_petsc(
            [0, size], [0, size], self.gf_model.tangent_matrix()
        )

        for _, brick in self.bricks.items():
            # Disable again all bricks previously enabled
            if brick.get_dt() and brick.get_linear():
                brick.disable_id_bricks(self.gf_model)

    def assemble_stiffness(self):
        """This function performs the assembly of the bricks dt=False and linear=True and set the PETSc.Mat attribute `stiffness`"""

        for _, brick in self.bricks.items():
            # Enable the bricks that are non-dynamical and linear
            if not brick.get_dt() and brick.get_linear():
                brick.enable_id_bricks(self.gf_model)

        self.gf_model.assembly(option="build_matrix")
        size = self.gf_model.nbdof()
        self.stiffness = extract_gmm_to_petsc(
            [0, size], [0, size], self.gf_model.tangent_matrix()
        )

        for _, brick in self.bricks.items():
            # Disable again all bricks previously enabled
            if not brick.get_dt() and brick.get_linear():
                brick.disable_id_bricks(self.gf_model)

    def assemble_rhs(self):
        """This function performs the assembly of the rhs position='source' and set the PETSc.Vec attribute `rhs`"""

        # Remark: I do not understand why enable_all_bricks give absurd results
        # Maybe due to mass matrices!!!
        for _, brick in self.bricks.items():
            # Enable the bricks that are in 'source' position
            if brick.get_position() == "source":
                brick.enable_id_bricks(self.gf_model)
            # And the non-linear ones (and not dt of course)
            if not brick.get_dt() and not brick.get_linear():
                brick.enable_id_bricks(self.gf_model)

        self.gf_model.assembly(option="build_rhs")
        self.rhs = PETSc.Vec().createWithArray(self.gf_model.rhs())

        for _, brick in self.bricks.items():
            # Disable again all bricks previously enabled
            if brick.get_position() == "source":
                brick.disable_id_bricks(self.gf_model)
            # And the non-linear ones (and not dt of course)
            if not brick.get_dt() and not brick.get_linear():
                brick.disable_id_bricks(self.gf_model)

    def assemble_nl_mass(self):
        """This function performs the assembly of the bricks dt=True and linear=False and set the PETSc.Mat attribute `nl_mass`"""

        for _, brick in self.bricks.items():
            # Enable the bricks that are dynamical and non-linear
            if brick.get_dt() and not brick.get_linear():
                brick.enable_id_bricks(self.gf_model)

        self.gf_model.assembly(option="build_matrix")
        size = self.gf_model.nbdof()
        self.nl_mass = extract_gmm_to_petsc(
            [0, size], [0, size], self.gf_model.tangent_matrix()
        )
        self.tangent_mass = self.mass + self.nl_mass
        self.tangent_mass.assemble()

        for _, brick in self.bricks.items():
            # Disable again all bricks previously enabled
            if brick.get_dt() and not brick.get_linear():
                brick.disable_id_bricks(self.gf_model)

    def assemble_nl_stiffness(self):
        """This function performs the assembly of the bricks dt=False and linear=False and set the PETSc.Mat attribute `nl_stiffness`"""

        for _, brick in self.bricks.items():
            # Enable the bricks that are non-dynamical and non-linear
            if not brick.get_dt() and not brick.get_linear():
                brick.enable_id_bricks(self.gf_model)

        self.gf_model.assembly(option="build_matrix")
        size = self.gf_model.nbdof()
        self.nl_stiffness = extract_gmm_to_petsc(
            [0, size], [0, size], self.gf_model.tangent_matrix()
        )
        self.tangent_stiffness = self.stiffness + self.nl_stiffness
        self.tangent_stiffness.assemble()

        for _, brick in self.bricks.items():
            # Disable again all bricks previously enabled
            if not brick.get_dt() and not brick.get_linear():
                brick.disable_id_bricks(self.gf_model)

    def disable_all_bricks(self):
        """This function disables all bricks in the `Model`"""
        for _, brick in self.bricks.items():
            brick.disable_id_bricks(self.gf_model)

    def enable_all_bricks(self):
        """This function enables all bricks in the `Model`"""
        for _, brick in self.bricks.items():
            brick.enable_id_bricks(self.gf_model)

    def IFunction(self, TS, t, z, zd, F):
        """
        IFunction for the time-resolution of the dpHs with PETSc TS fully implicit integrator

        :param TS: the PETSc TS object handling time-resolution
        :type TS: PETSc.TS
        :param t: time parameter
        :type t: float
        :param z: the state
        :type z: PETSc.Vec
        :param zd: the time-derivative of the state
        :type zd: PETSc.Vec
        :param F: the rhs vector
        :type F: PETSc.Vec
        """

        self.gf_model.set_time(t)

        self.assemble_nl_mass()
        self.assemble_nl_stiffness()

        self.tangent_mass.mult(zd, self.buffer)
        self.tangent_stiffness.multAdd(z, self.buffer, F)

        self.assemble_rhs()
        F.axpy(1, self.rhs)

    def IJacobian(self, TS, t, z, zd, sig, A, P):
        """
        IJacobian for the time-resolution of the dpHs with PETSc TS fully implicit integrator

        :param TS: the PETSc TS object handling time-resolution
        :type TS: PETSc.TS
        :param t: time parameter
        :type t: float
        :param z: the state
        :type z: PETSc.Vec
        :param zd: the time-derivative of the state
        :type zd: PETSc.Vec
        :param sig: a shift-parameter, depends on dt
        :type sig: float
        :param  A: the jacobian matrix
        :type A: PETSc.Mat
        :param  P: the jacobian matrix to use for pre-conditionning
        :type P: PETSc.Mat
        """

        self.gf_model.set_time(t)

        # Here improve to avoid re-build when time-step fails and is reduced ?
        # This would forbid or complexify time-varying parameter...
        # Do not modify for the moment, but keep in mind
        self.assemble_nl_mass()
        self.assemble_nl_stiffness()

        self.tangent_stiffness.copy(P)
        P.axpy(sig, self.tangent_mass)

        if A != P:
            print("Operator different from preconditioning")
            A.assemble()

    def get_cleared_TS_options(self):
        """
        To ensure a safe database for the PETSc TS environment
        """

        self.time_scheme = PETSc.Options()
        for key in self.time_scheme.getAll():
            self.time_scheme.delValue(key)

    def set_time_scheme(self, **kwargs):
        """
        Allows an easy setting of the PETSc TS environment

        :param \**kwargs: PETSc TS options and more (see examples)

        :return: set the `time_scheme` attribute as an options database for PETSc
        """

        for key, value in kwargs.items():
            self.time_scheme[key] = value

        if not self.time_scheme.hasName("ts_equation_type"):
            self.time_scheme[
                "ts_equation_type"
            ] = PETSc.TS.EquationType.DAE_IMPLICIT_INDEX2

        if not self.time_scheme.hasName("ts_type") and not self.time_scheme.hasName(
            "ts_ssp"
        ):
            self.time_scheme["ts_type"] = "cn"

        if not self.time_scheme.hasName("ksp_type"):
            self.time_scheme["ksp_type"] = "gmres"

        if not self.time_scheme.hasName("pc_type"):
            self.time_scheme["pc_type"] = "lu"

        if not self.time_scheme.hasName("pc_factor_mat_solver_type"):
            self.time_scheme["pc_factor_mat_solver_type"] = "superlu"

        if not self.time_scheme.hasName("t_0"):
            self.time_scheme["t_0"] = 0.0

        if not self.time_scheme.hasName("t_f"):
            self.time_scheme["t_f"] = 1.0

        if not self.time_scheme.hasName("dt"):
            self.time_scheme["dt"] = 0.01

        if not self.time_scheme.hasName("dt_save"):
            self.time_scheme["dt_save"] = 0.01

        if not self.time_scheme.hasName("ts_adapt_dt_max"):
            self.time_scheme["ts_adapt_dt_max"] = self.time_scheme["dt_save"]

        if not self.time_scheme.hasName("init_step"):
            self.time_scheme["init_step"] = True

        self.time_scheme["isset"] = True

    def exclude_algebraic_var_from_lte(self, TS):
        """
        Exclude the algebraic variable from the local error troncature in the time-resolution

        :param TS: the PETSc TS object handling time-resolution
        :type TS: PETSc.TS
        """

        atol = TS.getTolerances()[1]
        atol_v = atol * np.ones((self.gf_model.nbdof(),))
        id_alg = np.array([], dtype=int)
        # We add indices of all algebraic ports
        for name_port in self.ports.keys():
            if self.ports[name_port].get_algebraic():
                I = self.gf_model.interval_of_variable(self.ports[name_port].get_flow())
                id_alg = np.concatenate(
                    (id_alg, np.arange(I[0], I[0] + I[1], dtype=int))
                )
                # If it is not substituted, it also has the effort part
                if not self.ports[name_port].get_substituted():
                    I = self.gf_model.interval_of_variable(
                        self.ports[name_port].get_effort()
                    )
                    id_alg = np.concatenate(
                        (id_alg, np.arange(I[0], I[0] + I[1], dtype=int))
                    )
            # Else on dynamical ports
            else:
                # If costate are not substituted, we also add them
                if not self.ports[name_port].get_substituted():
                    I = self.gf_model.interval_of_variable(
                        self.ports[name_port].get_effort()
                    )
                    id_alg = np.concatenate(
                        (id_alg, np.arange(I[0], I[0] + I[1], dtype=int))
                    )
        id_alg.sort()
        atol_v[id_alg] = np.inf
        atol_v_petsc = PETSc.Vec().createWithArray(
            np.zeros(
                self.gf_model.nbdof(),
            ),
            comm=comm,
        )
        atol_v_petsc.setArray(atol_v)
        TS.setTolerances(atol=atol_v_petsc)

    def solve(self):
        """
        Perform the time-resolution of the dpHs thanks to PETSc TS

        The options database is setted in the `time_scheme` attribute

        :return:
            * fill the list solution['t'] with the saved times t
            * fill the list solution['z'] with the saved solutions at t as PETSc.Vec
        """

        assert self.time_scheme[
            "isset"
        ], "Time scheme must be defined before time-resolution"
        for name_variable in self.initial_value_setted.keys():
            assert self.initial_value_setted[name_variable], (
                "Initial value of",
                name_variable,
                "must be defined before time-resolution",
            )
        for name_control in self.controls.keys():
            assert self.controls[name_control]["isset"], (
                "Control",
                name_control,
                "must be defined before time-resolution",
            )

        # Initialize time scheme options (without override)
        self.set_time_scheme()

        # Load time scheme options to PETSc
        self.time_scheme = PETSc.Options()

        # Initialization
        self.disable_all_bricks()  # enabling is done in assembly accurately
        self.gf_model.set_time(float(self.time_scheme["t_0"]))
        self.gf_model.first_iter()

        self.assemble_mass()
        self.assemble_nl_mass()

        self.assemble_stiffness()
        self.assemble_nl_stiffness()

        self.J = self.tangent_stiffness.duplicate()
        self.J.assemble()

        self.assemble_rhs()

        self.F = self.tangent_mass.createVecRight()
        self.buffer = self.tangent_mass.createVecRight()

        self.ts_start = time.time()
        if (
            self.time_scheme["init_step"] == "true"
        ):  # Because PETSc DB store everything as str
            self.init_step()

        # TS
        TS = PETSc.TS().create(comm=comm)
        monitor = lambda TS, i, t, z: self.monitor(
            TS,
            i,
            t,
            z,
            dt_save=float(self.time_scheme["dt_save"]),
            t_0=float(self.time_scheme["t_0"]),
        )
        TS.setMonitor(monitor)
        TS.setIFunction(self.IFunction, self.F)
        TS.setIJacobian(self.IJacobian, self.J)
        TS.setTime(float(self.time_scheme["t_0"]))
        TS.setMaxTime(float(self.time_scheme["t_f"]))
        TS.setTimeStep(float(self.time_scheme["dt"]))
        TS.setMaxSNESFailures(-1)
        TS.setExactFinalTime(PETSc.TS.ExactFinalTime.INTERPOLATE)
        # TS.setExactFinalTime(PETSc.TS.ExactFinalTime.MATCHSTEP)
        TS.setFromOptions()
        self.exclude_algebraic_var_from_lte(TS)
        TS.solve(PETSc.Vec().createWithArray(self.gf_model.from_variables(), comm=comm))

        print(f"Elapsed time: {time.time() - self.ts_start:1.4g}s")
        print(
            f"Steps: {TS.getStepNumber()} ({TS.getStepRejections()} rejected, {TS.getSNESFailures()} Nonlinear solver failures)"
        )
        print(
            f"Nonlinear iterations: {TS.getSNESIterations()}, Linear iterations: {TS.getKSPIterations()}"
        )
        TS.reset()

        TS.destroy()

        self.enable_all_bricks()

        self.solve_done = True

    def init_step(self):
        """
        Perform a first initial step with a pseudo bdf

        It needs set_time['init_step']=True to be called.
        Can be setted up with set_time_scheme(init_step=True)

        :return: a consistent initial value in the `Model`
        """

        # TS
        TS = PETSc.TS().create(comm=comm)
        monitor = lambda TS, i, t, z: self.monitor(
            TS,
            i,
            t,
            z,
            dt_save=1.0,
            t_0=float(self.time_scheme["t_0"]),
            initial_step=True,
        )
        TS.setMonitor(monitor)
        TS.setIFunction(self.IFunction, self.F)
        TS.setIJacobian(self.IJacobian, self.J)
        TS.setTime(float(self.time_scheme["t_0"]))
        TS.setMaxTime(float(self.time_scheme["t_0"]) + 1e-6)
        TS.setTimeStep(1e-6)
        TS.setMaxSNESFailures(-1)
        TS.setMaxSteps(2.0)
        TS.setExactFinalTime(PETSc.TS.ExactFinalTime.MATCHSTEP)

        # Save options for later use and avoid overidden them with the initial step
        saved_options = dict()
        if self.time_scheme.hasName("ts_type"):
            saved_options["ts_type"] = self.time_scheme["ts_type"]
            self.time_scheme.delValue("ts_type")
        if self.time_scheme.hasName("ts_ssp"):
            saved_options["ts_ssp"] = self.time_scheme["ts_ssp"]
            self.time_scheme.delValue("ts_ssp")
        self.time_scheme["ts_type"] = "pseudo"
        self.time_scheme["ts_pseudo_increment"] = 1.1
        self.time_scheme["ts_pseudo_fatol"] = 1e-6
        self.time_scheme["ts_pseudo_frtol"] = 1e-9

        TS.setFromOptions()
        self.exclude_algebraic_var_from_lte(TS)

        print(
            "Perform an initial step using pseudo bdf scheme for initial value consistency"
        )
        TS.solve(PETSc.Vec().createWithArray(self.gf_model.from_variables(), comm=comm))
        print(f"Initialisation done in {time.time() - self.ts_start:1.4g}s")
        TS.reset()
        TS.destroy()

        # Delete options for initial step
        self.time_scheme.delValue("ts_type")
        self.time_scheme.delValue("ts_pseudo_increment")
        self.time_scheme.delValue("ts_pseudo_fatol")
        self.time_scheme.delValue("ts_pseudo_frtol")

        # Re-load previously saved options
        if "ts_type" in saved_options.keys():
            self.time_scheme["ts_type"] = saved_options["ts_type"]
        if "ts_ssp" in saved_options.keys():
            self.time_scheme["ts_ssp"] = saved_options["ts_ssp"]

    def monitor(self, TS, i, t, z, dt_save=1, t_0=0.0, initial_step=False):
        """
        Monitor to use during iterations of time-integration at each successful time step

        :param TS: the PETSc TS object handling time-resolution
        :type TS: PETSc.TS
        :param i: the iteration in the time-resolution
        :type i: int
        :param t: time parameter
        :type t: float
        :param z: the state
        :type z: PETSc.Vec
        :param dt_save: save the solution each dt_save s (different from the time-step `dt` used for resolution)
        :type dt_save: float
        :param t_0: the initial time
        :type t_0: float
        :param initial_step: `True` if this is the initial consistency step (default=`False`)
        :type initial_step: bool
        """

        if not initial_step:
            if self.solution["t"]:
                next_saved_at_t = self.solution["t"][-1] + dt_save
            else:
                next_saved_at_t = t_0
            if (next_saved_at_t - t <= 0.1 * TS.getTimeStep()) or (i == -1):
                sys.stdout.write(
                    f"\ri={i:8d} t={t:8g} * ({int(time.time() - self.ts_start)}s)         \n"
                )
                sys.stdout.flush()
                self.solution["t"].append(t)
                self.solution["z"].append(z.copy())
            else:
                sys.stdout.write(
                    f"\ri={i:8d} t={t:8g}   ({int(time.time() - self.ts_start)}s)         "
                )
                sys.stdout.flush()

        self.gf_model.to_variables(z.array)  # Update the state in the getfem `Model`
        self.gf_model.next_iter()  # Says to getfem that we go to the next iteration, not sure if needed

    def compute_Hamiltonian(self):
        """
        Compute each `term` constituting the Hamiltonian

        :return: fill the Hamiltonian[term]['values'] with a list of values at times t
        """

        assert (
            self.solve_done
        ), "System has not been solved yet, Hamiltonian can not be computed"

        print("Start computing the Hamiltonian")
        start = time.time()
        for t in range(len(self.solution["t"])):
            self.gf_model.to_variables(self.solution["z"][t])
            for term in range(len(self.hamiltonian)):
                term_value_at_t = 0.0
                for region in self.hamiltonian[term].get_regions():
                    term_value_at_t += gf.asm(
                        "generic",
                        self.domain._int_method[self.hamiltonian[term].get_mesh_id()],
                        0,
                        self.hamiltonian[term].get_expression(),
                        region,
                        self.gf_model,
                    )
                self.hamiltonian[term].set_value(term_value_at_t)

        self.hamiltonian.set_is_computed()
        logging.debug(f"Hamiltonian has been computed in {str(time.time() - start)} s")

    def plot_Hamiltonian(self, with_powers=True):
        """
        Plot each term constituting the Hamiltonian and the Hamiltonian

        May include the power terms

        :param with_powers: if `True` (default), the plot wil also contains the power of each algebraic ports
        :type with_powers: bool

        :return: a matplotlib figure
        """

        if not self.hamiltonian.get_is_computed():
            self.compute_Hamiltonian()

        t = np.array(self.solution["t"])
        fig = plt.figure(figsize=[8, 5])
        ax = fig.add_subplot(111)
        HamTot = np.zeros(t.size)

        for term in range(len(self.hamiltonian)):
            values = np.array(self.hamiltonian[term].get_values())
            HamTot += values
            ax.plot(t, values, label=self.hamiltonian[term].get_description())

        if len(self.hamiltonian) > 1:
            ax.plot(t, HamTot, label=self.hamiltonian.get_name())

        if with_powers:
            self.plot_powers(ax, HamTot=HamTot)

        ax.legend()
        ax.grid(axis="both")
        ax.set_xlabel("time t")
        ax.set_ylabel("Hamiltonian terms")
        ax.set_title("Evolution of Hamiltonian terms")
        plt.show()

    def compute_powers(self):
        """
        Compute each power associated to each algebraic port if it is not substituted (because of the parameter-dependency if it is)
        """

        assert (
            self.solve_done
        ), "System has not been solved yet, powers can not be computed"

        print("Start computing the powers (substituted ports are not automated)")
        start = time.time()
        for name_port in self.ports.keys():
            if (
                self.ports[name_port].get_algebraic()
                and not self.ports[name_port].get_substituted()
            ):
                if self.ports[name_port].get_region() == None:
                    region = -1
                else:
                    region = self.ports[name_port].get_region()

                if self.ports[name_port].get_kind() == "scalar-field":
                    times = "*"
                elif self.ports[name_port].get_kind() == "vector-field":
                    times = "."
                elif self.ports[name_port].get_kind() == "tensor-field":
                    times = ":"
                else:
                    raise ValueError(
                        "Unknown kind",
                        self.ports[name_port].get_kind(),
                        "for control port",
                    )

                # Control ports needs a minus for a better plot
                if name_port in self.controls.keys():
                    minus = "-"
                else:
                    minus = ""

                form = (
                    minus
                    + self.ports[name_port].get_flow()
                    + times
                    + self.ports[name_port].get_effort()
                )

                self.powers[name_port] = []
                for t in range(len(self.solution["t"])):
                    self.gf_model.to_variables(self.solution["z"][t])

                    power_value_at_t = gf.asm(
                        "generic",
                        self.domain._int_method[self.ports[name_port].get_mesh_id()],
                        0,
                        form,
                        region,
                        self.gf_model,
                    )
                    self.powers[name_port].append(power_value_at_t)

        print("Powers have been computed in", time.time() - start, "s")

        self.powers_computed = True

    def plot_powers(self, ax=None, HamTot=None):
        """
        Plot each power associated to each algebraic port

        The time integration for visual comparison with the Hamiltonian is done using a midpoint method

        If HamTot is provided, a `Balance` showing structure-preserving property is shown: must be constant on the plot

        :param ax: the gca of matplotlib when this function is called from plot_Hamiltonian()
        :type ax: matplotlib axes instance or None
        :param HamTot: the values of the Hamiltonian over time
        :type HamTot: numpy array

        :return: a matplotlib figure or terms in a matplotlib figure
        """

        if not self.powers_computed:
            self.compute_powers()

        need_to_set_figure = False

        t = np.array(self.solution["t"])
        if ax is None:
            fig = plt.figure(figsize=[8, 5])
            ax = fig.add_subplot(111)
            need_to_set_figure = True

        if HamTot is not None:
            SP_balance = HamTot.copy()

        for name_port in self.powers.keys():
            values = np.array(self.powers[name_port])
            power = np.zeros(t.size)
            # Integration in time
            for k in range(1, t.size):
                power[k] = power[k - 1] + 0.5 * (t[k] - t[k - 1]) * (
                    values[k - 1] + values[k]
                )
            ax.plot(t, power, label=name_port)
            if HamTot is not None:
                if name_port in self.controls.keys():
                    SP_balance -= power
                else:
                    SP_balance += power

        if HamTot is not None:
            # Check if the balance makes sense: should not have algebraic and substituted port
            check_makes_sense = True
            for name_port in self.ports.keys():
                if (
                    self.ports[name_port].get_algebraic()
                    and self.ports[name_port].get_substituted()
                ):
                    check_makes_sense = False
            if check_makes_sense:
                ax.plot(t, SP_balance, "--", label="Balance")

        if need_to_set_figure:
            ax.legend()
            ax.grid(axis="both")
            ax.set_xlabel("time t")
            ax.set_ylabel("Power terms")
            ax.set_title("Evolution of power terms")
            plt.show()

    def spy_Dirac(self, t=None, state=None):
        """
        !TO DO: improve a lot with position of bricks and no constitutive relations!!!
        """

        if not t == None:
            self.gf_model.set_time(t)
        else:
            self.gf_model.set_time(0.0)

        if not state == None:
            self.gf_model.to_variables(state)
        else:
            self.gf_model.to_variables(
                PETSc.Vec().createWithArray(
                    np.zeros(
                        self.gf_model.nbdof(),
                    ),
                    comm=comm,
                )
            )

        size = self.gf_model.nbdof()
        M = 0
        J = 0

        # Avoid missing a brick!
        self.enable_all_bricks()

        for _, brick in self.bricks.items():
            if brick.get_position() == "flow":
                for region in brick.get_regions():
                    matrix = gf.asm_generic(
                        self.domain._int_method[brick.get_mesh_id()],
                        2,
                        brick.get_form(),
                        region,
                        self.gf_model,
                    )
                    M += convert_PETSc_to_scipy(
                        extract_gmm_to_petsc([0, size], [0, size], matrix)
                    )
            if brick.get_position() == "effort":
                for region in brick.get_regions():
                    matrix = gf.asm_generic(
                        self.domain._int_method[brick.get_mesh_id()],
                        2,
                        brick.get_form(),
                        region,
                        self.gf_model,
                    )
                    J += convert_PETSc_to_scipy(
                        extract_gmm_to_petsc([0, size], [0, size], matrix)
                    )

        pl.spy(M, markersize=0.2)
        pl.show()
        pl.spy(J, markersize=0.2)
        pl.show()

    def export_matrices(self, t=None, state=None, path=None, to="matlab"):
        """
        TODO:
        """

        if path == None:
            path = set_default_path()

    def export_to_pv(self, name_variable, path=None, t="All"):
        """
        Export the solution to .vtu file(s) (for ParaView), with associated .pvd if t='All'

        :param name_variable: the variable to export
        :type name_variable: str
        :param path: the path of the output file (default: in the `outputs` folder of scrimp)
        :type path: str
        :param t: the time values of extraction (default: `All` the stored times), `All`, `Init`, `Final`
        :type t: numpy array
        """

        assert self.solve_done, "The system has not been solved yet"

        if not path:
            path = os.path.join(set_default_path(), name_variable)
            if not os.path.exists(path):
                os.makedirs(path)

        # Get the dofs of name_variable
        I = self.gf_model.interval_of_variable(name_variable)
        J = range(I[0], I[0] + I[1])

        if t == "All":
            sys.stdout.write(
                "Export all time values of " + name_variable + " is starting...\n"
            )
            with open(os.path.join(path, name_variable + ".pvd"), "w") as pvd_file:
                pvd_file.write('<?xml version="1.0"?>\n')
                pvd_file.write('<VTKFile type="Collection" version="0.1">\n')
                pvd_file.write("  <Collection>\n")
                for k in range(len(self.solution["t"])):
                    z = self.solution["z"][k][J]
                    FEM = self.gf_model.mesh_fem_of_variable(name_variable)
                    FEM.export_to_vtu(
                        os.path.join(path, name_variable + f"_{k:05d}.vtu"),
                        "ascii",
                        FEM,
                        z,
                        name_variable,
                    )
                    pvd_file.write(
                        '    <DataSet timestep="'
                        + str(self.solution["t"][k])
                        + '" file="'
                        + name_variable
                        + f"_{k:05d}.vtu"
                        + '" />\n'
                    )
                    sys.stdout.write("\rExport time %f" % self.solution["t"][k])
                    sys.stdout.flush()
                pvd_file.write("  </Collection>\n")
                pvd_file.write("</VTKFile>")
            pvd_file.close()
            sys.stdout.write("\rExport is done              \n")
        elif t == "Init":
            print("Export initial value of", name_variable)
            z = self.solution["z"][0][J]
            FEM = self.gf_model.mesh_fem_of_variable(name_variable)
            FEM.export_to_vtu(
                os.path.join(path, name_variable + "_init.vtu"),
                "ascii",
                FEM,
                z,
                name_variable,
            )
        elif t == "Final":
            print("Export final value of", name_variable)
            z = self.solution["z"][-1][J]
            FEM = self.gf_model.mesh_fem_of_variable(name_variable)
            FEM.export_to_vtu(
                os.path.join(path, name_variable + "_final.vtu"),
                "ascii",
                FEM,
                z,
                name_variable,
            )
        else:
            raise ValueError("t must be `All`, `Init` or `Final`. Unknown value:", t)

    def display(self, verbose=2):
        """
        A method giving infos about the dpHs

        !TO DO: improve presentation.

        :param verbose: the level of verbosity
        :type verbose: int
        """

        if verbose > 0:
            self.gf_model.variable_list()
            self.gf_model.brick_list()

        if verbose > 1:
            self.domain.display()

            print(self.states)
            print(self.costates)
            print(self.ports)
            print(self.bricks)
            print(self.controls)

            for term in self.hamiltonian.get_terms():
                print(term)

            for name_port in self.ports.keys():
                print(self.ports[name_port])
