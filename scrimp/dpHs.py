# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2022 Ghislain Haine
#
# See the LICENSE file in the root directory for license information.
#
# github: https://github.com/g-haine/scrimp

"""
- file:             dpHs.py
- authors:          Ghislain Haine, Florian Monteghetti
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

from scrimp.utils.linalg import extract_gmm_to_petsc, convert_PETSc_to_scipy
from scrimp import set_default_path


class dpHs:
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

        self.domain = (
            domain()
        )  #: The `domain` of a dpHs is an object that handle mesh(es) and dict of regions with getfem indices (for each mesh), useful to define `bricks` (i.e. forms) in the getfem syntax

        self.get_cleared_TS_options()  #: Clear and init `time_scheme` member, which embed PETSc TS options database

        self.states = dict()  #: The dict of `states`, store many infos for display()
        self.costates = (
            dict()
        )  #: The dict of `costates`, store many infos for display()
        self.ports = dict()  #: The dict of `ports`, store many infos for display()
        self.bricks = (
            dict()
        )  #: The dict of `bricks`, associating getfem `bricks` to petsc matrices obtained by the PFEM, store many infos for display()
        self.controls = (
            dict()
        )  #: The dict of `controls`, collecting information about control ports. Also appear in `ports` entries
        self.Hamiltonian = (
            list()
        )  #: The `Hamiltonian` of a dpHs is a list of dict containing several useful information for each term
        self.Hamiltonian_name = "Hamiltonian"  #: For plot purpose
        self.Hamiltonian_computed = (
            False  #: Top check if the Hamiltonian terms have been computed
        )
        self.powers = (
            dict()
        )  #: Store the computations of the powers f(t)*e(t) on each algebraic port
        self.powers_computed = False  #: To check if the powers have been computed
        self.tangent_mass = PETSc.Mat().create(
            comm=comm
        )  #: Tangent (non-linear + linear) mass matrix of the system in PETSc CSR format
        self.nl_mass = PETSc.Mat().create(
            comm=comm
        )  #: Non-linear mass matrix of the system in PETSc CSR format
        self.mass = PETSc.Mat().create(
            comm=comm
        )  #: Linear mass matrix of the system in PETSc CSR format
        self.tangent_stiffness = PETSc.Mat().create(
            comm=comm
        )  #: Tangent (non-linear + linear) stiffness matrix of the system in PETSc CSR format
        self.nl_stiffness = PETSc.Mat().create(
            comm=comm
        )  #: Non-linear stiffness matrix of the system in PETSc CSR format
        self.stiffness = PETSc.Mat().create(
            comm=comm
        )  #: Linear stiffness matrix of the system in PETSc CSR format
        self.rhs = PETSc.Vec().create(comm=comm)  #: rhs of the system in PETSc Vec
        self.initial_value_setted = (
            dict()
        )  #: To check if the initial values have been setted before time-resolution
        self.time_scheme[
            "isset"
        ] = False  #: To check if the PETSc TS time-integration parameters have been setted before time-integration
        self.ts_start = 0  #: For monitoring time in TS resolution
        self.solution = dict()  #: Will contain both time t and solution z
        self.solution["t"] = list()  #: Time t where the solution have been saved
        self.solution["z"] = list()  #: Solution z at time t
        self.solve_done = False  #: To check if the system has been solved
        self.F = PETSc.Vec().create(comm=comm)  #: A PETSc Vec for residual computation
        self.J = PETSc.Mat().create(comm=comm)  #: A PETSc Mat for Jacobian computation
        self.buffer = PETSc.Vec().create(
            comm=comm
        )  #: A PETSc Vec buffering computation

        self.gf_model = gf.Model(
            basis_field
        )  #: A getfem `Model` object that is use as core for the dpHs
        self.gf_model.add_initialized_data(
            "t", 0.0, sizes=1
        )  #: Says to getfem that the `Model` is time-dependent

        # HACK: Macros x, y and z are not available for source term otherwise (why?)
        # Relative problem: source term does not accept numpy functions for the moment
        self.gf_model.add_macro("x", "X(1)")
        self.gf_model.add_macro("y", "X(2)")
        self.gf_model.add_macro("z", "X(3)")

        print("A model with", basis_field, "unknowns has been initialized")

    def set_domain(self, name, parameters, built_in=True):
        """
        Define the object 'domain' of the pHs,
        either manually or using built_in geometries

        !TO DO: If not built_in, given from a script 'name.py' or a .geo file with args in the dict 'parameters'
        should be able to handle several meshes e.g. for interconnections, hence the list type

        :param name: id of the domain, either for built in, or user-defined auxiliary script
        :type name: str
        :param parameters: parameters for the construction, either for built in, or user-defined auxiliary script
        :type parameters: dict
        :param built_in: boolean to swith between built-in or user-defined domain
        :type built_in: bool

        :return: Construct the `domain` attribute of the dpHs
        """

        if built_in:
            import scrimp.utils.mesh

            mesh_function = getattr(scrimp.utils.mesh, name)
            gf_mesh = mesh_function(parameters, terminal=0)
            self.domain.mesh = [gf_mesh[0]]
            self.domain.dim = [gf_mesh[1]]
            self.domain.subdomains = [gf_mesh[2]]
            self.domain.boundaries = [gf_mesh[3]]

            self.domain.isset = True
            self.domain.set_mim_auto()

            if self.domain.isset:
                print("Domain has been setted")
                self.domain.display()

    def add_state(self, name, description, kind, region=None, mesh_id=0):
        """
        Add a variable to the dict `states` of the dpHs

        !TO DO: handling of `tensor-field` state

        :param name: the name of the state
        :type name: str
        :param description: a physically motivated description (e.g. `linear momentum`)
        :type description: str
        :param kind: the unknown type, must be `scalar-field` or `vector-field`
        :type kind: str
        :param region: the index of the region in mesh_id where the variable belong, useful with multi-domains mesh for interconnected dpHs
        :type region: int
        :param mesh_id: a state has to be associated to a unique mesh of the 'domain' attribute (overlap needs the definition of two different states with interface interaction)
        :type mesh_id: int

        :return: append the new state to the dict `states` of the dpHs
        """

        self.states[name] = {
            "description": description,
            "kind": kind,
            "region": region,
            "mesh_id": mesh_id,
            "costate": None,
            "port": None,
        }

        where = ""
        if region is not None:
            where = ", in region numbered " + str(region)
        print(
            "A state variable",
            name,
            ", describing '",
            description,
            "', has been initialized as a",
            kind,
            "on mesh",
            mesh_id,
            where,
        )

    def add_costate(self, name, description, state, substituted=False):
        """
        Add a variable to the dict `costates` of the dpHs

        A `state` is mandatory. The kind of field of the co-state is that of the state

        A `port` named as the associated state is created

        :param name: the name of the costate
        :type name: str
        :param description: a physically motivated description (e.g. `velocity`)
        :type description: str
        :param state: the name of the state
        :type state: str
        :param substituted: if 'True' (default: `False`) the constitutive relations are substituted into the dynamic
        :type substituted: bool

        :return:
            * append the new costate to the dict `costates` of the dpHs
            * append the non-algebraic port (state, costate) to the dict `ports` of the dpHs
        """

        assert state in self.states.keys(), (
            "State",
            state,
            "must be added before its co-state",
        )

        self.costates[name] = {
            "description": description,
            "kind": self.states[state]["kind"],
            "region": self.states[state]["region"],
            "mesh_id": self.states[state]["mesh_id"],
            "state": state,
            "port": None,
        }

        self.states[state]["costate"] = name

        self.add_port(
            state,
            state,
            name,
            self.costates[name]["kind"],
            self.states[state]["mesh_id"],
            algebraic=False,
            substituted=substituted,
            region=self.states[state]["region"],
        )

        self.states[state]["port"] = state
        self.costates[name]["port"] = state

        where = ""
        if self.states[state]["region"] is not None:
            where = ", in region numbered " + str(self.states[state]["region"])
        print(
            "A co-state variable",
            name,
            ", describing '",
            description,
            "', associated to state '",
            state,
            "' has been initialized as a",
            self.costates[name]["kind"],
            "on mesh",
            self.costates[name]["mesh_id"],
            where,
        )
        print(
            "The constitutive relations between the state",
            state,
            "and the co-state",
            name,
            "will"
            + (not substituted) * " not"
            + " be substituted for the resolution: variable",
            name,
            "will" + substituted * " not" + " be considered as an unknown",
        )

    def add_port(
        self,
        name,
        flow,
        effort,
        kind,
        mesh_id=0,
        algebraic=True,
        substituted=False,
        region=None,
    ):
        """
        Add a `port` object to the dpHs

        :param name: the name of the port
        :type name: str
        :param flow: the name id of the flow
        :type flow: str
        :param effort: the name id of the effort
        :type effort: str
        :param kind: the kind of the flow and effort variable (e.g. `scalar-field`)
        :type kind: str
        :param mesh_id: the id of the mesh where the variables belong to
        :type mesh_id: int
        :param algebraic: if `False` (default: `True`), the flow variable will be derivated in time at resolution
        :type algebraic: bool
        :param substituted: if `True`, the constitutive relation is substituted into the mass matrix of the flow, there is only one unknown in the getfem model
        :type substituted: bool
        :param region: the int identifying the region in mesh_id where the port belong
        :type region: int

        :return: create a `port` object and append it to the dict `ports` of the dpHs
        """

        self.ports[name] = port(
            name, flow, effort, kind, mesh_id, algebraic, substituted, region
        )

        where = ""
        if self.ports[name].region is not None:
            where = ", in region numbered " + str(self.ports[name].region)
        print(
            "A port '",
            name,
            "' (",
            flow,
            ",",
            effort,
            ") of",
            kind,
            "type has been added on mesh",
            mesh_id,
            where,
        )

    def add_FEM(self, name_port, order, FEM="CG"):
        """
        Define a FEM (Finite Element Method) for the variables associated to a `port` of the dpHs

        This FEM is a member of the `port`, but it is linked to the getfem `Model` at this stage

        The `kind` member of the `port` is used to deduce the dimension of the FE (`scalar-field` = 1, `vector-field` = dimension of the `mesh_id`-th mesh in `domain`)

        !TO DO: handle `tensor-field`

        :param name_port: the name of the `port` to discretize
        :type name_port: str
        :param order: the order of the FE
        :type order: int
        :param FEM: the FE to use, default `CG` for the classical Lagrange element, see port.set_FEM() for more details
        :type FEM: str

        :return: define a FEM for the flow and effort variables in the `port` name_port, and add them as fem_variable to the getfem `Model`
        """

        assert self.domain.isset, "Domain must be setted before adding FEM"

        assert name_port in self.ports.keys(), ("Port", name_port, "does not exist")

        # Select the right dimension
        if self.ports[name_port].kind == "scalar-field":
            dim = 1
        elif self.ports[name_port].kind == "vector-field":
            dim = self.domain.dim[self.ports[name_port].mesh_id]
        else:
            raise ValueError("Unknown kind of variables", self.ports[name_port]["kind"])

        self.ports[name_port].set_FEM(
            self.domain.mesh[self.ports[name_port].mesh_id], dim, order, FEM
        )

        # If region is not None, the variable is restricted to the region of mesh_id. Useful for boundary ports or interconnected dpHs on the same mesh
        if self.ports[name_port].region is not None:
            self.gf_model.add_filtered_fem_variable(
                self.ports[name_port].flow,
                self.ports[name_port].FEM,
                self.ports[name_port].region,
            )
        else:
            self.gf_model.add_fem_variable(
                self.ports[name_port].flow, self.ports[name_port].FEM
            )

        # If the port is algebraic and not substituted, we also add the effort to the list of variables
        if self.ports[name_port].algebraic and not self.ports[name_port].substituted:
            if self.ports[name_port].region is not None:
                self.gf_model.add_filtered_fem_variable(
                    self.ports[name_port].effort,
                    self.ports[name_port].FEM,
                    self.ports[name_port].region,
                )
            else:
                self.gf_model.add_fem_variable(
                    self.ports[name_port].effort, self.ports[name_port].FEM
                )

        # If the port is dynamic and the co-state is not substituted, the co-state must be added as an unknown variable
        if self.ports[name_port].effort in self.costates.keys():
            assert not self.ports[name_port].algebraic, (
                "Port",
                name_port,
                "should not be algebraic because",
                self.ports[name_port].effort,
                "is a co-state",
            )
            if not self.ports[name_port].substituted:
                if self.ports[name_port].region is not None:
                    self.gf_model.add_filtered_fem_variable(
                        self.ports[name_port].effort,
                        self.ports[name_port].FEM,
                        self.ports[name_port].region,
                    )
                else:
                    self.gf_model.add_fem_variable(
                        self.ports[name_port].effort, self.ports[name_port].FEM
                    )

        # If the port is dynamic, the state needs initialization for time-resolution
        if not self.ports[name_port].algebraic:
            self.initial_value_setted[self.ports[name_port].flow] = False

    def add_parameter(self, name, description, kind, expression, name_port):
        """
        Define a time-independent possibly space-varying parameter (x, y and z are the space variables to use) associated to a `port` of the dpHs

        This parameter is a member of the `port`, but it is added to the dpHs at this stage

        If the FEM of the port has already been setted, the parameter is initialized in this FEM without call needed

        :param name: the name of the parameter
        :type name: str
        :param description: a physically motivated description (e.g. `mass density`)
        :type description: str
        :param kind: the type of parameter, must be `scalar-field`, `vector-field` or `tensor-field`
        :type kind: str
        :param expression: the expression in getfem syntax (e.g. `x*y`, `[x, y]`, etc.)
        :type expression: str
        :param name_port: the name id of the `port` where the parameter belong to
        :type name_port: str

        :return: define a parameter to the `port` name_port, and initialized it if the FEM of the port is already defined
        """

        self.ports[name_port].add_parameter(name, description, kind, expression)

        if self.ports[name_port].isset:
            self.init_parameter(name, name_port)

    def init_parameter(self, name, name_port):
        """
        Initialize the parameter name in the FEM of the `port` of the dpHs

        A parameter is a member of the `port` where it belongs, but it is initialized from
        the parent dpHs of the `port` at this stage, where it is added to the getfem `Model` object

        :param name: the name of the parameter as defined with add_parameter()
        :type parameter: str
        :param name_port: the name of the `port` where the parameter belongs
        :type name_port: str

        :return: evaluate the parameter's expression in the FEM of the `port` name_port and add it to the getfem `Model`
        """

        evaluation = self.ports[name_port].init_parameter(
            name, self.ports[name_port].parameters[name]["expression"]
        )

        sizes = None
        if self.ports[name_port].parameters[name]["kind"] == "tensor-field":
            sizes = evaluation.shape[0]

        self.gf_model.add_initialized_fem_data(
            name, self.ports[name_port].FEM, evaluation, sizes=sizes
        )
        print("Parameter", name, "has been initialized with the FEM of port", name_port)

    def set_initial_value(self, name_variable, expression):
        """
        Set the initial value of the variable `name_variable` of the dpHs from an expression

        :param name_variable: the name of the variable to set
        :type name_variable: str
        :param expression: the expression of the function to use
        :type expression: str

        :return: set variable `name_variable` of the getfem Model of the dpHs to `expression`, evaluated on the FEM of the port variable
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
                name_variable == self.ports[name_port].flow
                or name_variable == self.ports[name_port].effort
            ):
                if self.ports[name_port].region == None:
                    initial_value = evaluation
                else:
                    nb_dofs_total = self.ports[name_port].FEM.nbdof()
                    dofs_on_region = self.ports[name_port].FEM.basic_dof_on_region(
                        self.ports[name_port].region
                    )
                    qdim = self.ports[name_port].FEM.qdim()
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
        print("Variable", name_variable, "has been initialized with:", expression)

    def set_from_vector(self, name_variable, x):
        """
        Set the value of the variable `name_variable` in the getfem `Model` from a numpy vector

        :param name_variable: the name of the variable to set
        :type name_variable: str
        :param x: the vector of values
        :type x: numpy array

        :return: set variable `name_variable` of the getfem `Model` of the dpHs to `x`
        """

        self.gf_model.set_variable(name_variable, x)

    def add_brick(
        self,
        name,
        form,
        regions,
        linear=True,
        dt=False,
        position="constitutive",
        mesh_id=0,
    ):
        """
        Add a `brick` in the getfem `Model` thanks to a form in GWFL getfem language

        The form may be non-linear

        :param name: the name of the brick, will be used mainly for plotting purpose
        :type name: str
        :param form: the form in GWFL getfem language
        :type form: str
        :param regions: the regions of mesh_id where the form applies
        :type regions: list(int)
        :param linear: default=True, parameter to help easy identification of linear bricks
        :type linear: bool
        :param dt: default=False, parameter to help easy identification of matrices applied to time-derivative of a variable (e.g. mass matrices)
        :type dt: bool
        :param position: default='constitutive', parameter to help easy identification of "where" is the form: in the Dirac structure ('flow' side or 'effort' side), or in the 'constitutive' relations. This serves for both the time-resolution and plotting purposes.
        :type position: str
        :param mesh_id: default=0, the id of the mesh where the form applies
        :type mesh: int

        :return: add bricks to the getfem `Model`
        """

        self.bricks[name] = {
            "id_bricks": [],
            "form": form,
            "mesh_id": mesh_id,
            "regions": regions,
            "linear": linear,
            "dt": dt,
            "position": position,
        }

        # Flows are on the left-hand side => need a minus for fully implicit formulation in time-resolution
        if position == "flow":
            form = "-(" + form + ")"

        if linear:
            for region in regions:
                id_brick = self.gf_model.add_linear_term(
                    self.domain.int_method[mesh_id], form, region
                )
                self.bricks[name]["id_bricks"].append(id_brick)
                print(
                    "Linear form '",
                    self.bricks[name]["form"],
                    "' has been added as",
                    position,
                    "relation on region",
                    region,
                    "of mesh",
                    mesh_id,
                )
        else:
            for region in regions:
                id_brick = self.gf_model.add_nonlinear_term(
                    self.domain.int_method[mesh_id], form, region
                )
                self.bricks[name]["id_bricks"].append(id_brick)
                print(
                    "Non-linear form '",
                    self.bricks[name]["form"],
                    "' has been added as",
                    position,
                    "relation on region",
                    region,
                    "of mesh",
                    mesh_id,
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
        """
        Add a control `port` to the dpHs

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
            name,
            flow,
            effort,
            kind,
            mesh_id,
            algebraic=True,
            substituted=False,
            region=region,
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
        """
        Apply a source term `expression` to the control port `name`

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
            "M_" + u,
            mass_form,
            [self.controls[name]["region"]],
            linear=True,
            dt=False,
            position="constitutive",
            mesh_id=self.controls[name]["mesh_id"],
        )

        # Construct the form
        expression_form = "-(" + expression + ")" + times + "Test_" + u
        # Add the source brick
        source_id = self.gf_model.add_source_term(
            self.domain.int_method[self.controls[name]["mesh_id"]],
            expression_form,
            self.controls[name]["region"],
        )

        self.bricks[u + "_source"] = {
            "id_bricks": [source_id],
            "form": expression_form,
            "mesh_id": self.controls[name]["mesh_id"],
            "regions": [self.controls[name]["region"]],
            "linear": False,
            "dt": False,
            "position": "source",
        }

        print(
            "Control function has been setted to",
            u,
            "=",
            expression,
            "on region",
            self.controls[name]["region"],
            "of mesh",
            self.controls[name]["mesh_id"],
        )
        self.controls[name]["isset"] = True

    def assemble_mass(self):
        """
        Perform the assembly of the bricks dt=True and linear=True and set the PETSc.Mat attribute `mass`
        """

        for name in self.bricks.keys():
            # Enable the bricks that are dynamical and linear
            if self.bricks[name]["dt"] and self.bricks[name]["linear"]:
                for k in self.bricks[name]["id_bricks"]:
                    self.gf_model.enable_bricks(k)

        self.gf_model.assembly(option="build_matrix")
        size = self.gf_model.nbdof()
        self.mass = extract_gmm_to_petsc(
            [0, size], [0, size], self.gf_model.tangent_matrix()
        )

        for name in self.bricks.keys():
            # Disable again all bricks previously enabled
            if self.bricks[name]["dt"] and self.bricks[name]["linear"]:
                for k in self.bricks[name]["id_bricks"]:
                    self.gf_model.disable_bricks(k)

    def assemble_stiffness(self):
        """
        Perform the assembly of the bricks dt=False and linear=True and set the PETSc.Mat attribute `stiffness`
        """

        for name in self.bricks.keys():
            # Enable the bricks that are non-dynamical and linear
            if not self.bricks[name]["dt"] and self.bricks[name]["linear"]:
                for k in self.bricks[name]["id_bricks"]:
                    self.gf_model.enable_bricks(k)

        self.gf_model.assembly(option="build_matrix")
        size = self.gf_model.nbdof()
        self.stiffness = extract_gmm_to_petsc(
            [0, size], [0, size], self.gf_model.tangent_matrix()
        )

        for name in self.bricks.keys():
            # Disable again all bricks previously enabled
            if not self.bricks[name]["dt"] and self.bricks[name]["linear"]:
                for k in self.bricks[name]["id_bricks"]:
                    self.gf_model.disable_bricks(k)

    def assemble_rhs(self):
        """
        Perform the assembly of the rhs position='source' and set the PETSc.Vec attribute `rhs`
        """

        # Remark: I do not understand why enable_all_bricks give absurd results
        # Maybe due to mass matrices!!!
        for name in self.bricks.keys():
            # Enable the bricks that are in 'source' position
            if self.bricks[name]["position"] == "source":
                for k in self.bricks[name]["id_bricks"]:
                    self.gf_model.enable_bricks(k)
            # And the non-linear ones (and not dt of course)
            if not self.bricks[name]["dt"] and not self.bricks[name]["linear"]:
                for k in self.bricks[name]["id_bricks"]:
                    self.gf_model.enable_bricks(k)

        self.gf_model.assembly(option="build_rhs")
        self.rhs = PETSc.Vec().createWithArray(self.gf_model.rhs())

        for name in self.bricks.keys():
            # Disable again all bricks previously enabled
            if self.bricks[name]["position"] == "source":
                for k in self.bricks[name]["id_bricks"]:
                    self.gf_model.disable_bricks(k)
            # And the non-linear ones (and not dt of course)
            if not self.bricks[name]["dt"] and not self.bricks[name]["linear"]:
                for k in self.bricks[name]["id_bricks"]:
                    self.gf_model.disable_bricks(k)

    def assemble_nl_mass(self):
        """
        Perform the assembly of the bricks dt=True and linear=False and set the PETSc.Mat attribute `nl_mass`
        """

        for name in self.bricks.keys():
            # Enable the bricks that are dynamical and non-linear
            if self.bricks[name]["dt"] and not self.bricks[name]["linear"]:
                for k in self.bricks[name]["id_bricks"]:
                    self.gf_model.enable_bricks(k)

        self.gf_model.assembly(option="build_matrix")
        size = self.gf_model.nbdof()
        self.nl_mass = extract_gmm_to_petsc(
            [0, size], [0, size], self.gf_model.tangent_matrix()
        )
        self.tangent_mass = self.mass + self.nl_mass
        self.tangent_mass.assemble()

        for name in self.bricks.keys():
            # Disable again all bricks previously enabled
            if self.bricks[name]["dt"] and not self.bricks[name]["linear"]:
                for k in self.bricks[name]["id_bricks"]:
                    self.gf_model.disable_bricks(k)

    def assemble_nl_stiffness(self):
        """
        Perform the assembly of the bricks dt=False and linear=False and set the PETSc.Mat attribute `nl_stiffness`
        """

        for name in self.bricks.keys():
            # Enable the bricks that are non-dynamical and non-linear
            if not self.bricks[name]["dt"] and not self.bricks[name]["linear"]:
                for k in self.bricks[name]["id_bricks"]:
                    self.gf_model.enable_bricks(k)

        self.gf_model.assembly(option="build_matrix")
        size = self.gf_model.nbdof()
        self.nl_stiffness = extract_gmm_to_petsc(
            [0, size], [0, size], self.gf_model.tangent_matrix()
        )
        self.tangent_stiffness = self.stiffness + self.nl_stiffness
        self.tangent_stiffness.assemble()

        for name in self.bricks.keys():
            # Disable again all bricks previously enabled
            if not self.bricks[name]["dt"] and not self.bricks[name]["linear"]:
                for k in self.bricks[name]["id_bricks"]:
                    self.gf_model.disable_bricks(k)

    def disable_all_bricks(self):
        """
        Disable all bricks in the `Model`
        """
        for name in self.bricks.keys():
            for k in self.bricks[name]["id_bricks"]:
                self.gf_model.disable_bricks(k)

    def enable_all_bricks(self):
        """
        Enable all bricks in the `Model`
        """
        for name in self.bricks.keys():
            for k in self.bricks[name]["id_bricks"]:
                self.gf_model.enable_bricks(k)

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
            if self.ports[name_port].algebraic:
                I = self.gf_model.interval_of_variable(self.ports[name_port].flow)
                id_alg = np.concatenate(
                    (id_alg, np.arange(I[0], I[0] + I[1], dtype=int))
                )
                # If it is not substituted, it also has the effort part
                if not self.ports[name_port].substituted:
                    I = self.gf_model.interval_of_variable(self.ports[name_port].effort)
                    id_alg = np.concatenate(
                        (id_alg, np.arange(I[0], I[0] + I[1], dtype=int))
                    )
            # Else on dynamical ports
            else:
                # If costate are not substituted, we also add them
                if not self.ports[name_port].substituted:
                    I = self.gf_model.interval_of_variable(self.ports[name_port].effort)
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

        print(f"Elapsed time: {time.time()-self.ts_start:1.4g}s")
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
        print(f"Initialisation done in {time.time()-self.ts_start:1.4g}s")
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
                    f"\ri={i:8d} t={t:8g} * ({int(time.time()-self.ts_start)}s)         \n"
                )
                sys.stdout.flush()
                self.solution["t"].append(t)
                self.solution["z"].append(z.copy())
            else:
                sys.stdout.write(
                    f"\ri={i:8d} t={t:8g}   ({int(time.time()-self.ts_start)}s)         "
                )
                sys.stdout.flush()

        self.gf_model.to_variables(z.array)  # Update the state in the getfem `Model`
        self.gf_model.next_iter()  # Says to getfem that we go to the next iteration, not sure if needed

    def set_Hamiltonian_term(self, description, expression, regions, mesh_id=0):
        """
        Defines a term of the Hamiltonian functional

        :param description: the name or description of the term (e.g. 'Kinetic energy')
        :type description: str
        :param expression: the formula, using the `Model` variables, defining the term. Parameters are allowed (e.g. '0.5*q.T.q')
        :type expression: str
        :param regions: the region ids of the mesh mesh_id where the expression has to be evaluated
        :type regions: list(int)
        :param mesh_id: the mesh id of the mesh where the regions belong
        :type mesh_id: int
        """

        term = dict()
        term["description"] = description
        term["expression"] = expression
        term["regions"] = regions
        term["mesh_id"] = mesh_id
        term["values"] = list()

        self.Hamiltonian.append(term)

    def set_Hamiltonian_name(self, name):
        """
        Set the `name` of the Hamiltonian. For plotting purposes.

        :param name: the name of the Hamiltonian, sum of all Hamiltonian terms
        :type name: str
        """

        self.Hamiltonian_name = name

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
            for term in range(len(self.Hamiltonian)):
                term_value_at_t = 0.0
                for region in self.Hamiltonian[term]["regions"]:
                    term_value_at_t += gf.asm(
                        "generic",
                        self.domain.int_method[self.Hamiltonian[term]["mesh_id"]],
                        0,
                        self.Hamiltonian[term]["expression"],
                        region,
                        self.gf_model,
                    )
                self.Hamiltonian[term]["values"].append(term_value_at_t)

        print("Hamiltonian has been computed in", time.time() - start, "s")

        self.Hamiltonian_computed = True

    def plot_Hamiltonian(self, with_powers=True):
        """
        Plot each term constituting the Hamiltonian and the Hamiltonian

        May include the power terms

        :param with_powers: if `True` (default), the plot wil also contains the power of each algebraic ports
        :type with_powers: bool

        :return: a matplotlib figure
        """

        if not self.Hamiltonian_computed:
            self.compute_Hamiltonian()

        t = np.array(self.solution["t"])
        fig = plt.figure(figsize=[8, 5])
        ax = fig.add_subplot(111)
        HamTot = np.zeros(t.size)

        for term in range(len(self.Hamiltonian)):
            values = np.array(self.Hamiltonian[term]["values"])
            HamTot += values
            ax.plot(t, values, label=self.Hamiltonian[term]["description"])

        if len(self.Hamiltonian) > 1:
            ax.plot(t, HamTot, label=self.Hamiltonian_name)

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
                self.ports[name_port].algebraic
                and not self.ports[name_port].substituted
            ):
                if self.ports[name_port].region == None:
                    region = -1
                else:
                    region = self.ports[name_port].region

                if self.ports[name_port].kind == "scalar-field":
                    times = "*"
                elif self.ports[name_port].kind == "vector-field":
                    times = "."
                elif self.ports[name_port].kind == "tensor-field":
                    times = ":"
                else:
                    raise ValueError(
                        "Unknown kind", self.ports[name_port].kind, "for control port"
                    )

                # Control ports needs a minus for a better plot
                if name_port in self.controls.keys():
                    minus = "-"
                else:
                    minus = ""

                form = (
                    minus
                    + self.ports[name_port].flow
                    + times
                    + self.ports[name_port].effort
                )

                self.powers[name_port] = []
                for t in range(len(self.solution["t"])):
                    self.gf_model.to_variables(self.solution["z"][t])

                    power_value_at_t = gf.asm(
                        "generic",
                        self.domain.int_method[self.ports[name_port].mesh_id],
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
                    self.ports[name_port].algebraic
                    and self.ports[name_port].substituted
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

        for name in self.bricks.keys():
            brick = self.bricks[name]
            if brick["position"] == "flow":
                for region in brick["regions"]:
                    matrix = gf.asm_generic(
                        self.domain.int_method[brick["mesh_id"]],
                        2,
                        brick["form"],
                        region,
                        self.gf_model,
                    )
                    M += convert_PETSc_to_scipy(
                        extract_gmm_to_petsc([0, size], [0, size], matrix)
                    )
            if brick["position"] == "effort":
                for region in brick["regions"]:
                    matrix = gf.asm_generic(
                        self.domain.int_method[brick["mesh_id"]],
                        2,
                        brick["form"],
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
        !TO DO:
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

            for term in self.Hamiltonian.keys():
                print(self.Hamiltonian[term])

            for name_port in self.ports.keys():
                self.ports[name_port].display()


class domain:
    """
    A class handling meshes and indices for regions

    Lists are used to handle interconnection of pHs, allowing for several meshes in the dpHs.
    """

    def __init__(self):
        """
        Constructor of the `domain` member of a dpHs
        """

        self.isset = False  #: A boolean to check if the domain has been setted
        self.mesh = list()  #: A list of getfem Mesh objects
        self.subdomains = list(
            dict()
        )  #: A list of dict `str`: `int` listing the indices of region of dim n in the getfem mesh
        self.boundaries = list(
            dict()
        )  #: A list of dict `str`: `int` listing the indices of region of dim n-1 in the getfem mesh
        self.dim = list()  #: A list of the dimensions of the getfem meshes
        self.int_method = list()  #: A list of getfem integration method MeshIm objects

    def set_mim_auto(self):
        """
        Define the integration method to a default choice
        """

        for k in range(len(self.mesh)):
            if self.dim[k] == 1:
                self.int_method.append(
                    gf.MeshIm(self.mesh[k], gf.Integ("IM_GAUSS1D(5)"))
                )
            elif self.dim[k] == 2:
                self.int_method.append(
                    gf.MeshIm(self.mesh[k], gf.Integ("IM_TRIANGLE(7)"))
                )
            elif self.dim[k] == 3:
                self.int_method.append(
                    gf.MeshIm(self.mesh[k], gf.Integ("IM_TETRAHEDRON(8)"))
                )
            else:
                self.isset = False
                print(
                    "Integration method has to be setted manually on mesh",
                    k,
                    " of dimension",
                    self.dim[k],
                )

    def display(self):
        """
        A method giving infos about the domain
        """

        assert self.isset, "Domain has not been setted yet"

        print("Domain is setted and contains", len(self.mesh), "mesh:")
        for k in range(len(self.mesh)):
            print("=== on mesh", k, " of dim", self.dim[k])
            print("* Subdomains are:", self.subdomains[k])
            print("* Boundaries are:", self.boundaries[k])
            self.mesh[k].display()


class port:
    """
    A class to handle a `port` of a dpHs

    It is mainly constituted of a flow variable, an effort variable, and a FEM
    """

    def __init__(
        self, name, flow, effort, kind, mesh_id, algebraic, substituted, region
    ):
        """
        Constructor of a `port` of a dpHs

        !TO DO: handle `tensor-field` variable

        :param name: the name of the port
        :type name: str
        :param flow: the name of the flow variable
        :type flow: str
        :param effort: the name of the effort variable
        :type effort: str
        :param kind: the type of the variables, must be `scalar-field` or `vector-field`
        :type kind: str
        :param mesh_id: the id of the mesh where the variables belong
        :type mesh_id: int
        :param algebraic: if `False`, the flow variable will be derivated in time at resolution
        :type algebraic: bool
        :param substituted: if `True`, the constitutive relation is substituted and there is only a getfem variable for the effort
        :type substituted: bool
        :param region: the int identifying the region in mesh_id where the port belong, useful for boundary ports
        :type region: int
        """

        self.isset = False  #: A boolean to check if the domain has been setted
        self.name = name  #: The name of the `port`
        self.flow = flow  #: The name of the flow variable
        self.effort = effort  #: The name of the effort variable
        self.kind = kind  #: The type of the variables (e.g. `scalar-field`)
        self.mesh_id = mesh_id  #: The id of the mesh where the variables belong
        self.algebraic = algebraic  #: If `True`, the equation associated to this port is algebraic, otherwise dynamic and the flow is derivated in time at resoltuion
        self.substituted = substituted  #: If `True, the getfem `Model` will only have an unknown variable for the effort: the constitutive relation is substituted into the mass matrix on the flow side
        self.parameters = (
            dict()
        )  #: A dict of parameters acting on the variables of the `port`
        self.FEM = None  #: A getfem MeshFem object to discretize the `port`
        self.region = region  #: If any, the int of the region of mesh_id where the flow/effort variables belong

    def set_FEM(self, mesh, dim, order, FEM):
        """
        Set the MeshFem getfem object defining the finite element method to use to discretize the port

        !TO DO: handle more FE

        :param mesh: the mesh where the FE are define
        :type mesh: getfem Mesh object
        :param dim: the dim of the FE !TO DO: improve this for automation
        :type dim: int
        :param order: the order of the FE
        :type order: int
        :param FEM: the FE to use
        :type FEM: str

        :return: set a MeshFem object on the mesh where the variables belong
        """

        if FEM == "CG":
            FEM_str = "FEM_PK(" + str(mesh.dim()) + "," + str(order) + ")"
        elif FEM == "DG":
            FEM_str = "FEM_PK_DISCONTINUOUS(" + str(mesh.dim()) + "," + str(order) + ")"
        else:
            raise ValueError(
                "Unknown FEM "
                + FEM
                + " for port "
                + self.name
                + "\nUse the gf_model `Model` attribute to set it directly"
            )

        self.FEM = gf.MeshFem(mesh, dim)
        self.FEM.set_fem(gf.Fem(FEM_str))

        self.isset = True
        print(FEM_str, "has been setted for port", self.name)
        self.FEM.display()

    def add_parameter(self, name, description, kind, expression):
        """
        Add a parameter acting on the variables of the `port`

        This parameter is given as a string expression, and may be space-varying (variable x, y and z)

        It may be a scalar-, vector- or tensor-field (e.g. `x*y', `[x, y]`, `[[1, 0],[0, y]]`, etc.)

        :param name: the name of the parameter
        :type name: str
        :param description: a physically motivated description (e.g. `mass density`)
        :type description: str
        :param kind: the type of the parameter (e.g. `scalar-field`)
        :type kind: str
        :param expression: the expression in getfem syntax
        :type expression: str

        :return: append a parameter to the dict `parameters` of the port
        """

        self.parameters[name] = {
            "description": description,
            "kind": kind,
            "expression": expression,
            "name_port": self.name,
        }
        print(
            "A parameter",
            name,
            ", describing '",
            description,
            "', of type",
            kind,
            "has been added to port '",
            self.name,
            "'",
        )

    def init_parameter(self, name, expression):
        """
        Set the parameter by initialization in the FE basis

        This method should not be called explicitly. The init_parameter() method of the parent `dpHs` should be called instead to ensure the parameter to be added to the getfem `Model` object

        :param name: the name of the parameter to initialize
        :type name: str
        :param expression: the expression in getfem syntax
        :type expression: str

        :return: the evaluation of the parameter in the FEM of the `port`
        :rtype: numpy array
        """

        assert self.isset, (
            "A FEM must be setted for port '",
            self.name,
            "' before initialization",
        )

        assert name in self.parameters.keys(), (
            "Parameter",
            name,
            "must be added before intialization",
        )

        evaluation = self.FEM.eval(expression, globals(), locals())
        print(
            "Parameter",
            name,
            "has been evaluated with the FEM of port '",
            self.parameters[name]["name_port"],
            "', with expression:",
            expression,
        )

        return evaluation

    def display(self):
        """
        A method giving infos about the `port`
        """

        assert self.isset, (
            "Port",
            self.name,
            "has not been setted yet (a FEM is probably missing)",
        )

        print(
            self.name,
            self.flow,
            self.effort,
            self.kind,
            self.mesh_id,
            self.algebraic,
            self.parameters,
        )
        self.FEM.display()
