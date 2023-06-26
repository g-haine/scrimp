# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2023 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""
- file:             dphs.py
- authors:          Giuseppe Ferraro, Ghislain Haine, Florian Monteghetti
- date:             22 nov. 2022
- brief:            class for distributed port-Hamiltonian system
"""

import getfem as gf
import logging
import matplotlib.pyplot as plt
import matplotlib.pylab as pl
import numpy as np
import os
import petsc4py
import sys
import time
import gc

petsc4py.init(sys.argv)
from petsc4py import PETSc
comm = PETSc.COMM_WORLD
rank = comm.getRank()

from scrimp.domain import Domain
from scrimp.state import State
from scrimp.costate import CoState
from scrimp.port import Parameter, Port
from scrimp.fem import FEM
from scrimp.control import Control_Port
from scrimp.brick import Brick
from scrimp.hamiltonian import Hamiltonian

from scrimp.utils.linalg import extract_gmm_to_petsc, extract_gmm_to_scipy

module_path = os.path.join(__file__[:-15], "outputs")

class DPHS:
    """A generic class handling distributed pHs using the GetFEM tools

    This is a wrapper in order to simplify the coding process

    Access to fine tunings is preserved as much as possible
    """

    def __init__(self, basis_field="real"):
        """Constructor of the `distributed port-Hamiltonian system` dphs class.

        Args:
            basis_field (str): basis field for unknowns (must be `real` or `complex`)
        """
                
        #: The `domain` of a dphs is an object that handle mesh(es) and dict of regions with getfem indices (for each mesh), useful to define `bricks` (i.e. forms) in the getfem syntax
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
        #: The `Hamiltonian` of a dphs is a list of dict containing several useful information for each term
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
        #: To check if the initial values have been set before time-resolution
        self.initial_value_set = dict()
        #: To check if the PETSc TS time-integration parameters have been set before time-integration
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
        #: To stop TS integration by keeping already computed timesteps if one step fails
        self.stop_TS = False
        #: A PETSc Vec for residual computation
        self.F = PETSc.Vec().create(comm=comm)
        #: A PETSc Mat for Jacobian computation
        self.J = PETSc.Mat().create(comm=comm)
        #: A PETSc Vec buffering computation
        self.buffer = PETSc.Mat().create(comm=comm)
        #: A getfem `Model` object that is use as core for the dphs
        self.gf_model = gf.Model(basis_field)
        #: Says to getfem that the `Model` is time-dependent
        self.gf_model.add_initialized_data("t", 0.0, sizes=1)

        # HACK: Macros x, y and z are not available for source term otherwise (why?)
        # Relative problem: source term does not accept numpy functions for the moment (why?)
        self.gf_model.add_macro("x", "X(1)")
        self.gf_model.add_macro("y", "X(2)")
        self.gf_model.add_macro("z", "X(3)")
    
        if rank==0:
            logging.info(
                f"A model with {basis_field} unknowns has been initialized"
            )

    def set_domain(self, domain: Domain):
        """This function sets a domain for the dphs.
        
        TODO: If not built_in, given from a script 'name.py' or a .geo file with args in the dict 'parameters' should be able to handle several meshes e.g. for interconnections, hence the list type

        Args:
            name (str): id of the domain, either for built in, or user-defined auxiliary script
            parameters (dict): parameters for the construction, either for built in, or user-defined auxiliary script
        """
        
        self.domain = domain
        if rank==0:
            logging.info(
                f"domain: {domain.get_name()} has been set"
            )

    def add_state(self, state: State):
        """This functions adds a state to state dict of the dphs.

        Args:
            state (State): the state
        """
        
        self.states[state.get_name()] = state
        
        if rank==0:
            logging.info(
                f"state: {state.get_name()} has been added"
            )

    def add_costate(self, costate: CoState):
        """This function adds a costate to the costate dict of the dphs, then defines and adds the dynamical port gathering the couple (State, CoState).

        Args:
            costate (CoState): the costate
        """
        
        # Set the `costate` in its `state`
        state = costate.get_state()
        state.set_costate(costate)
        
        # Add the `costate` to the list
        self.costates[costate.get_name()] = costate

        # Define the `port` gathering the `state` and the `costate`
        port = Port(
            state.get_name(),
            state.get_name(),
            costate.get_name(),
            costate.get_kind(),
            state.get_mesh_id(),
            algebraic=False,
            dissipative=False,
            substituted=costate.get_substituted(),
            region=state.get_region(),
        )
        # Add the `port` the the dphs
        self.add_port(port)

        # Set the `port` in the `state` and the `costate`
        state.set_port(port)
        costate.set_port(port)

        if rank==0:
            logging.info(
                f"costate: {costate.get_name()} has been added to state: {costate.get_state().get_name()}"
            )
            logging.info(
                f"state: {state.get_name()} has new costate: {state.get_costate().get_name()}"
            )

    def add_port(self, port: Port):
        """This function adds a port to the port dict of the dphs

        Args:
            port (Port): the port
        """
        
        self.ports[port.get_name()] = port
        
        if rank==0:
            logging.info(
                f"port: {port.get_name()} has been added"
            )

    def add_FEM(self, fem: FEM):
        """This function adds a FEM (Finite Element Method) for the variables associated to a port of the dphs

        TODO: handle `tensor-field`

        Args:
            fem (FEM): the FEM to use
        """

        # Check if the `domain` is set
        try:
            assert self.domain.get_isSet()
        except AssertionError as err:
            logging.error(
                "Domain must be set before adding FEM"
            )
            raise err

        # Check if the `port` is in the dict of the dphs
        name_port = fem.get_name()
        try:
            assert name_port in self.ports.keys()
        except AssertionError as err:
            logging.error(
                f"Port {name_port} does not exist"
            )
            raise err

        # Select the dimension, according to the `kind` of the `port`
        port = self.ports[name_port]
        if port.get_kind() == "scalar-field":
            dim = 1
        elif port.get_kind() == "vector-field":
            dim = self.domain.get_dim()[port.get_mesh_id()]
        else:
            logging.error(
                f"Unknown kind of variables {port.get_kind()}"
            )
            raise ValueError

        # Set the dimension in the `fem`
        fem.set_dim(dim)
        
        # Set the mesh of the `domain`, according to the `mesh_ID` of the `port`
        fem.set_mesh(self.domain.get_mesh()[port.get_mesh_id()])

        # Apply the fem to the `port`
        port.set_fem(fem)

        # If `region` is not None, the variable is restricted to the `region` of mesh_id. Useful for boundary ports or interconnected dphs on the same mesh.
        # Add the FEM to the `flow variable` of the `port` in getfem
        if port.get_region() is not None:
            self.gf_model.add_filtered_fem_variable(
                port.get_flow(),
                port.get_fem(),
                port.get_region(),
            )
        else:
            self.gf_model.add_fem_variable(port.get_flow(), port.get_fem())

        # If the port is algebraic AND not substituted, we also add the FEM to the `effort variable` of the `port` in getfem
        if port.get_algebraic() and not port.get_substituted():
            if port.get_region() is not None:
                self.gf_model.add_filtered_fem_variable(
                    port.get_effort(),
                    port.get_fem(),
                    port.get_region(),
                )
            else:
                self.gf_model.add_fem_variable(port.get_effort(), port.get_fem())

        # If the port is dynamic and the co-state is not substituted, the co-state must be added as an unknown variable
        if port.get_effort() in self.costates.keys():
            try:
                assert not port.get_algebraic()
            except AssertionError as err:
                logging.error(
                    f"Port {name_port} should not be algebraic because {port.get_effort()} is a co-state"
                )
                raise err
            if not port.get_substituted():
                if port.get_region() is not None:
                    self.gf_model.add_filtered_fem_variable(
                        port.get_effort(),
                        port.get_fem(),
                        port.get_region(),
                    )
                else:
                    self.gf_model.add_fem_variable(
                        port.get_effort(),
                        port.get_fem(),
                    )

        # If the port is dynamic, the state needs initialization for time-resolution
        if not port.get_algebraic():
            self.initial_value_set[port.get_flow()] = False

    def add_parameter(self, parameter: Parameter):
        """This function adds a time-independent (possibly space-varying parameter: `x`, `y` and `z` are the space variables) associated to a port of the dphs

        Args:
            parameter (Parameter): the parameter
        """
        
        name_port = parameter.get_name_port()
        self.ports[name_port].add_parameter(parameter)
        
        if rank==0:
            logging.info(
                f"{parameter.get_name()} has been added to port: {parameter.get_name_port()}"
            )
        
        # Initialization of parameters
        try:
            assert self.ports[name_port].get_is_set()
            self.init_parameter(parameter.get_name(), name_port)
        except AssertionError:
            logging.warning(
                f"{parameter.get_name()} has not been initialized because the port {parameter.get_name_port()} does not have a FEM yet"
            )

    def init_parameter(self, name, name_port):
        """This function initializes the parameter name in the FEM of the port name_port of the dphs and adds it to the getfem `Model`

        Args:
            name (str): the name of the parameter as defined with add_parameter()
            name_port (str): the name of the `port` where the parameter belongs
        """
        
        # Check if the `port` exists
        try:
            assert name_port in self.ports.keys()
        except AssertionError as err:
            logging.error(
                f"The port {name_port} does not exist"
            )
            raise err
        
        # Check if the FEM of the `port` has been set
        try:
            assert self.ports[name_port].get_is_set()
        except AssertionError as err:
            logging.error(
                "The FEM of the port must be set before initializing the parameter"
            )
            raise err
        
        # Evaluate the expression in the FEM of the `port`
        evaluation = self.ports[name_port].init_parameter(
            name, self.ports[name_port].get_parameter(name).get_expression()
        )
        
        if rank==0:
            logging.info(
                f"{name} has been set to {self.ports[name_port].get_parameter(name).get_expression()} in port: {name_port}"
            )

        # Check the size of the parameter according to its kind
        kind = self.ports[name_port].get_parameter(name).get_kind()
        if kind == "tensor-field":
            sizes = evaluation.shape[0]
        elif kind == "scalar-field":
            sizes = 1
        else:
            logging.error(
                f"Unknown kind of parameter {kind}"
            )
            raise ValueError
        
        # Add the parameter in getfem
        self.gf_model.add_initialized_fem_data(
            name, self.ports[name_port].get_fem(), evaluation, sizes=sizes
        )
        
        if rank==0:
            logging.info(
                f"{name} has been initialized with the FEM of port: {name_port}"
            )

    def set_initial_value(self, name_variable, expression):
        """This function sets the initial value of the variable `name_variable` of the dphs from an expression

        Args:
            name_variable (str): the name of the variable to set
            expression (str): the expression of the function to use
        """

        if name_variable in self.initial_value_set.keys():
            if self.initial_value_set[name_variable]:
                logging.warning(
                    f"defining again the initial value of {name_variable}"
                )
        else:
            logging.error(
                f"Variable {name_variable} does not need initialization"
            )
            raise ValueError

        evaluation = self.gf_model.mesh_fem_of_variable(name_variable).eval(
            expression, globals(), locals()
        )

        initial_value = None
        for _, port in self.ports.items():
            if name_variable == port.get_flow() or name_variable == port.get_effort():
                if port.get_region() == None:
                    initial_value = evaluation
                else:
                    nb_dofs_total = port.get_fem().nbdof()
                    dofs_on_region = port.get_fem().basic_dof_on_region(
                        port.get_region()
                    )
                    qdim = port.get_fem().qdim()
                    size = dofs_on_region.shape[0]
                    shape = (qdim, int(size / qdim))
                    evaluation_long = np.reshape(evaluation, (nb_dofs_total,))
                    initial_value_long = evaluation_long[dofs_on_region]
                    initial_value = np.reshape(initial_value_long, shape)
                continue

        try:
            assert not initial_value.all() == None
        except AssertionError as err:
            logging.error(
                f"{name_variable} can not be found in ports"
            )
            raise err

        self.set_from_vector(name_variable, initial_value)
        self.initial_value_set[name_variable] = True
        
        if rank==0:
            logging.info(
                f"{name_variable} has been initialized with: {expression}"
            )

    def set_from_vector(self, name_variable, x):
        """This function sets the value of the variable `name_variable` in the getfem `Model` from a numpy vector

        Args:
            name_variable (str): the name of the variable to set
            x (numpy array): the vector of values
        """

        self.gf_model.set_variable(name_variable, x)
        
        if rank==0:
            logging.info(
                f"{name_variable} has been set"
            )
            logging.debug(
                f"{name_variable} new value is: {x}"
            )

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
            if brick.get_linear() and not brick.get_position() == "source":
                id_brick = self.gf_model.add_linear_term(
                    self.domain._int_method[mesh_id], form, region
                )
                s = "Linear form "

            elif brick.get_position() == "source":
                id_brick = self.gf_model.add_source_term_brick(
                    self.domain._int_method[mesh_id], name_brick[:-7], form, region
                )
                s = "Source form "

            else:
                id_brick = self.gf_model.add_nonlinear_term(
                    self.domain._int_method[mesh_id], form, region
                )
                s = "Non-linear form "

            if rank==0:
                logging.info(
                    f"{s} '{self.bricks[name_brick].get_form()}' has been added as: {position} relation on region: {region} of mesh: {mesh_id}"
                )

            self.bricks[name_brick].add_id_brick_to_list(id_brick)
            
            if rank==0:
                logging.debug(
                    f"brick IDs: {id_brick} has been added to brick: {name_brick}"
                )

    def add_control_port(self, control_port: Control_Port):
        """This function adds a control `port` to the dphs

        Args:
            control_port (Control_Port): the control port.
        """

        self.controls[control_port.get_name()] = control_port
        self.add_port(control_port)

    def set_control(self, name, expression):
        """This function applies a source term `expression` to the control port `name`

        Args:
            name (str): the name of the port
            expression (str): the expression of the source term
        """

        kind = self.controls[name].get_kind()
        if kind == "scalar-field":
            times = "*"
        elif kind == "vector-field":
            times = "."
        elif kind == "tensor-field":
            times = ":"
        else:
            logging.error(
                f"Unknown kind of control {kind}"
            )
            raise ValueError
        
        # form of the mass matrix for the control variable (flow side type => adds a minus)
        u = self.controls[name].get_name_control()
        mass_form = "-" + u + times + "Test_" + u

        self.add_brick(
            Brick(
                "M_" + u,
                mass_form,
                [self.controls[name].get_region()],
                linear=True,
                dt=False,
                position="constitutive",
                mesh_id=self.controls[name].get_mesh_id(),
            )
        )

        # Add the source brick
        self.add_brick(
            Brick(
                u + "_source",
                expression,
                [self.controls[name].get_region()],
                linear=False,
                dt=False,
                position="source",
                mesh_id=self.controls[name].get_mesh_id(),
            )
        )

        self.controls[name].get_is_set()
        if rank==0:
            logging.info(
                f"Control port: {self.controls[name].get_name()} has been set to: {u} = {expression} on region: {self.controls[name].get_region()} of mesh: {self.controls[name].get_mesh_id()}"
            )

    def assemble_mass(self):
        """This function performs the assembly of the bricks dt=True and linear=True and set the PETSc.Mat attribute `mass`"""

        for _, brick in self.bricks.items():
            # Enable the bricks that are dynamical and linear
            if brick.get_dt() and brick.get_linear():
                brick.enable_id_bricks(self.gf_model)

        self.gf_model.assembly(option="build_matrix")
        rows = self.mass.getOwnershipRange()
        cols = self.mass.getOwnershipRangeColumn()
        extract_gmm_to_petsc(
            [rows[0], rows[1]], [cols[0], cols[1]], self.gf_model.tangent_matrix(), self.mass
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
        rows = self.stiffness.getOwnershipRange()
        cols = self.stiffness.getOwnershipRangeColumn()
        extract_gmm_to_petsc(
            [rows[0], rows[1]], [cols[0], cols[1]], self.gf_model.tangent_matrix(), self.stiffness
        )
        
        for _, brick in self.bricks.items():
            # Disable again all bricks previously enabled
            if not brick.get_dt() and brick.get_linear():
                brick.disable_id_bricks(self.gf_model)

    def assemble_rhs(self):
        """This function performs the assembly of the rhs and set the PETSc.Vec attribute `rhs`"""

        for _, brick in self.bricks.items():
            # Enable the bricks that are in 'source' position
            if brick.get_position() == "source" or brick.get_explicit():
                brick.enable_id_bricks(self.gf_model)
            # And the non-linear ones (and not dt of course)
            if not brick.get_dt() and not brick.get_linear():
                brick.enable_id_bricks(self.gf_model)

        self.gf_model.assembly(option="build_rhs")
        size = self.gf_model.nbdof()
        self.rhs.setValues(range(size), self.gf_model.rhs(),
                           addv=PETSc.InsertMode.INSERT_VALUES)
        self.rhs.assemble()
        
        for _, brick in self.bricks.items():
            # Disable again all bricks previously enabled
            if brick.get_position() == "source" or brick.get_explicit():
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
        rows = self.nl_mass.getOwnershipRange()
        cols = self.nl_mass.getOwnershipRangeColumn()
        extract_gmm_to_petsc(
            [rows[0], rows[1]], [cols[0], cols[1]], self.gf_model.tangent_matrix(), self.nl_mass
        )
        self.tangent_mass = self.mass.copy()
        self.tangent_mass.axpy(1,self.nl_mass)

        for _, brick in self.bricks.items():
            # Disable again all bricks previously enabled
            if brick.get_dt() and not brick.get_linear():
                brick.disable_id_bricks(self.gf_model)

    def assemble_nl_stiffness(self):
        """This function performs the assembly of the bricks dt=False and linear=False and set the PETSc.Mat attribute `nl_stiffness`"""

        for _, brick in self.bricks.items():
            # Enable the bricks that are non-dynamical and non-linear
            if not brick.get_dt() and not brick.get_linear() and not brick.get_explicit():
                brick.enable_id_bricks(self.gf_model)

        self.gf_model.assembly(option="build_matrix")
        rows = self.nl_stiffness.getOwnershipRange()
        cols = self.nl_stiffness.getOwnershipRangeColumn()
        extract_gmm_to_petsc(
            [rows[0], rows[1]], [cols[0], cols[1]], self.gf_model.tangent_matrix(), self.nl_stiffness
        )
        self.tangent_stiffness = self.stiffness.copy()
        self.tangent_stiffness.axpy(1,self.nl_stiffness)

        for _, brick in self.bricks.items():
            # Disable again all bricks previously enabled
            if not brick.get_dt() and not brick.get_linear() and not brick.get_explicit():
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
        IFunction for the time-resolution of the dphs with PETSc TS fully implicit integrator

        Args:
            TS (PETSc TS): the PETSc TS object handling time-resolution
            t (float): time parameter
            z (PETSc Vec): the state
            zd (PETSc Vec): the time-derivative of the state
            F (PETSc Vec): the rhs vector
        """

        self.assemble_nl_mass()
        self.assemble_nl_stiffness()
        
        self.tangent_mass.mult(zd, self.buffer)
        self.tangent_stiffness.multAdd(z, self.buffer, F)

        self.assemble_rhs()
        F.axpy(1, self.rhs)

    def IJacobian(self, TS, t, z, zd, sig, A, P):
        """
        IJacobian for the time-resolution of the dphs with PETSc TS fully implicit integrator

        Args:
            TS (PETSc TS): the PETSc TS object handling time-resolution
            t (float): time parameter
            z (PETSc Vec): the state
            zd (PETSc Vec): the time-derivative of the state
            sig (float): a shift-parameter, depends on dt
            A (PETSc Mat): the jacobian matrix
            P:(PETSc Mat) the jacobian matrix to use for pre-conditionning
        """

        self.assemble_nl_mass()
        self.assemble_nl_stiffness()

        self.tangent_stiffness.copy(P)
        P.axpy(sig, self.tangent_mass)

        if A != P:
            if rank==0:
                logging.info(
                    "Operator different from preconditioning"
                )
            A.assemble()
    
    def allocate_memory(self):
        """Pre-allocate memory for matrices and vectors"""
        
        self.mass.setSizes(self.gf_model.nbdof())
        self.mass.setOption(PETSc.Mat.Option.FORCE_DIAGONAL_ENTRIES, True)
        self.mass.setType("aij")
        self.mass.setUp()
        self.nl_mass.setSizes(self.gf_model.nbdof())
        self.nl_mass.setOption(PETSc.Mat.Option.FORCE_DIAGONAL_ENTRIES, True)
        self.nl_mass.setType("aij")
        self.nl_mass.setUp()
        self.tangent_mass.setSizes(self.gf_model.nbdof())
        self.tangent_mass.setOption(PETSc.Mat.Option.FORCE_DIAGONAL_ENTRIES, True)
        self.tangent_mass.setType("aij")
        self.tangent_mass.setUp()
        self.stiffness.setSizes(self.gf_model.nbdof())
        self.stiffness.setOption(PETSc.Mat.Option.FORCE_DIAGONAL_ENTRIES, True)
        self.stiffness.setType("aij")
        self.stiffness.setUp()
        self.nl_stiffness.setSizes(self.gf_model.nbdof())
        self.nl_stiffness.setOption(PETSc.Mat.Option.FORCE_DIAGONAL_ENTRIES, True)
        self.nl_stiffness.setType("aij")
        self.nl_stiffness.setUp()
        self.tangent_stiffness.setSizes(self.gf_model.nbdof())
        self.tangent_stiffness.setOption(PETSc.Mat.Option.FORCE_DIAGONAL_ENTRIES, True)
        self.tangent_stiffness.setType("aij")
        self.tangent_stiffness.setUp()
        
        self.rhs = self.nl_stiffness.createVecRight()
        self.F = self.nl_stiffness.createVecRight()
        self.buffer = self.nl_stiffness.createVecRight()
        
    def get_cleared_TS_options(self):
        """To ensure a safe database for the PETSc TS environment"""

        self.time_scheme = PETSc.Options()
        for key in self.time_scheme.getAll():
            self.time_scheme.delValue(key)

    def set_time_scheme(self, **kwargs):
        """Allows an easy setting of the PETSc TS environment

        Args:
            \**kwargs: PETSc TS options and more (see examples)
        """

        for key, value in kwargs.items():
            self.time_scheme[key] = value
        
        if not self.time_scheme.hasName("ts_equation_type"):
            self.time_scheme["ts_equation_type"] = PETSc.TS.EquationType.DAE_IMPLICIT_INDEX2
        
        if not self.time_scheme.hasName("ts_type") and not self.time_scheme.hasName("ts_ssp"):
            self.time_scheme["ts_type"] = "bdf"
            self.time_scheme["ts_bdf_order"] = 2
            
        if not self.time_scheme.hasName("ksp_type"):
            self.time_scheme["ksp_type"] = "gmres"
            
        if not self.time_scheme.hasName("pc_type"):
            self.time_scheme["pc_type"] = "lu"
            
        if not self.time_scheme.hasName("pc_factor_mat_solver_type"):
            self.time_scheme["pc_factor_mat_solver_type"] = "mumps"
            
        if not self.time_scheme.hasName("t_0"):
            self.time_scheme["t_0"] = 0.
            
        if not self.time_scheme.hasName("t_f"):
            self.time_scheme["t_f"] = 1.
            
        if not self.time_scheme.hasName("dt"):
            self.time_scheme["dt"] = 0.01
            
        if not self.time_scheme.hasName("dt_save"):
            self.time_scheme["dt_save"] = 0.01
            
        if not self.time_scheme.hasName("ts_adapt_dt_min"):
            self.time_scheme["ts_adapt_dt_min"] = 0.0001
            
        if not self.time_scheme.hasName("adapt_dt_max"):
            self.time_scheme["ts_adapt_dt_max"] = self.time_scheme["dt_save"]
            
        if not self.time_scheme.hasName("ts_max_snes_failures"):
            self.time_scheme["ts_max_snes_failures"] = -1
            
        if not self.time_scheme.hasName("ts_max_reject"):
            self.time_scheme["max_reject"] = -1
        
        if not self.time_scheme.hasName("init_step"):
            self.time_scheme["init_step"] = True
        
        self.time_scheme["isset"] = True

    def exclude_algebraic_var_from_lte(self, TS):
        """Exclude the algebraic variable from the local error troncature in the time-resolution

        Args:
            TS (PETSc TS): the PETSc TS object handling time-resolution
        """

        atol = TS.getTolerances()[1]
        atol_v = atol * np.ones((self.gf_model.nbdof(),))
        id_alg = np.array([], dtype=int)
        # We add indices of all algebraic ports
        for _, port in self.ports.items():
            if port.get_algebraic():
                I = self.gf_model.interval_of_variable(port.get_flow())
                id_alg = np.concatenate(
                    (id_alg, np.arange(I[0], I[0] + I[1], dtype=int))
                )
                # If it is not substituted, it also has the effort part
                if not port.get_substituted():
                    I = self.gf_model.interval_of_variable(port.get_effort())
                    id_alg = np.concatenate(
                        (id_alg, np.arange(I[0], I[0] + I[1], dtype=int))
                    )
            # Else on dynamical ports
            else:
                # If costate are not substituted, we also add them
                if not port.get_substituted():
                    I = self.gf_model.interval_of_variable(port.get_effort())
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
        
    def event(self, TS, t, z, fvalue):
        """Check if the time step is not too small"""
        
        fvalue[0] = TS.getTimeStep() - float(self.time_scheme['ts_adapt_dt_min'])

    def postevent(self, TS, event, t, z, forward):
        """If the time step is too small, ask for the end of the simulation"""
        
        self.stop_TS = True
        if rank==0:
            logging.info(
                "\n *** AUTOMATIC STOP -- LAST ITERATION *** \n"
            )
            logging.info(
                f"dt reached the limit {float(self.time_scheme['ts_adapt_dt_min']):8g}.\n"
            )
        
    def solve(self):
        """Perform the time-resolution of the dphs thanks to PETSc TS

        The options database is set in the `time_scheme` attribute
        """

        assert self.time_scheme[
            "isset"
        ], "Time scheme must be defined before time-resolution"
        for name_variable in self.initial_value_set.keys():
            try:
                assert self.initial_value_set[name_variable]
            except AssertionError as err:
                logging.error(
                    f"Initial value of {name_variable} must be defined before time-resolution"
                )
                raise err
        for name_control in self.controls.keys():
            try:
                assert self.controls[name_control].get_is_set()
            except AssertionError as err:
                logging.error(
                    f"Control {name_control} must be defined before time-resolution"
                )
                raise err

        # Initialize time scheme options (without override)
        self.set_time_scheme()

        # Load time scheme options to PETSc
        self.time_scheme = PETSc.Options()
        self.allocate_memory()

        # Initialization
        self.disable_all_bricks()  # enabling is done in assembly accurately
        self.gf_model.set_time(float(self.time_scheme["t_0"]))
        self.gf_model.first_iter()

        # Assembly
        self.assemble_mass()
        self.assemble_nl_mass()

        self.assemble_stiffness()
        self.assemble_nl_stiffness()

        self.assemble_rhs()

        # PETSc TS Jacobian and residual
        self.J = self.tangent_stiffness.duplicate()
        self.J.assemble()

        self.F = self.rhs.duplicate()
        self.F.assemble()
        
        self.buffer = self.rhs.duplicate()
        self.buffer.assemble()

        self.ts_start = time.perf_counter()
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
        TS.setEventHandler([0], [True], self.event, self.postevent)
        TS.setEventTolerances(1e-16, vtol=[1e-19])
        TS.setIFunction(self.IFunction, self.F)
        TS.setIJacobian(self.IJacobian, self.J)
        TS.setTime(float(self.time_scheme["t_0"]))
        TS.setMaxTime(float(self.time_scheme["t_f"]))
        TS.setTimeStep(float(self.time_scheme["dt"]))
        TS.setExactFinalTime(PETSc.TS.ExactFinalTime.INTERPOLATE)
        TS.setFromOptions()
        self.exclude_algebraic_var_from_lte(TS)
        
        InitVec = self.tangent_stiffness.createVecRight()
        InitVec.setValues(range(self.gf_model.nbdof()), self.gf_model.from_variables(),
                          addv=PETSc.InsertMode.INSERT_VALUES)
        InitVec.assemble()
        TS.solve(InitVec)

        if rank==0:
            logging.info(
                f"Elapsed time: {time.perf_counter() - self.ts_start:1.4g}s"
            )
            logging.info(
                f"Steps: {TS.getStepNumber()} ({TS.getStepRejections()} rejected, {TS.getSNESFailures()} Nonlinear solver failures)"
            )
            logging.info(
                f"Nonlinear iterations: {TS.getSNESIterations()}, Linear iterations: {TS.getKSPIterations()}"
            )
        TS.reset()

        TS.destroy()
        gc.collect()

        self.enable_all_bricks()
        self.solve_done = True

    def init_step(self):
        """Perform a first initial step with a pseudo bdf

        It needs set_time_scheme(init_step=True)
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
        TS.setTimeStep(float(self.time_scheme['dt'])/100.)
        TS.setMaxSteps(2)
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

        TS.setFromOptions()
        self.exclude_algebraic_var_from_lte(TS)
        
        if rank==0:
            logging.info(
                "Perform an initial step using pseudo bdf scheme for initial value consistency"
            )
        
        InitVec = self.tangent_stiffness.createVecRight()
        InitVec.setValues(range(self.gf_model.nbdof()), self.gf_model.from_variables(),
                          addv=PETSc.InsertMode.INSERT_VALUES)
        InitVec.assemble()
        TS.solve(InitVec)
        
        if rank==0:
            logging.info(
                f"Initialisation done in {time.perf_counter() - self.ts_start:1.4g}s"
            )
            
        TS.reset()
        TS.destroy()

        # Delete options for initial step
        self.time_scheme.delValue("ts_type")

        # Re-load previously saved options
        if "ts_type" in saved_options.keys():
            self.time_scheme["ts_type"] = saved_options["ts_type"]
        if "ts_ssp" in saved_options.keys():
            self.time_scheme["ts_ssp"] = saved_options["ts_ssp"]

    def monitor(self, TS, i, t, z, dt_save=1, t_0=0.0, initial_step=False):
        """Monitor to use during iterations of time-integration at each successful time step

        Args:
            TS (PETSc TS): the PETSc TS object handling time-resolution
            i (int): the iteration in the time-resolution
            t (float): time parameter
            z (PETSc Vec): the state
            dt_save (float): save the solution each dt_save s (different from the time-step `dt` used for resolution)
            t_0 (float): the initial time
            initial_step (bool): `True` if this is the initial consistency step (default=`False`)
        """

        if not initial_step:
            if self.solution['t']:
                next_saved_at_t = self.solution['t'][-1] + dt_save
            else:
                next_saved_at_t = t_0
            if (next_saved_at_t - t <= 0) or (i==-1) or self.stop_TS:
                if rank==0:
                    sys.stdout.write(
                        f"\ri={i:8d} t={t:8g} * ({int(time.perf_counter()-self.ts_start)}s)   dt={float(TS.getTimeStep()):8g}        \n"
                    )
                    sys.stdout.flush()
                self.solution["t"].append(t)
                self.solution["z"].append(z.copy())
            else:
                if rank==0:
                    sys.stdout.write(
                        f"\ri={i:8d} t={t:8g}   ({int(time.perf_counter()-self.ts_start)}s)   dt={float(TS.getTimeStep()):8g}        "
                    )
                    sys.stdout.flush()
        
        self.gf_model.set_time(TS.getTime()) # Update the time t in the getfem `Model`
        self.gf_model.to_variables(TS.getSolution().array) # Update the state in the getfem `Model`
        self.gf_model.next_iter() # Says to getfem that we go to the next iteration, not sure if needed
        
        PETSc.garbage_cleanup() # solution saved in getfem, cleaning is safe

    def compute_Hamiltonian(self):
        """Compute each `term` constituting the Hamiltonian"""

        try:
            assert self.solve_done
        except AssertionError as err:
            logging.error(
                "System has not been solved yet, Hamiltonian can not be computed"
            )
            raise err
        
        if rank==0:
            logging.info(
                "Start computing the Hamiltonian"
            )
        start = time.perf_counter()
        for t, _ in enumerate(self.solution["t"]):
            self.gf_model.to_variables(self.solution["z"][t])
            for _, term in enumerate(self.hamiltonian):
                term_value_at_t = 0.0
                for region in term.get_regions():
                    term_value_at_t += gf.asm(
                        "generic",
                        self.domain._int_method[term.get_mesh_id()],
                        0,
                        term.get_expression(),
                        region,
                        self.gf_model,
                    )
                term.set_value(term_value_at_t)

        self.hamiltonian.set_is_computed()
        
        if rank==0:
            logging.info(
                f"Hamiltonian has been computed in {str(time.perf_counter() - start)} s"
            )

    def plot_Hamiltonian(self, with_powers=True, save_figure=False, filename="hamiltonian.png"):
        """Plot each term constituting the Hamiltonian and the Hamiltonian

        May include the `power terms`, i.e. the sum over [t_0, t_f] of the flow/effort product of algebraic ports

        Args:
            with_powers (bool): if `True` (default), the plot will also contains the power of each algebraic ports
            save_figure (bool): if 'True' (defaults: False), save the plot
            filename (str): the name of the file where the plot is saved (defaults: `hamiltonian.png`)
        """

        if not self.hamiltonian.get_is_computed():
            self.compute_Hamiltonian()

        t = np.array(self.solution["t"])
        HamTot = np.zeros(t.size)
        if rank==0:
            fig = plt.figure(figsize=[8, 5])
            ax = fig.add_subplot(111)

        for _, term in enumerate(self.hamiltonian):
            values = np.array(term.get_values())
            HamTot += values
            if rank==0:
                ax.plot(t, values, label=term.get_description())

        if len(self.hamiltonian) > 1 and rank==0:
            ax.plot(t, HamTot, label=self.hamiltonian.get_name())

        if with_powers:
            self.plot_powers(ax, HamTot=HamTot)

        if rank==0:
            ax.legend()
            ax.grid(axis="both")
            ax.set_xlabel("time t")
            ax.set_ylabel("Hamiltonian terms")
            ax.set_title("Evolution of Hamiltonian terms")
            if save_figure:
                plt.savefig(os.path.join(module_path, "png", filename), dpi=300)
            plt.show()

    def compute_powers(self):
        """Compute each power associated to each algebraic port if it is not substituted (because of the parameter-dependency if it is)"""

        try:
            assert self.solve_done
        except AssertionError as err:
            logging.error(
                "System has not been solved yet, powers can not be computed"
            )
            raise err

        if rank==0:
            logging.info(
                "Start computing the powers (substituted ports are not automated)"
            )
        start = time.perf_counter()
        for _, port in self.ports.items():
            if port.get_algebraic() and not port.get_substituted():
                if port.get_region() == None:
                    region = -1
                else:
                    region = port.get_region()

                if port.get_kind() == "scalar-field":
                    times = "*"
                elif port.get_kind() == "vector-field":
                    times = "."
                elif port.get_kind() == "tensor-field":
                    times = ":"
                else:
                    logging.error(
                        f"Unknown kind {port.get_kind()} for control port"
                    )
                    raise ValueError

                form = port.get_flow() + times + port.get_effort()

                self.powers[port.get_name()] = []
                for t, _ in enumerate(self.solution["t"]):
                    self.gf_model.to_variables(self.solution["z"][t])

                    power_value_at_t = gf.asm(
                        "generic",
                        self.domain._int_method[port.get_mesh_id()],
                        0,
                        form,
                        region,
                        self.gf_model,
                    )
                    self.powers[port.get_name()].append(power_value_at_t)

        self.powers_computed = True
        
        if rank==0:
            logging.info(
                f"Powers have been computed in {time.perf_counter() - start} s"
            )

    def plot_powers(self, ax=None, HamTot=None):
        """Plot each power associated to each algebraic port

        The time integration for visual comparison with the Hamiltonian is done using a midpoint method

        If HamTot is provided, a `Balance` showing structure-preserving property is shown: must be constant on the plot

        ax (Matplotlib axis): the gca of matplotlib when this function is called from plot_Hamiltonian()
        HamTot (Numpy array): the values of the Hamiltonian over time
        """

        if not self.powers_computed:
            self.compute_powers()

        need_to_set_figure = False

        t = np.array(self.solution["t"])
        if ax is None and rank==0:
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
                power[k] = power[k - 1] + 0.5 * (t[k] - t[k - 1]) * (values[k - 1] + values[k])
            if rank==0:
                ax.plot(t, power, label=name_port)
            if HamTot is not None:
                if not self.ports[name_port].get_dissipative():
                    SP_balance -= power
                else:
                    SP_balance += power

        if HamTot is not None:
            # Check if the balance makes sense: should not have a port which is both algebraic and substituted
            check_makes_sense = True
            for _, port in self.ports.items():
                if port.get_algebraic() and port.get_substituted():
                    check_makes_sense = False
            if check_makes_sense and rank==0:
                ax.plot(t, SP_balance, "--", label="Balance")

        if need_to_set_figure and rank==0:
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
        
        # TODO: improve a lot!!!

        if not t == None:
            self.gf_model.set_time(t)
        else:
            self.gf_model.set_time(0.0)

        if not state == None:
            self.gf_model.to_variables(state)
        else:
            self.gf_model.to_variables(np.zeros((1,self.gf_model.nbdof())))

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
                    M += extract_gmm_to_scipy([0, size], [0, size], matrix)
            if brick.get_position() == "effort":
                for region in brick.get_regions():
                    matrix = gf.asm_generic(
                        self.domain._int_method[brick.get_mesh_id()],
                        2,
                        brick.get_form(),
                        region,
                        self.gf_model,
                    )
                    J += extract_gmm_to_scipy([0, size], [0, size], matrix)

        pl.spy(M, markersize=0.2)
        pl.show()
        pl.spy(J, markersize=0.2)
        pl.show()

    def export_matrices(self, t=None, state=None, path=None, to="matlab"):
        """
        TODO:
        """
        
        # TODO

        if path == None:
            path = module_path

    def export_to_pv(self, name_variable, path=None, t="All"):
        """Export the solution to .vtu file(s) (for ParaView), with associated .pvd if t='All'

        name_variable (str): the variable to export
        path (str): the path of the output file (default: in the `outputs` folder of scrimp)
        t (Numpy array): the time values of extraction (default: `All` the stored times), `All`, `Init`, `Final`
        """

        try:
            assert self.solve_done
        except AssertionError as err:
            logging.error(
                "The system has not been solved yet"
            )
            raise err

        if path is not None:
            path = os.path.join(path, name_variable)
        else:
            path = os.path.join(module_path, name_variable)

        if not os.path.exists(path):
            os.makedirs(path)

        # Get the dofs of name_variable
        I = self.gf_model.interval_of_variable(name_variable)
        J = range(I[0], I[0] + I[1])

        if t == "All":
            if rank==0:
                logging.info(
                    f"Export all time values of {name_variable} is starting...\n"
                )
            with open(os.path.join(path, name_variable + ".pvd"), "w") as pvd_file:
                pvd_file.write('<?xml version="1.0"?>\n')
                pvd_file.write('<VTKFile type="Collection" version="0.1">\n')
                pvd_file.write("  <Collection>\n")
                for k, _ in enumerate(self.solution["t"]):
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
                    if rank==0:
                        sys.stdout.write("\rExport time %f" % self.solution["t"][k])
                        sys.stdout.flush()
                pvd_file.write("  </Collection>\n")
                pvd_file.write("</VTKFile>")
            pvd_file.close()
            if rank==0:
                logging.info(
                    "\rExport is done              \n"
                )
        elif t == "Init":
            if rank==0:
                logging.info(
                    f"Export initial value of {name_variable}"
                )
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
            if rank==0:
                logging.info(
                    f"Export final value of {name_variable}"
                )
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
            logging.error(
                f"t must be `All`, `Init` or `Final`. Unknown value: {t}"
            )
            raise ValueError

    def display(self, verbose=2):
        """A method giving infos about the dphs

        verbose (int): the level of verbosity (defaults: 2)
        """

        # TODO: improve presentation.

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

            for _, port in self.ports.items():
                print(port)
