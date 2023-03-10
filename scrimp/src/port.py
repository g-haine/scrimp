import getfem as gf


class Parameter:
    """This class describes the Parameter for a Port."""

    def __init__(
        self, name: str, description: str, kind: str, expression: str, name_port: str
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
        mesh_id: int,
        algebraic: bool,
        substituted: bool,
        region: int,
    ):
        """Constructor of a `port` of a discrete port Hmiltonian system (dpHs).

        Args:
            name (str): the name of the port
            flow (str): the name of the flow variable
            effort (str): the name of the effort variable
            kind (str): the type of the variables, must be `scalar-field` or `vector-field`
            mesh_id (int): the id of the mesh where the variables belong
            algebraic (bool): if `False`, the flow variable will be derivated in time at resolution
            substituted (bool): if `True`, the constitutive relation is substituted and there is only a getfem variable for the effort
            region (int): the int identifying the region in mesh_id where the port belong, useful for boundary ports
        """

        self.__isSet = False  #: A boolean to check if the domain has been setted
        self.__name = name  #: The name of the `port`
        self.__flow = flow  #: The name of the flow variable
        self.__effort = effort  #: The name of the effort variable
        self.__kind = kind  #: The type of the variables (e.g. `scalar-field`)
        self.__mesh_id = mesh_id  #: The id of the mesh where the variables belong
        self.__algebraic = algebraic  #: If `True`, the equation associated to this port is algebraic, otherwise dynamic and the flow is derivated in time at resoltuion
        self.__substituted = substituted  #: If `True, the getfem `Model` will only have an unknown variable for the effort: the constitutive relation is substituted into the mass matrix on the flow side
        self.__parameters = (
            []
        )  #: A list of parameters acting on the variables of the `port`
        self.__fem = None  #: A getfem Meshfem object to discretize the `port`
        self.__region = region  #: If any, the int of the region of mesh_id where the flow/effort variables belong

    def get_is_set(self) -> bool:
        """This funcion gets the boolean value that indicates wether the port is set or not.


        Returns:
            bool: value that indicates if the port is set.
        """
        return self.__isSet

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

    def set_fem(self, mesh, dim: int, order: int, fem: str):
        """This function sets the Meshfem getfem object defining the finite element method to use to discretize the port.


        Args:
            mesh (Mesh): the mesh where the FE are define
            dim (int): the dim of the FE !TO DO: improve this for automation
            order (int): the order of the FE
            fem (str): the FE to use

        Raises:
            ValueError: Unknown fem for port.
        """
        # TO DO: handle more FE

        if fem == "CG":
            fem_str = "FEM_PK(" + str(mesh.dim()) + "," + str(order) + ")"
        elif fem == "DG":
            fem_str = "FEM_PK_DISCONTINUOUS(" + str(mesh.dim()) + "," + str(order) + ")"
        else:
            raise ValueError(
                "Unknown fem "
                + fem
                + " for port "
                + self.__name
                + "\nUse the gf_model `Model` attribute to set it directly"
            )

        self.__fem = gf.MeshFem(mesh, dim)
        self.__fem.set_fem(gf.Fem(fem_str))

        self.__isSet = True
        print(fem_str, "has been setted for port", self.__name)
        self.__fem.display()

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
        print(f"Parameter with name: {name} does not exit!")

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
            # assert False, f"Insertion parameter not valid expected {type(Parameter)} got {type(parameter)}"
            return False

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
                    "A fem must be setted for port '",
                    self.__name,
                    "' before initialization",
                )
                evaluation = self.__fem.eval(expression, globals(), locals())
                print(
                    "Parameter",
                    name,
                    "has been evaluated with the fem of port '",
                    p.get_name_port(),
                    "', with expression:",
                    expression,
                )
                return evaluation
        assert False, ("Parameter", name, "must be added before intialization")

    def __str__(self) -> str:
        """This function displays all the info of the port.

        Returns:
            str: detailed info of the port
        """

        assert self.__isSet, (
            "Port",
            self.__name,
            "has not been setted yet (a fem is probably missing)",
        )
        self.__fem.display()

        return f"{self.__name}, {self.__flow}, {self.__effort}, {self.__kind}, {str(self.__mesh_id)}, {str(self.__algebraic)}, {self.__parameters}"


if __name__ == "__main__":
    p = Port("name", "flow", "effort", "kind", 0, True, True, 1)
    param = Parameter(
        "name_para", "description", "kind_param", "expression", p.get_name()
    )
    p.add_parameter(p)
    print(p)

# class Port:
#     """
#     A class to handle a `port` of a dpHs
#
#     It is mainly constituted of a flow variable, an effort variable, and a fem
#     """
#
#     def __init__(self, name, flow, effort, kind, mesh_id, algebraic, substituted, region):
#         """
#         Constructor of a `port` of a dpHs
#
#         !TO DO: handle `tensor-field` variable
#
#         :param name: the name of the port
#         :type name: str
#         :param flow: the name of the flow variable
#         :type flow: str
#         :param effort: the name of the effort variable
#         :type effort: str
#         :param kind: the type of the variables, must be `scalar-field` or `vector-field`
#         :type kind: str
#         :param mesh_id: the id of the mesh where the variables belong
#         :type mesh_id: int
#         :param algebraic: if `False`, the flow variable will be derivated in time at resolution
#         :type algebraic: bool
#         :param substituted: if `True`, the constitutive relation is substituted and there is only a getfem variable for the effort
#         :type substituted: bool
#         :param region: the int identifying the region in mesh_id where the port belong, useful for boundary ports
#         :type region: int
#         """
#
#         self.__isSet = False  #: A boolean to check if the domain has been setted
#         self.__name = name  #: The name of the `port`
#         self.__flow = flow  #: The name of the flow variable
#         self.__effort = effort  #: The name of the effort variable
#         self.__kind = kind  #: The type of the variables (e.g. `scalar-field`)
#         self.__mesh_id = mesh_id  #: The id of the mesh where the variables belong
#         self.__algebraic = algebraic  #: If `True`, the equation associated to this port is algebraic, otherwise dynamic and the flow is derivated in time at resoltuion
#         self.__substituted = substituted  #: If `True, the getfem `Model` will only have an unknown variable for the effort: the constitutive relation is substituted into the mass matrix on the flow side
#         self.__parameters = dict()  #: A dict of parameters acting on the variables of the `port`
#         self.__fem = None  #: A getfem Meshfem object to discretize the `port`
#         self.__region = region  #: If any, the int of the region of mesh_id where the flow/effort variables belong
#
#     def set_fem(self, mesh, dim, order, fem):
#         """
#         Set the Meshfem getfem object defining the finite element method to use to discretize the port
#
#         !TO DO: handle more FE
#
#         :param mesh: the mesh where the FE are define
#         :type mesh: getfem Mesh object
#         :param dim: the dim of the FE !TO DO: improve this for automation
#         :type dim: int
#         :param order: the order of the FE
#         :type order: int
#         :param fem: the FE to use
#         :type fem: str
#
#         :return: set a Meshfem object on the mesh where the variables belong
#         """
#
#         if fem == 'CG':
#             fem_str = 'fem_PK(' + str(mesh.dim()) + ',' + str(order) + ')'
#         elif fem == 'DG':
#             fem_str = 'fem_PK_DISCONTINUOUS(' + str(mesh.dim()) + ',' + str(order) + ')'
#         else:
#             raise ValueError(
#                 'Unknown fem ' + fem + ' for port ' + self.__name + '\nUse the gf_model `Model` attribute to set it directly')
#
#         self.__fem = gf.Meshfem(mesh, dim)
#         self.__fem.set_fem(gf.fem(fem_str))
#
#         self.__isSet = True
#         print(fem_str, 'has been setted for port', self.__name)
#         self.__fem.display()
#
#     def add_parameter(self, name, description, kind, expression):
#         """
#         Add a parameter acting on the variables of the `port`
#
#         This parameter is given as a string expression, and may be space-varying (variable x, y and z)
#
#         It may be a scalar-, vector- or tensor-field (e.g. `x*y', `[x, y]`, `[[1, 0],[0, y]]`, etc.)
#
#         :param name: the name of the parameter
#         :type name: str
#         :param description: a physically motivated description (e.g. `mass density`)
#         :type description: str
#         :param kind: the type of the parameter (e.g. `scalar-field`)
#         :type kind: str
#         :param expression: the expression in getfem syntax
#         :type expression: str
#
#         :return: append a parameter to the dict `parameters` of the port
#         """
#
#         self.__parameters[name] = {'description': description,
#                                  'kind': kind,
#                                  'expression': expression,
#                                  'name_port': self.__name
#                                  }
#         print('A parameter', name, ', describing \'', description, '\', of type', kind, 'has been added to port \'',
#               self.__name, '\'')
#
#     def init_parameter(self, name, expression):
#         """
#         Set the parameter by initialization in the FE basis
#
#         This method should not be called explicitly. The init_parameter() method of the parent `dpHs` should be called instead to ensure the parameter to be added to the getfem `Model` object
#
#         :param name: the name of the parameter to initialize
#         :type name: str
#         :param expression: the expression in getfem syntax
#         :type expression: str
#
#         :return: the evaluation of the parameter in the fem of the `port`
#         :rtype: numpy array
#         """
#
#         assert self.__isSet, ('A fem must be setted for port \'', self.__name, '\' before initialization')
#
#         assert name in self.__parameters.keys(), ('Parameter', name, 'must be added before intialization')
#
#         evaluation = self.__fem.eval(expression, globals(), locals())
#         print('Parameter', name, 'has been evaluated with the fem of port \'', self.__parameters[name]['name_port'],
#               '\', with expression:', expression)
#
#         return evaluation
#
#     def display(self):
#         """
#         A method giving infos about the `port`
#         """
#
#         assert self.__isSet, ('Port', self.__name, 'has not been setted yet (a fem is probably missing)')
#
#         print(self.__name, self.__flow, self.__effort, self.__kind, self.__mesh_id, self.__algebraic, self.__parameters)
#         self.__fem.display()
