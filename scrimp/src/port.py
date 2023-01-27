class Port:
    """
    A class to handle a `port` of a dpHs

    It is mainly constituted of a flow variable, an effort variable, and a FEM
    """

    def __init__(self, name, flow, effort, kind, mesh_id, algebraic, substituted, region):
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
        self.parameters = dict()  #: A dict of parameters acting on the variables of the `port`
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

        if FEM == 'CG':
            FEM_str = 'FEM_PK(' + str(mesh.dim()) + ',' + str(order) + ')'
        elif FEM == 'DG':
            FEM_str = 'FEM_PK_DISCONTINUOUS(' + str(mesh.dim()) + ',' + str(order) + ')'
        else:
            raise ValueError(
                'Unknown FEM ' + FEM + ' for port ' + self.name + '\nUse the gf_model `Model` attribute to set it directly')

        self.FEM = gf.MeshFem(mesh, dim)
        self.FEM.set_fem(gf.Fem(FEM_str))

        self.isset = True
        print(FEM_str, 'has been setted for port', self.name)
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

        self.parameters[name] = {'description': description,
                                 'kind': kind,
                                 'expression': expression,
                                 'name_port': self.name
                                 }
        print('A parameter', name, ', describing \'', description, '\', of type', kind, 'has been added to port \'',
              self.name, '\'')

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

        assert self.isset, ('A FEM must be setted for port \'', self.name, '\' before initialization')

        assert name in self.parameters.keys(), ('Parameter', name, 'must be added before intialization')

        evaluation = self.FEM.eval(expression, globals(), locals())
        print('Parameter', name, 'has been evaluated with the FEM of port \'', self.parameters[name]['name_port'],
              '\', with expression:', expression)

        return evaluation

    def display(self):
        """
        A method giving infos about the `port`
        """

        assert self.isset, ('Port', self.name, 'has not been setted yet (a FEM is probably missing)')

        print(self.name, self.flow, self.effort, self.kind, self.mesh_id, self.algebraic, self.parameters)
        self.FEM.display()
