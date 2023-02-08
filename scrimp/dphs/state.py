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

    def get_name(self):
        return self._name

    def get_description(self):
        return self._description

    def get_kind(self):
        return self._kind

    def get_region(self):
        return self._region

    def set_costate(self, costate):
        from .costate import CoState

        if self.get_costate() is None and costate is not None:
            assert isinstance(costate, CoState)
            self._costate = costate
        else:
            print("A costate is already present for this state")

    def get_costate(self):
        return self._costate

    def set_port(self, port):
        from .port import Port

        if self.get_port() is None and port is not None:
            assert isinstance(port, Port)
            self._port = port
        else:
            print("A port is already present for this state")

    def get_port(self):
        return self._port

    def get_mesh_id(self):
        return self._mesh_id


# def add_state(self, name, description, kind, region=None, mesh_id=0):
#     """
#     Add a variable to the dict `states` of the dpHs
#
#     !TO DO: handling of `tensor-field` state
#
#     :param name: the name of the state
#     :type name: str
#     :param description: a physically motivated description (e.g. `linear momentum`)
#     :type description: str
#     :param kind: the unknown type, must be `scalar-field` or `vector-field`
#     :type kind: str
#     :param region: the index of the region in mesh_id where the variable belong, useful with multi-domains mesh for interconnected dpHs
#     :type region: int
#     :param mesh_id: a state has to be associated to a unique mesh of the 'domain' attribute (overlap needs the definition of two different states with interface interaction)
#     :type mesh_id: int
#
#     :return: append the new state to the dict `states` of the dpHs
#     """
#
#     self.states[name] = {'description': description,
#                          'kind': kind,
#                          'region': region,
#                          'mesh_id': mesh_id,
#                          'costate': None,
#                          'port': None}
#
#     where = ''
#     if region is not None:
#         where = ', in region numbered ' + str(region)
#     print('A state variable', name, ', describing \'', description, '\', has been initialized as a', kind, 'on mesh',
#           mesh_id, where)
