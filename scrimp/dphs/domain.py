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
            dict())  #: A list of dict `str`: `int` listing the indices of region of dim n in the getfem mesh
        self.boundaries = list(
            dict())  #: A list of dict `str`: `int` listing the indices of region of dim n-1 in the getfem mesh
        self.dim = list()  #: A list of the dimensions of the getfem meshes
        self.int_method = list()  #: A list of getfem integration method MeshIm objects

    def set_mim_auto(self):
        """
        Define the integration method to a default choice
        """

        for k in range(len(self.mesh)):
            if self.dim[k] == 1:
                self.int_method.append(gf.MeshIm(self.mesh[k],
                                                 gf.Integ('IM_GAUSS1D(5)')))
            elif self.dim[k] == 2:
                self.int_method.append(gf.MeshIm(self.mesh[k],
                                                 gf.Integ('IM_TRIANGLE(7)')))
            elif self.dim[k] == 3:
                self.int_method.append(gf.MeshIm(self.mesh[k],
                                                 gf.Integ('IM_TETRAHEDRON(8)')))
            else:
                self.isset = False
                print('Integration method has to be setted manually on mesh', k, ' of dimension', self.dim[k])

    def display(self):
        """
        A method giving infos about the domain
        """

        assert self.isset, ('Domain has not been setted yet')

        print('Domain is setted and contains', len(self.mesh), 'mesh:')
        for k in range(len(self.mesh)):
            print('=== on mesh', k, ' of dim', self.dim[k])
            print('* Subdomains are:', self.subdomains[k])
            print('* Boundaries are:', self.boundaries[k])
            self.mesh[k].display()