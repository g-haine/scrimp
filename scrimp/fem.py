class FEM:
    """This class defines what is a FEM object in SCRIMP."""

    def __init__(self, name, order, FEM="DG") -> None:
        self.__name = name
        self.__order = order
        self.__type = FEM
        self.__mesh = None
        self.__dim = None

    def get_name(self) -> str:
        """This function gets the name of the FEM.

        Returns:
            str: name of the FEM
        """
        return self.__name

    def get_order(self) -> int:
        """This function gets the order for the FEM.

        Returns:
            int: dim of the flow FEM
        """
        return self.__order

    def get_type(self) -> str:
        """This function gets the tyoe of the FEM.

        Returns:
            str: tyoe of the FEM
        """
        return self.__type

    def get_mesh(self):
        """This function gets the mesh of the FEM

        Returns:
            mesh: the mesh of FEM
        """
        return self.__mesh

    def get_dim(self) -> int:
        """This function gets the dimension of the FEM.

        Returns:
            int: the dimension of FEM
        """
        return self.__dim

    def set_mesh(self, mesh):
        """This function sets the Meshfem getfem object FEM object of scrimp.


        Args:
            mesh (Mesh): the mesh where the FE are define
        """
        self.__mesh = mesh

    def set_dim(self, dim: int):
        """This function sets the dimension for the FEM.

        Args:
            dim (int): the dimension fro the FEM
        """
        self.__dim = dim

    def __str__(self) -> str:
        """This function displays all the info of the FEM.

        Returns:
            str: detailed info of the FEM
        """
        return (
            f"{self.__name}, {self.__order}, {self.__type}, {self.__mesh}, {self.__dim}"
        )

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
