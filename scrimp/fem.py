# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""
- file:             fem.py
- authors:          Giuseppe Ferraro, Ghislain Haine
- date:             31 may 2023
- brief:            class for fem object
"""

import petsc4py
import sys

petsc4py.init(sys.argv)
from petsc4py import PETSc

comm = PETSc.COMM_WORLD
rank = comm.getRank()

import getfem as gf
import logging


class FEM:
    """Define a finite element for a SCRIMP state or port.

    Families: CG, DG, RT, BDM, NED1 (Nedelec first kind),
    NED2 (Nedelec second kind). The longer names NEDELEC1 and
    NEDELEC2 are accepted aliases. Here 'order' is the polynomial
    index k, **not** the Nedelec family number. NED2 supports only
    k=1,2 in the pinned SCRIMP GetFEM 5.5 patch.
    Passing a negative order with an explicit GetFEM descriptor
    remains supported for expert/custom elements.
    """

    def __init__(self, name, order, FEM="CG") -> None:
        self.__name = name
        self.__order = order
        self.__type = FEM
        self.__fem = None
        self.__mesh = None
        self.__dim = None
        self.__isSet = False

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
            str: type of the FEM
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

    def get_isSet(self) -> bool:
        """This function gets the flag to know if the FEM are set in getfem.

        Returns:
            bool: the flag to assert setting in getfem
        """

        return self.__isSet

    def get_fem(self):
        """This function returns the the FEM of getfem.

        Returns:
            gf.MeshFem: the FEM
        """

        return self.__fem

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

    def set_fem(self):
        """This function sets the Meshfem getfem object defining the finite element method to use to discretize the port.
        """

        try:
            assert self.get_mesh() is not None
        except AssertionError as err:
            if rank == 0:
                logging.error("A mesh must be set before adding a FEM.")
            raise err

        known = True
        if self.get_type() == "CG":
            fem_str = (
                "FEM_PK("
                + str(self.get_mesh().dim())
                + ","
                + str(self.get_order())
                + ")"
            )
        elif self.get_type() == "DG":
            fem_str = (
                "FEM_PK_DISCONTINUOUS("
                + str(self.get_mesh().dim())
                + ","
                + str(self.get_order())
                + ")"
            )
        elif self.get_type() == "RT":
            fem_str = (
                "FEM_RTK("
                + str(self.get_mesh().dim())
                + ","
                + str(self.get_order())
                + ")"
            )
            
        elif self.get_type() == "BDM":
            fem_str = (
                "FEM_BDMK("
                + str(self.get_mesh().dim())
                + ","
                + str(self.get_order())
                + ")"
            )
        elif self.get_type() in ("NED1", "NEDELEC1"):
            # SCRIMP GetFEM 5.5 patch 0003. k=0 delegates to GetFEM's
            # native FEM_NEDELEC(N), while k>=1 uses the patched family.
            dim, k = self.get_mesh().dim(), self.get_order()
            if dim not in (2, 3) or not isinstance(k, int) or k < 0:
                raise ValueError(
                    "NED1/NEDELEC1 expects a triangle/tetrahedron mesh "
                    "in dimension 2/3 and an integer order k>=0"
                )
            fem_str = f"FEM_NEDELEC1({dim},{k})"
        elif self.get_type() in ("NED2", "NEDELEC2"):
            # SCRIMP GetFEM 5.5 patch 0004 is *only* implemented for k=1,2.
            dim, k = self.get_mesh().dim(), self.get_order()
            if dim not in (2, 3) or not isinstance(k, int) or k not in (1, 2):
                raise ValueError(
                    "NED2/NEDELEC2 expects a triangle/tetrahedron mesh "
                    "in dimension 2/3 and an integer order k=1 or k=2"
                )
            fem_str = f"FEM_NEDELEC2({dim},{k})"
        elif self.get_order() < 0: # Hack to access directly to every FEM in GetFEM
            fem_str = self.get_type()
        else:
            logging.warning(
                f"Unknown fem {self.get_type()} in SCRIMP. \nUse the gf_model `Model` attribute to set it directly."
            )
            known = False

        if known:
            try:
                self.__fem = gf.MeshFem(self.get_mesh(), self.get_dim())
                self.__fem.set_fem(gf.Fem(fem_str))
            except AssertionError as err:
                if self.get_type() in ("NED1", "NEDELEC1", "NED2", "NEDELEC2"):
                    # These explicit H(curl) requests MUST NOT silently
                    # become componentwise CG1 when a GetFEM patch is absent.
                    raise RuntimeError(
                        f"Required intrinsic {self.get_type()} order "
                        f"{self.get_order()} unavailable in GetFEM: {fem_str}"
                    ) from err
                if rank == 0:
                    logging.error(f"Unable to set {self.get_type()}, use classical Lagrange of order 1 as default instead.")
                self.__order = 1
                self.__type = "CG"
                fem_str = (
                    "FEM_PK("
                    + str(self.get_mesh().dim())
                    + ","
                    + str(self.get_order())
                    + ")"
                )
                self.__fem = gf.MeshFem(self.get_mesh(), self.get_dim())
                self.__fem.set_fem(gf.Fem(fem_str))

            self.__isSet = True
            if rank == 0:
                logging.info(f"{fem_str} has been set for port {self.__name}")
                self.__fem.display()  # GetFEM infos
