# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2023 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""
- file:             domain.py
- authors:          Giuseppe Ferraro, Ghislain Haine
- date:             31 may 2023
- brief:            class for domain object
"""

import petsc4py
import sys

petsc4py.init(sys.argv)
from petsc4py import PETSc

comm = PETSc.COMM_WORLD
rank = comm.getRank()

import scrimp.utils.mesh
import getfem as gf
import logging

class Domain:
    """A class handling meshes and indices for regions

    Lists are used to handle interconnection of pHs, allowing for several meshes in the dpHs.
    """

    def __init__(self, name: str, parameters: dict):
        """Constructor of the `domain` member of a dpHs"""
        
        self._name = name
        self._isSet = False  #: A boolean to check if the domain has been setted
        self._mesh = list()  #: A list of getfem Mesh objects
        self._subdomains = list(
            dict()
        )  #: A list of dict `str`: `int` listing the indices of region of dim n in the getfem mesh
        self._boundaries = list(
            dict()
        )  #: A list of dict `str`: `int` listing the indices of region of dim n-1 in the getfem mesh
        self._dim = list()  #: A list of the dimensions of the getfem meshes
        self._int_method = list()  #: A list of getfem integration method MeshIm objects

        built_in_methods = dir(scrimp.utils.mesh)

        if name in built_in_methods:
            mesh_function = getattr(scrimp.utils.mesh, name)
            gf_mesh = mesh_function(parameters, terminal=0)
            self._mesh.append(gf_mesh[0])
            self._dim.append(gf_mesh[1])
            self._subdomains = [gf_mesh[2]]
            self._boundaries = [gf_mesh[3]]
            self._isSet = True

            self.set_mim_auto()

            if self._isSet and rank==0:
                logging.info(
                    "Domain has been setted"
                )
                self.display()

        # TODO: If not built_in, given from a script 'name.py' or a .geo file with args in the dict 'parameters'
        # should be able to handle several meshes e.g. for interconnections, hence the list type
        else:
            pass

    def set_mim_auto(self):
        """Define the integration method to a default choice"""

        # TODO: Allow for user definition
        
        for k in range(len(self._mesh)):
            if self._dim[k] == 1:
                self._int_method.append(
                    gf.MeshIm(self._mesh[k], gf.Integ("IM_GAUSS1D(5)"))
                )
            elif self._dim[k] == 2:
                self._int_method.append(
                    gf.MeshIm(self._mesh[k], gf.Integ("IM_TRIANGLE(7)"))
                )
            elif self._dim[k] == 3:
                self._int_method.append(
                    gf.MeshIm(self._mesh[k], gf.Integ("IM_TETRAHEDRON(8)"))
                )
            else:
                self._isSet = False
                if rank==0:
                    logging.warning(
                        f"Integration method has to be setted manually on mesh {k} of dimension {self._dim[k]}"
                    )

    def display(self):
        """A method giving infos about the domain"""

        try:
            assert self._isSet
        except AssertionError as err:
            logging.error(
                "Domain has not been setted yet"
            )
            raise err
        
        if rank==0:
            logging.info(
                f"Domain is setted and contains {len(self._mesh)} mesh(es):"
            )
            for k in range(len(self._mesh)):
                logging.info(
                    f"=== on mesh {k} of dim {self._dim[k]}"
                )
                logging.info(
                    f"* Subdomains are: {self._subdomains[k]}"
                )
                logging.info(
                    f"* Boundaries are: {self._boundaries[k]}"
                )
                
                self._mesh[k].display() # GetFEM infos

    def get_name(self) -> str:
        """This function get the name of the domain.

        Returns:
            str: name of the domain.
        """
        
        return self._name

    def get_mesh(self) -> list:
        """This function gets the list of mesh for the domain.

        Returns:
            list: list of mesh for the domain
        """
        
        return self._mesh

    def get_dim(self) -> list:
        """This function gets the list of dimensions for the domain.

        Returns:
            list: list of dimensions for the domain
        """
        
        return self._dim

    def get_subdomains(self) -> list:
        """This function gets the list of subdomains for the domain.

        Returns:
            list: list of the subdomains for the domain
        """
        
        return self._subdomains

    def get_boundaries(self) -> list:
        """This function gets the list of the bounderies for the domain.

        Returns:
            list: list of the boundaries for the domain
        """
        
        return self._boundaries

    def get_isSet(self) -> bool:
        """This function gets the boolean vale indicating wether if a Mesh has been set for the domain or not.

        Returns:
            bool: boolean indicating if a Mesh has been set
        """
        
        return self._isSet
