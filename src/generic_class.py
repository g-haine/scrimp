# -*- coding: utf-8 -*-
"""
Simplified generic class

@author: Ghislain Haine
"""

import numpy as np # Array manipulations
import matplotlib.pyplot as plt # To plot figures
import fenics as fe # Finite elements library
import multiphenics as mpfe # Allow easier use of FEniCS for Mixed Finite Elements spaces
import os

# Can be found in src/modules directory of SCRIMP
import gmsh_utils
import meshio_utils
import fenics_utils
import multiphenics_utils

fe.SubSystemsManager.init_petsc()
from petsc4py import PETSc
import PETSc_utils

# Parent class for both systems and interconnections
class _generic():
    
    def __init__(self):
        """
        Constructor
        """
        
        # Problem variables
        self.dim = 0
        

        # Mesh variables
        self.dmesh = None
        self.meshtag = None
        self.name_subdomains = None
        self.name_domain = None
        self.tag_domain = None
        self.name_boundaries = None
        self.tags_boundaries = None
        self.domain_restriction = None
        self.boundaries_restriction = dict()
        
        
        # Class variables
        self.mesh_done = False
        self.domain_done = False
        
        
    def set_dim(self,dim):
    
        self.dim = dim
        
    def mesh_with_gmsh(self,geofile,parameters=None,refinement=0):
    
        assert self.dim>0, \
            'Dimension must be a positive integer'
        
        assert os.path.exists(geofile), \
            f'{geofile} does not exist'
        
        # Mesh generation
        # /!\ : Dolfin (FEniCS) only fully supports first-order meshes
        gmshfile = geofile[0: -4]+'.msh'
        gmsh_utils.generate_mesh_cli(geofile,gmshfile,self.dim,refinement=refinement,log=1,\
                                     parameters=parameters,order=1,gmshformat=2,binary=True)
                                     
        # Print mesh informations
        gmsh_utils.print_summary(gmshfile)
        
        # Check .msh file
        meshio_utils.print_diagnostic_gmsh(gmshfile)
        
        # Convert to XDMF file
        xdmfiles = meshio_utils.convert_gmsh_to_XDMF(gmshfile,prune_z=True) # prune_z prunes the third coordinate
        
        # Load in Dolfin format (FEniCS)
        self.dmesh = fenics_utils.DolfinMesh.init_from_xdmf(xdmfiles['tetra'*(self.dim==3)+'triangle'*(self.dim==2)],\
                                                            xdmfiles['triangle'*(self.dim==3)+'line'*(self.dim==2)])
        
        # Get the physical tags from .geo file
        self.meshtag = meshio_utils.read_phys_gmsh(xdmfiles['gmsh_physical_entities'])
        
        # Get volumic/surfacic physical entities name
        self.name_subdomains = list(self.meshtag[self.dim-1].keys())
                
        # Get surfacic/lineic physical entities name and tags
        self.tags_boundaries = self.meshtag[self.dim-2]
        self.name_boundaries = list(self.tags_boundaries.keys())
        
        self.mesh_done = True
        
        
    def __call__(self):
        print('Generic class of SCRIMP.')
        print(40*'-', '\n')
        
        
    def set_dim(self,dim):
    
        self.dim = dim
        
        
    def mesh_with_gmsh(self,geofile,parameters=None,refinement=0):
    
        assert self.dim>0, \
            'Dimension must be a positive integer'
        
        assert os.path.exists(geofile), \
            f'{geofile} does not exist'
        
        # Mesh generation
        # /!\ : Dolfin (FEniCS) only fully supports first-order meshes
        gmshfile = geofile[0: -4]+'.msh'
        gmsh_utils.generate_mesh_cli(geofile,gmshfile,self.dim,refinement=refinement,log=1,\
                                     parameters=parameters,order=1,gmshformat=2,binary=True)
                                     
        # Print mesh informations
        gmsh_utils.print_summary(gmshfile)
        
        # Check .msh file
        meshio_utils.print_diagnostic_gmsh(gmshfile)
        
        # Convert to XDMF file
        xdmfiles = meshio_utils.convert_gmsh_to_XDMF(gmshfile,prune_z=True) # prune_z prunes the third coordinate
        
        # Load in Dolfin format (FEniCS)
        self.dmesh = fenics_utils.DolfinMesh.init_from_xdmf(xdmfiles['tetra'*(self.dim==3)+'triangle'*(self.dim==2)],\
                                                            xdmfiles['triangle'*(self.dim==3)+'line'*(self.dim==2)])
        
        # Get the physical tags from .geo file
        self.meshtag = meshio_utils.read_phys_gmsh(xdmfiles['gmsh_physical_entities'])
        
        # Get volumic/surfacic physical entities name
        self.name_subdomains = list(self.meshtag[self.dim-1].keys())
                
        # Get surfacic/lineic physical entities name and tags
        self.tags_boundaries = self.meshtag[self.dim-2]
        self.name_boundaries = list(self.tags_boundaries.keys())
        
        self.mesh_done = True
        
    def set_domain(self,subdomain):
    
        assert self.mesh_done, \
            'Please mesh before setting domain'
    
        assert subdomain in self.name_subdomains, \
            f'{subdomain} is not in {self.name_subdomains}'
        
        self.tag_domain = self.meshtag[self.dim-1][subdomain]
        
        self.name_domain = subdomain
        
        self.domain_restriction = multiphenics_utils.build_MeshRestriction_from_tags_subdomains(self.dmesh.mesh,\
                                                                                                self.dmesh.domains,self.tag_domain)
            
        for gamma in self.name_boundaries:
            self.boundaries_restriction[gamma] = multiphenics_utils.build_MeshRestriction_from_tags(self.dmesh.mesh,\
                                                                                                    self.dmesh.boundaries,\
                                                                                                    self.tags_boundaries[gamma])
        
        self.domain_done = True
        
