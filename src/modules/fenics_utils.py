#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Utility functions related to fenics.

@author: Florian Monteghetti
"""
import os
import copy
import fenics as fe
import dolfin
import meshio


############## Mesh handling

def mesh_load_xdmf(xdmfile,xdmfile_facets):
    """
    Load mesh from XDMF mesh. Dolfin is restricted to first-order tetrahedral
    meshes.

    Parameters
    ----------
    xdmfile : str
        XDMF mesh containing: nodes coordinates, triangle/tetra vertices and
        tags.
        
    xdmfile_facets : str
        XDMF mesh containing: nodes coordinates, line/triangle vertices and
        tags.

    Returns
    -------
    Identical to fenics_utils.mesh_load_dolfinXML().

    """
    
        # Read XDMF meshes
    msh = meshio.read(xdmfile)
    msh_facets = meshio.read(xdmfile_facets)
    
        # Test cell types and order
    def test_mixed_mesh(msh,filename):
        if len(msh.cells_dict)>1:
            raise ValueError(f"Dolfin does not support mixed meshes.\n\t Mesh: {filename}\n\tCell type: {msh.cells_dict.keys()}")    
    test_mixed_mesh(msh,xdmfile)
    test_mixed_mesh(msh_facets,xdmfile_facets)
    if ((msh.cells[0].type != 'triangle') and (msh.cells[0].type != 'tetra')):
        raise ValueError(f"Dolfin only supports first-order triangle(2D)/tetrahedral(3D) meshes.\n\t Mesh: {xdmfile}\n\t Cell type: {msh.cells[0].type}")        
    if ((msh_facets.cells[0].type != 'line') and (msh_facets.cells[0].type != 'triangle')):
        raise ValueError(f"Dolfin only supports triangle(2D)/tetrahedral(3D) as facets.\n\t Mesh: {xdmfile_facets}\n\t Cell type: {msh_facets.cells[0].type}")
    
        # Build dolfin mesh and mesh function over triangle/tetra
        # From https://fenicsproject.discourse.group/t/transitioning-from-mesh-xml-to-mesh-xdmf-from-dolfin-convert-to-meshio/412/157
        # (March 2021)
    dolfin_mesh = dolfin.Mesh()
    with dolfin.XDMFFile(xdmfile) as infile:
            # Dolfin mesh
        print(f"Building dolfin mesh from {xdmfile}...")
        infile.read(dolfin_mesh)
        dim = dolfin_mesh.topology().dim() # Mesh dimension
        N_triangle = dolfin_mesh.cells().shape[0]
        N_node = dolfin_mesh.coordinates().shape[0]
        order = dolfin_mesh.ufl_coordinate_element().degree()
        print(f"Done. Dolfin mesh properties:\n\tDimension: {dim}\n\tOrder: {order}\n\tTetra/triangle: {N_triangle}\n\tNodes: {N_node}")
            # Build mesh function over triangle/tetra
        print(f"Reading triangle/tetra tags from {xdmfile}...")
        mvc = dolfin.MeshValueCollection("size_t", dolfin_mesh, dim)
        infile.read(mvc, "physical_entities_tag")
        domains = dolfin.cpp.mesh.MeshFunctionSizet(dolfin_mesh, mvc)
        print(f"Done. Tags found: {domains.size()}")
        # Build mesh function over line/triangle
    print(f"Reading facets tags from {xdmfile_facets}...")
    mvc = dolfin.MeshValueCollection("size_t", dolfin_mesh, dim-1)
    with dolfin.XDMFFile(xdmfile_facets) as infile:
        infile.read(mvc, "physical_entities_tag")
    boundaries = dolfin.cpp.mesh.MeshFunctionSizet(dolfin_mesh, mvc)  
    print(f"Done. Tags found: {boundaries.size()}")
        # volume measure
    dx = dolfin.Measure("dx",domain=dolfin_mesh,subdomain_data=domains)
        # facet exterior measure
    ds = dolfin.Measure("ds",domain=dolfin_mesh,subdomain_data=boundaries)
        # facet interior measure
    dS = dolfin.Measure("dS",domain=dolfin_mesh,subdomain_data=boundaries)
    return (dolfin_mesh,domains,boundaries,dx,ds,dS)

def mesh_load_dolfinXML(xmlfile):
    """
    Load mesh from DOLFIN-XML file.    

    Parameters
    ----------
    xmlfile : str
        Mesh in DOLFIN-XML format.

    Returns
    -------
    mesh : dolfin.cpp.mesh.Mesh
        Mesh description.
        
    domains: dolfin.cpp.mesh.MeshFunctionSizet
        Function defined over each cell that indicates the domain to which
    the cell belong, expressed as a size_t.
        domains.array().shape=(mesh.num_cells(),)
        
    boundaries: dolfin.cpp.mesh.MeshFunctionSizet
        Function defined over each edge, indicating its surface/edge domain. 
        boundaries.array().shape=(mesh.num_edges(),)
        
    dx: ufl.measure.Measure
        Domain measure, associated with 'domains'.
        
    ds: ufl.measure.Measure
        Exterior measure, associated with 'boundaries'.
    
    dS: ufl.measure.Measure
        Interior measure, associated with 'boundaries'.

    """    
    __check_read_access(xmlfile)
    (fname,ext)=os.path.splitext(xmlfile)
    mesh = fe.Mesh(xmlfile)
    # load subdomains, boundaries    
    fname_domain=fname+"_physical_region.xml"
    domains=None
    if os.path.isfile(fname_domain):
        __check_read_access(fname_domain)
        domains=fe.MeshFunction("size_t",mesh,fname_domain)
    fname_boundaries=fname+"_facet_region.xml"
    boundaries=None
    if os.path.isfile(fname_boundaries):
        __check_read_access(fname_boundaries)
        boundaries=fe.MeshFunction("size_t",mesh,fname_boundaries)
        # surface measure
    dx = fe.Measure("dx",domain=mesh,subdomain_data=domains)
        # exterior measure
    ds = fe.Measure("ds",domain=mesh,subdomain_data=boundaries)
        # interior measure
    dS = fe.Measure("dS",domain=mesh,subdomain_data=boundaries)
    return (mesh,domains,boundaries,dx,ds,dS)


class DolfinMesh():

    def __init__(self,mesh,domains,boundaries,dx,ds,dS):
        self.mesh=mesh
        self.domains=domains
        self.boundaries=boundaries
        self.dx=dx
        self.ds=ds
        self.dS=dS

    @classmethod
    def init_from_gmsh(cls,xmlfile):
        (mesh,domains,boundaries,dx,ds,dS) = mesh_load_dolfinXML(xmlfile)
        return cls(mesh,domains,boundaries,dx,ds,dS)
    
    @classmethod
    def init_from_xdmf(cls,xdmfile,xdmfile_facets):
        (mesh,domains,boundaries,dx,ds,dS) = mesh_load_xdmf(xdmfile,xdmfile_facets)
        return cls(mesh,domains,boundaries,dx,ds,dS)

##############

def export_xdmf(xdmfile,sol_names,sol_spaces,t,y):
    """ Export time-domain solutions to XDMF file. All timesteps are stored
    in one file.
    TODO: Make it work with subspaces.
    Inputs
    ------
    xdmfile: str

    sol_names: list of str (length N_sol)
        Name of each exported function.
    
    sol_spaces: list of dolfin FunctionSpace (length N_sol)
        Function space of each exported function.
                        
    t,y: list of float, PETSc.Vec (length N)
        Computed solutions.
            
    """
    
    output = fe.XDMFFile(xdmfile)
    output.parameters["flush_output"] = True
    output.parameters["rewrite_function_mesh"] = False
    output.parameters["functions_share_mesh"] = True
    print(f"Export to {xdmfile}...")
    for i in range(len(y)):
        print(f"Step t={t[i]}")
        for j in range(len(sol_names)):
            # ui = fe.Function(Vi,fe.PETScVector(eigvec[i_plot]))
            fun = fe.Function(sol_spaces[j],fe.PETScVector(y[i]),name=sol_names[j])
            output.write(fun, t[i])
    print("Done.")

class PeriodicBoundary(fe.SubDomain):
    """
    Class that describes a list of periodic boundary conditions. It is
    designed to be compatible with dolfin.function.functionspace.FunctionSpace.

    Example
    --------
        # Create and initialize class
    pbc = fenics_utils.PeriodicBoundary()
    pbc.init()
        # Append periodic boundary conditions
    pbc.append(is_src1,is_dst1,map_src1_to_dst1)
    pbc.append(is_src2,is_dst2,map_src2_to_dst2)    
        # Build function space
    V = FunctionSpace(mesh, 'P', 1,constrained_domain=pbc)

    """
    
    def init(self):
        """
        Initialize class members.
        
        Remark
        ------
        This has to be called manually. Using __init__ would trigger an error
        in dolfin.

        Returns
        -------
        None.

        """
        self.is_src_lst = []
        self.is_dst_lst = []
        self.map_lst = []
        self.N = 0

    
    def append(self,is_src,is_dst,map_src_to_dst):
        """
        Append a periodic boundary condition. A periodic boundary condition
        is desribed by a source characteristic function, a destination
        characteristic function, and a map from source to destination.

        Parameters
        ----------
        is_src : function x -> bool
            is_src(x) is true if x belongs to the source boundary.
            
        is_dst : function y -> bool
            is_dst(y) is true if y belongs to the destination boundary.
            
        map_src_to_dst : function (x,y) -> void
            Map x (source) to y (destination). Only y is modified.

        Returns
        -------
        None.

        """
        self.is_src_lst.append(copy.deepcopy(is_src))
        self.is_dst_lst.append(copy.deepcopy(is_dst))
        self.map_lst.append(copy.deepcopy(map_src_to_dst))
        self.N += 1
        
    def inside(self,y,on_boundary):
        """
        Identify destination boundary.

        Parameters
        ----------
        y : real vector
            Coordinates.
            
        on_boundary : bool
            True if y belongs to a mesh boundary.

        Returns
        -------
        is_on_dest: bool
            True if y belongs to at least one of the destination boundaries.

        """
        is_inside = False
        for is_dst in self.is_dst_lst:
            is_inside = is_inside or is_dst(y)
        return is_inside and on_boundary
    
    def map(self,x,y):
        """
        Map x (destination) to y (source).

        Parameters
        ----------
        x,y : real vector
            Coordinates.

        Returns
        -------
        None.

        """
        y[0]=x[0]; y[1]=x[1] # dummy value
        for i in range(self.N): # for each periodic boundary conditions
            if (self.is_src_lst[i](x)): # x belongs to source boundary i?
                self.map_lst[i](x,y)
                
def assemble_GEP(A,a,B,b,bcs,diag_A=1e2,diag_B=1e-2):
    """
    Assemble the generalized eigenvalue problem
        a(u,v) = lambda*b(u,v) with Dirichlet boundary conditions 'bc'
    as
        A*U = lambda*B*U.

    Parameters
    ----------
    A,B : dolfin.cpp.la.PETScMatrix
        Output matrices, must be created by the calling code.
    a,b  : ufl.form.Form
        Input bilinear form.
    bcs : list of dolfin.fem.dirichletbc.DirichletBC
        Dirichlet boundary conditions
    diag_A, diag_B : float, optional
        Diagonal penalization to enforce Dirichlet boundary condition, see remark.

    Returns
    -------
    None.
    
    Remark
    ------
    Dirichlet boundary conditions are enforced by modifying A and B as follows:
        A[i,j]=A[j,i]=0, B[i,j]=B[j,i]=0, A_[i,i]=diag_A, B_[i,i]=diag_B,
    for i DoF on Dirichlet boundary. This induces a spurious eigenvalue
    lambda=diag_A/diag_B, whose multiplicity is the number of Dirichlet DoF.
    The associated spurious eigenfunctions are localized at each Dirichlet DoF.

    """
    fe.assemble(a,tensor=A)
    fe.assemble(b,tensor=B)
        # Enforce homogeneous Dirichlet boundary conditions
#    vec = fe.Function(V).vector()
    vec = fe.PETScVector(); vec.init(A.mat().size[0])
    for bc in bcs:
        bc.zero_columns(A,vec,diag_A)
        bc.zero_columns(B,vec,diag_B)
    
        
def create_GEP(A,B,opts={'solver': "arpack"}):
    """
    Create eigensolver for generalized eigenvalue problem A*x=lambda*B*x.

    Parameters
    ----------
    A,B: dolfin.cpp.la.PETScMatrix
        Matrices for generalized eigenvalue problem.
    opts : dict, optional
        Options. List available in dolfin cpp API/SLEPcEigenSolver.

    Returns
    -------
    EPS : dolfin.cpp.la.SLEPcEigenSolver
        Fenics wrapper around SLEPc eigensolver.
        
    Remark
    ------
    The Fenics wrapper around SLEPc is feature-limited. A more flexible
    alternative is to use SLEPc directly.

    """
    EPS = fe.SLEPcEigenSolver(A,B)
    for key in opts.keys():
        EPS.parameters[key]=opts[key]
    return EPS

def __check_read_access(file):
    file = open(file, "r"); file.close()
