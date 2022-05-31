#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Utility functions for multiphenics.

@author: Florian Monteghetti
"""
import numpy as np
import fenics as fe
import multiphenics as mpfe
from petsc4py import PETSc
import time
import PETSc_utils

def get_function(y,sol_idx,sol_name,W):
    """ Get dolfin.Function from BlockFunctionSpace.
    
    Inputs
    -------
    y: PETSc.Vec
    
    sol_idx, sol_name: int, str
    
    W: multiphenics BlockFunctionSpace
    """
    ul = mpfe.BlockFunction(W,mpfe.la.BlockPETScVector(y),name=sol_name, label=sol_name)
    u = mpfe.block_split(ul)[sol_idx]
    return u

def export_xdmf(xdmfile,sol_idxs,sol_names,t,y,W):
    """ Export time-domain solutions to XDMF file. All timesteps are stored
    in one file.
    
    Inputs
    ------
    xdmfile: str
    
    sol_idxs: list of int
        Indices of solutions to export.
        
    sol_names: list of str
        Names of exported solutions in xdmf file.
                
    t,y: list of float, PETSc.Vec
        Computed solutions.
        
    W: multiphenics BlockFunctionSpace
    
    """
    
    for j in range(len(sol_idxs)):
        xdmffile = xdmfile+sol_names[j]+".pvd"
        output = fe.XDMFFile(xdmffile)
        output.parameters["flush_output"] = True
        output.parameters["rewrite_function_mesh"] = False
        output.parameters["functions_share_mesh"] = True
        print(f"Export to {xdmffile}...")
        for i in range(len(y)):
            print(f"Step t={t[i]}")
            f = get_function(y[i], sol_idxs[j],sol_names[j],W)
            fun = fe.Function(W.sub(sol_idxs[j]),f.vector(),name=sol_names[j])
            output.write(fun, t[i])
    print("Done.")

def export_pvd(pvdfile,sol_idxs,sol_names,t,y,W):
    """ Export one time-domain solution to pvd file. Timesteps are stored in
    separate files.
    
    Inputs
    ------
    pvdfile: str
    
    sol_idx: int
        Index of solution to export.
    
    sol_name: str
        Name of exported solution in pvd file.
        
    t,y: list of float, PETSc.Vec
        Computed solutions.
        
    W: multiphenics BlockFunctionSpace
    
    """

    for j in range(len(sol_idxs)):
        pvdf = fe.File(pvdfile+sol_names[j]+".pvd")
        print(f"Export to {pvdf}...")
        for i in range(len(y)):
            print(f"Step t={t[i]}")
            f = get_function(y[i], sol_idxs[j],sol_names[j],W)
            fun = fe.Function(W.sub(sol_idxs[j]),f.vector(),name=sol_names[j])
            pvdf << (fun, t[i])
    print("Done.")

def build_MeshRestriction_from_tags(mesh,boundaries,target_tags):
    """
    Build a MeshRestriction object that describes the vertices, edges, facets,
    and cells that belong to the target boundaries.
    
    Parameters
    ----------
    mesh : dolfin.cpp.mesh.Mesh
        
    boundaries : dolfin.cpp.mesh.MeshFunctionSizet
        Tag of every edge/facets, size N_edge/N_facets.
        
    target_tags : list(int)
        Tags of target boundaries (2D mesh: 1D, 3D mesh: 2D).

    Returns
    -------
    list(dolfin.cpp.mesh.MeshFunctionBool)
        Functionally similar to MeshRestriction object. Each MeshFunctionBool
        is True where the entity belong to one of the target boundaries.
        2D: vertices, edges, cells.
        3D: vertices, edges, facets, cells.

    """
    dim = mesh.topology().dim()
        # -- Cells
    mf_cells = fe.MeshFunction("bool", mesh, 2 if dim==2 else 3)
    mf_cells.set_all(False)
        # -- Facets
    if dim==3:
        mf_facets = fe.MeshFunction("bool", mesh, 2)
            # True if facet on target boundary, False elsewhere
        target_facets = np.zeros(mesh.num_facets(),dtype='bool')
        for tag in target_tags:
            target_facets += boundaries.array()==tag
        mf_facets.set_values(target_facets)
        # -- Edge
    mf_edges = fe.MeshFunction("bool", mesh, 1)
        # True if edge on target boundary, False elsewhere
    target_edges = np.zeros(mesh.num_edges(),dtype='bool')
    if dim==2:            
        for tag in target_tags:
            target_edges += boundaries.array()==tag
    else:
            # Facet to edges connectivity (N_facets,3)
        facet_2_edges = np.array([edge.index() for facet in fe.facets(mesh) for edge in fe.edges(facet)])
        facet_2_edges = np.reshape(facet_2_edges,(mesh.num_facets(),3))
        target_edges[facet_2_edges[target_facets].flatten()] = True
    mf_edges.set_values(target_edges)
        # -- Vertex
    mf_vertices = fe.MeshFunction("bool", mesh, 0)
         # True if vertex on target boundary, False elsewhere
    target_vertices = np.zeros((mesh.num_vertices(),1),dtype='bool')
    if dim==2:
            # Edge to vertices connectivity (N_edge,2)
        edge_2_vertices = np.array([vert.index() for edge in fe.edges(mesh) for vert in fe.vertices(edge)])
        edge_2_vertices=np.reshape(edge_2_vertices,(mesh.num_edges(),2))
        target_vertices[edge_2_vertices[target_edges].flatten()] = True
    else:
            # Facet to vertices connectivity (N_facet,3)
        facet_2_vertices = np.array([vert.index() for facet in fe.facets(mesh) for vert in fe.vertices(facet)])
        facet_2_vertices = np.reshape(facet_2_vertices,(mesh.num_facets(),3))
        target_vertices[facet_2_vertices[target_facets].flatten()] = True
    mf_vertices.set_values(target_vertices)
        # -- Mesh restriction
    if dim==2:
        meshr = [mf_vertices,mf_edges,mf_cells]
    else:
        meshr = [mf_vertices,mf_edges,mf_facets,mf_cells]
    return meshr

def build_MeshRestriction_from_tags_subdomains(mesh,domains,target_tags):
    """
    Identical to build_MeshRestriction_from_tags, except that target_tags 
    identifies target subdomains (2D mesh: 2D, 3D mesh: 3D).
    """
    dim = mesh.topology().dim()
        # -- Cells
    mf_cells = fe.MeshFunction("bool", mesh, 2 if dim==2 else 3)
            # True if cell in target subdomain, False elsewhere
    target_cells = np.zeros(mesh.num_cells(),dtype='bool')
    for tag in target_tags:
        target_cells += domains.array()==tag            
    mf_cells.set_values(target_cells)
        # -- Facets
    if dim==3:
        mf_facets = fe.MeshFunction("bool", mesh, 2)
            # True if facet is on boundary of a cell in target subdomain, False elsewhere
        target_facets = np.zeros((mesh.num_facets(),),dtype='bool')
        cell_2_facets = np.array([cell.entities(2) for cell in fe.cells(mesh)])
        target_facets[cell_2_facets[target_cells].flatten()] = True
        mf_facets.set_values(target_facets)
        # -- Edges
    mf_edges = fe.MeshFunction("bool", mesh, 1)
            # True if edge is boundary of a cell in target subdomain, False elsewhere
    target_edges = np.zeros((mesh.num_edges(),),dtype='bool')
            # Cell to edges connectivity (N_cell,3/6)
    cell_2_edges = np.array([cell.entities(1) for cell in fe.cells(mesh)])
    target_edges[cell_2_edges[target_cells].flatten()] = True
    mf_edges.set_values(target_edges)
        # -- Vertex
    mf_vertices = fe.MeshFunction("bool", mesh, 0)
            # True if vertex on target subdomain, False elsewhere
    target_vertices = np.zeros((mesh.num_vertices(),1),dtype='bool')
            # Cell to vertices connectivity (N_cell,3/4)
    cell_2_vertices = mesh.cells()
    target_vertices[cell_2_vertices[target_cells].flatten()] = True
    mf_vertices.set_values(target_vertices)
        # -- Mesh restriction
    if dim==2:
        meshr = [mf_vertices,mf_edges,mf_cells]
    else:
        meshr = [mf_vertices,mf_edges,mf_facets,mf_cells]
    return meshr


class PDAE_linear(PETSc_utils.Fully_implicit_DAE):
    """ Linear PDAE F(t,y,yd) = 0 assembled with multiphenics and suitable for
    integration with PETSc TS fully implicit schemes.
    
    Description
    -----------
    This class describes a PDAE written as    
    
        F(t,y,yd) = 0 
        
    with residual vector given by
    
        F(t,y,yd) = J_yd * yd + J_y * y + F_src(t),
        
    where the jacobian matrices J_yd and J_y are constant.
    
    The PDAE is defined through vectors and matrices of (bi)linear forms, for
    assembly with multiphenics.
    
    Methods are provided to enable time integration with any fully implicit
    PETSc TS scheme. To reduce the computational cost, the jacobian matrices
    are assembled only once at the start of the integration. Results are stored
    in history.
    
    """ 
    def __init__(self,F_src_form,J_form,N,name,idx_alg=None):
        """
        Constructor.
        
        Parameters
        ----------
        F_src_form : function t-> list(ufl.form.Form)
            Time-dependent residual vector F_src(t).
            
        J_form : function returning tuple(list(list(ufl.form.Form)))
            Constant jacobian matrices (J_yd,J_y).
            
        N : int
            Size of unknown vector y.
            
        name : str
            Name.
            
        idx_alg : list(int) (Optional)
            Indices of algebraic variables.

        """
        super().__init__(N,name,idx_alg=idx_alg)
        self.F_src_form = F_src_form
        self.J_form = J_form
           
    def init(self):
        """ Initialize. Call before new time integration. """
            # History
        self.init_history()
             # Jacobian matrix setup
        self.init_jacobian()
            # Residual Vector setup
        self.init_residual()        

    def assemble_res_src(self,Fsrc,t):
        """ Assemble F_src(t). """
        mpfe.block_assemble(self.F_src_form(t),
                            block_tensor=mpfe.la.BlockPETScVector(Fsrc))
        
    def init_residual(self):
        """ Allocate vectors required by time integration. """
        self.F = self.J_yd.createVecRight()
        self.Fsrc = self.F.duplicate()
        self.Ftmp = self.J_yd.createVecRight()
            
    def init_jacobian(self):
        """ Allocate matrices required by time integration. """
        self.J_yd = fe.PETScMatrix().mat() # dF/dyd
        self.J_y = fe.PETScMatrix().mat() # dF/dy
        self.assemble_jacobian(self.J_yd,self.J_y)
        self.J = self.J_y.duplicate() # c*dF/dyd + dF/dy
        
    def assemble_jacobian(self,J_yd,J_y):
        """ Assemble jacobian matrices J_yd and J_y. """
        (J_yd_fun,J_y_fun) = self.J_form()

        mpfe.block_assemble(J_yd_fun,mpfe.la.BlockPETScMatrix(J_yd),keep_diagonal=True)
        mpfe.block_assemble(J_y_fun,
            mpfe.la.BlockPETScMatrix(J_y),keep_diagonal=True)

    def IFunction(self, ts, t, y, yd, F):
        """ Evaluate residual vector F(t,y,yd) = J_yd*yd + J_y*y + F_src(t)."""
        # print(f"IFunction: t={t}")
        # F = self.J_yd * yd + self.J_y * y + F_src
        self.J_yd.mult(yd,self.Ftmp)
        self.J_y.multAdd(y,self.Ftmp,F)
        self.assemble_res_src(self.Fsrc,t)
        F.axpy(1,self.Fsrc)
        
    def IJacobian(self, ts, t, y, yd, c, Amat, Pmat):
        """ Evaluate jacobian matrix J = c*J_yd + J_y."""
        # Pmat = c*J_yd + J_y
        self.J_y.copy(Pmat)
        Pmat.axpy(c,self.J_yd)
        if Amat != Pmat:
            print("Operator different from preconditioning")
            Amat.assemble()
            
            
class PDAE_nonlinear(PETSc_utils.Fully_implicit_DAE):
    """ Similar to PDAE_linear. The difference is that the residual is 
        F(t,y,yd) = J_yd * yd + J_y * y + F_nl(y) + F_src(t)."""
        
    def __init__(self,F_src_form,J_form,F_nl_form,J_nl_form,N,name,idx_alg=None,J_nl_skip=0):
        """
        Constructor.
        
        Parameters
        ----------
        F_src_form : function t-> list(ufl.form.Form)
            Time-dependent residual vector F_src(t).
            
        J_form : function returning tuple(list(list(ufl.form.Form)))
            Constant jacobian matrices (J_yd,J_y).
            
        F_nl_form : function PETSc.Vec -> list(ufl.form.Form)
            Non-linear part of the residual vector.
            
        J_nl_form : function ETSc.Vec -> tuple(list(list(ufl.form.Form)))
            Jacobian matrix associated with F_nl_form.
            
        N : int
            Size of unknown vector y.
            
        name : str
            Name.
            
        idx_alg : list(int) (Optional)
            Indices of algebraic variables.

        """
        super().__init__(N,name,idx_alg=idx_alg)
        self.F_src_form = F_src_form
        self.J_form = J_form
        self.F_nl_form = F_nl_form
        self.J_nl_form = J_nl_form
        self.J_nl_skip = J_nl_skip
        self.J_nl_last = self.J_nl_skip
        self.J_nl_assembly_count = 0 # for monitoring
        
        
    def init(self):
        """ Initialize. Call before new time integration. """
            # History
        self.init_history()
             # Jacobian matrix setup
        self.init_jacobian()
            # Residual Vector setup
        self.init_residual()        

    def assemble_res_src(self,Fsrc,t):
        """ Assemble F_src(t). """
        mpfe.block_assemble(self.F_src_form(t),
                            block_tensor=mpfe.la.BlockPETScVector(Fsrc))
        
    def init_residual(self):
        """ Allocate vectors required by time integration. """
        self.F = self.J_yd.createVecRight()
        self.Fsrc = self.F.duplicate()
        self.Ftmp = self.J_yd.createVecRight()
            
    def init_jacobian(self):
        """ Allocate matrices required by time integration. """
        self.J_yd = fe.PETScMatrix().mat() # dF/dyd
        self.J_y = fe.PETScMatrix().mat() # dF/dy
        self.assemble_jacobian(self.J_yd,self.J_y)
        self.J = self.J_y.duplicate() # c*dF/dyd + dF/dy
        self.J_nl = fe.PETScMatrix().mat()
        
    def assemble_jacobian(self,J_yd,J_y):
        """ Assemble jacobian matrices J_yd and J_y. """
        (J_yd_fun,J_y_fun) = self.J_form()

        mpfe.block_assemble(J_yd_fun,mpfe.la.BlockPETScMatrix(J_yd),keep_diagonal=True)
        mpfe.block_assemble(J_y_fun,
            mpfe.la.BlockPETScMatrix(J_y),keep_diagonal=True)

    def IFunction(self, ts, t, y, yd, F):
        """ Evaluate residual vector F(t,y,yd) = J_yd*yd + J_y*y + F_src(t)."""
        # print(f"IFunction: t={t}")
        # F = self.J_yd * yd + self.J_y * y + F_src
        self.J_yd.mult(yd,self.Ftmp)
        self.J_y.multAdd(y,self.Ftmp,F)
        self.assemble_res_src(self.Fsrc,t)
        F.axpy(1,self.Fsrc)
        mpfe.block_assemble(self.F_nl_form(y),block_tensor=mpfe.la.BlockPETScVector(self.Ftmp))
        F.axpy(1,self.Ftmp)
        
    def IJacobian(self, ts, t, y, yd, c, Amat, Pmat):
        """ Evaluate jacobian matrix J = c*J_yd + J_y."""
        # Pmat = c*J_yd + J_y
        self.J_y.copy(Pmat)
        Pmat.axpy(c,self.J_yd)
        if (self.J_nl_last == self.J_nl_skip):
            # print("Assenbly J_nl")
            mpfe.block_assemble(self.J_nl_form(y),mpfe.la.BlockPETScMatrix(self.J_nl))
            self.J_nl_last = 0
            self.J_nl_assembly_count += 1
        else:
            # print("skipped")
            self.J_nl_last += 1        
        Pmat.axpy(1.0,self.J_nl)
        if Amat != Pmat:
            print("Operator different from preconditioning")
            Amat.assemble()

        
        
        
