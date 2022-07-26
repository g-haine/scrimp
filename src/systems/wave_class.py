# -*- coding: utf-8 -*-
"""
Simplified class for the wave equation

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

# Can be found in src directory of SCRIMP
import generic_class

fe.SubSystemsManager.init_petsc()
from petsc4py import PETSc
import PETSc_utils

class wave(generic_class._generic):
    
    def __init__(self,interconnection=None):
        """
        Constructor
        """
        super().__init__()
        
        
        # Problem variables
        self.name_pHs = 'Wave equation'
        self.type_pHs = 'Distributed'
        self.mass_density = None
        self.Young_modulus = None
        self.Young_modulus_inv = None
        self.formulation = None
        self.causality = None
        self.Neumann_boundaries = None
        self.Dirichlet_boundaries = None
        self.interconnection = None
        
        
        # FEM variables
        self.type_FEM = None
        self.order_FEM = None
        self.W = None
        self.dofs = None
        self.fun_test = None
        self.fun_trial = None
        self.matrices = None
        self.descriptor = None
        
        
        # TS variables
        self.t_init = 0.
        self.dt = 0.01
        self.dt_export = 5.*self.dt
        self.t_final = 1.
        self.ts_scheme = 'cn'
        self.ts_order = None
        self.ksp_solver = 'preonly'
        self.ksp_pc = 'lu'
        self.Jacobian = None
        self.residual = None
        
        
        # Solution variables
        self.t = None
        self.z = None
        self.Hamiltonian = None
        self.kinetic = None
        self.potential = None
        self.supplied = None
        
        
        # Class variables
        self.formulation_done = False
        self.causality_done = False
        self.FEM_done = False
        self.mass_density_done = False
        self.Young_modulus_done = False
        self.matrices_done = False
        self.descriptor_done = False
        self.time_interval_done = False
        
        
        # Set mandatory variables
        self._set_interconnection(self.interconnection)
        
        
    def __call__(self):
        print('Class wave of SCRIMP')
        print(40*'-', '\n')
        
        
    def _set_interconnection(self,interconnection=None):
    
        self.interconnection = interconnection
        
        
    def mesh_with_gmsh(self,geofile,parameters=None,refinement=0):
        
        assert self.interconnection is None, \
            'Mesh is handled by interconnection object if it exists.'
        
        super().mesh_with_gmsh(geofile,parameters=parameters,refinement=refinement)
        
        
    def set_formulation(self,formulation):
    
        assert formulation in ['div', 'grad'], \
            f'unknown \'{formulation}\' formulation, choose \'div\' or \'grad\''
        
        self.formulation = formulation
        
        self.formulation_done = True
        
    def set_causality(self,causality):
    
        assert self.mesh_done, \
            'Please mesh before setting causality'
        
        for gamma in self.name_boundaries:
            assert gamma in list(causality.keys()), \
                f'Please add a condition on {gamma}: \'Neumann\', \'Dirichlet\' or \'None\''
        
        for gamma in list(causality.keys()):
            assert gamma in self.name_boundaries, \
                f'Unkown boundary {gamma}'
        
        self.causality = causality
        self.Neumann_boundaries = [k for k, v in self.causality.items() if v == 'Neumann']
        self.Dirichlet_boundaries = [k for k, v in self.causality.items() if v == 'Dirichlet']
        
        assert (len(self.Neumann_boundaries)>0 or len(self.Neumann_boundaries)>0), \
            'At least one boundary must have a \'Neumann\' or \'Dirichlet\' condition'
        
        self.causality_done = True
        
    def set_FEM(self,FEM):
    
        assert self.domain_done, \
            'Please set domain before setting FE families'
    
        assert self.causality_done, \
            'Please set causality before setting FE families'
        
        type_FEM = [val[0] for val in FEM]
        order_FEM = [val[1] for val in FEM]
        nb_Neumann_boundaries = len(self.Neumann_boundaries)
        nb_Dirichlet_boundaries = len(self.Dirichlet_boundaries)
        
        assert len(type_FEM)==(nb_Neumann_boundaries+nb_Dirichlet_boundaries+2), \
            'The number of FE must be equal to 2 (co-energy variables) plus the number of boundaries with condition different of \'None\''
            
        # Care must be taken: for Galerkin FE, VectorFunctionSpace must be used instead of FunctionSpace
        if ((type_FEM[0] == 'CG') or (type_FEM[0] == 'DG')):
            V_q = fe.VectorFunctionSpace(self.dmesh.mesh, type_FEM[0], order_FEM[0])
        else:
            V_q = fe.FunctionSpace(self.dmesh.mesh, type_FEM[0], order_FEM[0])
        V_p = fe.FunctionSpace(self.dmesh.mesh, type_FEM[1], order_FEM[1])
        spaces = [V_q, V_p]
        
        spaces = [V_q, V_p]
        restrictions = [self.domain_restriction, self.domain_restriction]
        
        for k in range(nb_Neumann_boundaries):
            # 2 times for input AND output
            spaces.append(fe.FunctionSpace(self.dmesh.mesh, type_FEM[k], order_FEM[k]))
            restrictions.append(self.boundaries_restriction[self.Neumann_boundaries[k]])
            spaces.append(fe.FunctionSpace(self.dmesh.mesh, type_FEM[k], order_FEM[k]))
            restrictions.append(self.boundaries_restriction[self.Neumann_boundaries[k]])
        for k in range(nb_Dirichlet_boundaries):
            # 2 times for input AND output
            spaces.append(fe.FunctionSpace(self.dmesh.mesh, type_FEM[k], order_FEM[k]))
            restrictions.append(self.boundaries_restriction[self.Dirichlet_boundaries[k]])
            spaces.append(fe.FunctionSpace(self.dmesh.mesh, type_FEM[k], order_FEM[k]))
            restrictions.append(self.boundaries_restriction[self.Dirichlet_boundaries[k]])
            
        self.W = mpfe.BlockFunctionSpace(spaces, restrict=restrictions)

        self.fun_test = mpfe.block_split(mpfe.BlockTestFunction(self.W))
        self.fun_trial = mpfe.block_split(mpfe.BlockTrialFunction(self.W))
        
        self.dofs = []
        for k in range(len(spaces)):
            self.dofs.append(len(self.W.block_dofmap().block_owned_dofs__local_numbering(k)))
        
        self.type_FEM = type_FEM
        self.order_FEM = order_FEM
        
        self.FEM_done = True
        
        
    def set_mass_density(self,rho_expression,degree=2):
    
        assert self.FEM_done, \
            'Please set FE families before mass density'
    
        self.mass_density = fe.Expression(rho_expression, degree=degree)
        
        self.mass_density_done = True
        
        
    def set_Young_modulus(self,T_expression,degree=2):
    
        assert self.FEM_done, \
            'Please set FE families before Young\'s modulus'
    
        self.Young_modulus = fe.Expression(T_expression, degree=degree)
    
        self.Young_modulus_inv = fe.inv(self.Young_modulus)
    
        self.Young_modulus_done = True
        
        
    def construct_matrices(self):
    
        assert self.formulation_done, \
            'Please set formulation: \'div\' or \'grad\''
    
        assert self.mass_density_done, \
            'Please set mass density'
    
        assert self.Young_modulus_done, \
            'Please set Young\'s modulus'
        
        u = self.fun_trial
        v = self.fun_test
        
        mass_matrices = [fe.dot(u[0]*fe.dot(self.Young_modulus_inv,v[0]))*self.dmesh.dx]
        mass_matrices.append(u[1]/self.mass_density*v[1]*self.dmesh.dx)
        for k in range(2*(len(self.Neumann_boundaries)+len(self.Dirichlet_boundaries))): # 2 times for input AND output
            mass_matrices.append(u[k+2]*v[k+2]*self.dmesh.ds)
        
        if self.formulation=='grad':
            Dt = fe.dot(fe.grad(v[1]),u[0])*self.dmesh.dx
            D = fe.dot(fe.grad(u[1]),v[0])*self.dmesh.dx
            
            
        elif self.formulation=='div':
            D = fe.div(u[0])*v[1]*self.dmesh.dx
            Dt = fe.div(v[0])*u[1]*self.dmesh.dx
            
            
        else:
            raise ValueError(f'unknown \'{self.formulation}\' formulation. Please set formulation to \'div\' or \'grad\'')
            
    
        return
        
        
    def construct_descriptor(self):
    
        return
        
        
    def set_Jacobian(self):
    
        return
        
        
    def set_control(self):
    
        return
        
        
    def set_residual(self):
    
        return
        
        
    def set_initial_data(self):
    
        return
        
        
    def TS_integration(self):
    
        return
        
        
    def Hamiltonian(self):
    
        return
        
        
    def plot_Hamiltonian(self):
        
        return
        
        
    def plot_state_at(self):
        
        return
        
        
    def export_to_pvd(self):
        
        return
        
        
    def export_to_mat(self):

        from scipy.io import savemat
        
        return
        
        
