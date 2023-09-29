# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2023 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""
- file:             utils/linalg.py
- authors:          Ghislain Haine, Florian Monteghetti
- date:             29 nov. 2022
- brief:            linear algebra functions
"""

import sys
import getfem as gf
import numpy as np
import petsc4py
petsc4py.init(sys.argv)
from petsc4py import PETSc

comm = PETSc.COMM_WORLD
import scipy.sparse as sp

def extract_gmm_to_petsc(I, J, M, B, comm=comm):
    """Extract a sub-matrix A from M, on interval I, J
    
    Args:
        I (Numpy array): line interval [begin, lenght]
        J (Numpy array): column interval [begin, lenght]
        M (SPMat GetFEM): matrix from which to extract the submatrix
        comm (MPI_Comm): MPI communicator
    
    Returns:
        PETSc.Mat: matrix with value M(I,J) in CSR format
    """
    
    R_I = range(I[0],I[1]) # Range of lines extraction
    R_J = range(J[0],J[1]) # Range of columns extraction
    
    # Because CSR is CSC of the transpose
    # and CSR is not available in getfem
    # I and J are switched
    # See if it can be optimised
    A = gf.Spmat('empty', I[1], J[1]) # Pre-allocation
    A.assign(R_I, R_J, gf.Spmat('copy', M, R_I, R_J))
    A.transpose()
    A.to_csc()
    
    A_ind = A.csc_ind() 
    indrow = A_ind[0] 
    indcol = A_ind[1] 
    data = A.csc_val()
    del A
    
    B.setValuesLocalCSR(indrow,indcol,data,addv=PETSc.InsertMode.INSERT_VALUES)
    B.assemble()
    
    return B
    
def extract_gmm_to_scipy(I,J,M):
    """Extract a sub-matrix A from M, on interval I, J
    
    Args:
        I (Numpy array): line interval [begin, lenght]
        J (Numpy array): column interval [begin, lenght]
        M (SPMat GetFEM): matrix from which to extract the submatrix
    
    Returns:
        PETSc.Mat: matrix with value M(I,J) in CSR format
    Returns:
        scipy.sparse.csr.csr_matrix: the matrix A in scipy.sparse.csr.csr_matrix format
    """
    
    R_I = range(I[0],I[1]) # Range of lines extraction
    R_J = range(J[0],J[1]) # Range of columns extraction
    
    # Because CSR is CSC of the transpose
    # and CSR is not available in getfem
    # I and J are switched
    # See if it can be optimised
    A = gf.Spmat('empty', I[1], J[1]) # Pre-allocation
    A.assign(R_I, R_J, gf.Spmat('copy', M, R_I, R_J))
    A.transpose()
    A.to_csc()
    
    A_ind = A.csc_ind() 
    indrow = A_ind[0] 
    indcol = A_ind[1] 
    data = A.csc_val()
    del A
    
    A_scipy = sp.csr_matrix((data, indcol, indrow), shape=(J[1],I[1]))
    
    return A_scipy
    
def convert_PETSc_to_scipy(A):
    """Convert from PETSc.Mat to scipy sparse (csr).
    
    Args:
        A (PETSc Mat): The matrix to convert

    Returns:
        scipy.sparse.csr.csr_matrix: the matrix A in scipy.sparse.csr.csr_matrix format
    """
    
    (indrow,indices,val) = A.getValuesCSR()
    A_scipy = sp.csr_matrix((val, indices, indrow), shape=A.size)
    
    return A_scipy
