# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2022 Ghislain Haine
#
# See the LICENSE file in the root directory for license information.
#
# github: https://github.com/g-haine/scrimp

"""
- file:             utils/linalg.py
- authors:          Ghislain Haine, Florian Monteghetti
- date:             29 nov. 2022
- last modified:    06 dec. 2022
- brief:            linear algebra functions
"""

import getfem as gf
from petsc4py import PETSc
import scipy.sparse as sp

def extract_gmm_to_petsc(I, J, M, comm=None):
    """
    Extract a sub-matrix A from M, on interval I, J
    
    :param I: line interval [begin, lenght]
    :type I: numpy array
    :param J: column interval [begin, lenght]
    :type J: numpy array
    :param M: matrix from which to extract the submatrix
    :type M: Spmat matrix (gmm format, default in getfem)
    :param comm: MPI communicator
    :type comm: MPI_Comm object or None
    
    :return: a PETSc.Mat matrix with value M(I,J) in CSR format
    """
    
    R_I = range(I[0],I[0]+I[1]) # Range of lines extraction
    R_J = range(J[0],J[0]+J[1]) # Range of columns extraction
    
    A = gf.Spmat('empty', I[1], J[1]) # Pre-allocation
    A.assign(range(I[1]), range(J[1]), gf.Spmat('copy', M, R_I, R_J))
    
    A.transpose() # Because CSR is CSC of the transpose
    A.to_csc() # and CSR is not available in getfem
    
    A_ind = A.csc_ind()
    indptr = A_ind[0]
    indices = A_ind[1]
    data = A.csc_val()
    
    B = PETSc.Mat().createAIJ(size=[I[1],J[1]], csr=(indptr,indices,data), comm=comm)
    
    B.assemble()
    
    return B
    
def convert_PETSc_to_scipy(A):
    """
    Convert from PETSc.Mat to scipy sparse (csr).
    
    @author: Florian Monteghetti
    
    :param A: The matrix to convert
    :type A: PETSc.Mat object

    :return: the matrix A in scipy.sparse.csr.csr_matrix format
    """
    
    (indptr,indices,val) = A.getValuesCSR()
    A_scipy = sp.csr_matrix((val, indices, indptr), shape=A.size)
    
    return A_scipy
