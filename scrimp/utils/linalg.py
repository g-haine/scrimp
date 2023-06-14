# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2023 Ghislain Haine
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

import sys
import getfem as gf
import petsc4py
petsc4py.init(sys.argv)
from petsc4py import PETSc
# // attempt gives strange results
# comm = PETSc.COMM_WORLD
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
    
    # Because CSR is CSC of the transpose
    # and CSR is not available in getfem
    # See if it can be optimised
    A.transpose()
    A.to_csc()
    A_ind = A.csc_ind() 
    indrow = A_ind[0] 
    indcol = A_ind[1] 
    data = A.csc_val()
    
    B = PETSc.Mat()
    # // attempt gives strange results
    # B.create(comm)
    # B.setSizes(A.size())
    # B.setType('aij')
    # B.setUp()
    
    B.createAIJ(size=A.size(), csr=(indrow,indcol,data))
    B.setOption(PETSc.Mat.Option.FORCE_DIAGONAL_ENTRIES, True) # Mandatory for some ksp solvers
    B.setUp()
    
    # // attempt gives strange results
    # Istart, Iend = B.getOwnershipRange()
    # for i in range(Istart,Iend):
    #     row_start = indrow[i]
    #     row_end = indrow[i+1]
    #     B.setValues(i,
    #                 indcol[row_start:row_end],
    #                 data[row_start:row_end],
    #                 addv=PETSc.InsertMode.INSERT_VALUES)
    
    # B.setOption(PETSc.Mat.Option.FORCE_DIAGONAL_ENTRIES, True)
    
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
    
    (indrow,indices,val) = A.getValuesCSR()
    A_scipy = sp.csr_matrix((val, indices, indrow), shape=A.size)
    
    return A_scipy
