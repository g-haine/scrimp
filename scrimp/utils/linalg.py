# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2025 ISAE-SUPAERO -- GNU GPLv3
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
PETSc.Options().setValue('info', None)

comm = PETSc.COMM_WORLD
import scipy.sparse as sp


def convert_gmm_to_petsc(M, B, comm=comm):
    """Convert a GetFEM matrix M to a PETSc one B

    Args:
        M (SPMat GetFEM): matrix to transfer
        B (PETSc.Mat): matrix to fill M with
        comm (MPI_Comm): MPI communicator

    Returns:
        None
    """

    A = gf.Spmat("copy", M)
    # Because CSR is CSC of the transpose
    # and CSR is not available in getfem
    A.transpose()
    row_ptr, col_idx = A.csc_ind()
    values = A.csc_val()
    
    B.setOption(PETSc.Mat.Option.NEW_NONZERO_ALLOCATION_ERR,False)
    B.setValuesLocalCSR(row_ptr, col_idx, values, addv=PETSc.InsertMode.INSERT_VALUES)
    B.setOption(PETSc.Mat.Option.NEW_NONZERO_ALLOCATION_ERR,True)

    # Assemble the PETSc matrix in parallel
    B.assemble()

def extract_gmm_to_scipy(I, J, M):
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

    R_I = range(I[0], I[1])  # Range of lines extraction
    R_J = range(J[0], J[1])  # Range of columns extraction

    # Because CSR is CSC of the transpose
    # and CSR is not available in getfem
    # I and J are switched
    # See if it can be optimised
    A = gf.Spmat("empty", I[1], J[1])  # Pre-allocation
    A.assign(R_I, R_J, gf.Spmat("copy", M, R_I, R_J))
    A.transpose()
    A.to_csc()

    A_ind = A.csc_ind()
    indrow = A_ind[0]
    indcol = A_ind[1]
    data = A.csc_val()
    del A

    A_scipy = sp.csr_matrix((data, indcol, indrow), shape=(J[1], I[1]))

    return A_scipy


def convert_PETSc_to_scipy(A):
    """Convert from PETSc.Mat to scipy sparse (csr).

    Args:
        A (PETSc Mat): The matrix to convert

    Returns:
        scipy.sparse.csr.csr_matrix: the matrix A in scipy.sparse.csr.csr_matrix format
    """

    (indrow, indices, val) = A.getValuesCSR()
    A_scipy = sp.csr_matrix((val, indices, indrow), shape=A.size)

    return A_scipy
