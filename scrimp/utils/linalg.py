# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
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
import os
import getfem as gf
import numpy as np
import petsc4py
from slepc4py import SLEPc

petsc4py.init(sys.argv)
from petsc4py import PETSc
PETSc.Options().setValue('info', None)

comm = PETSc.COMM_WORLD
import scipy.sparse as sp

def convert_gmm_to_petsc(M, B, comm=comm):
    """Accumulate a distributed GetFEM model matrix into a PETSc matrix.

    With GetFEM MPI level 2, model tangent matrices contain partial
    contributions on each rank, but use GLOBAL degree-of-freedom indices.
    Contributions to shared entries from different ranks must be added.

    Do not use PETSc's setValuesLocalCSR here: no GetFEM-to-PETSc
    local-to-global mapping has been installed and the GetFEM row indices
    are not the rows locally owned by PETSc.

    All ranks must call this routine. The target matrix retains its PETSc
    ownership distribution; PETSc routes off-process values at assembly.
    """
    A = gf.Spmat("copy", M)
    # GetFEM provides CSC for the transpose: turn it into the CSR of M.
    A.transpose()
    row_ptr, col_idx = A.csc_ind()
    values = A.csc_val()

    nrows, ncols = B.getSize()
    if len(row_ptr) != nrows + 1:
        raise ValueError(
            f"GetFEM/PETSc row count mismatch on rank {comm.getRank()}: "
            f"{len(row_ptr) - 1} vs PETSc {nrows}"
        )
    if len(col_idx) != len(values) or int(row_ptr[-1]) != len(values):
        raise ValueError("Invalid GetFEM CSR data (column/value count)")
    if len(col_idx) and (min(col_idx) < 0 or max(col_idx) >= ncols):
        raise ValueError("GetFEM CSR has columns outside the PETSc matrix")

    # This routine also reassembles nonlinear matrices many times during
    # Newton iterations. ADD_VALUES without clearing would accumulate stale
    # Jacobians from prior iterates.
    trace = os.environ.get("SCRIMP_MPI_TRACE") == "1"
    if trace:
        print(
            f"[MPI {comm.getRank()}/{comm.getSize()}] "
            f"GetFEM matrix -> PETSc: N={nrows}, nnz_local={len(values)}",
            file=sys.stderr, flush=True,
        )
    B.zeroEntries()
    B.setOption(PETSc.Mat.Option.NEW_NONZERO_ALLOCATION_ERR, False)
    try:
        if len(values):
            # IJV accepts global row/column numbering. Unlike setValuesCSR,
            # a global row map handles rows owned by other PETSc processes.
            B.setValuesIJV(
                row_ptr,
                col_idx,
                values,
                addv=PETSc.InsertMode.ADD_VALUES,
                rowmap=np.arange(nrows, dtype=PETSc.IntType),
            )
        if trace:
            print(f"[MPI {comm.getRank()}] before PETSc Mat.assemble", 
                  file=sys.stderr, flush=True)
        B.assemble()
        if trace:
            print(f"[MPI {comm.getRank()}] after PETSc Mat.assemble",
                  file=sys.stderr, flush=True)
    finally:
        B.setOption(PETSc.Mat.Option.NEW_NONZERO_ALLOCATION_ERR, True)


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

def monitor_EPS_short(EPS, it, nconv, eig, err,it_skip):
    """
    Concise monitor for EPS.solve().

    Parameters
    ----------
    eps : slepc4py.SLEPc.EPS
        Eigenvalue Problem Solver class.

    it : int
        Current iteration number.

    nconv : int
        Number of converged eigenvalue.

    eig : list
        Eigenvalues

    err : list
        Computed errors.

    it_skip : int
        Iteration skip.

    Returns
    -------
    None.
    """

    if (it==1):
        print('******************************')
        print('***  SLEPc Iterations...   ***')
        print('******************************')
        EPS_print_short(EPS)
        print("Iter. | Conv. | Max. error")
        print(f"{it:5d} | {nconv:5d} | {max(err):1.1e}")
    elif not it%it_skip:
        print(f"{it} | {nconv} | {max(err):1.1e}")

def EPS_print_short(EPS):
    print(f"Problem dimension: {EPS.getOperators()[0].size[0]}")
    print(f"Solution method: '{EPS.getType()}' with '{EPS.getST().getType()}'")
    nev,ncv = EPS.getDimensions()[0:2]
    print( f"Number of requested eigenvalues: {nev}")
    if ncv>0:
        print(f'Subspace dimension: {ncv}')
        
    tol, maxit = EPS.getTolerances()
    print( f"Stopping condition: tol={tol}, maxit={maxit}")

def EPS_print_results(EPS):
    print()
    print("******************************")
    print("*** SLEPc Solution Results ***")
    print("******************************")           
    its = EPS.getIterationNumber()
    print(f"Iteration number: {its}")
    nconv = EPS.getConverged()
    print( f"Converged eigenpairs: {nconv}")

    if nconv > 0:
        # Create the results vectors
        vr, vi = EPS.getOperators()[0].createVecs()
        print()
        print("Converged eigval.  Error ")
        print("----------------- -------")
        
        for i in range(nconv):
            k = EPS.getEigenpair(i, vr, vi)
            error = EPS.computeError(i)
            if k.imag != 0.0:
                print(f" {k.real:2.2e} + {k.imag:2.2e}j {error:1.1e}")
            else:
                print(f" {k.real:2.2e}         {error:1.1e}")
        print()

def solve_GEP_shiftinvert(A,B,
                          problem_type=SLEPc.EPS.ProblemType.GNHEP,
                          solver=SLEPc.EPS.Type.KRYLOVSCHUR,
                          which=SLEPc.EPS.Which.TARGET_MAGNITUDE,
                          nev=10,tol=1e-6,max_it=10,
                          target=0.0,shift=0.0,defl_space=[],
                          transform_hermitian=False):
    """
    Solve generalized eigenvalue problem A=lambda*B using shift-and-invert
    as spectral transform method.

    Parameters
    ----------
    A,B : PETSc.Mat
    problem_type : SLEPc.EPS.ProblemType, optional
    solver : SLEPc.EPS.Type., optional

    nev : int, optional
        Number of requested eigenvalues.

    tol : float, optional
        Tolerance.

    max_it : int, optional
        Maximum number of iterations.

    target : float, optional
        Target eigenvalue. Also used for sorting.

    shift : float, optional
        Shift 'sigma' used in shift-and-invert.

    defl_space: list of petsc4py.PETSc.Vec, optional
        Deflation space.
    
    transform_hermitian : bool
        Whether to transform the problem into a Hermitian generalized EVP.

    Returns
    -------
    eigval : list of complex
        Converged eigenvalues.

    eigvec_r : list of PETSc.Vec
        Converged eigenvector (real_part)

    eigvec_i : list of PETSc.Vec
        Converged eigenvectors (imag_part)

    """
    
    # Transformation H = iJ if requested
    if transform_hermitian:
        H = PETSc.Mat().createAIJ(size=A.getSize(), comm=comm)
        H.setFromOptions(); H.setUp()
        for i in range(A.getSize()[0]):
            cols, vals = A.getRow(i)
            for j, v in zip(cols, vals):
                H[i, j] = 1j * v
        H.assemble()
        A = H
        A.setOption(PETSc.Mat.Option.HERMITIAN, True)
        B.setOption(PETSc.Mat.Option.HERMITIAN, True)
        B.setOption(PETSc.Mat.Option.SPD, True)
        problem_type = SLEPc.EPS.ProblemType.GHEP

    # Build an Eigenvalue Problem Solver object
    EPS = SLEPc.EPS(); EPS.create(comm=comm)
    EPS.setOperators(A,B)

    # (G)HEP = (Generalized) Hermitian
    # (G)NHEP = (Generalized) Non-Hermitian
    EPS.setProblemType(problem_type)

    # set the number of eigenvalues requested
    EPS.setDimensions(nev=nev)

    # Set solver
    EPS.setType(solver)

    # set eigenvalues of interest
    EPS.setWhichEigenpairs(which)
    EPS.setTarget(target) # sorting

    # set tolerance and max iterations
    EPS.setTolerances(tol=tol,max_it=max_it)    

    # deflation space
    EPS.setDeflationSpace(defl_space)

    # Set up shift-and-invert

    # Only work if 'whichEigenpairs' is 'TARGET_XX'
    ST=EPS.getST()
    ST.setType(SLEPc.ST.Type.SINVERT)
    ST.setShift(shift)        
    EPS.setST(ST)

    # set monitor
    it_skip=1
    EPS.setMonitor(lambda eps, it, nconv, eig, err :
                   monitor_EPS_short(eps, it, nconv, eig, err, it_skip))

    # parse command line options
    EPS.setFromOptions()
    
    # Display all options (including those of ST object)
    EPS.view()            
    EPS.solve()
    EPS_print_results(EPS)

    return EPS

def EPS_get_spectrum(EPS, transform_hermitian=False):

    A = EPS.getOperators()[0]

    # Get results in lists
    eigval = [EPS.getEigenvalue(i) for i in range(EPS.getConverged())]
    eigvec_r = list(); eigvec_i = list()
    
    vr = A.createVecRight(); vi = A.createVecRight()

    for i in range(EPS.getConverged()):
        EPS.getEigenvector(i,vr,vi)
        eigvec_r.append(vr.copy())
        eigvec_i.append(vi.copy())
    
    # Transformation skew-symmétrique : λ = i * ν
    if transform_hermitian:
        eigval = [1j*v for v in eigval]

    # Sort by increasing real parts
    idx=np.argsort(np.real(np.array(eigval)),axis=0)
    eigval = [eigval[i] for i in idx]
    eigvec_r = [eigvec_r[i] for i in idx]
    eigvec_i = [eigvec_i[i] for i in idx]
    
    return (eigval,eigvec_r,eigvec_i)
