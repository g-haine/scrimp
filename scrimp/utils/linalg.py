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
import logging
from typing import Callable, Iterable, Optional, Sequence, Tuple

import getfem as gf
import numpy as np
import petsc4py
from slepc4py import SLEPc

petsc4py.init(sys.argv)
from petsc4py import PETSc

PETSc.Options().setValue('info', None)

comm = PETSc.COMM_WORLD
import scipy.sparse as sp

logger = logging.getLogger(__name__)


def _ensure_int_array(array: Sequence[int]) -> np.ndarray:
    """Return a contiguous numpy array compatible with PETSc indices."""

    if isinstance(array, np.ndarray) and array.dtype == PETSc.IntType:
        return array
    return np.asarray(array, dtype=PETSc.IntType)


def _ensure_scalar_array(array: Sequence[float]) -> np.ndarray:
    """Return a contiguous numpy array compatible with PETSc scalar storage."""

    if isinstance(array, np.ndarray) and array.dtype == PETSc.ScalarType:
        return array
    return np.asarray(array, dtype=PETSc.ScalarType)


def _infer_shape_from_spmat(spmat: gf.Spmat, indptr: Sequence[int], indices: Sequence[int]) -> Tuple[int, int]:
    """Infer the matrix shape using GetFEM meta-data when available."""

    for attr in ("size", "sizes"):
        getter = getattr(spmat, attr, None)
        if getter is None:
            continue
        try:
            shape = getter()
        except TypeError:
            shape = getter
        if isinstance(shape, (tuple, list)) and len(shape) == 2:
            return int(shape[0]), int(shape[1])

    n_rows = len(indptr) - 1
    if len(indices):
        n_cols = int(np.max(indices)) + 1
    else:
        n_cols = n_rows
    return int(n_rows), int(n_cols)


def extract_getfem_csr(M: gf.Spmat) -> Tuple[Tuple[int, int], np.ndarray, np.ndarray, np.ndarray]:
    """Extract the CSR representation from a GetFEM sparse matrix."""

    A = gf.Spmat("copy", M)
    A.transpose()
    A.to_csc()

    indptr, indices = A.csc_ind()
    values = A.csc_val()

    indptr_arr = _ensure_int_array(indptr)
    indices_arr = _ensure_int_array(indices)
    values_arr = _ensure_scalar_array(values)

    shape = _infer_shape_from_spmat(A, indptr_arr, indices_arr)

    return shape, indptr_arr, indices_arr, values_arr


def create_petsc_aij_from_csr(
    shape: Tuple[int, int],
    indptr: Sequence[int],
    indices: Sequence[int],
    values: Sequence[float],
    *,
    comm: PETSc.Comm = comm,
    mat: Optional[PETSc.Mat] = None,
    preferred_types: Optional[Iterable[str]] = None,
) -> PETSc.Mat:
    """Create or update a PETSc AIJ matrix from CSR data."""

    rowptr = _ensure_int_array(indptr)
    colidx = _ensure_int_array(indices)
    data = _ensure_scalar_array(values)

    if preferred_types is None:
        preferred_types = (PETSc.Mat.Type.AIJ,)

    if mat is None:
        mat = PETSc.Mat().create(comm=comm)
        mat.setSizes(shape)

        for mat_type in preferred_types:
            try:
                mat.setType(mat_type)
                break
            except PETSc.Error as exc:  # pragma: no cover - depends on PETSc build
                logger.warning(
                    "Failed to set PETSc matrix type '%s': %s. Falling back to next candidate.",
                    mat_type,
                    exc,
                )
        mat.setPreallocationCSR((rowptr, colidx))
    else:
        # Ensure sizes are up to date before setting new values
        mat.setSizes(shape)

    mat.setOption(PETSc.Mat.Option.NEW_NONZERO_ALLOCATION_ERR, False)
    mat.setValuesCSR(rowptr, colidx, data)
    mat.setOption(PETSc.Mat.Option.NEW_NONZERO_ALLOCATION_ERR, True)
    mat.assemble()

    return mat


def _collect_preferred_types(mat_type: Optional[str], use_gpu: bool) -> Tuple[str, ...]:
    if mat_type:
        return (mat_type,)

    if not use_gpu:
        return (PETSc.Mat.Type.AIJ,)

    gpu_candidates = []
    for candidate in ("AIJCUSPARSE", "AIJHIPSPARSE", "AIJVIENNACL"):
        value = getattr(PETSc.Mat.Type, candidate, None)
        if value is not None:
            gpu_candidates.append(value)

    gpu_candidates.append(PETSc.Mat.Type.AIJ)
    return tuple(gpu_candidates)


def convert_gmm_to_petsc(
    M: gf.Spmat,
    B: Optional[PETSc.Mat] = None,
    *,
    comm: PETSc.Comm = comm,
    use_gpu: bool = False,
    mat_type: Optional[str] = None,
) -> PETSc.Mat:
    """Convert a GetFEM matrix ``M`` to a PETSc ``Mat``.

    Args:
        M (SPMat GetFEM): matrix to transfer
        B (PETSc.Mat, optional): matrix to fill with the data from ``M``. If ``None``
            a new matrix is created.
        comm (MPI_Comm): MPI communicator
        use_gpu (bool): request a GPU matrix implementation when available
        mat_type (str, optional): force a specific PETSc matrix type

    Returns:
        PETSc.Mat: populated PETSc matrix
    """

    shape, indptr, indices, values = extract_getfem_csr(M)

    preferred_types = _collect_preferred_types(mat_type, use_gpu)

    return create_petsc_aij_from_csr(
        shape,
        indptr,
        indices,
        values,
        comm=comm,
        mat=B,
        preferred_types=preferred_types,
    )

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


class _ShellContext:
    """Helper context for PETSc matrix-free shell operators."""

    def __init__(self, action: Callable[[PETSc.Vec, PETSc.Vec], None], diagonal: Optional[PETSc.Vec] = None):
        self.action = action
        self.diagonal = diagonal

    def mult(self, mat: PETSc.Mat, x: PETSc.Vec, y: PETSc.Vec) -> None:
        self.action(x, y)


def create_shell_matrix(
    shape: Tuple[int, int],
    action: Callable[[PETSc.Vec, PETSc.Vec], None],
    *,
    comm: PETSc.Comm = comm,
    diagonal: Optional[PETSc.Vec] = None,
) -> PETSc.Mat:
    """Create a PETSc matrix-free shell matrix with an optional diagonal."""

    shell = PETSc.Mat().create(comm=comm)
    shell.setSizes(shape)
    shell.setType(PETSc.Mat.Type.PYTHON)
    shell.setPythonContext(_ShellContext(action, diagonal))
    shell.setUp()

    if diagonal is not None:
        shell.setDiagonal(diagonal, addv=PETSc.InsertMode.INSERT_VALUES)

    return shell

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
