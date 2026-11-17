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
from contextlib import nullcontext
import getfem as gf
import numpy as np
import petsc4py
from slepc4py import SLEPc

petsc4py.init(sys.argv)
from petsc4py import PETSc
PETSc.Options().setValue('info', None)

comm = PETSc.COMM_WORLD
import scipy.sparse as sp

def convert_gmm_to_petsc(
    M, B, comm=comm, *, profiler=None, profile_label=None,
    distributed_coo_eligible=False,
):
    """Transfer GetFEM's globally numbered matrix into distributed PETSc.

    SCRIMP_GETFEM_MATRIX_MODE selects ASSEMBLY SEMANTICS:
    - replicated: insert owned PETSc rows once (INSERT_VALUES).
    - distributed: sum all global contributions (ADD_VALUES).

    SCRIMP_GETFEM_BRIDGE_REPLICATED controls only extraction in
    replicated mode: 'full' (default) or experimental 'owned_rows'.

    For distributed mode, SCRIMP_GETFEM_BRIDGE_DISTRIBUTED=csr
    selects the historical IJV/ADD_VALUES path (also called ijv).
    SCRIMP_GETFEM_BRIDGE_DISTRIBUTED=coo enables
    native PETSc COO assembly only at call sites that explicitly declare
    distributed_coo_eligible=True (currently the initial LINEAR mass and
    stiffness matrices). Other calls retain the IJV/ADD_VALUES path,
    including state-dependent Jacobians and dynamic sparsity.
    This is an opt-in experimental path requiring MPI physical validation;
    it preserves the PETSc Mat handle and explicit diagonal slots.

    Optional profiler/profile_label use SCRIMP's existing TimingProfiler
    for bridge subphases. No MPI collectives or matrix changes are added.
    The bridge total is already measured by the DPHS assembly methods.
    """
    def measure(step):
        if profiler is None or profile_label is None:
            return nullcontext()
        return profiler.measure(f"{profile_label}.{step}")

    mode = os.environ.get(
        "SCRIMP_GETFEM_MATRIX_MODE", "replicated"
    ).strip().lower()
    if mode not in ("replicated", "distributed"):
        raise ValueError(
            "SCRIMP_GETFEM_MATRIX_MODE must be replicated or distributed"
        )
    extraction = os.environ.get(
        "SCRIMP_GETFEM_BRIDGE_REPLICATED", "full"
    ).strip().lower()
    if extraction not in ("full", "owned_rows"):
        raise ValueError(
            "SCRIMP_GETFEM_BRIDGE_REPLICATED must be full or owned_rows"
        )

    distributed_method = os.environ.get(
        "SCRIMP_GETFEM_BRIDGE_DISTRIBUTED", "ijv"
    ).strip().lower()
    if distributed_method not in ("ijv", "csr", "coo"):
        raise ValueError(
            "SCRIMP_GETFEM_BRIDGE_DISTRIBUTED must be csr, coo or ijv "
            "(legacy alias for csr)"
        )
    use_coo = (
        mode == "distributed" and distributed_method == "coo"
        and distributed_coo_eligible
    )

    nrows, ncols = B.getSize()
    lo, hi = B.getOwnershipRange()
    local_extraction = mode == "replicated" and extraction == "owned_rows"
    if local_extraction:
        # GetFEM Spmat('copy', M, rows, columns) preserves the global
        # column space; selecting only local columns loses couplings.
        gf_shape = tuple(int(d) for d in M.size())
        if gf_shape != (nrows, ncols):
            raise ValueError(
                f"GetFEM/PETSc shape mismatch on rank {comm.getRank()}: "
                f"{gf_shape} vs {(nrows, ncols)}"
            )
        if hi > lo:
            with measure("getfem.copy"):
                A = gf.Spmat("copy", M, range(lo, hi), range(ncols))
            with measure("getfem.transpose"):
                A.transpose()
            with measure("getfem.csr_extract"):
                row_ptr, col_idx = A.csc_ind()
                values = A.csc_val()
        else:
            # Avoid relying on GetFEM's empty-submatrix convention.
            row_ptr = np.zeros(1, dtype=PETSc.IntType)
            col_idx = np.empty(0, dtype=PETSc.IntType)
            values = np.empty(0, dtype=PETSc.ScalarType)
        extracted_rows = hi - lo
    else:
        # Stable reference AND the mandatory distributed path.
        with measure("getfem.copy"):
            A = gf.Spmat("copy", M)
        with measure("getfem.transpose"):
            A.transpose()
        with measure("getfem.csr_extract"):
            row_ptr, col_idx = A.csc_ind()
            values = A.csc_val()
        extracted_rows = nrows

    with measure("csr_validate"):
        row_ptr = np.asarray(row_ptr)
        col_idx = np.asarray(col_idx)
        values = np.asarray(values)
        if len(row_ptr) != extracted_rows + 1:
            raise ValueError(
                f"GetFEM/PETSc row count mismatch on rank {comm.getRank()}: "
                f"{len(row_ptr)-1} vs expected {extracted_rows}"
            )
        if (not len(row_ptr) or int(row_ptr[0]) != 0
                or len(col_idx) != len(values)
                or int(row_ptr[-1]) != len(values)):
            raise ValueError("Invalid GetFEM CSR data (row pointer/count)")
        if len(col_idx) and (min(col_idx) < 0 or max(col_idx) >= ncols):
            raise ValueError("GetFEM CSR has columns outside the PETSc matrix")

    trace = os.environ.get("SCRIMP_MPI_TRACE") == "1"
    if trace:
        print(
            f"[MPI {comm.getRank()}/{comm.getSize()}] "
            f"GetFEM -> PETSc: global_rows={nrows}, "
            f"extracted_rows={extracted_rows}, "
            f"nnz_stored_extracted={len(values)}, "
            f"nnz_nonzero_extracted={np.count_nonzero(values)}, "
            f"abs_sum_extracted={np.sum(np.abs(values)):.9g}, "
            f"mode={mode}, extraction={extraction}",
            file=sys.stderr, flush=True,
        )

    if use_coo:
        # PETSc COO accepts off-process global rows and duplicate
        # coordinates. It sums GetFEM contributions with the same
        # semantics as the distributed IJV ADD_VALUES transfer.
        #
        # SCRIMP.init_matrix() initially creates structural diagonal
        # slots even when absent in GetFEM. Preserve those entries in
        # the COO pattern to keep exact stored CSR and MUMPS semantics.
        # This replaces the graph in the SAME PETSc Mat object.
        with measure("petsc.coo.prepare"):
            counts = np.diff(row_ptr).astype(np.intp, copy=False)
            rows = np.repeat(
                np.arange(nrows, dtype=PETSc.IntType), counts
            )
            diagonal = np.arange(lo, hi, dtype=PETSc.IntType)
            coo_i = np.ascontiguousarray(
                np.concatenate((rows, diagonal))
            )
            coo_j = np.ascontiguousarray(
                np.concatenate((col_idx, diagonal))
            )
            coo_v = np.ascontiguousarray(
                np.concatenate((
                    values,
                    np.zeros(hi - lo, dtype=PETSc.ScalarType),
                ))
            )
        with measure("petsc.coo.preallocate"):
            B.setPreallocationCOO(coo_i, coo_j)
        with measure("petsc.coo.set_values"):
            # INSERT overwrites all previous entries. Repeated COO
            # coordinates are still summed, even with INSERT_VALUES.
            # PETSc performs its own final assembly here.
            B.setValuesCOO(
                coo_v, addv=PETSc.InsertMode.INSERT_VALUES
            )
        return

    # Old values must not accumulate across Newton reassemblies.
    with measure("petsc.zero_entries"):
        B.zeroEntries()
    B.setOption(PETSc.Mat.Option.NEW_NONZERO_ALLOCATION_ERR, False)
    try:
        with measure("petsc.insert"):
            if mode == "replicated" and hi > lo:
                if local_extraction:
                    local_ptr, local_cols, local_values = (
                        row_ptr, col_idx, values
                    )
                else:
                    begin, end = int(row_ptr[lo]), int(row_ptr[hi])
                    local_ptr = row_ptr[lo:hi + 1] - begin
                    local_cols = col_idx[begin:end]
                    local_values = values[begin:end]
                B.setValuesIJV(
                    local_ptr, local_cols, local_values,
                    addv=PETSc.InsertMode.INSERT_VALUES,
                    rowmap=np.arange(lo, hi, dtype=PETSc.IntType),
                )
            elif mode == "distributed" and len(values):
                # Overlapping and off-process rows MUST be summed globally.
                B.setValuesIJV(
                    row_ptr, col_idx, values,
                    addv=PETSc.InsertMode.ADD_VALUES,
                    rowmap=np.arange(nrows, dtype=PETSc.IntType),
                )
        if trace:
            print(f"[MPI {comm.getRank()}] before PETSc Mat.assemble",
                  file=sys.stderr, flush=True)
        with measure("petsc.assemble"):
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
