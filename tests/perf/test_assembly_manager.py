import pytest

np = pytest.importorskip("numpy")
petsc4py = pytest.importorskip("petsc4py")
from petsc4py import PETSc

from scrimp.core.assembly import MatrixAssemblyManager, MatrixAssemblyOptions


def _csr_identity(n: int):
    indptr = np.arange(n + 1, dtype=PETSc.IntType)
    indices = np.arange(n, dtype=PETSc.IntType)
    values = np.ones(n, dtype=PETSc.ScalarType)
    return (n, n), indptr, indices, values


@pytest.mark.perf
def test_sparsity_cache_reuse():
    manager = MatrixAssemblyManager()
    shape, indptr, indices, values = _csr_identity(4)

    matrix = manager.assemble_from_csr("mass", shape, indptr, indices, values)
    assert matrix.getSize() == shape

    new_values = np.full_like(values, 2.0)
    matrix_again = manager.assemble_from_csr("mass", shape, indptr, indices, new_values)

    assert matrix_again is matrix
    _, _, assembled_values = matrix_again.getValuesCSR()
    np.testing.assert_allclose(assembled_values, new_values)

    stats = manager.stats
    assert stats.builds == 1
    assert stats.reuses == 1
    manager.close()


@pytest.mark.perf
def test_async_assembly_future():
    options = MatrixAssemblyOptions(asynchronous=True, max_workers=1)
    manager = MatrixAssemblyManager(options=options)
    shape, indptr, indices, values = _csr_identity(5)

    future = manager.assemble_from_csr(
        "stiffness", shape, indptr, indices, values, asynchronous=True
    )

    matrix = future.result()
    assert matrix is manager.get_matrix("stiffness")
    manager.close()


@pytest.mark.perf
def test_gpu_backend_preference():
    options = MatrixAssemblyOptions(use_gpu=True)
    manager = MatrixAssemblyManager(options=options)
    shape, indptr, indices, values = _csr_identity(3)

    matrix = manager.assemble_from_csr("gpu", shape, indptr, indices, values)
    assert matrix.getType() in tuple(manager.preferred_types())
    manager.close()
