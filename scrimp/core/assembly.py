"""Assembly helpers and caches for PETSc matrices."""

from __future__ import annotations

import threading
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from typing import Dict, Iterable, Optional, Tuple

import numpy as np
from petsc4py import PETSc

from scrimp.utils.linalg import create_petsc_aij_from_csr, extract_getfem_csr


@dataclass
class SparsityPattern:
    shape: Tuple[int, int]
    indptr: np.ndarray
    indices: np.ndarray

    def matches(self, shape: Tuple[int, int], indptr: np.ndarray, indices: np.ndarray) -> bool:
        if self.shape != shape:
            return False
        return np.array_equal(self.indptr, indptr) and np.array_equal(self.indices, indices)


@dataclass
class MatrixAssemblyOptions:
    use_gpu: bool = False
    mat_type: Optional[str] = None
    asynchronous: bool = False
    max_workers: int = 2


@dataclass
class AssemblyStatistics:
    builds: int = 0
    reuses: int = 0


class MatrixAssemblyManager:
    """Manage PETSc matrix assembly using cached sparsity patterns."""

    def __init__(
        self,
        *,
        comm: PETSc.Comm = PETSc.COMM_WORLD,
        options: Optional[MatrixAssemblyOptions] = None,
    ) -> None:
        self.comm = comm
        self.options = options or MatrixAssemblyOptions()
        self._patterns: Dict[str, SparsityPattern] = {}
        self._matrices: Dict[str, PETSc.Mat] = {}
        self._lock = threading.RLock()
        self._executor: Optional[ThreadPoolExecutor] = None
        self._executor_workers: Optional[int] = None
        self._stats = AssemblyStatistics()

    # ------------------------------------------------------------------
    # lifecycle helpers
    # ------------------------------------------------------------------
    def close(self) -> None:
        if self._executor is not None:
            self._executor.shutdown(wait=False)
            self._executor = None
            self._executor_workers = None

    # ------------------------------------------------------------------
    # option helpers
    # ------------------------------------------------------------------
    def update_options(self, **kwargs) -> None:
        for key, value in kwargs.items():
            if hasattr(self.options, key):
                setattr(self.options, key, value)
        if not self.options.asynchronous and self._executor is not None:
            self.close()
        elif self.options.asynchronous and self._executor is not None:
            if (
                self._executor_workers is not None
                and self.options.max_workers != self._executor_workers
            ):
                self.close()

    # ------------------------------------------------------------------
    # statistics
    # ------------------------------------------------------------------
    @property
    def stats(self) -> AssemblyStatistics:
        return AssemblyStatistics(self._stats.builds, self._stats.reuses)

    # ------------------------------------------------------------------
    # matrix retrieval
    # ------------------------------------------------------------------
    def get_matrix(self, key: str) -> Optional[PETSc.Mat]:
        with self._lock:
            return self._matrices.get(key)

    def register_matrix(self, key: str, matrix: PETSc.Mat) -> None:
        with self._lock:
            self._matrices[key] = matrix

    def _preferred_types(self) -> Iterable[str]:
        if self.options.mat_type:
            return (self.options.mat_type,)
        if self.options.use_gpu:
            return tuple(
                t
                for t in (
                    getattr(PETSc.Mat.Type, "AIJCUSPARSE", None),
                    getattr(PETSc.Mat.Type, "AIJHIPSPARSE", None),
                    getattr(PETSc.Mat.Type, "AIJVIENNACL", None),
                    PETSc.Mat.Type.AIJ,
                )
                if t is not None
            )
        return (PETSc.Mat.Type.AIJ,)

    def preferred_types(self) -> Iterable[str]:
        return self._preferred_types()

    # ------------------------------------------------------------------
    # assembly routines
    # ------------------------------------------------------------------
    def assemble_from_getfem(
        self,
        key: str,
        matrix,
        *,
        target: Optional[PETSc.Mat] = None,
        asynchronous: bool = False,
    ):
        shape, indptr, indices, values = extract_getfem_csr(matrix)
        return self.assemble_from_csr(
            key,
            shape,
            indptr,
            indices,
            values,
            target=target,
            asynchronous=asynchronous,
        )

    def assemble_from_csr(
        self,
        key: str,
        shape: Tuple[int, int],
        indptr: np.ndarray,
        indices: np.ndarray,
        values: np.ndarray,
        *,
        target: Optional[PETSc.Mat] = None,
        asynchronous: bool = False,
    ):
        if asynchronous and self.options.asynchronous:
            if self._executor is None:
                self._executor = ThreadPoolExecutor(max_workers=self.options.max_workers)
                self._executor_workers = self.options.max_workers
            return self._executor.submit(
                self._assemble_internal,
                key,
                shape,
                indptr.copy(),
                indices.copy(),
                values.copy(),
                target,
            )

        return self._assemble_internal(key, shape, indptr, indices, values, target)

    # ------------------------------------------------------------------
    # internal implementation
    # ------------------------------------------------------------------
    def _assemble_internal(
        self,
        key: str,
        shape: Tuple[int, int],
        indptr: np.ndarray,
        indices: np.ndarray,
        values: np.ndarray,
        target: Optional[PETSc.Mat],
    ) -> PETSc.Mat:
        with self._lock:
            pattern = self._patterns.get(key)
            preferred_types = self._preferred_types()

            if pattern is None or not pattern.matches(shape, indptr, indices):
                matrix = create_petsc_aij_from_csr(
                    shape,
                    indptr,
                    indices,
                    values,
                    comm=self.comm,
                    mat=target or self._matrices.get(key),
                    preferred_types=preferred_types,
                )
                self._patterns[key] = SparsityPattern(shape, indptr.copy(), indices.copy())
                self._matrices[key] = matrix
                self._stats.builds += 1
                return matrix

            matrix = target or self._matrices.get(key)
            if matrix is None:
                matrix = create_petsc_aij_from_csr(
                    shape,
                    indptr,
                    indices,
                    values,
                    comm=self.comm,
                    mat=None,
                    preferred_types=preferred_types,
                )
                self._matrices[key] = matrix
                self._patterns[key] = SparsityPattern(shape, indptr.copy(), indices.copy())
                self._stats.builds += 1
                return matrix

            matrix.setOption(PETSc.Mat.Option.NEW_NONZERO_ALLOCATION_ERR, False)
            matrix.setValuesCSR(indptr, indices, values)
            matrix.setOption(PETSc.Mat.Option.NEW_NONZERO_ALLOCATION_ERR, True)
            matrix.assemble()
            self._stats.reuses += 1
            return matrix


__all__ = [
    "AssemblyStatistics",
    "MatrixAssemblyManager",
    "MatrixAssemblyOptions",
    "SparsityPattern",
]
