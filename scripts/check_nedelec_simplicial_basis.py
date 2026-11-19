#!/usr/bin/env python3
"""Small exact-arithmetic reference for simplicial Nédélec first-kind spaces.

This is NOT a GetFEM replacement or an H(curl)-conformity test. It checks a
parameterized polynomial generating set and a symmetric topological nodal
interpolation scheme on a reference triangle/tetrahedron. These are the
mathematical building blocks for an eventual separate GetFEM source patch.

Conventions: first-kind subdegree k>=0, lowest-order Whitney element k=0.
DoF points: interior nodes of the simplex barycentric lattice of degree k+2.
Edge: tangential samples; triangular faces: two tangential directions;
volume: component samples. Only triangles and tetrahedra.
"""
from __future__ import annotations

import argparse
from itertools import product
from math import comb

import sympy as sy


def exponents(dim: int, degree: int, homogeneous: bool = False):
    for idx in product(range(degree + 1), repeat=dim):
        if (sum(idx) == degree if homogeneous else sum(idx) <= degree):
            yield idx


def poly_space(dim: int, k: int, *, first_kind: bool):
    """Return a polynomial vector basis as component/exponent coefficient dicts."""
    basis = []
    for alpha in exponents(dim, k):
        for component in range(dim):
            basis.append({(component, alpha): 1})
    if not first_kind:
        return basis
    for alpha in exponents(dim, k, homogeneous=True):
        # Standard Koszul generators; in 3D exclude yz terms divisible by x.
        a = list(alpha)
        xp = tuple(a[i] + (i == 0) for i in range(dim))
        yp = tuple(a[i] + (i == 1) for i in range(dim))
        if dim == 2:
            basis.append({(0, yp): -1, (1, xp): 1})
        else:
            zp = tuple(a[i] + (i == 2) for i in range(dim))
            basis.append({(0, yp): -1, (1, xp): 1})
            basis.append({(0, zp): -1, (2, xp): 1})
            if a[0] == 0:
                basis.append({(1, zp): -1, (2, yp): 1})
    return basis


def reference_vertex(dim: int, vertex: int) -> sy.Matrix:
    return sy.Matrix([int(vertex == j + 1) for j in range(dim)])


def nodal_functionals(dim: int, k: int):
    """Symmetric collection of Nédélec first-kind local nodal functionals."""
    lattice_degree = k + 2
    for bary in product(range(lattice_degree + 1), repeat=dim + 1):
        if sum(bary) != lattice_degree:
            continue
        present = [i for i, a in enumerate(bary) if a != 0]
        x = sy.Matrix([sy.Rational(bary[j + 1], lattice_degree) for j in range(dim)])
        if len(present) == 1:  # simplex vertex
            continue
        if len(present) == 2:  # open edge: one tangent
            yield x, reference_vertex(dim, present[1]) - reference_vertex(dim, present[0])
        elif dim == 3 and len(present) == 3:  # open triangular face: two tangents
            for vertex in present[1:]:
                yield x, reference_vertex(dim, vertex) - reference_vertex(dim, present[0])
        elif len(present) == dim + 1:  # cell interior
            for component in range(dim):
                yield x, sy.eye(dim)[:, component]
        else:
            raise RuntimeError(f"Unexpected simplex entity: {bary}")


def vector_evaluate(poly: dict, position: sy.Matrix, dim: int) -> sy.Matrix:
    value = sy.zeros(dim, 1)
    for (component, alpha), coefficient in poly.items():
        monomial = coefficient
        for x, exp in zip(position, alpha):
            monomial *= x**exp
        value[component] += monomial
    return value


def n1_dimension(dim: int, k: int) -> int:
    return ((k + 1) * (k + 3)) if dim == 2 else (k + 1) * (k + 3) * (k + 4) // 2


def n2_dimension(dim: int, k: int) -> int:
    return dim * comb(k + dim, dim)


def audit(dim: int, k: int, *, exact: bool) -> tuple[int, int, int | None]:
    basis = poly_space(dim, k, first_kind=True)
    functionals = list(nodal_functionals(dim, k))
    expected = n1_dimension(dim, k)
    assert len(basis) == expected, (dim, k, "generators", len(basis), expected)
    assert len(functionals) == expected, (dim, k, "functionals", len(functionals), expected)
    rank = None
    if exact:
        mat = sy.Matrix(
            [[t.dot(vector_evaluate(f, x, dim)) for f in basis]
             for x, t in functionals]
        )
        rank = mat.rank()
        assert rank == expected, (dim, k, "rank", rank, expected)
    assert len(poly_space(dim, k, first_kind=False)) == n2_dimension(dim, k)
    return len(basis), len(functionals), rank



def n2_functionals(dim: int, k: int):
    """Pointwise tangent DoFs for experimental nodal Nedelec-II (k=1,2)."""
    assert dim in (2, 3) and k in (1, 2)
    n = k + 2
    vertices = [reference_vertex(dim, v) for v in range(dim + 1)]
    for bary in product(range(n + 1), repeat=dim + 1):
        if sum(bary) != n:
            continue
        present = [v for v, b in enumerate(bary) if b != 0]
        x = sum((sy.Rational(bary[v], n)*vertices[v]
                 for v in range(dim + 1)), sy.zeros(dim, 1))
        if len(present) == 2:
            yield x, vertices[present[1]] - vertices[present[0]]
        elif len(present) == 3 and k == 2:
            dominant = max(present, key=lambda v: bary[v])
            centroid = sum((vertices[v] for v in present),
                           sy.zeros(dim, 1)) / 3
            yield x, vertices[dominant] - centroid


def audit_n2(dim: int, k: int) -> int:
    """Exact rank of the nodal reference interpolation matrix."""
    basis = poly_space(dim, k, first_kind=False)
    functionals = list(n2_functionals(dim, k))
    n = n2_dimension(dim, k)
    assert len(basis) == len(functionals) == n
    matrix = sy.Matrix([
        [t.dot(vector_evaluate(p, x, dim)) for p in basis]
        for x, t in functionals
    ])
    rank = matrix.rank()
    assert rank == n, (dim, k, rank, n)
    return rank


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-order", type=int, default=2,
                        help="N1 orders 0..k; symbolic rank is tested up to k=2")
    args = parser.parse_args()
    if not 0 <= args.max_order <= 5:
        parser.error("Please use --max-order between 0 and 5")
    print("Reference simplicial Nedelec (NOT compiled GetFEM):")
    for dim in (2, 3):
        for k in range(args.max_order + 1):
            count, nfunc, rank = audit(dim, k, exact=(k <= 2))
            status = f"exact rank={rank}" if rank is not None else "count only"
            print(f"N1 dim={dim} k={k}: basis={count}, functionals={nfunc}; {status}")
    print("Nedelec-II nodal reference ranks (k=1,2):")
    for dim in (2, 3):
        for k in (1, 2):
            rank = audit_n2(dim, k)
            print(f"N2 dim={dim} k={k}: exact rank={rank}/{n2_dimension(dim, k)}")
    print("Reference spaces only: no H(curl) conformity or MPI claim.")


if __name__ == "__main__":
    main()
