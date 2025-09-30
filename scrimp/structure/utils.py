"""Utility helpers bridging symbolic structures and legacy objects."""

from __future__ import annotations

from typing import Any, Iterable, List

import sympy as sp

from .core import ConstitutiveRelation, DiracStructure, LagrangeSubspace

__all__ = ["sympy_to_getfem", "ensure_term", "ensure_brick"]


def sympy_to_getfem(expr: sp.Expr | Any) -> str:
    """Convert a SymPy expression to a GetFEM-ready string.

    The conversion uses :func:`sympy.printing.ccode` and applies a handful of
    substitutions so that the resulting string can be embedded in the legacy
    SCRIMP GWFL expressions.
    """

    if isinstance(expr, str):
        return expr
    expression = sp.sympify(expr)
    from sympy.printing.ccode import ccode

    code = ccode(sp.simplify(expression))
    # GetFEM expects power through `pow`, but accepts the C syntax produced by
    # sympy.  We only strip superfluous whitespace to obtain compact strings.
    return code.replace("\n", " ")


def _term_from_subspace(subspace: LagrangeSubspace, *, description: str | None = None):
    from scrimp.hamiltonian import Term

    desc = description or subspace.describe()
    expression = sympy_to_getfem(subspace.discrete_pairing())
    term = Term(
        description=desc,
        expression=expression,
        regions=list(subspace.regions),
        mesh_id=subspace.mesh_id,
        structure=subspace,
    )
    return term


def _term_from_dirac_structure(structure: DiracStructure) -> List[Any]:
    terms = []
    for index, subspace in enumerate(structure.subspaces):
        desc = f"{structure.describe()}[{index}]"
        terms.append(_term_from_subspace(subspace, description=desc))
    return terms


def ensure_term(term_like: Any) -> List[Any]:
    """Return a list of :class:`~scrimp.hamiltonian.Term` instances.

    ``term_like`` can be a legacy :class:`~scrimp.hamiltonian.Term`, a symbolic
    :class:`LagrangeSubspace` or a :class:`DiracStructure`.  The helper eases the
    integration in :class:`~scrimp.hamiltonian.Hamiltonian` while preserving
    backward compatibility with string-based expressions.
    """

    from scrimp.hamiltonian import Term

    if isinstance(term_like, Term):
        return [term_like]
    if isinstance(term_like, LagrangeSubspace):
        return [_term_from_subspace(term_like)]
    if isinstance(term_like, DiracStructure):
        return _term_from_dirac_structure(term_like)
    if isinstance(term_like, Iterable):
        collected: List[Any] = []
        for item in term_like:
            collected.extend(ensure_term(item))
        return collected
    raise TypeError(
        "Unsupported term specification, expected Term, LagrangeSubspace or DiracStructure"
    )


def ensure_brick(brick_like: Any):
    """Return a :class:`~scrimp.brick.Brick` from symbolic data."""

    from scrimp.brick import Brick

    if isinstance(brick_like, Brick):
        return brick_like
    if isinstance(brick_like, ConstitutiveRelation):
        form = sympy_to_getfem(brick_like.relation.lhs - brick_like.relation.rhs)
        metadata = brick_like.metadata
        regions = list(metadata.get("regions", []))
        mesh_id = int(metadata.get("mesh_id", 0))
        name = metadata.get("name", brick_like.name)
        position = metadata.get("position", "constitutive")
        linear = bool(metadata.get("linear", True))
        dt = bool(metadata.get("dt", False))
        explicit = bool(metadata.get("explicit", False))
        brick = Brick(
            name=name,
            form=form,
            regions=regions,
            mesh_id=mesh_id,
            linear=linear,
            dt=dt,
            position=position,
            explicit=explicit,
            structure=brick_like,
        )
        return brick
    if isinstance(brick_like, DiracStructure):
        raise TypeError(
            "DiracStructure instances should be decomposed before building bricks"
        )
    raise TypeError("Unsupported brick specification")
