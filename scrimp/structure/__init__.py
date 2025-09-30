"""Symbolic structural helpers for SCRIMP."""

from .core import ConstitutiveRelation, DiracStructure, LagrangeSubspace
from .utils import (
    ensure_term,
    ensure_brick,
    sympy_to_getfem,
)

__all__ = [
    "LagrangeSubspace",
    "ConstitutiveRelation",
    "DiracStructure",
    "ensure_term",
    "ensure_brick",
    "sympy_to_getfem",
]
