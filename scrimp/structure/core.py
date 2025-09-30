"""Symbolic discrete structure definitions for SCRIMP."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Mapping, MutableMapping, Optional, Sequence

import sympy as sp

__all__ = [
    "LagrangeSubspace",
    "ConstitutiveRelation",
    "DiracStructure",
]


def _ensure_symbol(value: Any) -> sp.Symbol:
    """Return a :class:`sympy.Symbol` for ``value``."""
    if isinstance(value, sp.Symbol):
        return value
    if isinstance(value, str):
        return sp.Symbol(value)
    raise TypeError(f"Cannot create sympy.Symbol from {value!r}")


def _ensure_expr(value: Any) -> sp.Expr:
    """Return a :class:`sympy.Expr` for ``value``."""
    if isinstance(value, sp.Expr):
        return value
    return sp.sympify(value)


@dataclass(frozen=True)
class LagrangeSubspace:
    """A symbolic description of a discrete Lagrange subspace.

    The object stores a flow variable, an effort variable and a density term
    describing how the discrete pairing is computed.  Metadata captures the
    geometric information required to build GetFEM expressions (mesh id,
    region, quadrature order, ...).
    """

    flow: sp.Symbol | str
    effort: sp.Symbol | str
    pairing_density: Any = 1
    metadata: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        object.__setattr__(self, "flow", _ensure_symbol(self.flow))
        object.__setattr__(self, "effort", _ensure_symbol(self.effort))
        object.__setattr__(self, "pairing_density", _ensure_expr(self.pairing_density))
        object.__setattr__(self, "metadata", dict(self.metadata))

    @property
    def regions(self) -> Sequence[int]:
        """Region identifiers where the pairing is evaluated."""
        regions = self.metadata.get("regions", ())
        if isinstance(regions, (list, tuple)):
            return regions
        if regions is None:
            return ()
        return (regions,)

    @property
    def mesh_id(self) -> int:
        """Mesh identifier for the pairing."""
        return int(self.metadata.get("mesh_id", 0))

    @property
    def sign(self) -> int:
        """Sign associated with the pairing contribution."""
        return int(self.metadata.get("sign", 1))

    def discrete_pairing(self) -> sp.Expr:
        """Return the symbolic discrete pairing density."""
        measure = _ensure_expr(self.metadata.get("measure", 1))
        return sp.simplify(self.sign * self.flow * self.effort * self.pairing_density * measure)

    def describe(self) -> str:
        """Human readable description."""
        default = f"<{self.flow}, {self.effort}>"
        return str(self.metadata.get("description", default))


@dataclass(frozen=True)
class ConstitutiveRelation:
    """Symbolic constitutive relation between flow and effort variables."""

    name: str
    relation: sp.Equality | Sequence[Any] | Mapping[str, Any]
    metadata: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        relation = self._build_relation(self.relation)
        object.__setattr__(self, "relation", relation)
        object.__setattr__(self, "metadata", dict(self.metadata))

    @staticmethod
    def _build_relation(value: Any) -> sp.Equality:
        if isinstance(value, sp.Equality):
            return value
        if isinstance(value, Mapping):
            if len(value) != 1:
                raise ValueError("Mapping constitutive relation must contain a single entry")
            symbol, expression = next(iter(value.items()))
            return sp.Eq(_ensure_symbol(symbol), _ensure_expr(expression))
        if isinstance(value, Sequence) and len(value) == 2:
            lhs, rhs = value
            return sp.Eq(_ensure_expr(lhs), _ensure_expr(rhs))
        raise TypeError("Unsupported constitutive relation specification")

    @property
    def variables(self) -> List[sp.Symbol]:
        return list(sorted(self.relation.free_symbols, key=lambda s: s.name))

    def substitution_map(self, solve_for: Optional[sp.Symbol | str] = None) -> Dict[sp.Symbol, sp.Expr]:
        """Return a substitution map derived from the relation."""

        eq = self.relation
        symbol = solve_for or self.metadata.get("solve_for")
        if symbol is None:
            if isinstance(eq.lhs, sp.Symbol):
                symbol = eq.lhs
            elif isinstance(eq.rhs, sp.Symbol):
                symbol = eq.rhs
            else:
                raise ValueError(
                    "Unable to infer substitution target, please pass 'solve_for' metadata"
                )
        symbol = _ensure_symbol(symbol)
        solutions = sp.solve(eq, symbol, dict=True)
        if not solutions:
            raise ValueError(f"Could not solve constitutive relation for {symbol!s}")
        expression = sp.simplify(solutions[0][symbol])
        return {symbol: expression}

    def describe(self) -> str:
        return self.metadata.get("description", self.name)


@dataclass
class DiracStructure:
    """A collection of Lagrange subspaces and constitutive relations."""

    subspaces: Sequence[LagrangeSubspace]
    constitutive_relations: Sequence[ConstitutiveRelation] = field(default_factory=list)
    metadata: MutableMapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.subspaces = list(self.subspaces)
        self.constitutive_relations = list(self.constitutive_relations)

    def total_power(self) -> sp.Expr:
        return sp.simplify(sum(subspace.discrete_pairing() for subspace in self.subspaces))

    def validate_power_balance(self) -> sp.Expr:
        power = self.total_power()
        if sp.simplify(power) != 0:
            raise ValueError(
                "DiracStructure power balance check failed: total power does not vanish"
            )
        return power

    def substitutions(self) -> Dict[sp.Symbol, sp.Expr]:
        mapping: Dict[sp.Symbol, sp.Expr] = {}
        for relation in self.constitutive_relations:
            mapping.update(relation.substitution_map())
        return mapping

    def describe(self) -> str:
        return self.metadata.get("description", "DiracStructure")

    def regions(self) -> Sequence[int]:
        regions: List[int] = []
        for subspace in self.subspaces:
            regions.extend(subspace.regions)
        return regions

    def mesh_ids(self) -> Sequence[int]:
        return sorted({subspace.mesh_id for subspace in self.subspaces})
