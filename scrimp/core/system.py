"""Light-weight containers describing the structural components of a DPHS model."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, TYPE_CHECKING, Any

from scrimp.hamiltonian import Hamiltonian

if TYPE_CHECKING:  # pragma: no cover - imported only for type checking
    from scrimp.brick import Brick
    from scrimp.control import Control_Port
    from scrimp.costate import CoState
    from scrimp.domain import Domain
    from scrimp.port import Port
    from scrimp.state import State


@dataclass
class SystemTopology:
    """Store the meshes, bricks and domains making up a model."""

    domain: Optional["Domain"] = None
    bricks: Dict[str, "Brick"] = field(default_factory=dict)

    def set_domain(self, domain: "Domain") -> None:
        """Attach a domain to the topology."""

        self.domain = domain


@dataclass
class StateSpace:
    """Keep track of state and costate variables registered in the model."""

    states: Dict[str, "State"] = field(default_factory=dict)
    costates: Dict[str, "CoState"] = field(default_factory=dict)

    def register_state(self, state: "State") -> None:
        """Register a new state instance."""

        self.states[state.get_name()] = state

    def register_costate(self, costate: "CoState") -> None:
        """Register a new costate instance."""

        self.costates[costate.get_name()] = costate


@dataclass
class IORegistry:
    """Collect boundary ports, control ports and Hamiltonian terms."""

    ports: Dict[str, "Port"] = field(default_factory=dict)
    controls: Dict[str, "Control_Port"] = field(default_factory=dict)
    hamiltonian: Hamiltonian = field(default_factory=lambda: Hamiltonian("Hamiltonian"))
    powers_computed: bool = False

    def register_port(self, port: "Port") -> None:
        self.ports[port.get_name()] = port

    def register_control(self, control: "Control_Port") -> None:
        self.controls[control.get_name()] = control


@dataclass
class TimeIntegrator:
    """Handle PETSc TS option bookkeeping and solution storage."""

    options: Dict[str, Any] = field(default_factory=dict)
    initial_values: Dict[str, bool] = field(default_factory=dict)
    solution: Dict[str, List[Any]] = field(
        default_factory=lambda: {"t": [], "z": []}
    )
    stop: bool = False
    ts_start: float = 0.0

    def reset_options(self) -> None:
        """Clear the PETSc options database."""

        self.options = {}

    def reset_solution(self) -> None:
        """Clear previously stored solutions."""

        self.solution = {"t": [], "z": []}

