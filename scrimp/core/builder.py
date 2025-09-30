"""Declarative helpers to populate :mod:`scrimp.dphs` models."""

from __future__ import annotations

from typing import Callable, Optional

from scrimp.core.system import IORegistry, StateSpace, SystemTopology
from scrimp.costate import CoState
from scrimp.port import Port
from scrimp.state import State


class TopologyBuilder:
    """Fluent helper that wires states, costates and ports together."""

    def __init__(
        self,
        topology: Optional[SystemTopology] = None,
        state_space: Optional[StateSpace] = None,
        io_registry: Optional[IORegistry] = None,
        port_callback: Optional[Callable[[Port], None]] = None,
    ) -> None:
        self.topology = topology or SystemTopology()
        self.state_space = state_space or StateSpace()
        self.io_registry = io_registry or IORegistry()
        self._port_callback = port_callback

    def with_domain(self, domain) -> "TopologyBuilder":
        """Attach a domain to the system topology."""

        self.topology.set_domain(domain)
        return self

    def add_state(self, state: State) -> "TopologyBuilder":
        """Register a state variable."""

        self.state_space.register_state(state)
        return self

    def add_costate(self, costate: CoState) -> "TopologyBuilder":
        """Register a co-state and automatically generate its dynamical port."""

        state = costate.get_state()
        if state.get_costate() is None:
            state.set_costate(costate)
        self.state_space.register_costate(costate)

        port = Port(
            state.get_name(),
            state.get_name(),
            costate.get_name(),
            costate.get_kind(),
            state.get_mesh_id(),
            algebraic=False,
            dissipative=False,
            substituted=costate.get_substituted(),
            region=state.get_region(),
        )
        self.add_port(port)
        state.set_port(port)
        costate.set_port(port)
        return self

    def add_port(self, port: Port) -> "TopologyBuilder":
        """Register a (possibly external) port."""

        self.io_registry.register_port(port)
        if self._port_callback is not None:
            self._port_callback(port)
        return self

    def build(self) -> tuple[SystemTopology, StateSpace, IORegistry]:
        """Return the composed managers."""

        return self.topology, self.state_space, self.io_registry

