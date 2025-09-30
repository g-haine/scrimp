# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2025 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""
- file:             __init__.py
- author:           Giuseppe Ferraro, Ghislain Haine
- date:             23 jun. 2023
- brief:            functions to initialize SCRIMP
"""

import os
import petsc4py
import sys

try:  # pragma: no cover - gmsh is optional during testing environments
    import gmsh  # type: ignore
except ModuleNotFoundError:  # pragma: no cover
    gmsh = None

petsc4py.init(sys.argv)
from petsc4py import PETSc

comm = PETSc.COMM_WORLD
rank = comm.getRank()

import scrimp.utils.config
scrimp.utils.config.set_paths()
scrimp.utils.config.set_verbose(1)

from scrimp.core import IORegistry, StateSpace, SystemTopology, TimeIntegrator, TopologyBuilder
from scrimp.dphs import DPHS
from scrimp.domain import Domain
from scrimp.state import State
from scrimp.costate import CoState
from scrimp.port import Parameter, Port
from scrimp.fem import FEM
from scrimp.control import Control_Port
from scrimp.brick import Brick
from scrimp.hamiltonian import Term, Hamiltonian
from scrimp.io import schema_loader

__all__ = [
    "DPHS",
    "Domain",
    "State",
    "CoState",
    "Parameter",
    "Port",
    "FEM",
    "Control_Port",
    "Brick",
    "Term",
    "Hamiltonian",
    "IORegistry",
    "StateSpace",
    "SystemTopology",
    "TimeIntegrator",
    "TopologyBuilder",
]
