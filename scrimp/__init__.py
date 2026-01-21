# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
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

import importlib.util
import os
import sys
import types

def _install_stub(module_name: str, factory):
    if module_name not in sys.modules:
        sys.modules[module_name] = factory()


def _make_petsc_stub():
    module = types.ModuleType("petsc4py")

    def _init(*args, **kwargs):
        return None

    class _FakeComm:
        def getRank(self) -> int:
            return 0

    class _FakePETSc(types.SimpleNamespace):
        COMM_WORLD = _FakeComm()

    module.init = _init
    module.PETSc = _FakePETSc()
    sys.modules.setdefault("petsc4py.PETSc", module.PETSc)
    return module


def _make_getfem_stub():
    module = types.ModuleType("getfem")

    class _FakeModel:
        def __init__(self, *args, **kwargs):
            self.args = args
            self.kwargs = kwargs

    def _asm(*args, **kwargs):
        raise RuntimeError("getfem is not available in this environment")

    module.Model = _FakeModel
    module.asm = _asm
    module.util_trace_level = lambda level: None
    module.util_warning_level = lambda level: None
    return module


_gmsh_module = sys.modules.get("gmsh")
if _gmsh_module is not None:
    _gmsh_spec = getattr(_gmsh_module, "__spec__", None)
else:
    _gmsh_spec = importlib.util.find_spec("gmsh")
if _gmsh_spec is not None:
    import gmsh  # type: ignore[import-untyped]
else:
    gmsh = None  # type: ignore[assignment]

_petsc_module = sys.modules.get("petsc4py")
if _petsc_module is not None:
    _petsc_spec = getattr(_petsc_module, "__spec__", None)
else:
    _petsc_spec = importlib.util.find_spec("petsc4py")
if _petsc_spec is None:
    _install_stub("petsc4py", _make_petsc_stub)

import petsc4py  # type: ignore[import-untyped]

try:  # pragma: no cover - gmsh is optional during testing environments
    import gmsh  # type: ignore
except ModuleNotFoundError:  # pragma: no cover
    gmsh = None

petsc4py.init(sys.argv)
from petsc4py import PETSc  # type: ignore[import-untyped]

comm = PETSc.COMM_WORLD
rank = comm.getRank()

_getfem_module = sys.modules.get("getfem")
if _getfem_module is not None:
    _getfem_spec = getattr(_getfem_module, "__spec__", None)
else:
    _getfem_spec = importlib.util.find_spec("getfem")
if _getfem_spec is None:
    _install_stub("getfem", _make_getfem_stub)

import scrimp.utils.config
if _getfem_spec is not None and _petsc_spec is not None:
    scrimp.utils.config.set_paths()
    scrimp.utils.config.set_verbose(1)

_slepc_spec = importlib.util.find_spec("slepc4py")
_runtime_ready = _getfem_spec is not None and _petsc_spec is not None and _slepc_spec is not None

if _runtime_ready:
    from scrimp.dphs import DPHS
    from scrimp.domain import Domain
    from scrimp.state import State
    from scrimp.costate import CoState
    from scrimp.fem import FEM
    from scrimp.control import Control_Port
else:  # pragma: no cover - optional runtime features
    DPHS = Domain = State = CoState = FEM = Control_Port = None

from scrimp.dphs import DPHS
from scrimp.domain import Domain
from scrimp.state import State
from scrimp.costate import CoState
from scrimp.port import Parameter, Port
from scrimp.brick import Brick

from scrimp.hamiltonian import Term, Hamiltonian

