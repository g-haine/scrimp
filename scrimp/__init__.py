# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2023 ISAE-SUPAERO -- GNU GPLv3
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

module_path = __file__[:-18]

for path in [os.path.join("outputs", "log"), os.path.join("outputs", "png"), os.path.join("outputs", "mesh")]:
    composed_path = os.path.join(module_path, path)
    if not os.path.isdir(composed_path):
        os.makedirs(composed_path)

from scrimp.utils.config import set_verbose
set_verbose()

from scrimp.dphs import DPHS
from scrimp.domain import Domain
from scrimp.state import State
from scrimp.costate import CoState
from scrimp.port import Parameter, Port
from scrimp.fem import FEM
from scrimp.control import Control_Port
from scrimp.brick import Brick
from scrimp.hamiltonian import Term, Hamiltonian
