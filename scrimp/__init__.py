from scrimp.state import State
from scrimp.costate import CoState
from scrimp.domain import Domain
from scrimp.fem import FEM
from scrimp.brick import Brick
from scrimp.port import Port, Parameter
from scrimp.control import Control_Port
from scrimp.hamiltonian import Term, Hamiltonian
from scrimp.dphs import DPHS


import os

module_path = __file__[:-18]

for path in [os.path.join("outputs", "log"), os.path.join("outputs", "png")]:
    composed_path = os.path.join(module_path, path)
    if not os.path.isdir(composed_path):
        os.makedirs(composed_path)
