# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2022 Ghislain Haine
#
# See the LICENSE file in the root directory for license information.
#
# github: https://github.com/g-haine/scrimp

from scrimp import *
from scrimp.utils.mesh import set_verbose_gf
from itertools import zip_longest


def ciao_eq():
    set_verbose_gf(0)

    # Init the distributed port-Hamiltonian system
    ciao = DPHS("real")

    # Set the domain (using the built-in geometry `Interval`)
    ciao.set_domain(Domain("Interval",{"L":1,"h":2}))    ## Define the variables and their discretizations

    # Define State/s`)
    states = [
    State("s1","","scalar-field",None,0),
    State("s2","","scalar-field",None,0),
    State("s3","","scalar-field",None,0),
    State("s4","","scalar-field",None,0)
    ]

    # Plot the Hamiltonian with the power supplied at the boundary
    ciao.plot_Hamiltonian(save_figure=True)

    return ciao

if __name__ == "__main__":
    ciao = ciao_eq()