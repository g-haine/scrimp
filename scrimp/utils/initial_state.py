# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Synchronize accepted PETSc pseudo-initial states to GetFEM variables.

Keep the actual data transfer independent of heavy optional GetFEM/gmsh
imports so the exact operation is unit-testable on lightweight CI.
"""

import numpy as np


def publish_accepted_pseudo_state(scatter, accepted, global_vec, gf_model):
    """Gather the full accepted state and publish it to GetFEM.

    The main TS constructs its initial PETSc vector from the GetFEM
    variables, not from the pseudo-TS's accepted vector. Every MPI rank
    must publish the same full vector, regardless of linearity flags.
    """
    scatter.scatter(accepted, global_vec)
    values = np.asarray(
        global_vec.getArray(readonly=True).real,
        dtype=np.float64,
    ).copy()
    gf_model.to_variables(values)
