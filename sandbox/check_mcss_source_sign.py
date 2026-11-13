# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Isolate GetFEM RHS sign conventions used by MCSS2024 set_control().

The legacy DPHS.add_brick/set_control uses TWO different kinds of
GetFEM bricks: nonlinear weak forms and source_term_brick. The
production IFunction recently changed from '+ model.rhs()' to
'- model.rhs()'. If source bricks enter rhs with the usual POSITIVE
sign while nonlinear weak forms enter with NEGATIVE sign, the new
blanket negative sign silently reverses prescribed controls.

No SCRIMP state, PETSc TS, SNES, KSP or MUMPS is constructed.
Run: mpiexec -n 1 python -u sandbox/check_mcss_source_sign.py
"""

import json

import getfem as gf
import numpy as np

mesh = gf.Mesh("cartesian", np.linspace(0.0, 1.0, 17))
mf = gf.MeshFem(mesh, 1)
mf.set_fem(gf.Fem("FEM_PK(1,1)"))
mim = gf.MeshIm(mesh, gf.Integ("IM_GAUSS1D(5)"))

# Nonlinear weak form N(u) = int(u^2 * v) dx at u=1.
nonlinear = gf.Model("real")
nonlinear.add_fem_variable("u", mf)
nonlinear.add_nonlinear_term(mim, "u*u*Test_u")
nonlinear.set_variable("u", np.ones(mf.nbdof()))
nonlinear.assembly(option="build_rhs")
nonlinear_rhs = float(np.real(np.asarray(nonlinear.rhs()).sum()))

# MCSS2024-style control: -u + desired = 0.
# The exact solution is u=desired=1. Evaluate its rhs at u=0,
# where the linear brick produces no state-dependent residual.
source = gf.Model("real")
source.add_fem_variable("u", mf)
source.add_initialized_data("desired", np.array([1.0]))
source.add_linear_term(mim, "-u*Test_u")
source.add_source_term_brick(mim, "u", "desired")
source.set_variable("u", np.zeros(mf.nbdof()))
source.assembly(option="build_rhs")
source_rhs = float(np.real(np.asarray(source.rhs()).sum()))

# At u=1 the linear brick contributes to the tangent, but the
# source RHS still has positive integral and no state derivative.
source.assembly(option="build_matrix")
import getfem as _gf
raw = _gf.Spmat("copy", source.tangent_matrix())
raw.transpose()
ptr, cols = (np.asarray(a) for a in raw.csc_ind())
vals = np.asarray(raw.csc_val())
tangent_integral = float(np.real(vals.sum()))

print(json.dumps({
    "mpi_rank_agnostic_example": "unit interval",
    "nonlinear_weak_form": "u*u*Test_u, u=1",
    "nonlinear_rhs_integral": nonlinear_rhs,
    "mcss_source_form": "-u*Test_u + desired*Test_u, u=0, desired=1",
    "source_rhs_integral": source_rhs,
    "source_tangent_integral": tangent_integral,
    "old_IFunction_control_row_at_u0": source_rhs,
    "new_direct_IFunction_control_row_at_u0": -source_rhs,
    "expected_control_row_at_u0": 1.0,
    "source_sign_reversal_if_both_rhs_are_negated":
        bool(source_rhs > 0 and nonlinear_rhs < 0),
}, indent=2), flush=True)
