# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Hard dependency/ABI gate for the experimental GetFEM 5.5+NumPy 2 image.

Run in any active target environment, also under mpiexec -n 2/4.
Checks Python ABI imports, GetFEM MPI-2, complex PETSc, and native
MUMPS factorization on distributed PETSc matrices. No claims about
physical SWE parity or scalability follow from this smoke alone.
"""
import json
import numpy as np
import scipy
import getfem as gf
import gmsh
import mpi4py
from petsc4py import PETSc
from slepc4py import SLEPc


comm = PETSc.COMM_WORLD
p = comm.getSize()
rank = comm.getRank()
assert int(np.__version__.split(".")[0]) >= 2, np.__version__
assert int(gf.util_mpi_parallelism_level()) == 2
assert np.issubdtype(PETSc.ScalarType, np.complexfloating)

# Do more than just import NumPy/GetFEM: cross their binary C-API
# boundary with a concrete model, assembly and owned PETSc insertion.
mesh = gf.Mesh("cartesian", np.linspace(0.0, 1.0, 5))
mf = gf.MeshFem(mesh, 1)
mf.set_fem(gf.Fem("FEM_PK(1,1)"))
mim = gf.MeshIm(mesh, gf.Integ("IM_GAUSS1D(5)"))
model = gf.Model("real")
model.add_fem_variable("u", mf)
model.add_nonlinear_term(mim, "u*u*Test_u")
model.set_variable("u", np.ones(model.nbdof()))
model.assembly(option="build_rhs")
gf_rhs = np.asarray(model.rhs())
assert gf_rhs.size == model.nbdof()
assert np.isfinite(gf_rhs).all()

# Deliberately avoid mixing an MPI implementation in a system package
# with another inside conda. This LU uses PETSc's actual factor backend.
n = max(8, 2 * p)
a = PETSc.Mat().createAIJ(size=(n, n), nnz=1, comm=comm)
a.setUp()
diag = a.createVecRight()
diag.set(2.0)
a.setDiagonal(diag)
a.assemble()
b = a.createVecRight()
b.set(1.0)
x = b.duplicate()
ksp = PETSc.KSP().create(comm=comm)
ksp.setType("preonly")
ksp.getPC().setType("lu")
ksp.getPC().setFactorSolverType("mumps")
ksp.setOperators(a)
ksp.solve(b, x)
assert ksp.getConvergedReason() > 0, ksp.getConvergedReason()
assert abs(x.norm() - (n**0.5)/2.0) < 1e-8
if rank == 0:
    print(json.dumps({
        "mpi_size": p,
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "getfem_mpi_level": gf.util_mpi_parallelism_level(),
        "petsc": PETSc.Sys.getVersion(),
        "petsc_scalar_complex": True,
        "mumps": "distributed LU passed",
        "gmsh": "imported",
        "slepc": "imported",
        "mpi4py": mpi4py.__version__,
    }, indent=2), flush=True)
ksp.destroy()
x.destroy()
b.destroy()
diag.destroy()
a.destroy()
