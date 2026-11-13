# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Read-only MPI/ABI inventory for locally compiled GetFEM + Conda PETSc.

Run from the SCRIMP repository, inside the active "scrimp" environment:

    python -u sandbox/audit_mpi_environment.py
    mpiexec -n 2 python -u sandbox/audit_mpi_environment.py
    mpiexec -n 4 python -u sandbox/audit_mpi_environment.py

No packages or configuration files are modified. Rank 0 additionally
calls ldd on the imported ELF extensions and records matching libraries.
This diagnoses mixed system/Conda MPI or BLAS/Fortran runtime linkage.
"""

import importlib
import json
import os
import re
import shutil
import subprocess
import sys

MATCH = re.compile(
    r"(?:mpi|mpich|open-rte|open-pal|pmix|petsc|slepc|mumps|"
    r"blas|lapack|gfortran|stdc\+\+|scalapack|metis|parmetis|boost)",
    re.IGNORECASE,
)


def command(argv):
    try:
        completed = subprocess.run(
            argv, check=False, stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT, text=True, timeout=8,
        )
        return {"code": completed.returncode, "stdout": completed.stdout[:15000]}
    except (OSError, subprocess.TimeoutExpired) as error:
        return {"error": str(error)}


def module_location(mod):
    return os.path.realpath(getattr(mod, "__file__", ""))


def run():
    # Import the native interfaces in the same process to see actual
    # dependencies and ABI conflicts rather than only Conda metadata.
    from mpi4py import MPI
    from petsc4py import PETSc
    from slepc4py import SLEPc
    import getfem as gf

    rank = MPI.COMM_WORLD.Get_rank()
    size = MPI.COMM_WORLD.Get_size()

    names = {
        "mpi4py.MPI": MPI,
        "petsc4py.PETSc": PETSc,
        "slepc4py.SLEPc": SLEPc,
        "getfem._getfem": importlib.import_module("getfem._getfem"),
    }
    details = {
        "rank": rank,
        "mpi4py_size": size,
        "petsc_rank": PETSc.COMM_WORLD.getRank(),
        "petsc_size": PETSc.COMM_WORLD.getSize(),
        "getfem_parallelism_level": int(gf.util_mpi_parallelism_level()),
        "python": sys.version.splitlines()[0],
        "python_executable": sys.executable,
        "conda_prefix": os.environ.get("CONDA_PREFIX"),
        "mpi_library_version": MPI.Get_library_version().strip(),
        "modules": {key: module_location(val) for key, val in names.items()},
    }
    try:
        details["petsc_version"] = PETSc.Sys.getVersion()
    except Exception as error:
        details["petsc_version_error"] = str(error)
    try:
        details["slepc_version"] = SLEPc.Sys.getVersion()
    except Exception as error:
        details["slepc_version_error"] = str(error)

    if (PETSc.COMM_WORLD.getSize() != size
            or PETSc.COMM_WORLD.getRank() != rank):
        raise RuntimeError("PETSc and mpi4py disagree on MPI world rank/size")

    if rank == 0:
        details["environment"] = {
            key: os.environ.get(key)
            for key in ("PATH", "LD_LIBRARY_PATH", "CONDA_PREFIX",
                        "PETSC_DIR", "PETSC_ARCH", "SLEPC_DIR",
                        "OMPI_MCA_btl", "OMP_NUM_THREADS")
        }
        details["executables"] = {
            name: shutil.which(name)
            for name in ("mpiexec", "mpicc", "mpicxx", "mpif90", "ldd")
        }
        details["mpiexec_version"] = command(["mpiexec", "--version"])
        details["mpicc_show"] = command(["mpicc", "-show"])
        details["loaded_elf_dependencies"] = {}
        for name, mod in names.items():
            location = module_location(mod)
            if not location:
                continue
            report = command(["ldd", location])
            if "stdout" in report:
                report["relevant_libraries"] = [
                    row.strip()
                    for row in report["stdout"].splitlines()
                    if MATCH.search(row)
                ]
                del report["stdout"]
            details["loaded_elf_dependencies"][name] = report

    print(json.dumps(details, indent=2, default=str), flush=True)


if __name__ == "__main__":
    run()
