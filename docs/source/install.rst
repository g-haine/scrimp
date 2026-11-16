How to install
==============

This page describes the conventional environment as well as the
**experimental GetFEM 5.5 / NumPy 2 / MPI level 2** installation.
The modern stack is NOT yet declared validated: the Docker build,
NumPy ABI import, MPI algebra tests and physics/scalability benchmarks
must pass on real machines before it becomes the release default.

Stable user environment (Conda)
-------------------------------

1. Install `Conda <https://docs.conda.io/>`_ or
   `Micromamba <https://mamba.readthedocs.io/>`_.
2. Clone the repository: ``git clone https://github.com/g-haine/scrimp``.
3. Enter the folder: ``cd scrimp``.
4. Run ``conda env create -f scrimp.yml``.
5. Activate: ``conda activate scrimp``.
6. Install SCRIMP: ``python -m pip install -e .``.

The existing ``scrimp.yml`` deliberately retains ``numpy<2`` and
the packaged GetFEM build for established users. It is **not**
silently upgraded by the experimental instructions below.

Modern parallel Docker image (GetFEM 5.5)
-----------------------------------------

A versioned image recipe lives in ``docker/Dockerfile``, with
``docker/environment.yml``. It uses one *consistent* conda-forge
OpenMPI stack for GetFEM, PETSc, petsc4py, SLEPc and MUMPS, with a
complex PETSc scalar type. GetFEM 5.5 is built from the official
release source with ``--enable-paralevel=2``. NumPy 2 is **required**,
not presumed to work: builds fail if its ABI smoke test fails.

.. code-block:: bash

    docker build --progress=plain -f docker/Dockerfile \
        -t scrimp-mpi:gf55 .

    docker run --rm scrimp-mpi:gf55 \
        mpiexec -n 2 python sandbox/smoke_modern_stack.py

    docker run --rm scrimp-mpi:gf55 \
        mpiexec -n 2 python -m unittest -v tests.test_mpi_algebra

To mount the repository for development, from its root:

.. code-block:: bash

    docker run --rm -it -v "$PWD:/work" \
        --workdir /work scrimp-mpi:gf55 bash

These examples launch multiple MPI processes **inside one container
on one host**. Running between HPC nodes requires a compatible
MPI launcher/transport (for example Apptainer with a site-supported
MPI ABI) and dedicated cluster testing. Avoid mixing Ubuntu MPI,
conda MPI and host MPI shared libraries.

The GitHub Actions workflow ``.github/workflows/docker-mpi.yml``
builds and tests the image. **Publication on GHCR is opt-in**:
manually run `Build and validate modern GetFEM MPI image` using
`workflow_dispatch` with `publish=true`. No image is published
unless building and smoke tests have succeeded. Record the produced
``docker-resolved-spec.txt`` artifact and image digest before
shipping a release.

Local GetFEM 5.5 MPI compilation (Conda environment)
----------------------------------------------------

This approach enables a source build without Docker and keeps a
separate, isolated environment instead of modifying your existing
SCRIMP installation:

.. code-block:: bash

    # Ubuntu: GL runtime for the Gmsh Python wheel.
    sudo apt-get update
    sudo apt-get install -y libglu1-mesa libgl1

    # Create a new environment from the same solver-compatible recipe.
    conda env create -n scrimp55 -f docker/environment.yml
    conda activate scrimp55

    # Run from the repository root. With GetFEM 5.5 release tarball:
    GETFEM_VERSION=5.5 MAKE_JOBS=4 \
        bash scripts/build-getfem-mpi.sh

    python -m pip install --no-deps --no-build-isolation -e .
    python sandbox/smoke_modern_stack.py
    mpiexec -n 2 python sandbox/smoke_modern_stack.py
    mpiexec -n 2 python -m unittest -v tests.test_mpi_algebra

The build script uses the active ``CONDA_PREFIX``, the
``mpicc/mpicxx/mpifort`` wrappers and a local build directory
``.build/getfem``. It refuses to silently install a serial
GetFEM: it checks the runtime value
``getfem.util_mpi_parallelism_level() == 2``.

To build the latest **development** GetFEM, use
``GETFEM_SOURCE=git GETFEM_REF=<commit-sha>`` and keep the
specific upstream commit in the experiment manifest.
That is separate from the stable **5.5 release**. The tarball
checksum may be pinned through ``GETFEM_SHA256``; this
must be set when turning an experimental recipe into a
reproducible release.

**NumPy 2 compatibility** is currently a validation goal, not
an asserted upstream guarantee. The C extension must be compiled
against the NumPy 2 headers in the environment and tested at
runtime. If either compilation or tests fail, do not downgrade
NumPy silently inside the modern environment; preserve the
error logs for debugging. Keep ``scrimp.yml`` as a fallback.

2D acoustic-wave strong scaling
-------------------------------

The non-graphical ``sandbox/bench_acoustic_2d.py`` builds the
*existing* ``examples/wave.py`` model without modifying its
physical variational forms. It replaces only the spatial mesh
spacing, final integration horizon, numerical solver and output
frequency. It reports separately model construction/meshing,
complete solve (assembly + DAE initialization + PETSc TS),
peak rank RSS and physical state norms.

For **strong scaling**, hold ``--mesh``, ``--dt``
and ``--final-time`` strictly fixed:

.. code-block:: bash

    mkdir -p sandbox/outputs/scaling

    for n in 1 2 4; do
        SCRIMP_MUMPS_RESET_PC_ON_JACOBIAN=1 \
        SCRIMP_GETFEM_SOURCE_RHS_MODE=raw \
        mpiexec -n "$n" python -u sandbox/bench_acoustic_2d.py \
            --mesh 0.05 --dt 0.001 --final-time 0.02 \
            --save-dt 0.01 \
            --output "sandbox/outputs/scaling/wave_${n}ranks.json" \
            > "sandbox/outputs/scaling/wave_${n}ranks.log" 2>&1
    done

Report strong-scaling speedup
:math:`S_P=T_1/T_P` and efficiency :math:`E_P=S_P/P`
only for the **same mesh, time scheme, time horizon and
solver configuration**. Check final field norms and
all MPI convergence statuses first. Repeat runs and report
medians and variation: a single wall-clock number is not
a meaningful scaling study. Use PETSc ``-log_view`` to
profile internal operations and separate GetFEM assembly,
MUMPS factorization, python-to-PETSc data transfer,
mesh construction and per-step all-rank gathers.

Tests
-----

Run targeted non-graphical tests with:

.. code-block:: bash

    python -m unittest -v tests.test_lumped_ports
    mpiexec -n 1 python -m unittest -v tests.test_mpi_algebra
    mpiexec -n 2 python -m unittest -v tests.test_mpi_algebra
    mpiexec -n 4 python -m unittest -v tests.test_mpi_algebra

The graphical examples below require a Matplotlib GUI backend:

.. code-block:: bash

    python examples/mass_spring_damper.py
    python examples/nonlinear_pendulum.py

Code structure
--------------

**SCRIMP** is a Python package. The ``scrimp/__init__.py`` file
is its root. For instance, the linear-algebra utilities can be
imported as ``import scrimp.utils.linalg``.

Documentation
-------------

Full documentation is published at
`SCRIMP documentation <https://g-haine.github.io/scrimp/>`_.
It is generated with `Sphinx <https://www.sphinx-doc.org/>`_.
