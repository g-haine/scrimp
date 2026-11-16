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
OpenMPI stack for native PETSc, SLEPc and MUMPS, with a
complex PETSc scalar type. The petsc4py/slepc4py Python bindings
are built against that stack and the NumPy 2 headers from source. GetFEM 5.5 is built from the official
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
    sudo apt-get install -y binutils file libglu1-mesa libgl1

    # Create a new environment from the same solver-compatible recipe.
    conda env create -n scrimp55 -f docker/environment.yml
    conda activate scrimp55

    # Build matching NumPy-2 Python wrappers against the installed
    # complex native PETSc/SLEPc (conda petsc4py may require NumPy<2).
    bash scripts/install-petsc-python-bindings.sh

    # Compile the official GetFEM 5.5 sources in MPI level 2.
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

**Python 3.12 installation compatibility:** The GetFEM 5.5 tarball
still bundles an older Automake ``py-compile`` that imports
``imp``, which was removed from Python 3.12. The shared
``scripts/build-getfem-mpi.sh`` installation therefore overrides
**only** Automake's bytecode compiler
(``make install 'am__py_compile=:'``), and immediately performs
``compileall`` with the active Python interpreter. This retains
native build and installation, rather than downgrading Python
or installing a global ``imp`` shim.

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

Native GetFEM MPI assembly audit (important)
---------------------------------------------

Before interpreting any parallel simulation, compare the **raw
GetFEM outputs** at P=1, 2 and 4. The official GetFEM MPI
documentation says:

* Model RHS vectors are **globally assembled and replicated** on
  every rank. A complete vector is expected and by itself is **not**
  a reason to divide its values by P.
* Tangent matrices are normally **distributed contributions**;
  the GetFEM-to-PETSc bridge must ADD these across ranks.
* The original maintainer's GetFEM 5.4.2 build instead returned
  full replicated tangent matrices and P-fold nonlinear RHS terms;
  this is a measured local behavior, **not** an assumption for 5.5.

We explicitly audit native matrices, generic nonlinear RHS,
prescribed source RHS and their sum, separately, without invoking
SCRIMP's normalization, PETSc transfer or MUMPS:

.. code-block:: bash

    mkdir -p sandbox/outputs/getfem-assembly
    for n in 1 2 4; do
        mpiexec -n "$n" python -u \
            sandbox/diagnose_getfem_mpi_assembly.py \
            --output "sandbox/outputs/getfem-assembly/${n}.json"
    done

    python sandbox/analyze_getfem_mpi_assembly.py \
        sandbox/outputs/getfem-assembly/*.json

The analytical references on [0,1], with u=1 and source=2,
are -1 for the generic RHS, +2 for the source-only RHS,
+1 for the mixed RHS, and 1/2 for the linear/nonlinear
tangent actions on a vector of ones.

The analyzer reports whether the generic RHS or source RHS
are **physical** or **P-fold duplicated**, and whether tangent
matrices are **full and replicated** or **distributed by
contribution**. The GetFEM 5.5 Docker build was measured with
correct raw generic and source RHS and distributed matrix
contributions on P=2/4. The Docker image therefore sets:

.. code-block:: bash

    SCRIMP_GETFEM_MATRIX_MODE=distributed
    SCRIMP_GETFEM_RHS_MODE=auto
    SCRIMP_GETFEM_SOURCE_RHS_MODE=raw

The default for non-Docker environments remains
``SCRIMP_GETFEM_MATRIX_MODE=replicated``, preserving the behavior
of the previously validated *local* GetFEM 5.4.2 installation.
For a freshly compiled GetFEM **never** copy these flags blindly:
first measure the native assembly with the diagnostic above, then
run its analyzer with the matching ``--matrix-mode replicated``
or ``--matrix-mode distributed`` parameter. For this GetFEM 5.5
container the mandatory CI command is:

.. code-block:: bash

    python sandbox/analyze_getfem_mpi_assembly.py \
        --matrix-mode distributed --require-physical-rhs \
        sandbox/outputs/getfem-assembly/*.json

The RHS bridge calibrates generic contributions in
``SCRIMP_GETFEM_RHS_MODE=auto``, while
``SCRIMP_GETFEM_SOURCE_RHS_MODE=raw``
preserves unscaled GetFEM source terms that are assembled once.

**Do not change RHS normalization solely because GetFEM returns
full vectors or matrices.** First compare raw numerical values with
the known integrals and the P=1 reference. Preserve logs from
both GetFEM 5.4.2 and GetFEM 5.5 to distinguish an installation
issue from an upstream GetFEM convention.

GetFEM parallelism reference:
`GetFEM parallel documentation <https://getfem.org/userdoc/parallel.html>`_.

2D wave Hamiltonian and boundary-power balance (MPI)
----------------------------------------------------

Before timing the parallel solver, validate the **undamped** wave
example's port-Hamiltonian power identity using the same GetFEM
Hamiltonian integrands as ``examples/wave.py``. The prescribed
left-boundary input in the published example is nearly zero during
short smoke runs; a separate **nonzero boundary drive** prevents
vacuous energy-balance checks without changing the wave equations.

For the quadratic Hamiltonian and Crank--Nicolson, check the discrete
midpoint identity (all four control ports have
``dissipative=True``, meaning a positive ``flow*effort`` is
subtracted from energy in SCRIMP's existing plot convention):

.. math::

   H(z_{k})-H(z_{k-1}) + (t_k-t_{k-1})
   \sum_{\Gamma} \int_{\Gamma}
   f_{\Gamma}(z_{k-1/2}) e_{\Gamma}(z_{k-1/2}) \,d\Gamma = 0.

The residual is calculated **per time step** at midpoint *state and
time*, not inferred from Hamiltonian conservation with zero input.
It also checks each scalar GetFEM integral is the same on every
rank before accepting it. The full Hamiltonian, each boundary
power trace, the signed work and the numerical defect are exported.

Inside the validated GetFEM 5.5 container (from repository root):

.. code-block:: bash

    mkdir -p sandbox/outputs/energy
    for n in 1 2 4; do
        mpiexec --oversubscribe -n "$n" \
            python sandbox/check_acoustic_energy_mpi.py \
            --mesh 0.2 --dt 0.002 --final-time 0.02 \
            --drive-left '0.08*sin(20*t)' \
            --output "sandbox/outputs/energy/wave_${n}.json"
    done

    python sandbox/analyze_acoustic_power_balance.py \
        sandbox/outputs/energy/wave_*.json

Use ``--diagnose-only`` to archive **raw measurements**, without
asserting an energy tolerance, when investigating a previously
unvalidated GetFEM build; do not describe such a run as validated.
CI uses a strict numeric gate and compares the complete Hamiltonian
and boundary-power trajectories at MPI 1, 2 and 4. Hosted CI uses
oversubscribed cores and its timing data are **not** suitable for
strong-scaling claims.

This check covers the *linear undamped acoustic problem*. It does
not yet prove energy consistency for nonlinear Saint-Venant models
or for substituted/dissipative ports.

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

To summarize multiple comparable runs automatically:

.. code-block:: bash

    python sandbox/analyze_acoustic_scaling.py \
        sandbox/outputs/scaling/wave_*ranks.json

For more robust measurements, repeat each configuration several
times and provide all of the resulting JSON files to the analyzer.
It refuses mismatched meshes, end times, integrator settings or
significantly different physical field norms. The benchmark itself
exposes per-rank peak memory usage; for detailed PETSc timing,
append ``-log_view`` to the command line and capture its log.

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
