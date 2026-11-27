How to install
==============

**Docker is the recommended/default installation for SCRIMP v1.3.0.**
The published Linux/AMD64 image provides a precompiled environment
containing Python, GetFEM 5.5, PETSc/SLEPc (complex), MUMPS, OpenMPI and
SCRIMP itself. User simulations and their outputs live in an
**independent, host-owned project**, not inside the image or a SCRIMP
source checkout. A source checkout is **not** needed by a released-image
user.

.. important::

   **SCRIMP v1.3.0 (27 November 2026):** the stable Git tag and
   **public** ``ghcr.io/g-haine/scrimp:1.3.0`` Docker image are available.
   The stable image retains the accepted ``dist-rc4`` OCI index digest
   ``sha256:737f58f12475723c0924f25528a0794a767fe6051dbecbd80a194d3b15eb3f01``
   and was built from Git revision
   ``17ae5a1ce063ff5481812dc73f3d3fc88c153c70``.
   A fresh Ubuntu 24.04/AMD64 runner validated the independent user
   project/host CSV, MPI startup and native Heat--Wave Hamiltonian
   parity (`CI evidence
   <https://github.com/g-haine/scrimp/actions/runs/46030753404>`_).
   The standalone user ZIP is distributed separately through the
   `v1.3.0 GitHub Release
   <https://github.com/g-haine/scrimp/releases/tag/v1.3.0>`_.
   Windows/macOS/ARM64 remain best-effort.

Default: precompiled Docker and standalone project
-------------------------------------------------

1. Install Docker (on Windows, Docker Desktop with WSL2; on macOS,
   Docker Desktop). VS Code and Dev Containers are recommended for
   editing/debugging but **not required** for command-line execution.
2. Download the matching
   ``scrimp-project-template.zip`` asset from the SCRIMP v1.3.0
   GitHub Release and unpack it into your own project directory. The
   release ZIP embeds the **validated** image reference in
   ``.devcontainer/devcontainer.json``.
3. Pull the released Docker image:

   .. code-block:: bash

      docker pull ghcr.io/g-haine/scrimp:1.3.0

4. Open the project directory in VS Code and choose
   **Dev Containers: Reopen in Container**. In its terminal, run:

   .. code-block:: bash

      python simulations/first_simulation.py
      mpiexec -n 2 python simulations/mpi_check.py

   The oscillator is a real serial SCRIMP simulation and saves
   ``outputs/oscillator.csv`` **on the host**. The second command
   checks only MPI startup, not a distributed PDE power/physics result.

5. Without VS Code, execute the same project by mounting your
   independent workspace (Linux shell):

   .. code-block:: bash

      cd ~/MySCRIMPProject
      docker run --rm --user "$(id -u):$(id -g)" \
          -e HOME=/tmp -e XDG_CACHE_HOME=/tmp/.cache \
          -v "$PWD:/work:rw" -w /work \
          ghcr.io/g-haine/scrimp:1.3.0 \
          python simulations/first_simulation.py

   On Docker Desktop, permissions and mount semantics differ; consult
   the Windows/macOS guidance below. No numerical output belongs
   in a throwaway container filesystem.

The complete standalone workspace, release ZIP generation, explicit
image tags and F5 debugging are documented in :doc:`vscode`.
The distribution qualification contract (including exact source
revision, image metadata, host CSV and MPI tests) is documented
in ``docs/V130_DOCKER_FIRST_RELEASE_GATE.md``. All required
Linux/AMD64 distribution checks have passed; that does not
qualify other platforms.

Alternative: native Conda installation
--------------------------------------

The original Conda workflow is retained for current users and
advanced developers. It is **not** the v1.3.0 default.

1. Install `Conda <https://docs.conda.io/>`_ or
   `Micromamba <https://mamba.readthedocs.io/>`_.
2. Clone the repository: ``git clone https://github.com/g-haine/scrimp``.
3. Enter the folder: ``cd scrimp``.
4. Run ``conda env create -f scrimp.yml``.
5. Activate: ``conda activate scrimp``.
6. Install SCRIMP: ``python -m pip install -e .``.

The historical ``scrimp.yml`` deliberately retains its older
packaged GetFEM/NumPy dependencies for existing users; **do not**
interpret native Conda as a validated GetFEM 5.5/MPI/NumPy 2
replacement for the candidate Docker stack.

Modern parallel Docker image (GetFEM 5.5)
-----------------------------------------

The v1.3.0 **Docker build** recipe is in ``docker/Dockerfile``,
with dependencies in ``docker/environment.yml``. It uses one
consistent conda-forge OpenMPI stack and complex PETSc/SLEPc/MUMPS.
The PETSc/SLEPc Python bindings are compiled against the installed
native libraries and NumPy 2 headers. GetFEM 5.5 is compiled with
``--enable-paralevel=2``, including SCRIMP's approved simplex
BDM/Nédélec source patches.

Native dependencies are now compiled in **cache-friendly layers
before copying SCRIMP's Python sources**. Docker BuildKit retains
APT, Micromamba, pip and verified GetFEM tarball downloads.

The candidate build driver fetches the selected GitHub revision
by default into a temporary checkout and stages only audited tracked
sources in a *second* isolated context. It **never hands Docker the
live SCRIMP repository**, even if the user has tens of gigabytes
of simulations under sandbox/outputs. For intentionally local
library development, set `SCRIMP_DOCKER_SOURCE=local`; this still
copies only tracked allowlisted source files into the build context.
Neither mode includes user simulation data.

Use an **independently tagged candidate**, not the qualified existing
``scrimp-mpi:local`` image:

.. code-block:: bash

    cd ~/Simulations/scrimp
    git pull --ff-only
    bash scripts/build-v130-docker-candidate.sh
    bash scripts/smoke-v130-docker-candidate.sh

By default these target ``scrimp-mpi:v1.3.0-rc1``. The build script
writes ``sandbox/outputs/v130-docker-build.log`` and records the
built image ID; the smoke writes
``sandbox/outputs/v130-docker-smoke.log``, exports the exact
resolved Conda package spec and checks the GetFEM simplex descriptor
support, NumPy ABI, complex PETSc and distributed MUMPS at 1/2 ranks.
It explicitly rejects Nédélec-II order k>=3, which has **not**
been implemented.

The ``rc1`` script default is a **local development candidate tag**,
not the public image. To use the stable qualified environment,
pull ``ghcr.io/g-haine/scrimp:1.3.0`` instead. Do not overwrite a
validated local image while experimenting with source builds.

For interactive development with an **independent local
simulation project** (not the source checkout), see :doc:`vscode`.
For a manual command-line session from the user's simulation folder:

.. code-block:: bash

    cd ~/MySCRIMPProject
    sudo docker run --rm -it -v "$PWD:/work" \
        --workdir /work ghcr.io/g-haine/scrimp:1.3.0 bash

This binds the local workspace read-write, so simulations and outputs
persist on the host even when the container is removed. On Linux,
check UID/GID mapping before relying on direct docker-run commands;
the VS Code Dev Containers template enables UID matching.

Further instructions and the cache-validation protocol are in
``docs/V130_DOCKER_BUILD.md``. The stable image is publicly available as ``ghcr.io/g-haine/scrimp:1.3.0``.
Existing Sphinx/Conda instructions for the stable ``scrimp.yml``
environment remain unchanged.

The GitHub Actions workflow ``.github/workflows/docker-mpi.yml``
is still manual-only and separately checks MPI/native assembly,
algebra and representative wave physics. GHCR publication remains
opt-in and must follow actual image validation; the local candidate
build script never publishes.

Windows and macOS users: best-effort Docker/Dev Containers
------------------------------------------------------------

SCRIMP has users on Windows and macOS. These are **best-effort
installation workflows**, not independently qualified native
Windows/macOS GetFEM/PETSc/MUMPS configurations. The released
Linux/AMD64 image is ``ghcr.io/g-haine/scrimp:1.3.0``.
Check its platform before use; historical local tags like
``scrimp-mpi:v1.3.0-rc3`` are not published GHCR tags.

**Windows 10/11 (WSL2):** install Docker Desktop, enable its WSL2
backend and integration for your Linux distribution, and install
VS Code plus the Dev Containers extension. Open a **WSL2 Linux**
terminal and verify ``docker version``. Keep your SCRIMP user
project under the WSL2 Linux filesystem (rather than an NTFS
``/mnt/c`` mount) if you will create many numerical outputs.
From this terminal, after pulling a compatible image, run:

.. code-block:: bash

    cd ~/MySCRIMPProject
    docker images
    docker run --rm -v "$PWD:/work" -w /work \
        ghcr.io/g-haine/scrimp:1.3.0 \
        python simulations/first_simulation.py

These commands run a Linux container through WSL2; they do **not**
constitute a native Windows compilation. If Docker Desktop's
integration/permissions differ, use its documented WSL2 integration
settings. Use :doc:`vscode` for the independent project workflow.

**macOS (Intel/Apple Silicon):** install Docker Desktop for the
machine architecture and VS Code/Dev Containers. Once an image
exists locally, inspect its architecture using
``docker image inspect IMAGE:TAG``; a Linux AMD64-only image
will not become ARM64-native merely because Docker Desktop
can start it with emulation. On an Apple Silicon machine,
``--platform linux/amd64`` is a possible **experimental**
compatibility option for an AMD64 image, but its numerical
performance and GetFEM/PETSc/MUMPS MPI behavior remain unqualified:

.. code-block:: bash

    cd ~/MySCRIMPProject
    docker run --rm --platform linux/amd64 \
        -v "$PWD:/work" -w /work \
        ghcr.io/g-haine/scrimp:1.3.0 \
        python simulations/first_simulation.py

Do not infer scientific accuracy or MPI safety from a successful
container launch. Prefer a separately built and validated ARM64
image when one becomes available. Docker Desktop's bind-mount
permissions and filesystem performance differ from native Linux.

In **all** cases, begin with the bundled oscillator and a bounded
two-rank MPI start-up check, then progress to validated PDE models.
If your setup differs from the tested Linux baseline, record image
digest/platform, Python/PETSc/GetFEM versions and raw logs when
reporting a problem. Windows/macOS use is welcome, but it is
**not a mandatory v1.3.0 release-blocking acceptance gate**.

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

Native GetFEM MPI assembly timing verbosity
-------------------------------------------

GetFEM 5.5 prints ``Assembly time ...`` for **every model assembly**
on the master MPI rank, even when SCRIMP calls
``gf.util_trace_level(0)``. This comes from an unconditional
C++ ``std::cout`` in ``src/getfem_models.cc``, rather than
SCRIMP's Python logger.

For **new native GetFEM builds made with SCRIMP's build script**,
``scripts/patches/getfem-5.5/0000-mpi-assembly-timing.patch`` changes
only this diagnostic print, after extracting the upstream source and
before compilation. It is applied by the same ordered patch loop as
the BDM fix. Normal execution is quiet. To restore upstream-style
per-assembly timing information for diagnosis, set:

.. code-block:: bash

    SCRIMP_GETFEM_ASSEMBLY_TIMING=1 mpiexec -n 4 \
        python sandbox/smoke_modern_stack.py

The flag must be present in the environment **inside** the
container and on all MPI ranks. The source patch fails to apply if
GetFEM's expected code has changed: it must be reviewed rather than
silently forced. Neither mathematical expressions, assembly results,
nor MPI communication are changed. Avoid silencing all stdout/stderr
in Python: GetFEM's native errors and diagnostics must remain visible.

**An image built before this change still prints the messages.**
Rebuilding the native GetFEM layer is necessary to change that image;
defer it to the next planned image rebuild rather than recompiling
GetFEM solely to remove cosmetic output.

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

Optional SCRIMP phase profiling (without rebuilding GetFEM)
-----------------------------------------------------------

For v1.3.0 development, SCRIMP can measure its own GetFEM assembly
steps, GetFEM-to-PETSc transfers, DAE consistency initialization,
PETSc TS calls and implicit residual/Jacobian callbacks **per MPI rank**.
Profiling is disabled by default. It adds no MPI collectives inside
PETSc callbacks; reports are gathered explicitly after the solve.

Update the Python sources from the branch under development, then
reuse the **existing local Docker image** with a bind mount; no new
GetFEM compilation is necessary:

.. code-block:: bash

    git pull --ff-only
    mkdir -p sandbox/outputs/scaling-profile
    sudo docker run --rm \
        --cpuset-cpus="0,2,4,6" \
        -e SCRIMP_MUMPS_RESET_PC_ON_JACOBIAN=1 \
        -v "$PWD:/workspace/scrimp" \
        scrimp-mpi:local \
        mpiexec --map-by hwthread --bind-to hwthread -n 4 \
        python -u sandbox/bench_acoustic_2d.py \
        --mesh 0.08 --dt 0.002 --final-time 0.04 \
        --save-dt 0.01 --profile \
        --output /workspace/scrimp/sandbox/outputs/scaling-profile/wave_4.json

The container's non-root user needs write access to the chosen output
directory. Alternatively, write JSON under ``/tmp`` inside the
container and pipe/copy it to a local file. Without ``--profile``,
the original benchmark behaves as before.

In other DPHS scripts, set ``SCRIMP_PROFILE=1`` **before** constructing
the DPHS and call ``system.get_profile()`` on **every** MPI rank
after ``system.solve()``. The return value is JSON-serializable.

Reported wall durations use per-event MPI rank min/mean/max, and
record the number of calls. **Nested timers overlap**:
``main.petsc.ts_solve`` includes time spent in Python IFunction
and IJacobian callbacks, whose timings in turn include GetFEM and
PETSc bridge work. The derived ``main.petsc.non_callback_remainder``
subtracts callback durations on each rank; it includes factorization,
solve, monitoring and MPI waiting, and is **not pure MUMPS time**.
For the acoustic benchmark, add ``--petsc-log`` on a **separate
profiling run** to activate PETSc's native event logger programmatically
(``PETSc.Log.begin()`` before setup and ``PETSc.Log.view()`` after
the simulation, collectively on all ranks). This distinction matters:
DPHS currently clears PETSc's options database at construction,
so relying on the command-line ``-log_view`` alone can be misleading.
PETSc's printed summary contains native matrix, factorization and
solve events, including ``MatLUFactorSym``, ``MatLUFactorNum`` and
``MatSolve`` when exposed by the linked PETSc build.

With ``SCRIMP_MUMPS_RESET_PC_ON_JACOBIAN=1``, the profiler also
reports ``main.petsc.pc_reset`` (inclusive in
``main.petsc.callback.ijacobian``). Do not disable this stability
workaround simply to improve benchmark timings without verifying
sparsity lifecycle and physical parity.

These measurements do not alter matrices, RHS scaling, equations
or the Crank--Nicolson scheme.

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
