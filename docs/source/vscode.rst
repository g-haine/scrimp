First SCRIMP project with VS Code (v1.3.0 preview)
======================================================

.. important::

   This is the **target v1.3.0 experience**, being developed on
   the release/1.3.0-preparation branch. The public GHCR v1.3.0 image
   and the Dev Containers end-to-end qualification are **not yet released**.

What is on the user's computer?
-------------------------------

After installation the user has **three separate components**:

1. Local VS Code (Python and, for Docker, Dev Containers extensions).
2. The SCRIMP numerical environment (a Docker image containing SCRIMP,
   PETSc, SLEPc, GetFEM, Gmsh and Python; or a native Conda environment).
3. An **independent local project** holding their simulations, geometry
   inputs and persistent output files.

The project is independent from SCRIMP's source repository. It works
with both a precompiled and locally built Docker image; the same
directory also works with a native Conda installation.

The supplied project template creates this structure:

.. code-block:: text

   MySCRIMPProject/
     .devcontainer/devcontainer.json
     .vscode/launch.json
     .vscode/tasks.json
     .vscode/settings.json
     simulations/first_simulation.py
     simulations/mpi_check.py
     geometries/
     outputs/

The first example is a genuine serial port-Hamiltonian
mass-spring-damper simulation and saves outputs/oscillator.csv.
The second script checks MPI startup, **not physical MPI parity**.

Creating the independent workspace
----------------------------------

When SCRIMP v1.3.0 is released, download the standalone
scrimp-project-template.zip archive from GitHub Releases and extract
it as your own project. That release asset **is not published yet**.
At release packaging time the maintainer must embed the **validated**
GHCR tag in the ZIP's Dev Container metadata, for example:

.. code-block:: bash

   python3 scripts/package-scrimp-project-template.py \
       --image ghcr.io/g-haine/scrimp:1.3.0

This command only creates a local ZIP: it does not publish that
image or asset. Only run it with the GHCR tag after qualification.

During development, obtain the template from GitHub once:

.. code-block:: bash

   BOOTSTRAP="$(mktemp -d -t scrimp-starter.XXXXXXXX)"
   git clone --depth 1 --branch release/1.3.0-preparation \
       https://github.com/g-haine/scrimp.git "$BOOTSTRAP/source"
   python3 "$BOOTSTRAP/source/scripts/create-scrimp-project.py" \
       "$HOME/MySCRIMPProject"

The BOOTSTRAP directory is a separate **newly created** temporary
source checkout. The user's new project is in MySCRIMPProject and
can be retained indefinitely. The temporary source checkout can
be removed after successful project creation; never remove an
existing user directory to make room for it.

The generator refuses to overwrite any pre-existing directory and
never copies the source repository's sandbox, output data or Git history.
For future release packaging, the ZIP is built by:

.. code-block:: bash

   python3 scripts/package-scrimp-project-template.py

A. Precompiled Docker installation (planned release path)
---------------------------------------------------------

After the image has **passed validation and actually been published**,
the planned pull command is:

.. code-block:: bash

   docker pull ghcr.io/g-haine/scrimp:1.3.0

No compilation or SCRIMP source checkout is required. The image is
stored in Docker's image store, *not* in the user's simulations folder.
The published project template should select that image by default.
For the pre-release template, change the image field in
.devcontainer/devcontainer.json to this tag **only after it exists**.

B. Locally compiled Docker installation
---------------------------------------

A developer can rebuild the same candidate from a specific GitHub
revision. The official build driver fetches that revision into a
temporary checkout and stages only tracked, allowlisted source files
into an isolated Docker context:

.. code-block:: bash

   # In a temporary checkout of SCRIMP:
   bash scripts/build-v130-docker-candidate.sh
   bash scripts/smoke-v130-docker-candidate.sh

**Dev Containers image selection (Linux validation finding):**
Use an **explicit complete image reference** in
``.devcontainer/devcontainer.json``, including its tag:

.. code-block:: json

   "image": "scrimp-mpi:v1.3.0-rc1"

In the first VS Code GUI attempt, the template's expression
`${localEnv:SCRIMP_DOCKER_IMAGE:scrimp-mpi:v1.3.0-rc1}`
was incorrectly resolved to just `scrimp-mpi`. VS Code
then attempted to fetch `scrimp-mpi:latest` from Docker Hub and
failed. This was **not** a Docker permissions or SCRIMP runtime
failure. The starter now uses the literal candidate tag. For a
different built image, edit the ``image`` property or generate
the project with ``--image IMAGE:TAG``. Projects created before
the fix must update their **own copy** of devcontainer.json;
pulling the SCRIMP Git repository does not modify existing
independent projects.

The development tag is scrimp-mpi:v1.3.0-rc1. It **passed its
first native Docker build and 1/2-rank CLI smoke on 19 November 2026**,
including a newly scaffolded project, an actual serial oscillator
simulation and a CSV saved on the host. It is not yet release-qualified
or published. This tag is the *pre-release default* in
the project's .devcontainer configuration. The Docker build **does
not access any local simulations**, even from a developer checkout.
Explicit library development with uncommitted changes uses
SCRIMP_DOCKER_SOURCE=local (still using physical source staging).
A user only modifying simulations never rebuilds the image.

VS Code workflow for either Docker installation
------------------------------------------------

1. Install Docker and VS Code with Python and Dev Containers extensions.
   On Windows use Docker Desktop with WSL2; on macOS use Docker Desktop.
2. Open **MySCRIMPProject** in VS Code (not SCRIMP's source checkout).
3. Choose **Dev Containers: Reopen in Container**.
4. Check the selected interpreter is
   /opt/conda/envs/scrimp/bin/python.
5. Run the following in the VS Code integrated terminal:

.. code-block:: bash

   python simulations/first_simulation.py
   mpiexec -n 2 python simulations/mpi_check.py

The host project is bind-mounted at /work. The first command
writes MySCRIMPProject/outputs/oscillator.csv **on the host**.
Closing VS Code or removing the container does not delete that file.

The included VS Code Run and Debug configurations support **F5**
for a serial simulation or the active Python file; MPI multiprocess
debugging requires a separate setup. VS Code Tasks also provide
the first simulation and a two-rank MPI smoke.

SCRIMP initializes output directories relative to the working
directory (outputs/log, outputs/mesh, outputs/png, outputs/pv).
Linux Dev Containers attempts to map the mambauser UID/GID to that
of the host user. Test host output-directory permissions before
large runs; Windows and macOS use different Docker Desktop volume
semantics.

CPU/NUMA selection should be recorded for performance benchmarks.
Do not use MPI --oversubscribe as evidence of strong scaling: it is
appropriate only for small functional tests when the host has too
few CPUs.

Using an existing local sandbox directory
-----------------------------------------

Developers may already have many personal simulations under their
SCRIMP source checkout's sandbox/. **They do not need to move or copy
those files.** The Linux command-line runner binds the existing
sandbox directory at /work at runtime:

.. code-block:: bash

   cd ~/Simulations/scrimp
   SCRIMP_SIM_DIR="$PWD/sandbox" \
       bash scripts/run-v130-docker-simulation.sh wave_1D.py

This command is **not an image build**. The Python script is read from
the local sandbox; its results stay under sandbox/outputs/. Change the
script name to run another local simulation, or set SCRIMP_MPI_RANKS=2
for a script that is actually MPI-compatible. The runner is currently
a Linux/macOS shell development convenience, **not** the final
cross-platform user installer. It still requires a successfully
built/pulled SCRIMP Docker image.

C. Native Conda / source installation (without Docker)
------------------------------------------------------

Use the native build instructions in :doc:`install` for the intended
platform. Then activate the installed Conda environment:

.. code-block:: bash

   conda activate scrimp55  # for the native modern-stack example in install.rst
   cd ~/MySCRIMPProject
   python simulations/first_simulation.py
   mpiexec -n 2 python simulations/mpi_check.py

Open this **same project** in VS Code, select **Python: Select
Interpreter** and choose the Conda interpreter. Do not select
Reopen in Container. Native installation is distinct from compiling
a Docker image locally, which still targets a Linux container.

Known limitations / release gates
---------------------------------

- The candidate image has **passed its native Linux build and
  1/2-rank CLI smoke**. Actual VS Code Dev Containers GUI workflow,
  native compatibility and Windows/macOS/Apple Silicon remain
  unqualified; do not extrapolate the Linux smoke to these platforms.
- The public GHCR v1.3.0 tag and starter project ZIP are not
  published; both require an accepted release.
- The introductory oscillator is a serial model; use validated
  PDE scripts for MPI and Hamiltonian/power-balance qualification.
- The old >45 GB Docker-context regression is addressed by staging
  only approved Git-tracked source files in a separate temporary
  directory. User-owned mounted workspaces are never part of the
  image build.

See the installation guide for the numerical stack, the release
plan for remaining scientific acceptance tests and the developer
Docker build guide for cache and provenance details.
