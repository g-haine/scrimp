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

Python is installed, even if VS Code says "No Python found"
-----------------------------------------------------------

.. important::

   **Do not click "Install uv and Python"** when opening a SCRIMP
   Dev Container. The image **already contains Python**, installed
   in the named Micromamba environment ``scrimp`` along with
   SCRIMP, GetFEM, PETSc and SLEPc. Installing another Python with
   ``uv`` would create an unrelated interpreter and might leave
   VS Code running code outside the qualified numerical environment.

During the first Linux Dev Containers GUI test (19 November 2026),
VS Code displayed::

   No Python found. Would you like to install uv and use it to
   install Python?

This message comes from **VS Code's Python environment discovery /
installation UI**, not from SCRIMP or the Python runtime. The
subsequent integrated terminal showed ``(scrimp)``, and both
``python simulations/first_simulation.py`` and
``mpiexec -n 2 python simulations/mpi_check.py`` **passed**.
Therefore that specific prompt was a **false negative for the
already installed container interpreter**.

**When this prompt appears:**

1. Select **Don't ask again** (or **Cancel** to dismiss without
   installing anything). Do not select **Install uv and Python**.
2. Wait until the Dev Container has fully opened and the Python
   extension has completed initialization.
3. Open a Python file, press **Ctrl+Shift+P**, select
   **Python: Select Interpreter**, then choose **Enter interpreter
   path** and enter ``/opt/conda/envs/scrimp/bin/python`` if
   that interpreter is not already selected. Depending on VS Code's
   extension version, the command/labels may differ slightly.
4. Confirm the environment in the **Dev Container's integrated
   terminal**, not in a host terminal:

   .. code-block:: bash

      python -c "import sys, scrimp; print(sys.executable); print('SCRIMP OK')"
      python simulations/first_simulation.py

   Expected first line of the interpreter check::

      /opt/conda/envs/scrimp/bin/python

5. If the terminal test passes but **F5** still cannot find Python,
   check that the Microsoft **Python** and **Python Debugger**
   extensions are installed **in Dev Container: SCRIMP**, not
   only in the local VS Code window. Restart the Python extension /
   reload the VS Code window, and explicitly select the interpreter
   above. This is an editor integration problem, not an indication
   that the numerical image needs rebuilding.

The starter's ``.devcontainer/devcontainer.json`` already declares
``python.defaultInterpreterPath`` as the Conda interpreter
**inside the container**. VS Code's environment-discovery UI may
still present a generic first-run message before selection
completes. The native Conda installation path uses the
**actual local** interpreter instead; do not force the container
path into ``.vscode/settings.json`` shared by both modes.

Automatic interpreter discovery and Pylance imports (v1.3.0)
--------------------------------------------------------------

The first real Ubuntu VS Code test confirmed two *editor-only*
issues even though SCRIMP and debugpy executed correctly:

- The Python Environments extension repeatedly proposed
  **"No Python found / Install uv and Python"**.
- Pylance reported `Import "scrimp" could not be resolved`
  despite a successful `import scrimp` and F5 debugging.

The template now declares **three complementary container-scoped
settings** (all in
``.devcontainer/devcontainer.json`` under
``customizations.vscode.settings``):

.. code-block:: json

   {
     "python.defaultInterpreterPath": "/opt/conda/envs/scrimp/bin/python",
     "python-envs.globalSearchPaths": ["/opt/conda/envs"],
     "python.interpreter.infoVisibility": "always",
     "python.analysis.extraPaths": ["/workspace/scrimp"]
   }

The interpreter and discovery root reference SCRIMP's **already
installed Micromamba environment**. The source search path is the
parent directory of the installed editable
``/workspace/scrimp/scrimp/`` package. Pylance can then analyze
the Python source directly, without relying on how ``pip install
-e .`` implements editable-install import hooks. The configuration
does **not** disable or suppress the ``reportMissingImports``
diagnostic for other imports.

**Updating a project created before this change:** its files are
independent of SCRIMP's Git repository. Pulling SCRIMP does **not**
alter an existing ``MySCRIMPProject``. Open
``MySCRIMPProject/.devcontainer/devcontainer.json`` and add
the three *new* settings to the existing
``customizations.vscode.settings`` object, keeping the complete
Docker image tag. Select **Dev Containers: Rebuild and Reopen in
Container** to apply the updated container-specific settings.
This recreates the editor's development container from the existing
image; it **does not recompile GetFEM, PETSc, SLEPc or SCRIMP**.
The local ``simulations/`` and ``outputs/`` remain mounted
and persistent.

**Verification:** inside the re-opened container, the Python status
bar should identify the selected interpreter. Open a Python file,
wait for Pylance analysis to refresh, and confirm the yellow
``reportMissingImports`` underline on ``import scrimp``
disappears. The ``Python: Select Interpreter`` list should
offer the SCRIMP environment; ``Python Environments: Refresh All
Environment Managers`` can refresh a stale discovery cache.
If the interpreter remains missing, check the **Python
Environments** output panel and use **Python Environments: Run
Python Environment Tool (PET) in Terminal** to inspect discovery
before creating/installing anything. If the interpreter is found
but Pylance still complains, use **Python: Restart Language
Server** and inspect the **Pylance** output log.

These settings match Microsoft's Python Environments discovery and
Pylance ``extraPaths`` guidance, but **automatic popup
suppression is not yet verified on a new VS Code launch**. The
new verification is a separate GUI acceptance gate. Do not treat
a Python language-server error as a failure of the already
passing Python/MPI runtime.

Container-specific settings belong in Dev Container metadata,
**not** shared ``.vscode/settings.json``: users who select a
native Conda interpreter must not inherit hard-coded Linux Docker
paths.

The host project is bind-mounted at /work. The first command
writes MySCRIMPProject/outputs/oscillator.csv **on the host**.
Closing VS Code or removing the container does not delete that file.

The included VS Code Run and Debug configurations support **F5**
for a serial simulation or the active Python file; MPI multiprocess
debugging requires a separate setup. VS Code Tasks also provide
the first simulation and a two-rank MPI smoke.

**Verified on Ubuntu Dev Containers (20 November 2026):** the user
opened `simulations/first_simulation.py`, started the
`SCRIMP: first simulation` F5 configuration, and successfully
stepped with **F10** from the DPHS initialization to the
`p = S.State("p", "Momentum", "lumped")` instruction. The
debugger's **Locals** pane showed the live `system` (DPHS)
and `q` (State), and **Call Stack** contained `main()`.
The terminal showed execution via
`/opt/conda/envs/scrimp/bin/python` and the Python
`debugpy` launcher.

For a first debugging session:

1. Open `simulations/first_simulation.py` in VS Code.
2. Click in the gutter beside `system = S.DPHS("real")`
   to add a **red breakpoint**.
3. Press **F5** and select **SCRIMP: first simulation**.
4. When paused (yellow current-line arrow), press **F10**
   to step over one Python statement and expand `system`
   or `q` in the **Locals** pane to inspect their attributes.
5. Press **F5** to resume execution, or **Shift+F5** to
   stop the debugger. On **20 November 2026**, the user
   confirmed that execution resumed and the simulation
   finished normally after stepping and inspecting the state.

This is a **serial Python debugging** validation. Debugging
simultaneous MPI worker ranks, and breakpoints inside compiled
GetFEM/PETSc/MUMPS C/C++/Fortran code, are different workflows
and are **not** qualified by this Python test.

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

- The candidate image **passed its Linux native build** and the
  1/2-rank CLI smoke. A **real VS Code Dev Container session** on
  Ubuntu subsequently ran the 10-step oscillator, wrote the CSV
  to the host, and started SCRIMP on two MPI ranks.
  **Serial Python debugging has also passed** with F5,
  a breakpoint, F10 stepping, Locals, Call Stack and normal
  completion on resume. Source-level settings now address
  automatic interpreter discovery and editable-SCRIMP imports;
  their actual fresh-window GUI retest remains open. Native
  compatibility, Windows/macOS/Apple Silicon, MPI multi-process
  debugging and power-balance/physics qualification remain
  separate gates.
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
