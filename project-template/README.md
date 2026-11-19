# My SCRIMP simulations — VS Code starter project

This is **your own project**, independent of the SCRIMP source repository.
All files in `simulations/`, `geometries/` and `outputs/` are yours.
Changing a simulation **never rebuilds** the SCRIMP Docker image.

## What is installed, and where?

- **On the computer:** VS Code; Docker (if using containers); this
  project directory, with editable Python scripts and persistent output.
- **In the Docker image:** SCRIMP, Python, GetFEM, Gmsh, PETSc, SLEPc
  and the numerical runtime. The image does not contain your simulations.
- **With native Conda:** SCRIMP and dependencies instead live in your
  activated Conda environment; this project directory stays the same.

## Option A: precompiled Docker image (after v1.3.0 is published)

The future released image is expected to be
`ghcr.io/g-haine/scrimp:1.3.0` (**not yet published**).
The final GitHub Release project ZIP will need its Dev Container
image name set to that **verified published** tag at packaging time.
Once published, install Docker (Docker Desktop on Windows/macOS)
and pull the image:

```bash
docker pull ghcr.io/g-haine/scrimp:1.3.0
```

Then change `image` in `.devcontainer/devcontainer.json` to
`ghcr.io/g-haine/scrimp:1.3.0`.

## Option B: Docker image built locally

Build SCRIMP via the official candidate build procedure (from a
temporary checkout of the GitHub repository, **not this project**).
It produces `scrimp-mpi:v1.3.0-rc1` for development/validation.
The template's default Dev Container image is that local candidate
tag **until the public image is qualified and released**.

Alternatively, edit `.devcontainer/devcontainer.json` to refer to
another **already built** local SCRIMP image. Always use the **complete
image reference including its tag**, such as
`"image": "scrimp-mpi:v1.3.0-rc1"`. Avoid a nested
`${localEnv:...:fallback-with-colons}` expression here: a prior
Dev Containers build truncated the candidate tag to `scrimp-mpi`,
causing an erroneous Docker Hub lookup for `:latest`.
For a newly created project, the generator's `--image IMAGE:TAG`
option embeds the chosen tag explicitly. If you created a project
before this fix, update **its own** `.devcontainer/devcontainer.json`;
pulling newer SCRIMP sources will not alter your personal project.

## Work in VS Code with Docker (both options)

1. Install VS Code and the **Dev Containers** extension. Docker must
   be running. On Windows, use Docker Desktop with WSL2.
2. Open **this directory** in VS Code.
3. Choose **Dev Containers: Reopen in Container**.
4. Select the interpreter `/opt/conda/envs/scrimp/bin/python`
   if VS Code has not selected it automatically.
5. Open the integrated terminal and run:

```bash
python simulations/first_simulation.py
mpiexec -n 2 python simulations/mpi_check.py
```

### If VS Code says "No Python found" or offers "Install uv and Python"

**Do not install another Python.** SCRIMP's numerical image already
contains Python and all numerical libraries in the **`scrimp` Conda
environment**. On the first Ubuntu Dev Containers launch, VS Code
briefly offered to download `uv` and install Python even though
the container's Python worked: the simulation and two-rank MPI
check both passed in the integrated terminal.

Click **Don't ask again** or **Cancel**, then open a Python file and
run **Python: Select Interpreter** (`Ctrl+Shift+P`). If necessary,
choose **Enter interpreter path** and enter:

```text
/opt/conda/envs/scrimp/bin/python
```

Verify from the **container terminal** (which should start with
`(scrimp)`):

```bash
python -c "import sys, scrimp; print(sys.executable); print('SCRIMP OK')"
```

The printed interpreter should be
`/opt/conda/envs/scrimp/bin/python`. The starter
`.devcontainer/devcontainer.json` sets
`python.defaultInterpreterPath` to this value already, but
the editor may still show a generic discovery popup at first launch.
If **F5** fails after the terminal succeeds, check the Microsoft
Python and Python Debugger extensions are installed **inside the
Dev Container**, reload VS Code, and re-select that interpreter.
Do **not** rebuild SCRIMP or create an `uv` environment to
address this editor prompt.

The first command is a **real, short serial port-Hamiltonian
mass–spring–damper simulation**. It saves `outputs/oscillator.csv`.
The second only checks MPI startup. To test MPI physics use a
separately validated SCRIMP PDE example (not this serial oscillator).

The `.vscode/launch.json` configurations let you press **F5**
to debug the first example or the active Python file. The task menu
also contains an ordinary run task and a 2-rank MPI startup task.
Parallel Python debugging needs more configuration than ordinary F5.

**All output is persistent on your computer**: Dev Containers mounts
this directory read-write at `/work`. SCRIMP uses the current working
directory to create `outputs/log`, `outputs/mesh`, `outputs/png`
and `outputs/pv`. Container deletion, image rebuilds and VS Code
restarts do not remove those files.

On Linux, Dev Containers attempts to match the container user's
UID/GID to yours. Ownership and bind-mount permissions still need
qualification on Windows, macOS and multi-user Linux machines.

## Option C: local native installation (no Docker)

Follow SCRIMP's native Conda/source installation instructions in
the official installation guide. Activate the resulting environment:

```bash
conda activate scrimp55  # or your actual SCRIMP Conda env name
python simulations/first_simulation.py
```

In VS Code, use **Python: Select Interpreter** to select the native
Conda interpreter. **Do not choose Reopen in Container** for this
installation mode. The same project files and `outputs/` persist
on the host.

## Customize and retain your simulations

- Put independent simulations under `simulations/`.
- Put your Gmsh geometry and other inputs under `geometries/`.
- Keep outputs under `outputs/` or an explicit subdirectory of this
  workspace. Large outputs are ignored by the included `.gitignore`.
- Use an external volume for very large datasets if needed;
  the mounted project path is not copied when building an image.
- Version your own project using Git if desired. It is **not**
  a checkout of the SCRIMP library.

## Release status

This starter project is being prepared on the
`release/1.3.0-preparation` branch. The public `v1.3.0` image
and native cross-platform installation are **not yet qualified**.
The **Linux VS Code Dev Container terminal run** passed the actual
first SCRIMP simulation, host-persistent CSV and two-rank MPI
startup on 19 November 2026. Python debugging with **F5** and
cross-platform VS Code behavior still require separate acceptance.
