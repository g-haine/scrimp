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
another **already built** local SCRIMP image.

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
`release/1.3.0-preparation` branch. The public `v1.3.0` image,
native cross-platform installation and Dev Containers end-to-end
tests are **not yet qualified**. The starter configuration and
scientific example still require first-runtime acceptance.
