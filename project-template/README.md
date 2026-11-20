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

### Python environment and Pylance in Dev Containers

The SCRIMP runtime has already validated Python, SCRIMP imports,
MPI and F5/debugpy under Ubuntu. Its interpreter is:

```text
/opt/conda/envs/scrimp/bin/python
```

Pylance once reported a **false** `import scrimp` warning because
SCRIMP is installed editably from `/workspace/scrimp`. The
`python.analysis.extraPaths` setting fixed the warning in
the actual Ubuntu VS Code GUI.

A different problem persisted: **"Invalid Python interpreter is
selected"**, accompanied by repeated offers to install Python
via `uv`. Even after clearing the workspace interpreter choice,
the Python extension logged an **active interpreter at the correct
path** but said *"No Python environment is set for this
resource"*. The **Python Environments** extension was managing
discovery, but did not assign our Micromamba environment to /work.

The container template therefore **opts out of that environment
manager**, retaining the original Python extension, Pylance and
Python Debugger. Its Dev Container settings include:

```json
"python.defaultInterpreterPath": "/opt/conda/envs/scrimp/bin/python",
"python.useEnvironmentsExtension": false,
"python.interpreter.infoVisibility": "always",
"python.analysis.extraPaths": ["/workspace/scrimp"]
```

This is an **editor-only** change: Python, Conda, GetFEM, MPI
and the already validated SCRIMP image are untouched.
**Confirmed on Ubuntu:** after adding the missing
`"python.useEnvironmentsExtension": false` line to the
existing user's Dev Container settings, VS Code no longer
reported an invalid/missing interpreter. The editable SCRIMP
import warning had already disappeared with `extraPaths`.
A clean-room new-user test and other platforms remain
unqualified.

**If you created your workspace before this fix:** edit *its own*
`.devcontainer/devcontainer.json` and set
`"python.useEnvironmentsExtension": false` inside
`customizations.vscode.settings`. Remove
`"python-envs.globalSearchPaths"` if present. Retain the
explicit interpreter path, Pylance search path and full
Docker image tag. Select **Developer: Reload Window** or,
if the new remote settings are not picked up,
**Dev Containers: Rebuild and Reopen in Container**.
This does not recompile the SCRIMP image or erase your project.

If the old interpreter selection is still cached, use
**Python: Clear Workspace Interpreter Setting** and
**Python: Select Interpreter** to choose the existing
`/opt/conda/envs/scrimp/bin/python`. **Do not click
"Install uv and Python"**: it would create an unrelated
environment. Do not disable Pylance missing-import diagnostics
or modify the shared `.vscode/settings.json`, which must also
work for native Conda users.

The first command is a **real, short serial port-Hamiltonian
mass–spring–damper simulation**. It saves `outputs/oscillator.csv`.
The second only checks MPI startup. To test MPI physics use a
separately validated SCRIMP PDE example (not this serial oscillator).

The `.vscode/launch.json` configurations let you press **F5**
to debug the first example or the active Python file. The task menu
also contains an ordinary run task and a 2-rank MPI startup task.
**Verified on Ubuntu Dev Containers:** F5 starts `debugpy`; a red
breakpoint at `system = S.DPHS("real")` suspends execution;
**F10** steps through state initialization; the `system` and
`q` objects are visible in **Locals**, and `main()` appears
in **Call Stack**. Resume with **F5** or stop with
**Shift+F5**. The *resume-to-completion* check **passed
on 20 November 2026**. Parallel MPI Python debugging needs
a different multi-process configuration and has not been qualified.

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
startup on 19 November 2026. On **20 November**, Python
debugging with **F5/debugpy**, breakpoints, F10 stepping,
Locals, Call Stack and successful **resume-to-completion** passed.
The Ubuntu user also confirmed **interpreter recognition
without errors** after opting out of Python Environments,
and Pylance's editable-import diagnostic is resolved. Clean-room
onboarding, multi-rank debugging and cross-platform behavior
remain separate acceptance checks.
