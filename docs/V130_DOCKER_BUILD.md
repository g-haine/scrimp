# SCRIMP v1.3.0 — Docker candidate build and cache qualification

**Status (19 November 2026): candidate image `scrimp-mpi:v1.3.0-rc1` built successfully from Git revision `96961d647c2caf41dd77d2a7b5926e86e86c4b69`. GetFEM 5.5 native compilation, FEM constructors, MPI-2 GetFEM, complex PETSc/SLEPc 3.26.1, distributed MUMPS (1/2 MPI ranks), a real serial SCRIMP model with host-persistent CSV, and a two-rank user-project MPI startup all PASS. The full scientific physics/scaling and actual VS Code Dev Containers acceptance gates remain open; no GHCR release or qualified-image replacement has occurred.**
**Latest candidate status (20 November 2026): native stack and
warm-cache qualification of `scrimp-mpi:v1.3.0-rc2` at
`ae96cb7334d7c88cd2e0dbd4d9615af66f1f4a61` **PASS**.
The targeted ComplexWarning, host CSV and GetFEM/MPI 1/2-rank smoke
passed, and the second build reused all Dockerfile stages. See the
[formal rc2 acceptance report](V130_DOCKER_RC2_ACCEPTANCE.md).
The post-rc2 PETSc `-info` source fix is **not** included in rc2
and will require rc3 qualification; no GHCR release occurred.

This document is part of [v1.3.0 release point 2](V130_RELEASE_PLAN.md).
Always build/test a **new tag**. The existing qualified
`scrimp-mpi:local` image must not be retagged, overwritten, pruned
or deleted until the candidate has passed its gates.

## 1. Purpose and cache boundaries

The source `docker/Dockerfile` is a BuildKit Dockerfile based on
`mambaorg/micromamba:2.10.0-ubuntu24.04`. It now has these boundaries:

| Layer | Inputs controlling invalidation | Persistent BuildKit download cache |
| --- | --- | --- |
| OS/Ubuntu packages | Dockerfile APT declarations/base image | apt indexes and archives |
| Micromamba/NumPy/PETSc/SLEPc native packages | `docker/environment.yml` | `/opt/conda/pkgs` |
| petsc4py/slepc4py native wrappers | `scripts/install-petsc-python-bindings.sh`, env layer | pip wheel/source downloads |
| GetFEM 5.5 MPI core and Python extension | `scripts/build-getfem-mpi.sh`, approved top-level `*.patch` | verified upstream GetFEM tarball |
| SCRIMP Python runtime, examples and tests | explicit `COPY` of allowlisted sources only | none required for native libs |

Only **explicitly selected source paths** are copied after GetFEM
and PETSc/SLEPc are built, so editing `scrimp/`, `sandbox/` or examples
does **not** invalidate expensive native layers. Documentation and
simulation artifacts are *not* copied into the runtime image. A change
to an approved `0000`–`0004` patch intentionally invalidates the
GetFEM layer, but *not* the earlier Conda/PETSc layers. Edits to the
patch README or to `experimental/` do not invalidate native GetFEM.

GetFEM 5.5 source is downloaded into a separate BuildKit cache using
an atomic `*.part` download. Its SHA-256 is checked **before**
extraction. The candidate Dockerfile pins:

```text
0e042770b86152b09be4bfd0ee160b8d10d49c3f8a0b9b06d55f1ffaa1d30fc4
```

(The digest is for the **reviewed Savannah** `getfem-5.5.tar.gz` release archive;
a mismatch **fails the build**, removes the bad cached archive and
must be investigated rather than bypassed.) The native build directory
is **not** a BuildKit cache: its compilation state is kept in the
Docker layer itself, so patch changes cannot silently reuse stale
objects.

BuildKit cache mounts persist across builds using the same Docker
builder and are **not included in the resulting image**.
BuildKit caches can be cleared by a Docker build-cache prune; no
SCRIMP script prunes them automatically.

## 1a. Physical isolation of the Docker context (critical)

**19 November 2026 regression caught on the Linux workstation:**
the old `.dockerignore` with a deny-all `**` followed by
directory negations such as `!sandbox/` **did not protect against
recursive simulation data**. The host Python preflight counted
212 apparently allowed files (1.78 MiB), but BuildKit transferred
**more than 45.60 GB** before the build was interrupted with Ctrl+C.
This was a real build-context regression; the earlier source-side
audit provided a false sense of safety. There was **no native
numerical compilation outcome** in this run.

**Permanent fix: Docker never receives the working tree as its
build context.** The root `.dockerignore` now contains only
`**` (deny everything). Direct
`docker build ... .` is intentionally unsupported; it cannot
access the sources needed by the Dockerfile.

The supported build driver
`scripts/build-v130-docker-candidate.sh` fetches a **clean GitHub
revision by default** (override via `SCRIMP_GITHUB_REF`). For
deliberate uncommitted library development only, set
`SCRIMP_DOCKER_SOURCE=local`. In *both* modes it creates an **empty
second temporary Docker context outside the repository**, with a
cleanup trap. Its audited staging script:

1. Reads `git ls-files --cached -z` (tracked files only), and
   matches exact, slash-aware allowlist patterns for installation
   metadata, SCRIMP Python sources, examples, sandbox diagnostics,
   tests, build scripts and five approved GetFEM patches.
2. Rejects symlink source paths, missing essential tracked inputs,
   unexpectedly oversized files, and context totals over 64 MiB.
3. Copies only these files into the temporary directory and checks
   the staged byte total before invoking Docker.
4. Invokes `docker buildx build --load ... "$STAGED_CONTEXT"`,
   **not** `docker build ... .`; cleanup removes the staging
   directory after completion, failure or Ctrl+C.

This changes the protection from pattern-dependent Dockerignore
semantics to **physical separation**. Nothing resembling
`sandbox/outputs`, `.build`, simulation archives, plots, videos,
untracked local scripts or user data exists in the supplied
temporary directory. A separately generated tiny `.dockerignore`
inside that directory has no broad negations.

To inspect the tracked-source selection without Docker:

```bash
python3 scripts/audit-v130-docker-context.py
```

To materialize it independently for inspection (remove only
that explicitly created temporary directory afterwards):

```bash
CTX="$(mktemp -d -t scrimp-context-inspect.XXXXXXXX)"
python3 scripts/audit-v130-docker-context.py --stage-to "$CTX"
du -sh "$CTX"
find "$CTX" -type f | wc -l
rm -rf -- "$CTX"
```

The new driver performs staging automatically. At the previous
commit, the repository's Git-tracked allowlisted inputs totaled
approximately **1.5 MiB**. The first fixed BuildKit run must now
show a transfer of that order of magnitude, not gigabytes. If
it does not, stop immediately and investigate before compiling.

The manual GitHub Actions Docker build also uses physical staging,
to avoid bypassing this safety boundary.

## 1b. BuildKit cache ownership and CPU-only complex solver variants

**First clean build finding (19 November 2026, source bcab1507):**
the physically staged GitHub context contained only 189 files /
1.51 MiB; BuildKit transferred **1.60 MB**, closing the former
>45 GB regression. Micromamba created `scrimp`, and the NumPy-2 /
PETSc-3.26.1 / SLEPc-3.26.1 Python ABI checks **passed**.
The subsequent GetFEM release tarball download failed with curl
exit **23** (`Permission denied` on
`/tmp/scrimp-getfem-tarballs/getfem-5.5.tar.gz.part.1`).
The same build also reported that pip had disabled its cache and
micromamba could not open `/opt/conda/pkgs/pkgs.lock`.

**Cache repair:** Before switching to the non-root `mambauser`,
the Dockerfile now temporarily mounts **all three existing BuildKit
package caches as root**, reassigns their ownership to
`mambauser:mambauser`, and then explicitly checks the resulting
non-root writability before installation. This preserves downloaded
packages and does not prune caches or alter the host Conda
environment. This correction still needs a real Docker re-run.

**Avoid unwanted CUDA packages:** The first solve selected
`petsc-3.26.1-cuda13_complex` and
`slepc-3.26.1-cuda13_complex`, which pulled a **2 GB**
environment with substantial unused CUDA dependencies. The
general Linux CPU/MPI candidate now requires
`petsc=3.26.*=complex_*` and `slepc=3.26.*=complex_*`:
the leading `complex_` excludes `cuda12_complex_` and
`cuda13_complex_`. A Docker build stage checks the resolved
Conda metadata and fails if a GPU PETSc/SLEPc build or large
CUDA numerical runtime libraries sneak back in. A separate
GPU-specific environment may be provided in a future release;
the CPU candidate must work without CUDA.

**Cache reuse after the repair:** Changing `docker/environment.yml`
intentionally invalidates the Micromamba and later image layers,
so `petsc4py` and `slepc4py` must be built again for the
**new CPU libraries**. Generic downloads in the BuildKit caches
are preserved where their package hashes still match.
The previous CUDA binding builds are NOT ABI-compatible with the
new CPU variant and must not be reused.

**No full build, native GetFEM patch acceptance, MPI smoke or VS
Code runtime acceptance is claimed yet.** Continue with the normal
candidate builder after pulling the corrected branch.

## 1c. GetFEM 5.5 archive discrepancy — reviewed and resolved

**Source provenance (19 November 2026):** the Savannah download
`https://download-mirror.savannah.gnu.org/releases/getfem/stable/getfem-5.5.tar.gz`
is a valid gzip release archive of **16,911,362 bytes**, with SHA-256
`0e042770b86152b09be4bfd0ee160b8d10d49c3f8a0b9b06d55f1ffaa1d30fc4`.
It is not byte-identical to the earlier Fossies snapshot
(16,911,397 bytes, SHA-256
`f6c09fb2edd112b6abe04e0485534f5e965d6f8f9ad853308b7120bfc6722ab5`).
The pinned immutable upstream Git reference is
`v5.5` at `c242f7fe0989e9a90846c63eaccc6b495d774724`.

**Independent comparison performed on the developer's workstation:**
1,170 files common to the release archive and upstream Git tag
were compared byte-for-byte. Only two shared files differ:

- `CMakeLists.txt`: corrects a stale version string from **5.4.2**
  to **5.5.0**, unifies GMM/GetFEM version handling, and cleans
  `check_include_file_cxx` syntax. Our GetFEM build uses
  **Autotools**, not CMake.
- `interface/src/gf_util.cc`: moves the existing
  `#if defined(GMM_USES_MUMPS)` binary availability test into a
  `MUMPS_LINKED` macro defined through `#ifdef GMM_USES_MUMPS`.
  Both report 1 when the preprocessor symbol is defined and 0
  otherwise. No change to the MPI numerical algorithm was found.

The remaining **six absent native-looking files** are
`interface/src/gf_poly.cc`, three
`interface/src/octave/` headers, and
`src/getfem_enumeration_dof_para.cc` /
`src/getfem_interpolation_on_deformed_domains.cc`.
We checked the upstream Autotools `src/Makefile.am` and
`interface/src/Makefile.am`: these missing files are not
among the compiled source lists for our GetFEM MPI/Python
configuration. In particular, `getfem_enumeration_dof_para.cc`
is explicitly commented out. Their omission is therefore
acceptable for **this** build, not a claim that the archive
contains the entire Git source tree.

**All five SCRIMP GetFEM patches `0000`–`0004` applied
sequentially without errors** to a disposable extraction of the
Savannah archive. The new accepted exact digest is recorded in
`docker/Dockerfile`; the mandatory pre-extraction SHA-256 check
remains in place. Previous Conda/PETSc layers are unaffected
because the GetFEM digest build argument is declared **after**
those layers.

**Reproducible source audit** (separate from native compilation):

```bash
cd ~/Simulations/scrimp
git pull --ff-only
python3 scripts/audit-getfem55-release-source.py \
  2>&1 | tee sandbox/outputs/getfem55-release-source-audit.log
```

The updated audit requires this exact Savannah archive digest,
the immutable upstream tag, *exactly these two reviewed shared-file
changes*, *exactly these six reviewed source omissions*, and
clean application of all five patches. Any additional difference
or missing patch target **fails closed**. This is a reviewed
distribution variant, not a claim of bit-for-bit identity to
the Git tag or of an upstream cryptographic signature.

**Still required:** run the candidate Docker build with its
new strict digest, verify GetFEM MPI level 2, the complex ABI,
BDM/Nédélec FEM qualification, 1/2-rank SCRIMP/MUMPS numerical
smokes and the independent VS Code user-project output.
No release image has been promoted or published.

## 2. Build a candidate without disturbing the qualified image

From a temporary checkout of SCRIMP's preparation branch,
or from an existing developer checkout (the actual compiled sources
are fetched from **GitHub**, not your local worktree):

```bash
cd ~/Simulations/scrimp
git pull --ff-only
bash scripts/build-v130-docker-candidate.sh
```

For an exact revision, set `SCRIMP_GITHUB_REF` to an immutable
Git commit hash. The build logs the selected commit and adds it to
the image's OCI revision label. If the branch advances, the next
invocation may select a newer revision unless pinned. Only the
deliberate `SCRIMP_DOCKER_SOURCE=local` mode sees tracked edits
from the developer checkout, and those are still staged in a
separate context.

Default tag: `scrimp-mpi:v1.3.0-rc1`; default native Make jobs: 4.
Override deliberately, for example:

```bash
MAKE_JOBS=2 SCRIMP_DOCKER_CANDIDATE_IMAGE=scrimp-mpi:v1.3.0-rc2 \
  bash scripts/build-v130-docker-candidate.sh
```

The build driver **requires** `docker buildx build --load`;
it deliberately refuses to fall back to the legacy Docker builder.
It never pushes images and refuses `scrimp-mpi:local` and
`scrimp-mpi:latest` as output tags.
If the current user has permission to run Docker without `sudo`,
set `SCRIMP_DOCKER_NO_SUDO=1`.

Output is logged to `sandbox/outputs/v130-docker-build.log`.
On success, the script records the immutable Docker image ID and
`org.opencontainers.image.revision` label. Keep this log with the
release evidence. A build interrupted during native compilation may
restart that step but should **not redownload** a valid cached
tarball or cached pip/Conda packages.

**Validated candidate (19 Nov):** the full image build and native,
FEM/MPI 1/2-rank, independent user-workspace and persistent host-output
smokes **passed**. The historical cache-ownership and GetFEM source
checksum blockers are closed. The actual VS Code GUI/Dev Container
interaction, whole-science MPI parity, physical power balances,
cross-platform compatibility and warm-cache reproducibility remain
**separate acceptance gates**. This candidate has not been published
or promoted over `scrimp-mpi:local`.

### Buildx prerequisite (Ubuntu 22.04 / Docker from Ubuntu repositories)

The candidate Dockerfile uses BuildKit cache mounts and requires the
**Docker Buildx CLI plugin**. A plain `docker build` with
`DOCKER_BUILDKIT=1` is not a fallback when Buildx is unavailable.

On the currently tested host, Docker is installed as Ubuntu's
`docker.io` package (29.1.3-0ubuntu3~22.04.2), **not** Docker CE.
The correct matching package is `docker-buildx`:

```bash
sudo apt-get install docker-buildx
sudo docker buildx version
git pull --ff-only
bash scripts/build-v130-docker-candidate.sh
```

Do **not** use `docker-buildx-plugin` on this host: that is the
package name used by the upstream Docker CE repository, which is not
configured here. The build runner now checks `sudo docker buildx
version` before any image construction and fails with an actionable
message rather than starting a build that cannot work.

**Historical initial build blocker:** before the plugin was
installed, Docker reported `BuildKit is enabled but the buildx
component is missing or broken`. This was corrected by installing
Ubuntu's `docker-buildx`; the candidate image now builds.

## 2a. After installation: independent user projects in VS Code

**An image is not a simulation workspace.** SCRIMP's own `sandbox/`
belongs to the source repository; users instead create their own
project with an editable `simulations/` folder, `geometries/`,
VS Code `.devcontainer/` setup and writable, persistent `outputs/`.

The project starter is under `project-template/`; make a workspace
in a **new local directory** with:

```bash
python3 scripts/create-scrimp-project.py "$HOME/MySCRIMPProject"
```

The project template's Dev Container defaults to
`scrimp-mpi:v1.3.0-rc1` during development. With a published
GHCR image, point the Dev Container to the corresponding qualified
tag instead. For native Conda, select its Python interpreter in VS
Code and do not reopen the project in a container.

The first simulation saves `outputs/oscillator.csv` with
time, displacement, momentum and Hamiltonian, **in the user's
directory**. The MPI script is only a startup smoke. Neither script
causes a Docker rebuild. See
[the VS Code user guide](source/vscode.rst).

The candidate packaging command is
`python3 scripts/package-scrimp-project-template.py`, which creates
a **preview** ZIP pointing to the local candidate image.
**At release publication only**, prepare the ZIP with the actually
validated and published image tag:

```bash
python3 scripts/package-scrimp-project-template.py \
  --image ghcr.io/g-haine/scrimp:1.3.0
```

The command only creates a local ZIP and never publishes images or
release assets. A corresponding public GHCR image and release asset
are **not available yet**; they remain gated on the complete
release validation.

## 3. Run the candidate smoke and capture reproducible inputs

After the candidate image has built:

```bash
bash scripts/smoke-v130-docker-candidate.sh
```

This checks:

1. The GetFEM MPI parallelism level remains 2, PETSc remains
   complex, and RT, BDM, Nédélec-I/II known k=1/2 constructors have
   their expected local DoFs. Nédélec-II k>=3 is deliberately
   **rejected**, not a supported element.
2. The existing `sandbox/smoke_modern_stack.py` ABI and MUMPS solve
   pass inside the candidate with **1 and 2 MPI ranks**.
3. A **new independent host-side project** is created in a temporary
   directory. The image runs the actual mass–spring–damper SCRIMP
   example against this writable bind mount, and the smoke refuses
   to pass unless `outputs/oscillator.csv` is present **on the
   host**. It then checks an MPI startup in the same workspace.
   Neither test reads or writes the developer's real sandbox.
4. The Conda resolved **explicit package spec**, GetFEM top-level
   patch hashes, image ID and source Git revision are preserved.

The smoke log is
`sandbox/outputs/v130-docker-smoke.log`. The resolved package spec
and patch checksum list are exported to:
`sandbox/outputs/v130-docker-resolved-spec.txt` and
`sandbox/outputs/v130-getfem-patch-sha256.txt`.

**Scope limitation:** these are candidate-stack, ABI, user-workspace
and FEM constructor smokes. They do not replace the earlier full wave
1/2/4/6-rank strong-scaling qualification or certify physical Maxwell
power balance, which are distinct release workstreams.

## 3a. First complete native and user-workspace qualification — PASS

**Evidence supplied 19 November 2026 (build revision
`96961d647c2caf41dd77d2a7b5926e86e86c4b69`):**

- Buildx exported/loaded
  `scrimp-mpi:v1.3.0-rc1` with image ID
  `sha256:ab08d9beaa60a1685b7e8df1e69f5a085f5571ca46c2b7d8554f53ec42ace3e8`
  (the recorded ID is a manifest-list digest, not a statement
  about other architectures).
- Native GetFEM MPI level **2**; PETSc/SLEPc **3.26.1 complex**,
  NumPy **2.5.5**, SciPy **1.18.1**, mpi4py **4.1.2**.
- Local simplex descriptors/DoF counts pass for RT k=0 (2D/3D),
  BDM k=1 (2D/3D) and k=2 (3D), Nédélec I k=1 (2D/3D)
  and k=2 (3D), Nédélec II k=1/2 (2D/3D). Unsupported
  Nédélec-II k>=3 correctly raises an exception.
- PETSc distributed MUMPS factor/solve passes at 1 and 2 ranks.
- A **real** independently scaffolded user simulation solves the
  four-DoF mass–spring–damper system for 10 steps
  (`t=0..0.1`, `dt=0.01`), with zero rejections/failures.
  It saves 11 rows with `time,q,p,H` to a host-mounted
  `outputs/oscillator.csv`; the smoke checks its presence
  **outside** Docker.
- The same independent user project imports SCRIMP and starts on
  **two MPI ranks**; both print their rank and return success.

**Observed non-blocking diagnostics (not hidden):**
the real GetFEM model receives PETSc's complex trial vector.
The legacy `dtype=np.float64` cast in
`DPHS._sync_lumped_trial_state` produced a NumPy
`ComplexWarning` even though the real first-user simulation
completed. The source now **explicitly copies the real part**,
matching the existing TS monitor's behavior; this is a narrow
post-image change requiring an image rebuild to validate.
The MPI startup printed a PETSc option-database warning about
unused `-info`; its origin and need for suppression have not
yet been qualified. GetFEM also prints a trace during the
stack smoke; neither message affects this acceptance result.

**Subsequent VS Code Dev Container Ubuntu qualification (19 Nov):**
the user opened the independent project inside a real VS Code
Dev Container with `(scrimp)` active and `/work` mounted.
The integrated terminal ran the ten-step `first_simulation.py`
successfully and created `/work/outputs/oscillator.csv`.
`mpiexec -n 2 python simulations/mpi_check.py` also passed.
The candidate image therefore passes **both the CLI host-volume
smoke and the Linux VS Code integrated-terminal workflow**.

The editor separately displayed **"No Python found — Install uv
and Python"** despite the working container interpreter at
`/opt/conda/envs/scrimp/bin/python`. This generic VS Code
discovery/installer prompt is **not evidence that Python or
SCRIMP is absent**. The user should dismiss it and, if needed,
select the existing interpreter with **Python: Select Interpreter**.
No extra Python installation or image rebuild is appropriate;
see [the first-user guide](source/vscode.rst). The template's
Dev Container already declares the path.

**VS Code F5/debugpy qualification (20 Nov): PASS for serial Python.**
A real Ubuntu Dev Container debug session used the interpreter
`/opt/conda/envs/scrimp/bin/python` and started under the
`debugpy` launcher. The first SCRIMP simulation accepted a
breakpoint on the DPHS constructor line, **F10 / Step Over**
advanced through model initialization to line 31, and
`system` (DPHS) and `q` (State) were inspectable under
**Locals**. `main()` appeared in **Call Stack**. This confirms
actual VS Code breakpoint/step/inspection operation, beyond
the integrated-terminal simulation acceptance.
On **20 November**, the paused simulation was resumed with
F5 and completed normally. Multi-rank MPI debugging remains
untested. The first-use editor exposed repeated **"No Python
found / Install uv"** prompts and a false Pylance
`reportMissingImports` for editable `scrimp`, even though
the Python interpreter worked. The initially attempted Micromamba discovery path
(`python-envs.globalSearchPaths=/opt/conda/envs`) was
**insufficient**. The Python log knew the absolute interpreter,
but `ms-python.vscode-python-envs` reported no environment for
the workspace. The project template now explicitly selects
`/opt/conda/envs/scrimp/bin/python`, opts out of that
new environment manager using
`"python.useEnvironmentsExtension": false` (keeping Python,
Pylance and debugpy), and configures
`python.analysis.extraPaths=/workspace/scrimp` for
editable-package source analysis.

**GUI retest PASS (20 Nov):** the user reported that adding the
missing opt-out setting to their *existing personal project*
removed the interpreter error; the Pylance import warning had
already disappeared. The development container retained all
previous runtime, MPI and F5 functionality. This was an editor
configuration repair, **not** a numerical image rebuild.
A fresh-user onboarding run is still required before broad
cross-platform distribution claims.

**Qualification boundary:** these results are *smokes*, not proof of
Hamiltonian or power-balance preservation, parallel strong scaling,
full Maxwell physics, MPI multi-process debugging or
cross-platform VS Code onboarding. Follow those distinct roadmap
gates before publishing an image.

## 4. Verify the warm-cache behavior

With no source changes, a second
`bash scripts/build-v130-docker-candidate.sh` should show `CACHED`
on native build layers and finish without repeating them. A
SCRIMP-only Python/source change should rebuild the **final SCRIMP**
layer and retain the native PETSc and GetFEM layers. Changing an
approved GetFEM patch should rebuild GetFEM and later layers but
reuse the earlier PETSc/SLEPc build.

Record the successful repeated build in the release audit before
describing the Docker image as cache-qualified. The rebuilt candidate
should carry the corresponding Git revision label.

### 4a. Next-image qualification — exact source SHA, warm-cache evidence (20 Nov 2026)

The Ubuntu/VS Code workflow was accepted with source revision
`96961d647c2caf41dd77d2a7b5926e86e86c4b69`. The subsequent
`DPHS._sync_lumped_trial_state` real-part / `ComplexWarning` fix is
**not in that image**. Preserve `scrimp-mpi:local` and preferably also
the previously accepted `scrimp-mpi:v1.3.0-rc1` tag. Use an independent
`scrimp-mpi:v1.3.0-rc2` candidate; no GHCR publication.

On the validated Ubuntu host, with Docker accessible without sudo:

~~~bash
cd ~/Simulations/scrimp
git switch release/1.3.0-preparation
git pull --ff-only
mkdir -p sandbox/outputs
export SCRIMP_GITHUB_REF="$(git rev-parse HEAD)"
export SCRIMP_DOCKER_CANDIDATE_IMAGE=scrimp-mpi:v1.3.0-rc2
export SCRIMP_DOCKER_NO_SUDO=1

docker image inspect --format '{{.Id}}' scrimp-mpi:local \
  > sandbox/outputs/v130-qualified-local-before.txt
bash scripts/build-v130-docker-candidate.sh
cp sandbox/outputs/v130-docker-build.log \
   sandbox/outputs/v130-docker-build-first-rc2.log
bash scripts/smoke-v130-docker-candidate.sh
cp sandbox/outputs/v130-docker-smoke.log \
   sandbox/outputs/v130-docker-smoke-rc2.log
bash scripts/build-v130-docker-candidate.sh
cp sandbox/outputs/v130-docker-build.log \
   sandbox/outputs/v130-docker-build-warm-rc2.log

docker image inspect --format '{{.Id}}' scrimp-mpi:local \
  > sandbox/outputs/v130-qualified-local-after.txt
cmp sandbox/outputs/v130-qualified-local-before.txt \
    sandbox/outputs/v130-qualified-local-after.txt
docker image inspect --format \
  '{{.Id}} {{index .Config.Labels "org.opencontainers.image.revision"}}' \
  scrimp-mpi:v1.3.0-rc2
~~~

If the host requires sudo, omit `SCRIMP_DOCKER_NO_SUDO=1` and
prefix the explicit `docker image inspect` commands with `sudo`.
The build and smoke scripts will then use sudo themselves.
**Keep the Git SHA, `MAKE_JOBS`, Docker/Buildx builder, GetFEM
digest and Conda inputs unchanged** across the two builds. Never
use `--no-cache`, prune BuildKit caches or rebuild directly from
the live tree. The driver fetches a GitHub snapshot and constructs
a physically isolated, Git-tracked allowlisted context.

**Acceptance gates and recorded evidence:**

1. **Exact provenance / safety:** build stdout shows the chosen
   full Git SHA, with the OCI `revision` label matching it.
   BuildKit transfers a small context (previously ~1.6 MB).
   No untracked experiments, simulations or `sandbox/outputs`
   enter the context. Keep both complete build logs and image IDs.
2. **Native cache reuse:** the *second* build must report `CACHED`
   for the Conda solve/install, PETSc/SLEPc Python bindings and
   GetFEM 5.5 compile layers. No native rebuilding or fresh source
   download should occur. A previously interrupted native `RUN`
   may execute again on the first build; cache mounts should still
   preserve downloaded inputs. A miss in the warm run is a failed
   cache gate requiring diagnosis, **not** a reason to delete caches.
3. **Targeted numerical/warning regression:** the candidate smoke
   runs the real starter simulation under
   `warnings.simplefilter("error", np.exceptions.ComplexWarning)`.
   Reappearance of the warning is fatal. Host output must contain
   11 finite rows at t=0..0.1 and satisfy
   `H=5*q*q/2+p*p/(2*1.5)` at each saved point to numerical
   tolerance (`rel_tol=1e-10`, `abs_tol=1e-12`). This is an
   **internal energy formula** check, not a power-balance proof.
4. **Native / MPI stack:** maintain GetFEM MPI level 2,
   NumPy 2 and complex PETSc/SLEPc; validated RT/BDM/Nédélec
   I/II simplex k=1/2 descriptors and Nédélec-II k>=3 rejection;
   distributed MUMPS 1/2-rank solves, independent host-mounted
   simulation and 2-rank MPI startup.
5. **No regression to the qualified image:** the ID of
   `scrimp-mpi:local` is unchanged; neither `latest` nor GHCR
   receives this candidate. The existing F5/F10/Pylance/interpreter
   GUI acceptance remains valid for the *workflow*; it does not
   replace this new-image numerical check.
6. **Qualification boundary:** a warm-cache rebuild is *not*
   bitwise-clean-room reproducibility. Base-image digest/platform
   and exact Conda lock are still unpinned release requirements.
   Full physical DAE/power and Maxwell 3D acceptance are separate.

**PETSc unused `-info` option:** the earlier 2-rank starter
reported an unused-option warning. Do not suppress it without
identifying the source. Compare a bare petsc4py 2-rank process
against the SCRIMP starter, recording stderr, exit codes, `sys.argv`,
the PETSc option database and any `PETSC_OPTIONS` passed to the
container. For example:

~~~bash
docker run --rm scrimp-mpi:v1.3.0-rc2 \
  mpiexec --oversubscribe -n 2 python -c \
  'import sys, petsc4py; petsc4py.init(sys.argv); from petsc4py import PETSc; print("argv=",sys.argv,"options=",PETSc.Options().getAll(),flush=True)' \
  2>&1 | tee sandbox/outputs/v130-petsc-info-bare.log
~~~

Compare with the saved `v130-docker-smoke-rc2.log` and
`project-template/simulations/mpi_check.py`. This is a diagnostic,
not a demonstrated fix. Do not change normal PETSc option handling
until the cause is established.

**Source CI audit:** commit `170f26d1` failed in the initial
Python compilation step because it still named the *deleted*
`scripts/patch-getfem-mpi-timing.py`. The workflow and its test
were corrected to target the approved
`0000-mpi-assembly-timing.patch`. The new branch commit must
obtain its own green source CI run; that older failure is not a
native GetFEM or physics failure.


### 4b. Actual rc2 host results and next candidate (20 Nov 2026)

The requested exact-SHA rebuild and second identical build have now
been run natively on the Ubuntu developer host: **PASS** for physical
context isolation, native GetFEM/PETSc layer cache reuse and the
strictened runtime smoke. Read
[V130_DOCKER_RC2_ACCEPTANCE.md](V130_DOCKER_RC2_ACCEPTANCE.md)
for the full matrix, exact revision, numerical checks and image
hashes. The `-info` diagnostic investigation in 4a is **closed at
the source level**: `scrimp/utils/linalg.py` was explicitly
injecting the PETSc option at import. The post-rc2 source
fix removes that side effect, and the two-rank Docker smoke now
asserts `not PETSc.Options().hasName("info")`.

The rc2 runtime manifestation is **not fixed retroactively**:
build a separate `scrimp-mpi:v1.3.0-rc3` candidate from an exact
Git commit and run the new smoke, retaining rc1/rc2 and
`scrimp-mpi:local`:

~~~bash
cd ~/Simulations/scrimp
git switch release/1.3.0-preparation
git pull --ff-only
mkdir -p sandbox/outputs
export SCRIMP_GITHUB_REF="$(git rev-parse HEAD)"
export SCRIMP_DOCKER_CANDIDATE_IMAGE=scrimp-mpi:v1.3.0-rc3
# Set SCRIMP_DOCKER_NO_SUDO=1 only if Docker runs without sudo.
bash scripts/build-v130-docker-candidate.sh
bash scripts/smoke-v130-docker-candidate.sh
~~~

**Before the build**, save the `scrimp-mpi:local` image ID;
afterwards compare it with the new ID using the same Docker identity
(sudo or non-sudo) to establish tag preservation. Keep all logs.
No long physics/performance campaign is justified solely by
removing an import-time optional PETSc diagnostic.

**BuildKit manifest nuance:** the first and warm rc2 exports
share the exact same main OCI manifest digest
(`sha256:82a74c8e35fa49a1f447a97184c666745e8ead4d96390708c13522ac167d59a8`)
and image config
(`sha256:08d06d2b7317e4a2f10ddfd6d3fe08e4ae35714e47fda7bed62d7fd5d7368fce`).
Their BuildKit attestation manifests and resulting manifest-list IDs
**differ**, so the log proves runtime-cache reuse, *not* a stable
OCI manifest-list digest across builds. Do not claim bitwise image
reproducibility until attestation behavior and the clean-room
base-image + Conda lock are qualified.

**Packaging observation:** `pip` still reports `scrimp-1.2.0`
for the candidate, matching the current `setup.py`. The version
change is intentionally gated with the release packaging freeze.

### 4c. rc3 native `-info` correction — PASS; harness-only `-c` retest

The **20 November 2026 rc3 build and initial smoke** at immutable Git
revision `1d6d0c90115213298047b28237a7f38087563a35` have
**passed**: all native layers `CACHED`, correct GetFEM MPI level 2,
complex PETSc/MUMPS 1/2-rank execution, real oscillator, host CSV and
NumPy `ComplexWarning` gate. Both MPI import ranks confirmed no
implicit PETSc `-info`. Before/after `scrimp-mpi:local` image IDs
matched. Full evidence:
[V130_DOCKER_RC3_ACCEPTANCE.md](V130_DOCKER_RC3_ACCEPTANCE.md).

The final message `Option left: name:-c ... source: command line`
came from using `python -c` **inside our diagnostic check**, not
from SCRIMP's runtime code. The smoke script now creates a disposable
Python file under the isolated user test workspace, bind-mounts it at
runtime (never in the image context), and rejects any PETSc unused
option warning. This harness correction requires **no new Docker
image build**. On the Ubuntu host:

~~~bash
cd ~/Simulations/scrimp
git switch release/1.3.0-preparation
git pull --ff-only
SCRIMP_DOCKER_CANDIDATE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
SCRIMP_DOCKER_NO_SUDO=1 \
  bash scripts/smoke-v130-docker-candidate.sh
cp sandbox/outputs/v130-docker-smoke.log \
   sandbox/outputs/v130-docker-smoke-rc3-warning-free.log
cp sandbox/outputs/v130-petsc-options-smoke.log \
   sandbox/outputs/v130-petsc-options-rc3.log
~~~

Omit `SCRIMP_DOCKER_NO_SUDO=1` on hosts using `sudo docker`.
**Final diagnostic evidence received (20 November 2026): PASS.**
The developer executed the corrected host-side two-rank import
check against existing rc3. Both ranks printed
`PASS: PETSc import options clean`, followed by
`PASS: PETSc -info and -c absent, no unused-option warning.`
The excerpt confirms this targeted test; the previously supplied
full rc3 native/MPI smoke remains the supporting evidence for
other gates. **Ubuntu Docker native/cache qualification is closed**,
without rebuilding rc3. Exact evidence:
[V130_DOCKER_RC3_ACCEPTANCE.md](V130_DOCKER_RC3_ACCEPTANCE.md).

## 5. Reproducibility and work not done yet

- **Pinned today:** base image *version tag*, PETSc/SLEPc native
  `3.26.*` complex variant, Python major/minor, and GetFEM 5.5
  source SHA-256. The exact solved Conda packages are exported.
- **Still needed before publication:** choose and commit a
  platform-specific explicit Conda lock from the tested image;
  verify a fresh build from that lock; pin the actual base-image
  digest/platform; and decide on a downloadable prebuilt image
  versus rebuilding native dependencies.
- **Image size:** this initial cache-layer refactor still contains
  the installed native build toolchain and GetFEM worktree. A
  multi-stage production-runtime slimming pass can follow once
  the new build is stable; do not change the known-good numerical
  stack and Docker layout simultaneously.
- **Distribution:** Ubuntu/Debian launcher, Windows
  Docker Desktop/WSL2 instructions and macOS/Apple Silicon support
  and limitations are **subsequent parts of release point 2**.
  They are not yet validated by this Linux-host candidate.
- **Higher polynomial order:** the accepted FEM scope is k=1,2.
  For k>=3 extension work see
  [the deferred higher-order plan](V130_GETFEM_HIGHER_ORDER_TODO.md).

**Next action:** advance the remaining reproducible
distribution (base image digest, Conda lock and fresh build) and
scientific qualification gates (initial DAE consistency, complete
wave trajectories, port-Hamiltonian power balance and Maxwell 3D).
Preserve rc1/rc2/rc3 and `scrimp-mpi:local`; no rebuild or repeated
large-wave performance campaign is justified by the now-closed
PETSc import diagnostic.
