# SCRIMP v1.3.0 — Docker candidate build and cache qualification

**Status: build wiring prepared, native clean rebuild NOT YET EXECUTED.**
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
f6c09fb2edd112b6abe04e0485534f5e965d6f8f9ad853308b7120bfc6722ab5
```

(The digest is for the official `getfem-5.5.tar.gz` release archive;
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

**Warning:** a clean build of the new Dockerfile has not yet run on
the user's host. The first native build is the next qualification
action; no image has been promoted or released by committing these
files.

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

**Observed initial build blocker:** before the plugin was installed,
Docker reported `BuildKit is enabled but the buildx component is
missing or broken`. This occurs before the Dockerfile starts;
it is **not** evidence of an error in the native numerical libraries
or source patches. The new candidate image has not yet built.

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

The candidate packaging command for a future GitHub release asset
is `python3 scripts/package-scrimp-project-template.py`. The
source-only project ZIP is not published yet; no public image or
starter release artifact is claimed at this stage.

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

**Next action:** run the new candidate build and the two-rank smoke
on the Linux workstation, review both logs, then decide whether to
freeze the complete environment lock and continue with distribution.
