# SCRIMP v1.3.0 — Docker-first distribution and release gate

**Decision, 27 November 2026:** precompiled Docker is the **default
user installation** for v1.3.0. Locally built Docker remains the
developer/reproducibility alternative; native Conda remains an
advanced/legacy alternative. A final user needs only Docker, a
version-matched standalone project ZIP and optionally VS Code with Dev
Containers. Never require cloning the SCRIMP library into a simulation
workspace.

## Proven, not yet proven

- **Qualified earlier, with bounded scope:** local Ubuntu image
  `scrimp-mpi:v1.3.0-rc3`, GetFEM 5.5, PETSc/SLEPc complex,
  MUMPS, MPI 1/2/4/6 where explicitly tested, scientific model
  checks and the independent Dev Container preview.
- **Canonical examples:** full nine-case user-output campaign at
  `20d79d0d` passed with no generic output siblings.
- **Documentation:** `cf75f4f761c2854788e5ef519375b31e7f040509`
  passed the strict Sphinx HTML/LaTeX build (zero Sphinx warnings)
  and compiled a real 168-page PDF. This is **not** an image
  distribution qualification.
- **Still pending:** new distribution image built from final sources;
  its exact OCI revision and installed Python version; standalone
  release ZIP and host-volume permissions; a genuinely fresh-machine
  install; GUI/F5 acceptance with the new image; wheel build/metadata
  audit; registry publication and image digest verification; optional
  tests on Windows/macOS. No public v1.3.0 image, release tag or ZIP
  has been published.

## Candidate construction — do not replace rc3 or local images

Run from the SCRIMP source checkout on Ubuntu with working Docker
Buildx. Use a **new, unique local tag**, never the validated `rc3`
or `scrimp-mpi:local`. Build sources from an exact Git SHA (not a
moving branch) and keep the pre-existing Docker cache:

```bash
cd ~/Simulations/scrimp
git pull --ff-only
REV="$(git rev-parse HEAD)"
IMAGE="scrimp-mpi:v1.3.0-dist-rc1"

SCRIMP_DOCKER_NO_SUDO=1 \
SCRIMP_DOCKER_CANDIDATE_IMAGE="$IMAGE" \
SCRIMP_DOCKER_SOURCE=github \
SCRIMP_GITHUB_REF="$REV" \
bash scripts/build-v130-docker-candidate.sh

SCRIMP_DOCKER_NO_SUDO=1 \
SCRIMP_DOCKER_CANDIDATE_IMAGE="$IMAGE" \
bash scripts/smoke-v130-docker-candidate.sh
```

The build script physically stages only allowlisted Git-tracked files
into a fresh Docker context, excluding user simulations and untracked
artifacts. It never pushes, prunes or deletes images. Verify the
resolved package spec and approved patch SHA256 files from the native
smoke. The Dockerfile currently uses an editable source installation
**inside the image**; this is not a source checkout or dependency
installation on the user's host.

## Independent ZIP/image runtime acceptance

**This is a new gate, not yet run on a freshly built image.** It creates
the same standalone release ZIP used for distribution, extracts it
into an isolated test workspace, then mounts **only that project**
into the container. No source-tree volume is mounted. The exact
`org.opencontainers.image.revision` label must match the Git SHA,
the installed `scrimp` distribution version must equal `1.3.0`,
the ZIP must pin the tested image, the real oscillator must save a
host-persistent and physically consistent CSV, and two-rank MPI
startup must succeed.

```bash
SCRIMP_DOCKER_NO_SUDO=1 \
SCRIMP_V130_DISTRIBUTION_IMAGE="$IMAGE" \
SCRIMP_V130_EXPECT_REVISION="$REV" \
bash scripts/run-v130-docker-distribution-gate.sh
```

The script writes a new, non-destructive
`sandbox/outputs/v130-docker-distribution-*/` directory containing
`gate.log`, `summary.json`, the tested ZIP and the extracted
project/host CSV. A serial oscillator plus `mpi_check.py` is an
**installation and MPI startup gate**, not an MPI PDE physics audit.
Native GetFEM/PETSc/MUMPS physics remains separately covered by
`scripts/smoke-v130-docker-candidate.sh` and existing PDE tests.

## Wheel and fresh-machine checks (separate required gates)

Build a wheel from clean release sources, check its `METADATA`
name/version, payload and installation into an isolated Python venv
(without importing missing native libraries):

```bash
python3 scripts/audit-v130-wheel.py
```

The audit keeps the wheel and `summary.json` in a new
`sandbox/outputs/v130-wheel-audit-*/` directory. It proves packaging,
not that a host with no GetFEM can run SCRIMP. Run the
tests inside the **candidate image**, rather than treating a
host-side Conda installation as proof of Docker portability.

On an independent fresh Linux/Ubuntu Docker environment (no SCRIMP
source checkout and no cached image), pull the **tested, immutable
image digest** only after publishing an **explicitly identified
pre-release** image. Use the standalone ZIP, verify UID/GID
ownership, persistent CSV, GUI VS Code F5/F10 and MPI startup and
at least one bounded actual MPI PDE physical test. Record
platform, `docker image inspect` ID/digest, image Git revision,
PETSc/SLEPc/GetFEM versions, logs and hashes. Do not infer AMD64
compatibility for Apple Silicon without an explicit platform
decision and testing.

## Publication policy

1. Reconcile `setup.py` and Sphinx versions to v1.3.0 (prepared in
   this Docker-first branch; **not yet a published package**).
2. Require successful source tests, strict HTML/LaTeX/PDF,
   candidate build, native smoke, independent ZIP/image gate,
   wheel audit, and real fresh-install test on the intended Linux
   architecture. Keep previous numerical acceptance proofs scoped.
3. Once all gates pass, explicitly publish the verified Docker image
   to GHCR, check the digest by pulling it from a separate machine,
   then build `scrimp-project-template.zip` with
   `--image ghcr.io/g-haine/scrimp:1.3.0` **only after that reference
   actually exists and is qualified**. A digest-pinned image may
   additionally be used for exact reproducibility.
4. Update release notes and documentation to remove
   pre-release/pending caveats only where newly proven. Review open
   issues before tagging `v1.3.0`.

Never change `examples/published/`, silently switch mathematical
solvers, assume that Docker startup proves MPI physics, or treat
Windows/macOS as native qualified platforms.
