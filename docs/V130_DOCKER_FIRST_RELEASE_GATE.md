# SCRIMP v1.3.0 — Docker-first distribution and release gate

**Decision, 27 November 2026:** precompiled Docker is the **default
user installation** for v1.3.0. Locally built Docker remains the
developer/reproducibility alternative; native Conda remains an
advanced/legacy alternative. A final user needs only Docker, a
version-matched standalone project ZIP and optionally VS Code with Dev
Containers. Never require cloning the SCRIMP library into a simulation
workspace.

## Candidate qualification and remaining stable-release work

- **Qualified earlier, with bounded scope:** local Ubuntu image
  `scrimp-mpi:v1.3.0-rc3`, GetFEM 5.5, PETSc/SLEPc complex,
  MUMPS, MPI 1/2/4/6 where explicitly tested, scientific model
  checks and the independent Dev Container preview.
- **Canonical examples:** full nine-case user-output campaign at
  `20d79d0d` passed with no generic output siblings.
- **Documentation:** the initial PDF qualification at
  `cf75f4f761c2854788e5ef519375b31e7f040509` was superseded
  by the actual version-`1.3.0` build at
  `e81003b7ef32360c755a4db56b65d4093b744aea`: 8/8
  source tests, strict Sphinx HTML/LaTeX with **zero Sphinx
  warnings** and genuine `scrimp.pdf` (**168 pages**).
  The final pdfLaTeX log has no unresolved references but
  still contains three `hyperref` PDF-string warnings and
  small layout box warnings. These are not Sphinx warnings.
- **Now qualified, bounded to Ubuntu/amd64:** built and inspected
  `scrimp-mpi:v1.3.0-dist-rc2` (revision `91d7eb1a`); real native
  stack/MUMPS MPI1/2 smoke, standalone ten-file project ZIP with
  source-independent simulation and host-persistent CSV, wheel metadata
  audit, and the user-confirmed VS Code Dev Container F5/F10/F5
  debugging path. The Python source CI became green on the subsequent
  release-preparation branch, after standardizing 89 Python headers.
- **Now PASS on a genuinely new Ubuntu 24.04/AMD64 runner:**
  GHCR prerelease `ghcr.io/g-haine/scrimp:1.3.0-dist-rc3` was
  published and pulled **by OCI index digest**
  `sha256:0cbbce6ad42b23552aff9d12555342f3864b50522054fe039c161a7145b2732a`;
  the actual pulled local image ID was
  `sha256:c8307d382ed88f449d86bd09336ad75908a1fc58ba49d80baae2719c0d5de798`.
  Image source label = `58dbe187de308979938153f9adfcf94dac4df85a`.
  Independent ZIP, oscillator/host CSV, MPI2 startup and 10-step
  Heat--Wave MPI1/MPI2 native Hamiltonian parity **PASS**
  (relative L2 `4.275706381199414e-15`, tolerance `1e-7`).
  [GitHub Actions run 46025625273, attempt 2](
  https://github.com/g-haine/scrimp/actions/runs/46025625273) archives the
  actual logs and evidence. **This is not full-field equivalence,
  graphical VS Code testing on a new desktop, or Windows/macOS proof.**
- **Current docs/build progress:** v1.3.0 wheel metadata/build and
  isolated install PASS, Sphinx HTML/LaTeX/PDF PASS on
  `0cdb2ad8`. A separate source-distribution audit
  now builds the `sdist` from tracked files and reconstructs a
  byte-equivalent Python wheel; it has passed its CI step.
- **Still pending:** stable `ghcr.io/g-haine/scrimp:1.3.0` tag,
  final release project ZIP and GitHub tag/release, issue triage,
  and a deliberate stable publication decision. Windows/macOS
  remain best-effort.

## Ubuntu native and independent ZIP acceptance (2026-11-27)

**PASS, bounded to the exact local x86_64/Ubuntu environment**:
candidate image `scrimp-mpi:v1.3.0-dist-rc2` was built from
`91d7eb1a12329d54854f13cd60d7860c36b856a4`; OCI image ID
`sha256:94b5c725c213ee3d3737d26ad190b2c103deec1b9b573b3dca340939f7a6a270`
and its revision label both matched. The native smoke passed
GetFEM descriptors and MUMPS/PETSc/SLEPc imports under MPI1/2,
a genuine lumped oscillator (10 accepted BDF2 steps, no solver
failures), 11 persistent CSV samples and exact analytic
`H(q,p)=2.5*q²+p²/3` to floating-point round-off, two-rank MPI
startup, and import-only PETSc options diagnostics.

The independent distribution gate
`PASS_STANDALONE_ZIP_IMAGE_SERIAL_MPI_STARTUP` then PASSed on
the same image, extracting the tested ten-file
`scrimp-project-template.zip` and running it in a separate
user-only bind mount, with physically verified host-persistent
`outputs/oscillator.csv`. Evidence:
`sandbox/outputs/v130-docker-distribution-YQpvlaPX/`.
The package wheel `scrimp-1.3.0-py3-none-any.whl` was also
built, payload-audited and installed into an isolated metadata
environment, without claiming native scientific dependencies
are available outside Docker. The qualified *distribution*
image is not yet public or pushed to GHCR.

**PASS — Ubuntu VS Code GUI and debugpy using the independent ZIP:**
the user extracted the qualified distribution ZIP into a separate
project, opened the project in VS Code Dev Containers (a VS Code
close/reopen was needed once), and executed the original starter
simulation. Its host-persistent `oscillator.csv` contains 11 samples
at `t=0,0.01,...,0.10` with the analytical energy identity to
about `5.6e-17`, `H(0)=0.4`, `H(0.1)=0.39928159`.
The user **explicitly confirmed** the debugging sequence on
2026-11-27: **F5** stopped at a breakpoint, **F10** advanced
instruction by instruction, and **F5** resumed the simulation
through successful completion. This qualifies the tested local
Ubuntu Dev Container and independent ZIP at image `dist-rc2`;
it does **not** qualify an image freshly pulled from a registry,
Windows/macOS, or multi-process MPI debugging.

**Now PASS (additional evidence received 2026-11-27):**
the original image-internal Heat–Wave PDE, using GetFEM and MUMPS
on 3,154 mesh elements / 7,689 DoFs, ran 10 accepted CN steps
to `t=0.02` with 1 and 2 MPI ranks. Both computed 11
genuine FE Hamiltonians; the full trajectory relative
`L2` difference was `4.252240939683056e-15`,
versus the independently specified `1e-7` acceptance tolerance.
The gate explicitly mounted only its independent probe directory,
no SCRIMP sources, and returned
`PASS_NATIVE_HEAT_WAVE_MPI1_MPI2_HAMILTONIAN_PARITY`.
This is **energy observable parity**, not full-state FE equivalence
or a claim for other MPI sizes. Evidence:
`sandbox/outputs/v130-distribution-heat-wave-mpi-wb9SrBRd/`.

**Historical dist-rc2 status at the time of that earlier audit:**
a fresh GHCR pull was still pending. That gate was subsequently
completed on `dist-rc3` (see the updated evidence above). GUI
F5/F10 remains proven only with the local Ubuntu Dev Container;
the stable GHCR tag, final ZIP, release approval and GitHub
release are still outstanding.

## Final release-candidate acceptance without hard-coded provenance

After source freeze, the workflow
`.github/workflows/v130-ghcr-clean-install.yml` is triggered by the
**release-candidate Git tag** `v1.3.0-rcN`. It derives the exact
checkout SHA from that immutable tag, resolves the corresponding
`ghcr.io/g-haine/scrimp:1.3.0-dist-rcN` **registry OCI index digest**,
then pulls **by digest**, checks the embedded OCI source revision and
runs the standalone project ZIP/host-owned CSV/MPI and native Heat--Wave
MPI1/MPI2 physics gates on a fresh Ubuntu 24.04/AMD64 runner.

**Order matters:** publish and verify the candidate **image first**,
and only then push its **matching Git RC tag**. Pushing a tag before the
image is accessible makes the CI correctly fail before runtime testing.
The workflow is also manually dispatchable with an existing candidate
tag once it is present on the GitHub default branch. Never modify
the repository just to hard-code a new digest.

For the next *separately named* candidate, from an already checked-out
and finalized source commit (example: `rc4`):

```bash
set -euo pipefail
test -z "$(git status --porcelain --untracked-files=no)"
REV="$(git rev-parse HEAD)"
RC=4
IMAGE="scrimp-mpi:v1.3.0-dist-rc${RC}"
REMOTE="ghcr.io/g-haine/scrimp:1.3.0-dist-rc${RC}"

SCRIMP_DOCKER_NO_SUDO=1 \
SCRIMP_DOCKER_CANDIDATE_IMAGE="$IMAGE" \
SCRIMP_DOCKER_SOURCE=github \
SCRIMP_GITHUB_REF="$REV" \
bash scripts/build-v130-docker-candidate.sh

SCRIMP_DOCKER_NO_SUDO=1 \
SCRIMP_DOCKER_CANDIDATE_IMAGE="$IMAGE" \
bash scripts/smoke-v130-docker-candidate.sh

test "$(docker image inspect --format \
  '{{index .Config.Labels "org.opencontainers.image.revision"}}' \
  "$IMAGE")" = "$REV"

# Requires docker login ghcr.io with write:packages, or other
# authorized GHCR publication identity. Does NOT publish stable 1.3.0.
docker tag "$IMAGE" "$REMOTE"
docker push "$REMOTE"
docker buildx imagetools inspect "$REMOTE"

# Do NOT reuse or move an already published candidate Git tag.
git tag -a "v1.3.0-rc${RC}" "$REV" \
  -m "SCRIMP v1.3.0 release candidate ${RC}"
git push origin "v1.3.0-rc${RC}"
```

The final command starts the hosted acceptance without creating a
new source commit. Store the **registry digest** and the actual
`docker image inspect .Id` separately; these are distinct identifiers.
If the tag/image revision check or any numerical gate fails,
**do not** re-tag the image as the stable release. Correct the defect,
freeze a new commit and prepare a **new RC number** instead.

This acceptance is limited to the recorded Ubuntu/AMD64 runtime,
independent user project, host outputs and MPI1/MPI2 physics. It does
not test the SCRIMP graphical GUI, Windows, macOS, higher MPI sizes,
or arbitrary nonlinear DAE initialization.

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

**PASS on the local dist-rc2 Ubuntu/amd64 image (2026-11-27).** It creates
the same standalone release ZIP intended for distribution, extracts it
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

## Original-image MPI PDE parity (independent of source checkout)

A new bounded gate `scripts/run-v130-distribution-heat-wave-mpi-parity.sh`
executes the original, image-internal
`examples.heat_wave.heat_wave_eq(physics_only=True)` on one and
two ranks of **the same tagged Docker image**, with identical genuine
GetFEM weak forms, mesh, physical initialization, MUMPS and fixed
Crank--Nicolson steps. The *only* host mount is an isolated directory
containing the MPI probe and its output files: there is no SCRIMP
source or canonical-example bind mount. Rank zero writes the 11
native Hamiltonian samples (0 to 0.02 by 0.002), and a standard-library
checker compares the complete Hamiltonian trajectories to
`relative_l2 <= 1e-7`. A PASS proves only this bounded physical
observable's MPI1/MPI2 parity, **not** fieldwise or universal MPI
equivalence.

**PASS, native execution on Ubuntu/amd64, 2026-11-27, image
`dist-rc2`.** Full gate result
`PASS_NATIVE_HEAT_WAVE_MPI1_MPI2_HAMILTONIAN_PARITY`:
MPI1/MPI2, 10 accepted steps, 11 Hamiltonian samples,
relative L2 difference `4.252240939683056e-15` versus
`1e-7` tolerance. From this feature worktree, reproduce
the qualification using the image's **actual** OCI revision,
not the moving worktree HEAD:

```bash
SCRIMP_DOCKER_NO_SUDO=1 \
SCRIMP_V130_DISTRIBUTION_IMAGE=scrimp-mpi:v1.3.0-dist-rc2 \
SCRIMP_V130_EXPECT_REVISION=91d7eb1a12329d54854f13cd60d7860c36b856a4 \
bash scripts/run-v130-distribution-heat-wave-mpi-parity.sh
```

It creates non-destructive
`sandbox/outputs/v130-distribution-heat-wave-mpi-*/` evidence with
rank-specific `series.json`, logs, and a fail-closed `summary.json`.

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
