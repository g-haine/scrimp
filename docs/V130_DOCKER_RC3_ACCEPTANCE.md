# SCRIMP v1.3.0 — Docker rc3 targeted qualification (Ubuntu)

**Date:** 20 November 2026  
**Candidate:** `scrimp-mpi:v1.3.0-rc3`  
**Exact runtime revision:** `1d6d0c90115213298047b28237a7f38087563a35`  
**Status:** **PASS — targeted Ubuntu native/BuildKit and PETSc diagnostic qualification complete.**
Original rc3 smoke covered the native/MPI stack, host-persistent oscillator and
PETSc `-info` correction. The corrected host-side two-rank import diagnostic
was rerun against the **existing** rc3 image and is **warning-free**. This does
not imply full scientific or clean-room/bitwise release acceptance.

This report is based on the Ubuntu developer's complete terminal
output of the rc3 build and smoke. It is not evidence of a fresh
machine/locked-environment build or an independent full MPI physics run.

## Provenance, isolation, BuildKit reuse

- The build fetched GitHub at exactly
  `1d6d0c90115213298047b28237a7f38087563a35`, independently of
  the host's modified `sandbox/sloshing_tank.py`.
- Audited temporary context: **191 tracked source files / 1.54 MiB**,
  BuildKit transfer **1.62 MB**. Local simulations and
  `sandbox/outputs` were excluded.
- Existing native dependency layers were **CACHED**: APT, Conda,
  complex CPU PETSc/SLEPc bindings, GetFEM 5.5 and its five patches.
  Rebuilt only source-copy, SCRIMP installation, quick ABI checks and
  final image metadata.
- Candidate manifest-list digest as exported:
  `sha256:b7e15a628b538bc313117d008a56a2d6d41b5e9774330a584e920fc61dd92c3d`.
  Runtime image manifest:
  `sha256:0ab4d372455a4dd9c11870d6b370c033605f4c7c83602bb27c77413e55a9000c`.
- The developer compared `scrimp-mpi:local` ID before/after with
  `cmp`; there is no difference/error in the complete shell
  transcript. No image was published to GHCR.

## Native, numerical and integration smoke

| Gate | Observed |
| --- | --- |
| GetFEM 5.5 MPI parallelism level 2 | PASS |
| Complex PETSc/SLEPc 3.26.1, NumPy 2.5.5 and MPI bindings | PASS |
| Approved local RT/BDM/Nédélec-I/Nédélec-II DoF constructors | PASS |
| Nédélec-II k>=3 unsupported guard | PASS |
| Distributed MUMPS LU factor/solve, 1 and 2 ranks | PASS |
| Independent host-mounted user project | PASS |
| Real driven/damped oscillator, 10 accepted timesteps, zero failures | PASS |
| 11 samples `t=0..0.1` in host CSV, Hamiltonian formula | PASS |
| `np.exceptions.ComplexWarning` promoted to error | PASS |
| Independent project two-rank MPI startup and SCRIMP import | PASS |
| `PETSc.Options().hasName("info")` under two ranks | **False on both** |

The oscillator CSV checks the *algebraic* Hamiltonian
`H = 5*q*q/2 + p*p/(2*1.5)`, not the full port-Hamiltonian
power balance. Nédélec constructor counts do not qualify Maxwell
3D conformity or orientation.

## Final `-c` warning comes from the test harness, not SCRIMP

The original rc3 smoke used `mpiexec -n 2 python -c '...'` to check
that SCRIMP no longer injected `-info`. The check printed
`PASS: PETSc -info not forced by SCRIMP import` on **both** MPI ranks.
PETSc then reported

~~~text
Option left: name:-c (no value) source: command line
~~~

because `petsc4py.init(sys.argv)` also saw Python's `-c` token.
No code-based `-info` warning remained. This is a
**diagnostic harness artifact**, not a numerical failure.

The host-side `scripts/smoke-v130-docker-candidate.sh` was
subsequently fixed to create a short, temporary **real .py file**
in the independent test workspace, execute it at two ranks from a
runtime bind mount, check `-info` and `-c` are absent, and
**fail on any PETSc unused-option warning**. The file never enters
a Docker build context; no new image or native compilation is
required. Preserve the existing rc3 image and run only:

~~~bash
cd ~/Simulations/scrimp
git switch release/1.3.0-preparation
git pull --ff-only
SCRIMP_DOCKER_NO_SUDO=1 \
SCRIMP_DOCKER_CANDIDATE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
  bash scripts/smoke-v130-docker-candidate.sh
~~~

Omit `SCRIMP_DOCKER_NO_SUDO=1` if the host uses `sudo docker`.
Capture both `sandbox/outputs/v130-docker-smoke.log` and
`sandbox/outputs/v130-petsc-options-smoke.log` before any further
qualification. **Follow-up host evidence (20 November 2026) — PASS:** the developer
supplied the actual diagnostic output from the corrected test:

~~~text
=== PETSc optional diagnostics: no implicit -info under two MPI ranks ===
PASS: PETSc import options clean (rank 1/2)
PASS: PETSc import options clean (rank 2/2)
PASS: PETSc -info and -c absent, no unused-option warning.
~~~

This closes the final **PETSc option diagnostic** gate with both MPI
ranks, independently of the prior rc3 native smoke. The supplied
excerpt confirms the targeted diagnostic rather than reproducing
the entire second smoke log. The underlying rc3 image was **not
rebuilt**, because the correction only changed the host-side
test harness.

## Follow-up gates

- The existing `setup.py` still installs `scrimp-1.2.0` from a
  v1.3.0 candidate; defer version metadata freeze to packaging.
- Base-image digest/platform and Conda explicit lock must still be
  committed and clean-room-tested before claiming full
  reproducibility. Warm-cache reuse is already proven by rc2.
- Full source-to-source wave parity, initial DAE consistency, full
  discrete port-power identities, Maxwell 3D, memory and performance
  remain separate scientific/release qualifications.
- The local tracked edit `sandbox/sloshing_tank.py` was not built
  into rc3 and must not be overwritten by unrelated checkout cleanup.

**No GHCR publication, no mutation of `master`, and no promotion over
`scrimp-mpi:local`.**
