# SCRIMP v1.3.0 — Ubuntu Docker candidate rc2, native/cache acceptance

**Date:** 2026-11-20  
**Branch:** `release/1.3.0-preparation`  
**Image:** `scrimp-mpi:v1.3.0-rc2`  
**Exact image source:** `ae96cb7334d7c88cd2e0dbd4d9615af66f1f4a61`  
**Status:** **PASS — targeted native/runtime and BuildKit warm-cache qualification**, not general v1.3.0 scientific release acceptance.

This acceptance reflects **two real Docker Buildx builds and the native/user-project
smoke performed by the developer on Ubuntu**, as transcribed in the
20 November developer terminal report. No GitHub Actions remote Docker build
or other host is claimed.

## Observed image provenance and BuildKit behavior

- The supported builder fetched the exact GitHub commit above (not the
  working tree) and physically staged **190 allowlisted tracked files,
  1.54 MiB**. BuildKit transferred **1.62 MB**. The log explicitly
  excludes local simulations, `sandbox/outputs`, graphics and untracked data.
- **First build:** Conda install, petsc4py/slepc4py binding build, and patched
  GetFEM 5.5 compilation were all `CACHED`; only the application source
  COPY layers and editable install/smoke were re-executed.
- **Second build at the identical SHA:** every Dockerfile build stage,
  including the final editable SCRIMP install/smoke layer, reported `CACHED`.
  This validates warm-cache reuse, **not a clean-machine rebuild**.
- GetFEM uses the reviewed Savannah 5.5 source archive and the five
  versioned patches. Native GetFEM MPI parallelism level **2** and
  CPU complex PETSc/SLEPc **3.26.1**, NumPy **2.5.5**,
  SciPy **1.18.1**, mpi4py **4.1.2** are confirmed in smoke output.
- `scrimp-mpi:local` was never named as a build target, nor was an image
  pushed. **The supplied transcript does not include actual before/after
  `scrimp-mpi:local` IDs**, so tag-preservation by configuration is confirmed,
  but the requested independent ID comparison remains unrecorded.

### Why the two reported Docker IDs differ despite all layers being cached

| OCI component | First build | Warm build |
| --- | --- | --- |
| Application image manifest | `sha256:82a74c8e35fa49a1f447a97184c666745e8ead4d96390708c13522ac167d59a8` | **same** |
| Image configuration | `sha256:08d06d2b7317e4a2f10ddfd6d3fe08e4ae35714e47fda7bed62d7fd5d7368fce` | **same** |
| BuildKit attestation manifest | `sha256:cbb45b917244840e52d42e739aeb6813d7c19e8391ba92d4bb43aae2ad1a943a` | `sha256:7c60d30001e95dfdc57bbe2418ad32d3021d03279d7972f4b39e9e25e43643ef` |
| Reported manifest-list ID | `sha256:5ae22339eaa9483cf3ac1bf16012bc1fe15281e1e7d34e2290360843ba32e866` | `sha256:4eb62339124d78cc72baaf8baa695aea7b3552b30d5eb93cf87ee48c847328b9` |

The logs therefore support **identical main runtime manifest/configuration
and different BuildKit attestation data**. They do **not** establish a
bit-for-bit reproducible OCI manifest-list digest. Do not use
the manifest-list ID alone as a scientific solver regression metric.
Investigate attestation provenance/determinism separately when
preparing public immutable releases; do not modify the build
architecture prematurely.

## Native/runtime acceptance matrix

| Gate | Observed result |
| --- | --- |
| GetFEM 5.5 MPI-2 | PASS |
| Selected RT, BDM, Nédélec-I/II triangle/tetrahedron local DoF counts | PASS, restricted release orders |
| Explicitly unsupported Nédélec-II k>=3 | Correctly rejected |
| PETSc complex scalar and GetFEM/NumPy ABI | PASS |
| Distributed PETSc/MUMPS solves | PASS at 1/2 ranks |
| Independent user-project generation | PASS |
| 0D driven/damped oscillator under NumPy `ComplexWarning` promoted to error | 10 accepted steps, 0 rejections, 0 nonlinear failures; **PASS** |
| Host-persistent `outputs/oscillator.csv` | 11 finite samples over t=0..0.1; **PASS** |
| Algebraic CSV Hamiltonian `H=5*q²/2+p²/(2*1.5)` | PASS at rel=1e-10, abs=1e-12 |
| Independent user project MPI startup | PASS at 2 ranks |

The CSV Hamiltonian test checks algebraic export consistency, **not**
the discrete port-Hamiltonian power balance. The two-rank starter
checks imports and MPI startup, **not** distributed PDE equivalence.
The accepted Ubuntu VS Code F5/F10/interpreter/Pylance workflow belongs
to the separately documented **rc1** GUI acceptance; no fresh rc2 GUI
session was performed.

## Diagnostics found and subsequent source changes

**PETSc `-info` warning persists in rc2's two-rank starter:**

```text
WARNING! There are options you set that were not used!
Option left: name:-info (no value) source: code
```

Source inspection found the direct cause:
`scrimp/utils/linalg.py` unconditionally executed
`PETSc.Options().setValue('info', None)` during import.
This is not required for solver execution. A **post-rc2 source
change** removes the implicit diagnostic; a source regression and
a two-rank Docker smoke assertion now check that the option is not
forced. **This fix is not present in the accepted rc2 image.**
A separately tagged `rc3` built from a new immutable commit must
confirm the fix in native runtime before it is marked PASS.

**Packaging version:** the rc2 build pip output reports
`scrimp-1.2.0` (the still-current version in `setup.py`).
This is expected for a pre-release integration image but is a
**release metadata gate**: pin an appropriate v1.3.0 prerelease
or final version during packaging freeze, without silently treating
`rc2` as a published Python package version.

**Local source tree:** the developer's log showed a tracked
modification to `sandbox/sloshing_tank.py`; the isolated GitHub
checkout means this local change was not part of the image and
must not be discarded or staged accidentally.

## Remaining release gates

1. Run a new, precisely pinned **rc3** image smoke to qualify the
   source-only PETSc option fix; preserve rc1/rc2 and
   `scrimp-mpi:local`.
2. Verify and retain a before/after `scrimp-mpi:local` image ID
   comparison in the next run.
3. Freeze the base-image digest/platform and explicit tested
   Conda environment; audit BuildKit attestation repeatability.
4. Keep full wave MPI trajectory parity, DAE initial consistency,
   discrete Hamiltonian/port-power balances, Maxwell 3D, long-time
   memory, packaging and cross-platform qualification as **distinct**
   scientific/release acceptance gates.

**Conclusion:** Ubuntu Docker **native stack + BuildKit reuse +
independent starter candidate rc2** PASS. The complete SCRIMP
v1.3.0 release, the post-rc2 PETSc change and bitwise OCI list
reproducibility remain **unqualified**.
