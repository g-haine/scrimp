# SCRIMP v1.3.0 — release highlights (draft)

> **Pre-publication draft, 27 November 2026.** This document is a short,
> user-facing release summary, **not** a release announcement. The stable
> GitHub tag, GHCR `1.3.0` image and official project ZIP are **not yet
> published**. Publish these notes only after the stable release gates pass.
> Detailed change history and scoped scientific evidence remain in
> [CHANGELOG.md](../CHANGELOG.md) and
> [v1.3.0 release readiness](V130_RELEASE_READINESS_2026-11-27.md).

## TL;DR

**Faster. More reliable in MPI. More nonlinear models. Easier to install.**

- **Performance:** distributed GetFEM-to-PETSc sparse COO assembly, lower
  memory overhead and selective MUMPS symbolic/numeric factor reuse
  (`on_change`). A specifically measured 2-D acoustic configuration
  showed **up to ~48× end-to-end speedup vs v1.2.0**; this is a
  configuration-specific benchmark, **not** a general speedup guarantee.
- **MPI correctness:** corrected matrix/RHS assembly and distributed
  solve conventions, collective power computations and safer shared
  ParaView outputs. Native 1/2/4/6-rank campaigns and scoped full-field
  and Hamiltonian checks establish reproducibility in the tested models.
- **GetFEM 5.5 source patches:** five SCRIMP-maintained patches fix MPI
  assembly verbosity, tetrahedral BDM degrees of freedom and Python
  extension linking, and provide first-/second-kind simplicial
  Nédélec elements. These are **SCRIMP-side additions**, not claims
  that the patches have been merged upstream. Nédélec-II is supported
  only at polynomial degrees **k=1 and k=2**.
- **Nonlinear and DAE support:** correct state-dependent mass Jacobians,
  fully explicit nonlinear weak-form terms evaluated at the last
  accepted state, and consistent **linear** algebraic DAE initialization
  without changing prescribed differential initial data.
- **Physical observables:** safer accepted-time sampling, more reliable
  Hamiltonian and port-power diagnostics and explicit distinction
  between genuine every-step CN quadrature and approximate nodal
  integration for decimated CN/BDF/BDF2/ARKIMEX.
- **Examples:** canonical examples produce genuine GetFEM fields,
  plots, archives and, where appropriate, animations themselves,
  under `outputs/<system_name>/`. A clean-layout campaign checked
  eight serial examples and a ninth Heat–Wave MPI2 case.
- **Docker first:** the recommended user path is a prebuilt Linux/AMD64
  scientific environment (GetFEM 5.5, complex PETSc/SLEPc, MUMPS, MPI,
  Gmsh) plus a **standalone, host-owned project ZIP**. VS Code Dev
  Containers, standard Python editing/debugging and persistent local
  results are supported on the qualified Ubuntu setup.
  Native Conda and local Docker builds remain alternatives.

## What changed for existing users?

Most public modeling concepts and original published examples are kept
intact; the eight historical `examples/published/` source files were
not rewritten. There are, however, **observable default changes**:

- `init_step_ts_type="algebraic"` is the default for supported **linear**
  DAEs. Incompatible initial values or unsupported layouts now report
  errors rather than silently changing prescribed differential data.
  Historical `pseudo` and `bdf` strategies remain explicitly
  available, as does `init_step=False`.
- Explicit nonlinear `Brick(explicit=True)` expressions are now
  consistently held at the **last accepted time/state** during
  implicit Newton iterations. This is a deliberate lag, not midpoint
  or fully implicit evaluation.
- MUMPS `SCRIMP_MUMPS_RESET_POLICY=on_change` becomes the default.
  `always` and `never` remain configurable.
- Automatic conforming interface detection and rank-safe shared
  ParaView output are enabled in the documented supported cases;
  ambiguous interfaces or incomplete distributed-only exports
  fail explicitly. Legacy opt-outs are available in the documentation.
- Canonical example results now have a consistent per-model output
  layout, with no duplicate generic output directories.

## SCRIMP-maintained patches to GetFEM 5.5

The prebuilt SCRIMP Docker environment uses upstream **GetFEM 5.5
with five explicit, version-pinned source patches**. They are applied
in a reproducible sequence by `scripts/build-getfem-mpi.sh`, and
incompatible source revisions are rejected rather than silently
building an unpatched library.

| Patch | User-visible purpose |
| --- | --- |
| `0000-mpi-assembly-timing.patch` | Make GetFEM's MPI assembly-time printing opt-in via `SCRIMP_GETFEM_ASSEMBLY_TIMING=1`, avoiding noisy default logs. |
| `0001-bdm-tetrahedron-unisolvent-dofs.patch` | Correct face and interior degrees of freedom for tetrahedral Brezzi–Douglas–Marini elements; targeted native checks cover **k=1,2**. |
| `0002-python-interface-link-in-tree.patch` | Make the Python `_getfem` extension link against the freshly built GetFEM archives, not an unrelated installed library. |
| `0003-nedelec-first-kind-simplices.patch` | Add `FEM_NEDELEC1(N,K)` on triangles and tetrahedra. The **k=0** descriptor delegates to GetFEM's built-in lowest-order Nédélec; targeted higher-order tests cover **k=1,2**. |
| `0004-nedelec-second-kind-simplices.patch` | Add `FEM_NEDELEC2(N,K)` on triangles and tetrahedra for **k=1,2 only**; **k≥3 is not implemented**. |

SCRIMP exposes `NED1`/`NEDELEC1` and
`NED2`/`NEDELEC2` aliases through `S.FEM`; unsupported
descriptors fail explicitly instead of silently falling back to
another element family. These patches were evaluated with targeted
finite-element constructors, unisolvency/DoF, mass-matrix and
continuity checks. **That is not an exhaustive certification of all
orders, geometries or MPI physics**, nor a claim of upstream GetFEM
acceptance. See
[GetFEM 5.5 patch inventory](../scripts/patches/getfem-5.5/README.md)
and [higher-order follow-up](V130_GETFEM_HIGHER_ORDER_TODO.md).

## Scientific qualification — what the tests actually establish

- A local Ubuntu/AMD64 Docker runtime has undergone source, FE family,
  nonlinear Jacobian, DAE, wave, heat, Saint–Venant and 3-D Maxwell
  tests within their recorded scopes. Not every model is validated on
  all meshes, solvers, MPI counts or time integrators.
- On 27 November, a **fresh GitHub-hosted Ubuntu 24.04/AMD64 runner**
  downloaded the pre-release `dist-rc3` image **by immutable GHCR
  digest**, confirmed its OCI source revision, tested the independently
  generated user ZIP and persistent real oscillator output, and passed
  a genuine Heat–Wave PDE Hamiltonian trajectory parity test in MPI1/2.
  Hamiltonian relative L2 difference: **4.28e-15** (threshold 1e-7).
  Evidence: [GitHub Actions run 46025625273, attempt 2](
  https://github.com/g-haine/scrimp/actions/runs/46025625273).

## Installation (once the stable v1.3.0 release is published)

1. Download `scrimp-project-template.zip` from the matching
   [GitHub Release](https://github.com/g-haine/scrimp/releases).
2. Pull the released Docker image:
   `docker pull ghcr.io/g-haine/scrimp:1.3.0`
3. Open the extracted project in VS Code and select
   **Dev Containers: Reopen in Container**, or run with Docker from
   a terminal. Your `simulations/` and `outputs/` remain on the host.
4. Run `python simulations/first_simulation.py`, then
   `mpiexec -n 2 python simulations/mpi_check.py`.

**Do not use the commands above to infer the stable image already
exists.** They are final-user instructions conditional on publication.

## Known limits and planned follow-ups

- The default algebraic initializer does **not** solve arbitrary
  nonlinear/higher-index DAEs; general nonlinear fixed-differential-
  state initialization is a separate research and implementation task.
- BDF/BDF2, delayed explicit nonlinear terms, decimated output and
  arbitrary spatial discretizations do not automatically guarantee
  exact discrete energy conservation.
- The qualified pre-release distribution is **Linux/AMD64**. Windows,
  macOS, Apple Silicon, multi-process GUI debugging and MPI VTK
  without full replicas are not generally qualified.
- Modernization of SCRIMP's own graphical GUI (not VS Code),
  including its handling of supported and unsupported DAE
  initialization strategies, is a **post-v1.3.0** follow-up
  (tracked in [issue #32](https://github.com/g-haine/scrimp/issues/32)).

For full version-by-version corrections, reproducibility details,
numerical acceptance gates, and past diagnostic attempts, see
[CHANGELOG.md](../CHANGELOG.md) and the
[documentation](https://g-haine.github.io/scrimp/).
