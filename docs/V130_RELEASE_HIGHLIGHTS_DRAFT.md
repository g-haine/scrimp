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
- A separate original-state viscous 1-D Saint–Venant / translating-tank
  user project exercised nonlinear mass on 647 DoFs for up to 6,250
  accepted CN steps without nonlinear failures. Its physical energy
  balance and wall treatment still require investigation; it is
  **evidence of solver robustness, not a proof of exact passivity**.

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
- The original viscous sloshing-tank model's wall and energy audit, as
  well as modernization of SCRIMP's own graphical GUI (not VS Code),
  are explicitly **post-v1.3.0** follow-ups.

For full version-by-version corrections, reproducibility details,
numerical acceptance gates, and past diagnostic attempts, see
[CHANGELOG.md](../CHANGELOG.md) and the
[documentation](https://g-haine.github.io/scrimp/).
