# Changelog

Changes to SCRIMP are recorded here. See also the
[GitHub releases](https://github.com/g-haine/scrimp/releases).

## 1.3.0 — in development (not released)

### Performance and MPI
- Use exact `on_change` MUMPS factor reuse by default. Reuse the
  symbolic factor when the assembled PETSc CSR graph is unchanged and
  the numerical factor when the stored coefficients are identical on
  every MPI rank. Changes to the time-scheme shift or nonlinear
  tangents invalidate numerical reuse automatically.
- Prefer MUMPS parallel analysis `ICNTL(28)=2` for multi-rank LU/MUMPS
  runs, respecting user-provided PETSc options.
- Preserve the legacy `never`/`always` overrides and the physical
  model/discretization. No unconditional Jacobian freeze is introduced.
- Maintain the previously developed distributed GetFEM-to-PETSc COO
  bridge and compact-matrix-memory path.
- Add a release-to-current wave benchmark with native v1.2.0 serial
  reference, staged current PETSc profiling, complete state-vector
  comparisons and 1/2/4/6-rank strong scaling.
- Record the deferred MPI scalability, nonlinear-physics and
  experimental Schur results in `docs/V130_DEVELOPMENT_ROADMAP.md`.

### Release qualification outstanding
- Confirm one numerical LU (rather than ten) on the actual production
  `on_change` CN wave run; qualify invalidation when coefficients change.
- Execute the v1.2.0-vs-current serial and current MPI benchmark locally,
  validate complete states and discrete Hamiltonian/port-power balance,
  rerun the full test suite and build Sphinx HTML/LaTeX.
- Only then bump packaged/documentation versions and create a release tag.

## 1.2.0 — 2026-11-12

### Added

- **Global scalar (lumped) ports** through `DPHS.add_lumped_port()`:
  finite-dimensional state, co-state and algebraic variables with exactly
  one degree of freedom, independent of mesh resolution.
- Support for assembling coupled distributed PDE and global ODE weak-form
  contributions, with GetFEM variables synchronized to PETSc Newton
  trial iterates.
- Two new runnable, documented port-Hamiltonian ODE examples:
  `examples/mass_spring_damper.py` and
  `examples/nonlinear_pendulum.py`.
- `docs/source/lumped.rst` and worked Sphinx pages explaining the
  formulation, implementation and validation of the new examples.
- An **experimental** viscous shallow-water / moving-tank PDE–ODE
  prototype in `sandbox/sloshing_tank.py`, plus convergence diagnostics
  in `sandbox/tank_convergence.py`.

### Changed

- Corrected the nonlinear implicit residual assembly path for lumped
  systems having **linear time-derivative mass matrices**.
- Standardized the SCRIMP copyright/license header in tracked Python
  sources and the GUI code-generation template.
- Updated the test suite to use the currently supported public API,
  rather than removed symbolic/schema modules.
- Updated the package and documentation version to **1.2.0**.
- Adopted the standard setuptools PEP 517 build backend with
  `pyproject.toml`, retaining `setup.py` for package metadata.
- Used the root project README for the package's long description.

### Validation

- Python syntax check passed, excluding generated `outputs/` data.
- Full unit-test discovery: **57 tests passed**.
- Full default runs of both new finite-dimensional examples passed.
- Existing `examples/wave.py` completed: **501 steps**, no rejected
  steps and no nonlinear solver failures; Hamiltonian and boundary
  powers computed.
- Sphinx HTML and LaTeX built with `-W` and no warnings;
  the documentation PDF contains **100 pages**.
- The distribution wheel was built successfully. After the final
  README metadata correction, rebuild the wheel once more before
  tagging.

### Known limitations

- The sloshing tank remains **experimental**: wall impermeability and
  convergence of spatial fields have not yet been fully validated.
- The corrected nonlinear residual path does **not** address general
  state-dependent nonlinear mass matrices.
- PETSc may report an unused `-isset false` option; some real-valued
  runs emit complex-to-real conversion warnings.
- In the historical wave example, the last monitor entry was
  `t = 5.01` for `t_f = 5.0`; this endpoint discrepancy is not
  addressed in this release.
- Standalone examples have passed runtime checks; quantitative
  time-refinement and energy-balance studies for both new ODE
  examples remain future validation work.

## 1.1.1 — 2024-03-21

- Changed the output-directory behavior to use the current working
  directory.
- Corrected the sign of a source term and improved documentation.
- Added notebooks.

See the [v1.1.1 release notes](https://github.com/g-haine/scrimp/releases/tag/v1.1.1).
