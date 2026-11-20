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

### Consistent DAE initialization (new default candidate; native gate pending)

- Add a **distributed PETSc algebraic-only initialization** method for
  supported real, linear DAE models. It keeps the prescribed differential
  initial coefficients unchanged, solves only for algebraic efforts/ports
  with bounded LSQR refinements, and verifies the *actual*
  `IFunction` algebraic residual before starting the physical TS.
  Unsupported layouts or incompatible fixed initial data raise
  explicit errors; there is **no automatic pseudo-time fallback**.
- **Stage `init_step_ts_type="algebraic"` as the implicit default**
  on the development branch. This change itself has **not yet passed
  the independent native no-option migration gate and is not released**.
  `init_step_ts_type="pseudo"` and `init_step_ts_type="bdf"`
  remain available explicitly; `init_step=False` remains valid
  without an initialization strategy.
- Document the distinction between **exact algebraic projection**,
  historical **pseudo-time relaxation** (which can move the
  prescribed differential initial data), and **small negative-start
  BDF2** (which also changes them, even with a very small step).
  See `docs/source/initialization.rst`.
- Preserve the historical pseudo behavior explicitly in the
  original top-level examples, legacy notebooks, published Maxwell
  example, GUI-generated legacy scripts and selected sandbox workloads.
  Their physical initial data and equations are not modified.
- **Native core opt-in qualification PASSED** on two independent
  port-Hamiltonian linear acoustic models, conservative and dissipative,
  on MPI 1/2 with ten fixed CN steps. Native algebraic residuals
  were about `1e-13`, initial dynamic drift zero, and full
  signed Hamiltonian + nonzero boundary work (+ positive internal
  damping on the dissipative model) defects about `1e-15`.
  The incompatible original Gaussian was explicitly rejected on
  both MPI ranks. The historical pseudo trajectory matched its
  frozen reference; tiny negative-time BDF2 remained usable with
  measured nonzero initial drift.
- **Remaining release gates:** independently test the implicit
  (no-option) new default on both models using the frozen rc3 runtime,
  verify explicit pseudo/BDF2 regressions, complete broader example
  compatibility and exceptional MPI cleanup/memory audits, and
  rerun the full native source/test/doc/release qualification.
  Do not tag or rebuild the release image before these pass.

### Release qualification outstanding
- Confirm native numeric and symbolic LU call counts agree with the
  exact operator coefficient/graph changes detected collectively by
  `on_change`; qualify invalidation when coefficients change. The
  final fixed-CN MATCHSTEP may legitimately change the shift by
  floating-point roundoff (2 numeric LU for 10 steps observed).
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
