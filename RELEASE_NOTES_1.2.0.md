# SCRIMP v1.2.0

This release introduces **spatially global scalar (lumped) ports** for
finite-dimensional port-Hamiltonian systems and selected PDE–ODE
couplings, plus new examples, regression tests and documentation.

## Highlights

- **Global scalar ports:** `DPHS.add_lumped_port()` registers one
  GetFEM degree of freedom for a scalar state, co-state or algebraic
  variable, independent of the distributed mesh.
- **Coupled implicit systems:** GetFEM nonlinear assembly is
  synchronized with PETSc Newton trial iterates. The corrected residual
  applies to lumped systems with a **linear time-derivative mass
  matrix**, not arbitrary nonlinear mass matrices.
- **New examples:** a forced, damped mass–spring oscillator and a
  nonlinear pendulum with viscous damping and optional torque.
- **Documentation:** a guide to lumped ports, worked ODE examples,
  refreshed Sphinx pages and a 100-page reference PDF.
- **Project maintenance:** 57 passing unit tests, standardized Python
  file headers, updated test suite, version 1.2.0 metadata and
  PEP 517 packaging via `pyproject.toml`.

## Validated

- Full unit-test suite: **57 tests passed**.
- Both new examples ran successfully over their default full durations.
- Historical 2D wave example: **501 steps**, no rejected steps or
  nonlinear solver failures; Hamiltonian and boundary powers computed.
- HTML and LaTeX documentation compiled with Sphinx's strict `-W`
  setting; PDF built successfully.
- Wheel build succeeded under the PEP 517 backend.

## Known limitations

- `sandbox/sloshing_tank.py` is an **experimental prototype**.
  Its impermeable wall conditions and spatial-field convergence still
  require validation. It is *not* promoted as a supported physical
  benchmark in this release.
- General state-dependent nonlinear mass-matrix systems are outside
  the present residual correction.
- The historical wave example monitored `t=5.01` despite
  `t_f=5.0`. PETSc may report unused options and some complex-to-real
  conversion warnings.
- Quantitative time-refinement and strict power-balance experiments
  for the two new ODE examples remain additional validation targets.

See [CHANGELOG.md](https://github.com/g-haine/scrimp/blob/v1.2.0/CHANGELOG.md)
for the change history and the
[online documentation](https://g-haine.github.io/scrimp/).
