# SCRIMP 1.2.0 — Release preparation (DRAFT, 2026-11-12)

> **Not released.** Do not tag, publish, or merge until the gates below pass.
> This preparation branch is based on the **experimental**
> `feature/lumped-pde-ode-sloshing-tank` branch, not the default
> `master`. Package and Sphinx versions are intentionally unchanged
> (currently 1.1.1).

## Proposed user-visible changes

- **Finite-dimensional global scalar ports**: `State(kind="lumped")`
  and `DPHS.add_lumped_port()` register one GetFEM degree of freedom,
  enabling port-Hamiltonian ODEs and certain coupled PDE--ODE models.
- Improved implicit nonlinear residual evaluation for systems with
  global scalar ports and **linear** time-derivative mass matrices.
  The new path synchronizes GetFEM with each PETSc Newton trial.
- New runnable `examples/mass_spring_damper.py` demonstrating an
  externally forced, damped oscillator and the input-power balance.
- New runnable `examples/nonlinear_pendulum.py` demonstrating a
  nonlinear gravitational constitutive law, damping and optional torque.
- New Sphinx guide `docs/source/lumped.rst` and two worked example
  pages; refreshed example navigation and installation/tests.
- Experimental `sandbox/sloshing_tank.py` (from the feature branch),
  and `sandbox/tank_convergence.py` with fixed-step energy diagnostics.
  **Keep in the sandbox.** The boundary conditions require further
  physical and numerical validation.

## Compatibility

The historical `DPHS.IFunction` branch for models without lumped
variables has intentionally not been changed. The new nonlinear
residual path currently applies only to systems containing lumped
variables and having a linear time-derivative mass matrix.
Substituted dynamic lumped co-states are not supported.

The scalar port backend currently uses a reference GetFEM integration
domain for 0D weak forms; this is an implementation limitation, not a
physical discretization of the oscillator or pendulum.

## Release gates

- [ ] Run `python -m unittest -v tests.test_lumped_ports` in Conda SCRIMP.
- [ ] Run `python -m unittest -v tests.test_lumped_examples`.
- [ ] Run `python examples/mass_spring_damper.py` and verify energy/work
      balance and exact final time.
- [ ] Run `python examples/nonlinear_pendulum.py` and verify damped
      energy decreases for zero applied torque.
- [ ] Build the complete Sphinx docs, resolve warnings and missing imports.
- [ ] Re-run at least one existing PDE example for regression
      (e.g. `examples/wave.py`).
- [ ] Perform quantitative time-refinement tests on the two new ODEs.
- [ ] Decide whether the sloshing-tank boundary conditions and
      nonlocal coupling meet the criteria for inclusion as a
      supported feature, or explicitly retain experimental status.
- [ ] Fix or document PETSc complex-to-real casting warnings and
      numerical limitations.
- [ ] Choose release version, update `setup.py` and Sphinx
      `conf.py` together; create change log and tag **only after** review.

## Observations already supplied by a SCRIMP user (tank prototype)

For mesh sizes h=0.05 and 0.025, with BDF2 time steps
0.008/0.004/0.002, the maximum energy-balance defects were about
6.1e-5 / 1.2e-5 / 3.0e-6 and the observed temporal orders were
approximately 2.35, then 1.99. Mass conservation and global
velocity constraint were near machine precision. These are
**balance-defect convergence tests**, not proof of spatial-field
convergence or impermeable boundaries.

## Review strategy

Continue work in `release/1.2.0-preparation`, based on
`feature/lumped-pde-ode-sloshing-tank`. Keep the release PR a
draft until validation and version agreement. No publication action
is performed by this document.
