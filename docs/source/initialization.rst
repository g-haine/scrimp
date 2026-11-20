========================================
DAE initial conditions and initialization
========================================

.. important::

   On the v1.3.0 development branch, the new `algebraic` mode
   is **opt-in and restricted to qualifying linear systems**.
   The current default remains `pseudo`. It will be changed
   only after native tests of the production integration,
   compatibility checks of the examples and a separate
   default-promotion decision. Do not conflate the qualified
   *diagnostic prototype* with a qualified production feature.

Mathematical motivation
=======================

A discretized port-Hamiltonian system contains differential
variables (such as `q` and `p` in the acoustic wave) and
algebraic variables (such as coenergy and boundary-port fields).
The prescribed initial states alone need not satisfy the
weak constitutive and boundary constraints.

With :math:\`z=(z_d,z_a)\`, only the algebraic rows of
the actual PETSc residual need satisfy

.. math::

   F_a(t_0,z_{d,0},z_a,0) = 0.

The complete differential residual at zero time derivative
need not vanish: the differential equations determine
the initial time derivatives.

Three alternative strategies
============================

Algebraic projection — `init_step_ts_type="algebraic"`
------------------------------------------------------

For a compatible **linear** model, keep the prescribed
differential coefficients *exactly fixed*, and find a
correction only to the algebraic coefficients:

.. math::

   K_{aa}\delta_a = -F_a(t_0,z_{d,0},z_a,0).

The distributed PETSc implementation extracts the
algebraic block using ownership-aware index sets,
solves a least-squares problem with KSP LSQR,
and re-evaluates **the original SCRIMP IFunction**
after correction. PETSc KSP convergence alone
is not a sufficient DAE-consistency check.

The method uses a separate algebraic feasibility criterion
and iterative accuracy target, with a bounded number of
refinement passes. It creates **no pseudo-time TS**,
preserves every prescribed dynamic coefficient, and
**rejects** unsupported or incompatible initial data
instead of silently changing them.

For the v1.3.0 experimental opt-in:

.. code-block:: python

   wave.set_time_scheme(
       ts_type="cn",
       dt=0.002,
       init_step=True,
       init_step_ts_type="algebraic",
       init_algebraic_atol=1e-10,
       init_algebraic_rtol=1e-8,
       init_algebraic_accuracy=1e-11,
       init_algebraic_max_passes=3,
   )

The sample accuracy target `1e-11` is not a universal
physical tolerance; the appropriate precision depends on
the model and desired MPI/trajectory reproducibility.
Internal LSQR controls also include
`init_algebraic_ksp_rtol`, `init_algebraic_ksp_atol`
and PETSc's `-scrimp_dae_init_` options prefix.

Currently only real systems with linear mass/stiffness,
a complete nonoverlapping algebraic/differential layout,
zero algebraic mass rows, and no lumped or substituted
ports are admitted. Other models require additional
methods or a deliberate legacy selection.

Historical pseudo-time relaxation — `pseudo`
--------------------------------------------

The historically default SCRIMP initialization is:

.. code-block:: python

   wave.set_time_scheme(
       ts_type="cn",
       dt=0.002,
       init_step=True,
       init_step_ts_type="pseudo",
       init_step_nb_iter=1,
       # init_step_dt defaults to physical dt**2
   )

It builds an **additional PETSc pseudo-time TS**,
takes its configured tiny steps, and synchronizes the
accepted state to GetFEM before the actual physical
CN or BDF solve.

**Unlike an algebraic projection**, this artificial-time
integration can change the prescribed differential
coefficients as well as the algebraic unknowns.
The compatible native wave experiment measured
coefficient-vector changes
:math:\`\|\Delta q_0\|_2=1.95958812\cdot10^{-3}\`
and :math:\`\|\Delta p_0\|_2=1.30196443\cdot10^{-5}\`.
The resulting physical CN trajectory still satisfied
its signed Hamiltonian/port-power identity: conserving
power from a **modified** initial state does not
recover the original prescribed Cauchy problem.

Explicit `pseudo` is useful for reproducing old SCRIMP
trajectories and when its relaxation behavior is intended.
Reducing the pseudo-time step alone cannot guarantee
both algebraic consistency and exact state preservation.

Tiny-step BDF2 initialization — `bdf`
-------------------------------------

The alternative historical BDF mode can be configured
to start slightly **before** the physical initial time:

.. code-block:: python

   wave.set_time_scheme(
       ts_type="cn",
       dt=0.002,  # physical time-step
       init_step=True,
       init_step_ts_type="bdf",
       init_step_bdf_order=2,
       init_step_dt=1e-6,
       init_step_from_negative_time=True,
   )

With `init_step_from_negative_time=True`, the extra
initialization TS starts at :math:\`t_0-\Delta t_{\mathrm{init}}\`
and takes a fixed step to :math:\`t_0\`.
For PETSc 3.21, its order-two BDF implementation
internally bootstraps the first accepted step
with a half-step. A **new** physical TS, with
fresh history, then starts at :math:\`t_0\`.

The small `init_step_dt` is **not** the main
physical `dt`. It must remain positive and
numerically resolvable. This approach is still
a **time integration**, not a guarantee that
initial differential variables are preserved
or all DAE constraints are satisfied.
Roundoff, stiffness, bootstrap behavior and
the DAE index matter. A tiny BDF2 step is **not**
a universal nonlinear or high-index DAE
initialization method.

Without `init_step_from_negative_time=True`,
the historical BDF initialization starts at
:math:\`t_0\` and advances through its configured
initialization steps; it is a different procedure.

Quick comparison
================

.. list-table::
   :header-rows: 1
   :widths: 20 25 25 30

   * - Strategy
     - Preserves prescribed differential DOFs?
     - Extra initialization TS?
     - Current support
   * - `algebraic`
     - **Yes**, for qualified linear cases
     - **No**
     - Explicit opt-in; rejects unsupported/incompatible data
   * - `pseudo`
     - **Not guaranteed**
     - **Yes**
     - Historical default, reproduces legacy relaxation
   * - `bdf`, order 2 and tiny step
     - **Not guaranteed**
     - **Yes**
     - Negative-time small-step alternative; DAE-dependent

The main solver choice (`ts_type="cn"`,
`ts_type="bdf"`, etc.) remains independent
of the consistency initializer mode.

Scientific qualification and limits
===================================

The **independent diagnostic** was validated with
the original acoustic wave equations, the compatible
nonzero momentum `p0=x*(3-x)`, ten fixed CN
steps, and a nonzero left-boundary input.
At MPI 1/2 its algebraic residual was approximately
:math:\`10^{-13}\`, the maximum full-trajectory
MPI discrepancy about :math:\`9.22\cdot10^{-13}\`,
and the Hamiltonian + boundary work defect around
:math:\`6\cdot10^{-15}\`.
The **published Gaussian** was rejected by the
diagnostic at its specified numerical tolerance.

This is **not** validation for arbitrary nonlinear,
high-index, substituted or lumped-port models.
Although the algebraic PETSc matrix is distributed,
the current GetFEM state synchronization still
uses `Scatter.toAll` and replicates the full
initial state per rank. Memory scalability
requires independent qualification.

The detailed mathematical derivation, exact
measurements, provenance, earlier failed experiments
and production acceptance contract are retained in:

* `docs/V130_DAE_ALGEBRAIC_INITIALIZER_DESIGN.md`
* `docs/V130_DAE_INITIAL_VALIDATION.md`
* `docs/V130_RELEASE_PLAN.md`

A change of default from `pseudo` to `algebraic`
must be accompanied by an explicit audit of
examples and unsupported models, separate
native regression evidence and updated
release documentation.
