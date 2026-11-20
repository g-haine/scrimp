========================================
DAE initial conditions and initialization
========================================

.. important::

   **v1.3.0 development branch (default-migration candidate):**
   The default strategy is now `algebraic` **when `init_step=True`
   and no initialization strategy is explicitly supplied**.
   Its support is limited to qualifying real, linear DAE layouts.
   The true **implicit-default native gate PASSED on 20 November
   2026**: two independent compatible linear acoustic pH-DAEs
   (conservative and dissipative), MPI 1/2, full physical CN state
   and signed power checks. The original Gaussian is rejected
   before CN, and explicit pseudo/BDF2 still work. See section 27
   of `docs/V130_DAE_INITIAL_VALIDATION.md`.

   The method still **does not support arbitrary nonlinear,
   higher-index, substituted or lumped DAE layouts**, and
   v1.3.0 as a complete release is not yet qualified.
   Unsupported/incompatible conditions raise explicit errors;
   there is no silent relaxation with `pseudo`.
   Original examples, notebooks and published Maxwell explicitly
   retain `pseudo` for historical compatibility. Broader
   memory, MPI and example gates remain pending.

.. important::

   Before blaming a nonlinear initializer for a boundary mismatch,
   check the GetFEM GWFL **matrix orientation**. Unlike Python/NumPy,
   nested matrices in GWFL are constructed from **column vectors**.
   See :doc:`getfem_gwfl` for the official upstream reference,
   correct rotation syntax and native Saint-Venant evidence.
   Correct sign conventions alone do not establish DAE consistency.

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

For the **explicit algebraic mode**, or to override
its tolerances while using the default:

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

**Core native result (20 November 2026):**
the historical ``pseudo`` main physical CN
trajectory is identical to its independently
frozen baseline (maximum scaled difference **0**).
The nonzero dynamic drift values above therefore
represent the actual accepted historical method,
not a regression caused by the new opt-in code.

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

**Core native result (20 November 2026):** on the
same compatible polynomial acoustic wave, the
explicit negative-time BDF2 ``init_step_dt=1e-6``
produced coefficient-vector drift
:math:`\|\Delta q_0\|_2=1.9595912757338825\cdot10^{-5}`
and
:math:`\|\Delta p_0\|_2=8.680599894760661\cdot10^{-10}`.
The maximum defect in the signed physical
Hamiltonian plus boundary-work balance
was :math:`4.74270349\cdot10^{-15}`.
Although smaller than the measured pseudo-TS
drift, this is **nonzero**: the tiny BDF2 step
is not an exact fixed-differential-data projection.

The opt-in **production core** algebraic method
was also qualified on the same 1/2-MPI positive
wave, matching all eleven frozen CN states
to a maximum scaled error of ``0`` on MPI 1
and ``5.30974918e-18`` on MPI 2, while
preserving every prescribed dynamic coefficient.
See section 21 of
``docs/V130_DAE_INITIAL_VALIDATION.md``
for the complete experiment and provenance.

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
     - Default candidate for supported linear DAE; rejects unsupported/incompatible data
   * - `pseudo`
     - **Not guaranteed**
     - **Yes**
     - Explicit historical relaxation; reproduces legacy trajectories
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

**Current source migration (20 November 2026):**
The development branch now selects `algebraic` as
the implicit default. The qualified core opt-in
projection and a second dissipative model both
PASSED their native MPI 1/2 physics gates.
A dedicated frozen-image no-explicit-initializer
regression and the wider examples/memory audit
must pass **before this behavior is released**.
The published Maxwell example and historical
notebooks now explicitly retain `pseudo`.

With `init_step=False`, omitting the initializer
choice remains valid: SCRIMP performs no
initialization step. An explicit contradictory
request `init_step_ts_type="algebraic"` together
with `init_step=False` raises an error.

See `scripts/run-v130-dae-default-promotion.sh`
for the native migration acceptance test.

Mesh dependence of the algebraic residual criterion
=================================================

.. important::

   The present algebraic-initializer acceptance threshold is applied
   to the Euclidean norm of a **finite-element coefficient residual**.
   It is **not mesh independent**, even though its formula has
   no explicit mesh-size parameter. Do not interpret a lower raw
   residual on a finer mesh as a smaller physical constraint defect.

At a fixed mesh, the current diagnostic checks

.. math::

   \|F_{a,h}\|_2 \leq \varepsilon_{\rm abs}
   + \varepsilon_{\rm rel}\|F_{a,h}^{(0)}\|_2.

Here :math:`h` means the spatial mesh step, **not** the
fluid-height state named ``h`` in the shallow-water model.
The residual entries, number of FE basis functions, weak
integration weights, conditioning and units change on refinement.
For instance, a smooth volume load in two dimensions can yield
a Euclidean residual that scales like the mesh step even when
the underlying physical defect does not decrease.

A future mesh-aware acceptance test must distinguish:

* Algebraic convergence and native GetFEM/PETSc residuals at
  each fixed mesh, with solver accuracy checked separately.
* **Per-constraint physical residuals**, measured with
  suitable dual norms. For an appropriate Riesz/Gram
  matrix :math:`R_{a,h}`, one possible norm is

  .. math::

     \|F_{a,h}\|_{V_{a,h}\u0027}
       =\sqrt{F_{a,h}^{\mathsf T}R_{a,h}^{-1}F_{a,h}}.

  Volume and boundary blocks require their own Gram matrices
  and physically meaningful characteristic scales.
* Boundary compatibility measured through
  :math:`L^2(\Gamma_h)` velocity/traction defects, with
  appropriate nonzero reference scales.
* Refinement studies using consistent physical initial data,
  tracking these errors over several spatial mesh sizes.

A physically incompatible *discrete* initial condition must
not be accepted simply because its residual vector becomes
small under refinement. Conversely, finite precision or
a poorly scaled coefficient norm is not a proof of physical
incompatibility.

Finite-difference perturbations and roundoff
===========================================

The original rectangular shallow-water example has
fluid heights of order three, density :math:`\rho=1000`
and gravity :math:`g=10`, so the pressure/coenergy variable
is of order :math:`3\times10^4`. When a normalized
algebraic direction is perturbed with a very small
step :math:`\delta=10^{-6}`, subtracting nearly equal
assembled residual vectors may amplify roundoff.

The first native #28 shallow-water Jacobian check returned
a relative difference of :math:`5.53585698\times10^{-5}`
for this single small perturbation. This is a numerical
**failure at that perturbation**, not yet proof of
an incorrect assembled Jacobian.

The follow-up source-only diagnostic evaluates central
differences for several step sizes, maintains the
**original strict tolerance**, checks the independent
affine algebraic identity, and retains a JSON report
even on numerical rejection. Native evidence remains
pending. No production nonlinear initializer is
qualified by this diagnostic alone.
