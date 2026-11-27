Energy, powers and time discretization (v1.3.0)
==================================================

This chapter explains what the public SCRIMP post-processing API **actually
computes**. The meaning of an energy residual depends on the discrete
Hamiltonian, the port orientations and the **time integrator**. A successful
PETSc TS solve is not, by itself, a numerical proof of passivity.

Hamiltonian: physical value, not a time quadrature
--------------------------------------------------

After calling :meth:`scrimp.dphs.DPHS.solve`, the method
:meth:`scrimp.dphs.DPHS.compute_Hamiltonian` evaluates every Hamiltonian
term from the **genuine stored GetFEM state coefficients** at the corresponding
stored physical time. Expressions and spatial integration regions are those
declared by the user. They are evaluated using GetFEM's own weak-form
integration, not by interpolating nodal arrays or fitting a time polynomial.

In v1.3.0, the model time is set for **every** evaluation. Recomputing a
Hamiltonian replaces the previous history; it never appends a second copy
of the same trajectory. A subsequent solve invalidates the cached values,
and invalid/nonfinite intermediate values cannot be published as a
successfully computed trajectory. This makes the sampling consistent; it
**does not** imply that the time integrator conserves (or dissipates)
that Hamiltonian exactly.

The following calls are intended to work in a regular user script:

.. code-block:: python

   system.solve()
   system.compute_Hamiltonian()
   H = system.get_Hamiltonian()
   system.plot_Hamiltonian(
       with_powers=True,
       save_figure=True,
       filename="Hamiltonian.png",
       show=False,
   )

With ``system.set_output_name("my_system")``, the generated files include
``outputs/my_system/png/Hamiltonian.png``,
``outputs/my_system/svg/Hamiltonian.svg`` and
``outputs/my_system/npz/Hamiltonian.npz``. The NumPy file records the
actual time grid, total Hamiltonian, individual Hamiltonian terms and, when
requested, **raw oriented algebraic-port powers** and their sampling
conventions. All ranks participate in collective GetFEM work, while
shared figures and data are written only on rank 0.

Port powers: raw signs, not universal energy balance
----------------------------------------------------

:meth:`scrimp.dphs.DPHS.compute_powers` uses the ports' actual GetFEM
flow-effort weak products, preserving their declared boundary orientation
and dissipative status. Neither this method nor the integration helper
silently changes the physical sign of a control or an interface port.

For a physical model whose sign conventions yield

.. math::

   H(t_b)-H(t_a)
   = \int_{t_a}^{t_b} \bigl(P_{\rm boundary}(t)
     -D(t)\bigr)\,\mathrm dt,

the correct assignment of a port to boundary input, output, interconnection
or dissipation is **model-specific**. For example, two oriented faces of a
conforming heat-wave interface have equal exchanged power with opposite
contribution to the total energy balance. A generic port plot cannot
establish that identity without verifying the model's actual weak forms.

Crank--Nicolson: only consecutive *accepted* PETSc states
---------------------------------------------------------

For PETSc ``ts_type="cn"``, SCRIMP uses the ``cn_step`` power convention
**only if every accepted TS step has its two consecutive end states saved**.
The check uses the real accepted-step count, the saved times, and (when
time stepping is fixed) the declared ``dt``. For each step `k >= 1`,
the current implementation evaluates the power using the actual
GetFEM coefficient-space midpoint

.. math::

   z_{k-1/2}=\tfrac12\bigl(z_{k-1}+z_k\bigr),

with the existing SCRIMP convention for the model time, and accumulates
its signed contribution as

.. math::

   W_{\rm CN}(t_n)
   = \sum_{k=1}^{n}(t_k-t_{k-1})P_k^{\rm CN}.

The initial value ``P[0]`` is **instantaneous** and is not a previous
step's work. It does not enter the sum.

For a quadratic, time-independent Hamiltonian with the properly matched
discrete weak forms and constraints, this convention permits the familiar
*model-specific* exact CN energy-power identity, up to solver tolerances.
It is **not universally exact** for nonquadratic Hamiltonians,
time-dependent constitutive data, or arbitrarily prescribed input
time functions. In particular, one must not mistake a GetFEM evaluation
using the source's stored time convention for an independent proof that
all time-dependent forcings were integrated by a midpoint rule.

CN with ``save_step > dt`` does **not** provide the adjacent states
needed for each step power. The API therefore switches to
``instantaneous`` samples at each truly stored `(t_k,z_k)`
and uses the *approximate nodal trapezoidal integral* below.
It must **never** label differences between nonconsecutive saved
states as genuine CN step work.

BDF, BDF2, ARKIMEX and other PETSc TS types
-------------------------------------------

For non-CN schemes, ``compute_powers()`` samples each port
instantaneously at the **stored** PETSc solution snapshots.
It does not reconstruct derivative histories, hidden stage states,
stage weights, BDF startup steps or solver-internal quadratures.

The default cumulative work drawn by
:meth:`scrimp.dphs.DPHS.plot_Hamiltonian` is the nodal trapezoid

.. math::

   W_{\rm trap}(t_n)
   = \sum_{k=1}^{n}\frac{t_k-t_{k-1}}{2}
     \bigl(P(t_{k-1},z_{k-1})+P(t_k,z_k)\bigr).

This is **only a post-processing approximation** to the continuous-time
work integral. The approximation can differ from the method's internal
discrete work even when the solution trajectory converges correctly.
With sparse saved states, rapidly varying controls, sharp transients,
stiff or nonlinear dynamics, it can differ noticeably. No general
accuracy order for the energy balance follows from the name BDF2,
the TS solver order, or the trapezoidal sampling formula alone.

In particular:

* **BDF/BDF2:** the trapezoidal integral does not implement a BDF
  discrete chain rule, G-stability energy estimate, or the correct
  multistep method's internal work. A quantity displayed as
  ``Balance (legacy port signs; model-dependent)`` is a **diagnostic**,
  not an automatically conserved physical invariant.
* **ARKIMEX / other Runge--Kutta and IMEX methods:** the stored
  endpoints are not a substitute for the internal stage values and
  Butcher weights. The plotted nodal integral is not an exact
  solver-stage quadrature.
* **Adaptive stepping or decimated output:** use the true saved
  time differences; do not assume they equal the originally
  requested fixed ``dt``. If every accepted CN state cannot be
  confirmed, the method conservatively selects instantaneous sampling.

To establish an *exact* discrete energy or passivity assertion for a
particular model and time integrator, calculate the corresponding
discrete power balance using the actual method/stage data and verify
it independently. Do not re-label an approximate generic plot as
a numerical acceptance test.

How to inspect the sampling used in a run
------------------------------------------

After ``system.compute_powers()``, the original algebraic ports expose
``port.get_power_sampling()``, which is ``"cn_step"`` or
``"instantaneous"``. Inspect their raw ``port.get_power()`` time
series and the time grid ``system.solution["t"]``, or load the
``npz`` saved by ``plot_Hamiltonian()``.
In the figure, approximate cumulative work is explicitly described as
``nodal trapezoid (approximate)``.

The source of these conventions is
``scrimp/utils/power_sampling.py`` and
``scrimp/utils/power_integration.py``.
Negative tests are in
``tests/test_v130_user_output_power_hamiltonian.py``.
For a genuine two-sided coupled example that makes the oriented powers
and dissipation visible, see :doc:`examples/heat_wave`.
