Finite-dimensional port-Hamiltonian systems
=============================================

.. _lumped:

Global scalar ports (experimental)
----------------------------------

SCRIMP's :class:`scrimp.dphs.DPHS` can represent **global scalar
states** and co-states in addition to distributed finite-element fields.
This feature is currently under validation for the next release.

For instance, a mass--spring--damper has two energy states, displacement
:math:`q(t)` and momentum :math:`p(t)`, and no spatial dependence.
The corresponding energy is

.. math::

   H(q,p) = \tfrac{k}{2} q^2 + \tfrac{1}{2m} p^2.

There must be exactly **one degree of freedom** for each variable,
independent of the number of elements in any mesh used elsewhere.

.. code-block:: python

   import scrimp as S

   sys = S.DPHS("real")
   sys.set_domain(S.Domain("Interval", {"L": 1.0, "h": 1.0}))
   q = S.State("q", "Displacement", "lumped")
   sys.add_state(q)
   sys.add_costate(S.CoState("e_q", "Elastic effort", q))
   sys.add_lumped_port("q")
   sys.set_initial_value("q", 0.3)

A global variable is registered with GetFEM's constant-size
`Model.add_variable(name, 1)`, not with `add_fem_variable`.
The method :meth:`scrimp.dphs.DPHS.add_lumped_port` is used **instead
of** `add_FEM` for such ports.

Global algebraic variables
--------------------------

An algebraic variable such as the tank velocity :math:`v(t)` can also
be registered as a single degree of freedom:

.. code-block:: python

   sys.add_port(S.Port("velocity", "v", "unused", "lumped",
                       algebraic=True, substituted=True))
   sys.add_lumped_port("velocity")

A `substituted=True` algebraic port registers only the flow variable,
which can participate in nonlinear constitutive relations. The flow
has no time derivative and needs no initial value through
`set_initial_value`. Dynamic substituted co-states are not supported
by this initial implementation.

Why do bricks use a spatial domain?
-----------------------------------

SCRIMP's existing bricks are assembled using GetFEM weak forms. To
reuse those same operators for global variables, the example integrates
the scalar residual over a convenient **reference interval** of length
:math:`L`, dividing by its measure:

.. math::

   \int_0^L (\dot q-e_p)\,\frac{\mathrm dx}{L}=0
   \quad\Longleftrightarrow\quad
   \dot q=e_p.

In SCRIMP, these two contributions are represented as

.. code-block:: python

   S.Brick("M_q", "q*Test_q/L", [1], dt=True, position="flow")
   S.Brick("J_qp", "e_p*Test_q/L", [1], position="effort")

Here `L` is a GetFEM macro defined by the user. This reference
interval is an **implementation detail**, not a physical continuum
model for the oscillator. There is **no mesh refinement to perform**
for these 0D examples.

In a **coupled PDE--ODE** system, integration also enables a true
nonlocal constraint. In the experimental sloshing tank, for example,

.. math::

   M_T v + \int_0^L m(x,t)\,\mathrm dx - P = 0

is assembled by the single weak-form contribution

.. code-block:: python

   S.Brick("velocity_constraint",
           "(mom+(Mt*v-P)/L)*Test_v", [1])

where `mom` is a spatially distributed finite-element field.
Only one global row is assembled for `v`.

Solver and present limits
-------------------------

The extension updates GetFEM at each PETSc Newton trial iterate,
which is required for implicit nonlinear PDE--ODE coupling.

For systems with global scalar ports **and a linear time-derivative
mass matrix**, the fully implicit residual takes the form

.. math::

   F(t,z,\dot z)=M\dot z + K_{\mathrm{lin}}z
                - \mathrm{rhs}_{\mathrm{GetFEM}}(z,t).

The nonlinear tangent is used for the Jacobian, not multiplied by
the trial state when reconstructing the residual. This distinction
matters for external forcing and nonlinear constitutive laws.

.. warning::

   The global-port backend has been exercised on a sloshing-tank
   PDE--ODE prototype, but this does not constitute full validation.
   Boundary conditions, convergence of spatial fields and the
   original state-dependent nonlinear mass-matrix path remain open
   issues. Do not infer that all PDE--ODE couplings are supported.

Examples
--------

* :doc:`examples/mass_spring_damper` -- linear damped oscillator with a
  sinusoidal force.
* :doc:`examples/nonlinear_pendulum` -- nonlinear restoring torque
  :math:`mg\ell\sin(\theta)` and viscous damping.
* `sandbox/sloshing_tank.py` -- experimental coupling to one-dimensional
  viscous shallow-water dynamics; retained in `sandbox` until its
  wall conditions have been validated.
