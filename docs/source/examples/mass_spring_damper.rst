Mass--spring--damper
====================

.. _mass-spring-damper:

A two-state finite-dimensional port-Hamiltonian model
-----------------------------------------------------

Let :math:`q` be the displacement of a point mass, :math:`p=m\dot q`
its momentum, :math:`m>0` the mass, :math:`k>0` the stiffness and
:math:`c\geq0` the viscous damping coefficient. The applied force
is :math:`f(t)`. Then

.. math::

   \dot q = \frac{p}{m},\qquad
   \dot p = -kq-\frac{c}{m}p+f(t).

Use the energy variables :math:`z=(q,p)^\top` and Hamiltonian

.. math::

   H(z)=\frac12 kq^2+\frac{p^2}{2m}.

The co-energy variables are

.. math::

   e_q=\frac{\partial H}{\partial q}=kq,\qquad
   e_p=\frac{\partial H}{\partial p}=\frac{p}{m}.

The port-Hamiltonian representation is

.. math::

   \dot z =
   \underbrace{\begin{bmatrix}0&1\\-1&0\end{bmatrix}}_{J}
   \begin{bmatrix}e_q\\e_p\end{bmatrix}
   -\underbrace{\begin{bmatrix}0&0\\0&c\end{bmatrix}}_{R}
   \begin{bmatrix}e_q\\e_p\end{bmatrix}
   +\begin{bmatrix}0\\1\end{bmatrix}f(t),\qquad
   y=e_p.

Since :math:`J=-J^\top` and :math:`R\succeq0`, the continuous
power balance is

.. math::

   \dot H=-c\,e_p^2+f(t)e_p.

Discretization and code
-----------------------

This is an ODE, not a spatial PDE. The new `lumped` port kind
assigns **one** GetFEM DOF to each state and its co-state.
The reference interval :math:`(0,1)` in the script is only used
to assemble scalar GetFEM weak forms; it is not a physical spring mesh.

The implementation follows the standard SCRIMP workflow:

#. Create a `DPHS` and its reference `Domain`.
#. Add two `State` objects with `kind="lumped"`, and their `CoState` objects.
#. Call `add_lumped_port` for `q` and `p` (not `add_FEM`).
#. Add mass, interconnection, damping, forcing and constitutive `Brick` objects.
#. Initialize the two states **and their constitutive co-energies**,
   integrate with implicit PETSc BDF2, then plot.

The factor `/ell` in each brick compensates for the integration over
the reference interval length `ell`, so that the resulting scalar
equations are exactly those above. The time-dependent forcing is
assembled as a nonlinear/residual brick rather than as a fixed matrix
term; this is needed in SCRIMP's current GetFEM assembly path.

**Initial consistency is not a time-integration step.** The dynamical
initial values are prescribed as `q0=0.4`, `p0=0`. Since this
finite-dimensional model has explicit constitutive laws, the exact
algebraic initialization is simply

.. math::

   e_q(0)=kq_0=2,\qquad e_p(0)=p_0/m=0.

The example now assigns these two genuine scalar GetFEM co-energy
variables directly and disables pseudo-time relaxation. This preserves
**all** physical initial data instead of moving the momentum to
satisfy the constraints. It does not change the mass, interconnection,
damping, forcing or the subsequent BDF2 integrator. The generic SCRIMP
DAE initializers are unchanged.

To reproduce the historical pseudo-time trajectory for regression,
call `mass_spring_damper(initialization="legacy_pseudo")` explicitly.
The default is `initialization="exact"`.

The complete runnable source is:

.. literalinclude:: ../../../examples/mass_spring_damper.py
   :language: python
   :caption: examples/mass_spring_damper.py

Run from the repository root (in the SCRIMP Conda environment):

.. code-block:: console

   python examples/mass_spring_damper.py

**The user-facing example itself now saves its visualizations** in
`outputs/showcase/mass_spring_damper/figures/`:

- `mass_spring_trajectory.png/.svg`: original GetFEM displacement,
  momentum, velocity and applied force over the saved time grid.
- `mass_spring_phase_portrait.png/.svg`: phase path
  :math:`(q,p/m)`, colored by accepted PETSc time.
- `mass_spring_energy_and_power.png/.svg`: original Hamiltonian,
  viscous dissipation, applied-force power, both sides of the
  exact quadratic BDF2 identity **for steps** :math:`n\\geq2`,
  and its true discrete residual.
- `mass_spring_motion.gif`: **schematic** mass/spring motion,
  using selected original saved state coefficients. This does
  not construct or imply a finite-element spatial model.
  No interpolated or fabricated time states are included.

A machine-readable `manifest.json` and the original
`native_lumped_time_history.npz` are also written next to the
`figures` directory. The optional `plot=False` function argument
turns off figure creation without changing the physical solve.
The separate research diagnostics do **not** render these figures;
they only inspect the output produced by this example.

Verification
------------

The expected energy variation satisfies

.. math::

   H(t)-H(0)=\int_0^t\bigl[f(s)e_p(s)-c\,e_p(s)^2\bigr],\mathrm ds.

The continuous work integral is useful for understanding the physics,
but it should **not** be confused with an exact energy identity for
BDF2. This example uses the original quadratic Hamiltonian and
fixed-step BDF2, so every step :math:`n\geq2` also has the exact
*discrete* quadratic identity

.. math::

   \frac{
       H(z_n)-H(z_{n-1})
       +H(2z_n-z_{n-1})-H(2z_{n-1}-z_{n-2})
       +H(z_n-2z_{n-1}+z_{n-2})
   }{2\Delta t}
   +c\left(\frac{p_n}{m}\right)^2
   -f(t_n)\frac{p_n}{m}=0.

This identity is obtained by multiplying the **original**
semi-discrete oscillator equation by its current co-energy, not
by applying a trapezoidal rule to a saved power plot.
It covers BDF2 steps :math:`n\geq2`, not the internal PETSc
**BDF2 startup stages**, which require a separate analysis.

The original BDF2/pseudo short pilot at `dt=.01`, `t_f=.08`
ran successfully in MPI1 and MPI2 and confirmed the correct
differential equations and step-energy balance at residual
`4.13e-14`, with full MPI state parity `8.88e-16`.
**However, its physical-initial-state gate FAILED:** the pseudo
initialization performed one small `dt=0.0001` evolution and replaced
`p(0)=0` by `p(0)=-0.00019999666`, with a corresponding
`q(0)-0.4=-1.33331e-8`. That historical configuration is
**not qualified as a faithful initialization of the prescribed IVP**.

The corrected native qualification pilot checks the *same* physical
model, BDF2 time stepping and GetFEM weak forms, but starts from
exactly consistent algebraic co-energies, with no pseudo step.
It uses MUMPS in MPI 1 and 2 and measures the actual algebraic
relations
:math:`e_q=kq`, :math:`e_p=p/m`, both differentiated BDF2
equations, the exact step energy-power defect and the complete
MPI state trajectory. The independent reference is the
**forced two-state ODE**, integrated by a high-accuracy RK4
scheme, not a replacement SCRIMP model:

.. code-block:: console

   SCRIMP_DOCKER_NO_SUDO=1 \\
   SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \\
   bash scripts/run-v130-mass-spring-physical-pilot.sh

**The corrected exact-initialization short pilot has now PASSED**
its native GetFEM/PETSc MPI1/MPI2 scientific checks:
`q_0=.4,p_0=0,e_{q,0}=2,e_{p,0}=0` exactly,
maximum exact BDF2 step energy-power defect
`4.21e-14`, maximum complete DAE MPI1/2 discrepancy
`3.47e-18`, and an independent RK4 forced-ODE error
`6.00e-5` over the eight physical steps.

The remaining full-horizon **800-step** user-visualization and MPI
qualification is executed by

.. code-block:: console

   SCRIMP_DOCKER_NO_SUDO=1 \\
   SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \\
   bash scripts/run-v130-mass-spring-full-showcase.sh

The long-run script checks original program-generated PNG/SVG/GIF
on eight steps *before* running the two full MPI campaigns,
then checks every saved state and every step's exact BDF2
power identity. This long experiment is **staged, not yet
native-qualified**. A time-step refinement study and an
unforced, undamped regression (`f=0,c=0`) are distinct tests.

.. note::

   The historical pseudo-started example passed the smoke and full
   standalone execution tests, but changed the originally prescribed
   initial momentum. This motivated the exact-constitutive-data fix.
   Do not treat successful execution as proof of correct IVP data.
   The corrected default has subsequently passed short MPI1/MPI2
   scientific acceptance; the complete eight-second showcase
   remains pending until native host execution.
