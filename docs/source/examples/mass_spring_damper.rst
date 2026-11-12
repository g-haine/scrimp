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
#. Initialize the states, integrate with implicit PETSc BDF2, then plot.

The factor `/ell` in each brick compensates for the integration over
the reference interval length `ell`, so that the resulting scalar
equations are exactly those above. The time-dependent forcing is
assembled as a nonlinear/residual brick rather than as a fixed matrix
term; this is needed in SCRIMP's current GetFEM assembly path.

The complete runnable source is:

.. literalinclude:: ../../../examples/mass_spring_damper.py
   :language: python
   :caption: examples/mass_spring_damper.py

Run from the repository root (in the SCRIMP Conda environment):

.. code-block:: console

   python examples/mass_spring_damper.py

Verification
------------

The expected energy variation satisfies

.. math::

   H(t)-H(0)=\int_0^t\bigl[f(s)e_p(s)-c\,e_p(s)^2\bigr],\mathrm ds.

Compute the two integrals from the stored time histories to verify
the numerical energy balance. BDF2 is not an exactly energy-preserving
integrator: the measured defect should be examined under time-step
refinement. A free undamped oscillator (`f=0`, `c=0`) is a useful
additional regression case.

.. note::

   This example is part of the next-release preparation branch.
   The Python source has not been runtime-tested in the authoring
   environment; full GetFEM/PETSc testing is required before release.
