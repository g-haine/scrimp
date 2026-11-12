Nonlinear pendulum
==================

.. _nonlinear-pendulum:

A nonlinear port-Hamiltonian ODE
--------------------------------

Consider a point mass :math:`m` on a massless rod of length
:math:`\ell`, rotating in a vertical plane under gravity :math:`g`.
Let :math:`\theta` denote its angle, :math:`I=m\ell^2` the
moment of inertia and :math:`p=I\dot\theta` its angular momentum.
An optional external torque :math:`\tau(t)` and viscous damping
:math:`d\geq0` complete the model:

.. math::

   \dot\theta = \frac{p}{I},\qquad
   \dot p=-mg\ell\sin\theta-\frac{d}{I}p+\tau(t).

The energy is

.. math::

   H(\theta,p)=\frac{p^2}{2I}
                 +mg\ell(1-\cos\theta).

It follows that

.. math::

   e_\theta=mg\ell\sin\theta,\qquad
   e_p=\frac{p}{I}.

The canonical port-Hamiltonian representation uses

.. math::

   \begin{bmatrix}\dot\theta\\\dot p\end{bmatrix}
   =
   \begin{bmatrix}0&1\\-1&0\end{bmatrix}
   \begin{bmatrix}e_\theta\\e_p\end{bmatrix}
   -
   \begin{bmatrix}0&0\\0&d\end{bmatrix}
   \begin{bmatrix}e_\theta\\e_p\end{bmatrix}
   +
   \begin{bmatrix}0\\1\end{bmatrix}\tau(t).

The power balance reads

.. math::

   \dot H=-d\,e_p^2+\tau(t)e_p.

This example intentionally uses the full
:math:`\sin\theta` relation, **not** the small-angle approximation.

Implementation notes
--------------------

The states `theta` and `p` and their corresponding efforts are
global `lumped` variables with one degree of freedom each. The
nonlinear constitutive law is written with an explicitly nonlinear
GetFEM brick:

.. code-block:: python

   S.Brick("CR_theta_nonlinear",
           "gravity_torque*sin(theta)*Test_e_theta/ell",
           [1], linear=False)

The GetFEM reference interval is only a weak-form assembly support,
and the `/ell` normalization yields a single ODE residual.
All variables are spatially constant.

The complete runnable example is:

.. literalinclude:: ../../../examples/nonlinear_pendulum.py
   :language: python
   :caption: examples/nonlinear_pendulum.py

Run from the repository root:

.. code-block:: console

   python examples/nonlinear_pendulum.py

By default the example starts with a large angle
:math:`\theta(0)=1.2` radians and no external torque.
It plots the time trajectory, phase portrait and Hamiltonian.
For periodic forcing, call:

.. code-block:: python

   from examples.nonlinear_pendulum import nonlinear_pendulum
   pend = nonlinear_pendulum(torque_amplitude=0.1, plot=True)

Numerical validation
--------------------

With zero torque and positive damping, the continuous Hamiltonian is
non-increasing. With forcing, compare its change with the integrated
input power minus viscous losses. The observed numerical defect should
decrease under temporal refinement. Energy dissipation from BDF2 is
numerical in addition to the physical damping.

.. note::

   This example has passed the automated numerical smoke tests and a
   full default standalone run in the SCRIMP GetFEM/PETSc environment.
   Quantitative temporal-refinement studies remain a separate check.
