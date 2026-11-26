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

The example now initializes its **exact constitutive co-energies**
from the original nonlinear sine law, preserving the prescribed
initial differential data, without a pseudo-time evolution:

.. math::

   \theta_0=1.2,\qquad p_0=0,\qquad
   e_{\theta,0}=mg\ell\sin(1.2),\qquad e_{p,0}=0.

Only the initialization path has changed: the nonlinear
GetFEM weak formulation, constitutive bricks, damping,
forcing and original physical BDF2 scheme remain unchanged.
The historical one-step pseudo-time relaxation is still
available explicitly via
`nonlinear_pendulum(initialization="legacy_pseudo")`.
The corrected default has not been declared native-qualified
before the independent pilot evidence is received.

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
decrease under temporal refinement. BDF2 adds its own time-discretization
effects to the physical dynamics.

**Important distinction:** unlike the linear mass--spring example,
the pendulum Hamiltonian is **not quadratic**, so a telescoping
quadratic BDF2 identity expressed purely in differences of
`H(z_n)` is **not exact**. For the original discrete BDF2
equations at every physical step :math:`n\geq2`, the correct
test is the true co-energy-weighted differential identity

.. math::

   e_{\theta,n}\,
      \frac{3\theta_n-4\theta_{n-1}+\theta_{n-2}}{2\Delta t}
   +e_{p,n}\,
      \frac{3p_n-4p_{n-1}+p_{n-2}}{2\Delta t}
   +d\,e_{p,n}^2-\tau(t_n)e_{p,n}=0.

This is the power pairing of the actual **nonlinear residuals**,
not a fictitious discrete gradient formula for the Hamiltonian.
It does not cover PETSc's internal startup half-step.

An independent MPI1/MPI2 diagnostic with both original
zero-torque and externally forced `0.1*sin(2*pi*0.3*t)`
cases is prepared:

.. code-block:: console

   SCRIMP_DOCKER_NO_SUDO=1 \\
   SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \\
   bash scripts/run-v130-nonlinear-pendulum-physical-pilot.sh

The diagnostic checks all four true GetFEM lumped unknowns,
the exact initial nonlinear constitutive laws, the original
Hamiltonian, both BDF2 differential equations, their
co-energy power pairing, the independent **full sine**
RK4 reference, and complete MPI1/MPI2 parity.

.. note::

   The example has previously completed numerical smoke tests and
   an original standalone execution. That alone is not an
   independent physical qualification. The exact-initial-state
   revised pilot, full default horizon and reproducible example-native
   PNG/SVG/GIF remain **native pending** until their actual evidence
   has been reviewed.
