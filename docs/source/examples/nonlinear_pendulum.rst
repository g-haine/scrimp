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

The full sine GetFEM weak formulation, constitutive bricks,
damping and forcing are retained. The canonical example uses
**one** physical time scheme: PETSc `ts_type="cn"`, fixed time step,
without any pseudo-time initialization. The former BDF2 short
qualification was a useful diagnostic but **does not qualify CN**.
The CN native qualification remains pending.

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

Numerical validation — Crank–Nicolson
--------------------------------------

The canonical example now uses **PETSc Crank–Nicolson** only.
In PETSc 3.26, `ts_type="cn"` means the **endpoint trapezoidal**
method, not the implicit midpoint method. Define the endpoint
average :math:`\bar a_n=(a_{n-1}+a_n)/2`. At every step,
the expected equations are

.. math::

   \begin{aligned}
   \frac{\theta_n-\theta_{n-1}}{\Delta t}
      &= \bar e_{p,n},\\
   \frac{p_n-p_{n-1}}{\Delta t}
      &= -\bar e_{\theta,n}-d\,\bar e_{p,n}
         +\bar\tau_n.
   \end{aligned}

The corresponding **exact discrete co-energy power pairing**
tests the actual CN residual:

.. math::

   \bar e_{\theta,n}
      \frac{\theta_n-\theta_{n-1}}{\Delta t}
   +\bar e_{p,n}
      \frac{p_n-p_{n-1}}{\Delta t}
   +d\,\bar e_{p,n}^2
   -\bar\tau_n\,\bar e_{p,n}=0.

However, :math:`H=p^2/(2I)+mg\ell(1-\cos\theta)`
is nonquadratic. The endpoint-averaged nonlinear effort
is **not an exact discrete gradient**. Consequently,

.. math::

   \frac{H_n-H_{n-1}}{\Delta t}
   =\bar\tau_n\,\bar e_{p,n}
    -d\,\bar e_{p,n}^{2}
    +\Gamma_n,

where :math:`\Gamma_n` is the **measured CN nonlinear
chain-rule defect**

.. math::

   \Gamma_n=
      \frac{H_n-H_{n-1}
      -\bar e_{\theta,n}(\theta_n-\theta_{n-1})
      -\bar e_{p,n}(p_n-p_{n-1})}{\Delta t}.

The canonical `python examples/nonlinear_pendulum.py` directly
saves native trajectory, phase portrait, nonlinear energy/power
plots (PNG/SVG) and a **schematic** GIF, all from actually
accepted GetFEM states. Files are written to
`outputs/showcase/nonlinear_pendulum/free/figures/`.
The forced variant, called with `torque_amplitude=0.1`,
writes instead to `outputs/showcase/nonlinear_pendulum/forced/figures/`.
The plots explicitly distinguish true changes of :math:`H`,
CN discrete port work and dissipation, and the chain-rule
defect; no exact Hamiltonian conservation is implied.

The short *historical BDF2* experiment (26 November 2026,
source `2ae4bc7d`) passed MPI1/MPI2 with and without forcing.
**It is not a native qualification of CN**, which must be tested
again with the original GetFEM/PETSc implementation:

.. code-block:: console

   SCRIMP_DOCKER_NO_SUDO=1 \\
   SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \\
   bash scripts/run-v130-nonlinear-pendulum-physical-pilot.sh

The pilot checks both torque cases, full sine constitutive
data, all CN step equations and power pairings, independent
RK4 trajectories, and complete DAE MPI parity.

.. note::

   **CN native qualification is pending** until the above
   user-host calculation has been run and examined.
   The original default twelve-second horizon and
   example-native graphics require a further release
   acceptance run. Do not infer a convergence order or
   exact energy conservation from the earlier BDF2 pilot.
