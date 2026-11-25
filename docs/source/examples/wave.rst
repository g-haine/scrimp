The wave equation — 2A: conservative anisotropic membrane
=========================================================

.. _Wave_2D:

.. automodule:: examples.wave
   :members:
   :undoc-members:
   :show-inheritance:

Why this example?
-----------------

The :doc:`first tutorial <../started>` showed a **one-dimensional**
string and its two boundary ports. Our next model is genuinely
**two-dimensional**, heterogeneous and anisotropic, with **mixed
boundary control**: three boundary segments prescribe an outward
normal traction, while the fourth prescribes a *velocity*. We will
keep exactly the same materials, mesh and boundary convention when
we add an internal dissipative port in example **2B**. We will
then study the complementary co-energy formulation in example
**2C** (:doc:`wave_coenergy`).

This page describes the actual code in ``examples/wave.py`.
Its snippets are included **directly from that file**, avoiding
manual discrepancies between the tutorial and the program.
The validated source is the same for all PETSc solver backends.

The physical model
------------------

Let :math:`\Omega=(0,2)\times(0,1)` and let :math:`w(t,x,y)`
denote the transverse displacement of a membrane. We consider

.. math::

   \rho(x)\,\partial_{tt} w
   - \operatorname{div}(T(x,y)\nabla w)=0,

where the in-plane tension tensor and areal mass density are

.. math::

   T(x,y)=\begin{bmatrix}5+x&xy\\xy&2+y\end{bmatrix},
   \qquad \rho(x)=3-x.

Here :math:`T` is **not** a scalar Young modulus. It is a
symmetric, positive-definite anisotropic tensor; both :math:`T`
and :math:`\rho` vary spatially. The lower, right and upper
boundaries are :math:`\Gamma_N` (GetFEM regions 10, 11 and 12).
The left boundary is :math:`\Gamma_D` (region 13). Writing
:math:`n` for the **outward** unit normal, the controls are

.. math::

   n\cdot(T\nabla w)=u_N \quad\hbox{on }\Gamma_N,
   \qquad \partial_t w=u_D \quad\hbox{on }\Gamma_D.

The second condition prescribes *velocity*, **not** displacement.
This is sometimes called a velocity-Dirichlet port in the
port-Hamiltonian literature.

Energy variables, co-energy variables and power
-----------------------------------------------

We choose the energy states

.. math::

   q=\nabla w,\qquad p=\rho\,\partial_t w,

and their conjugate co-energies

.. math::

   e_q=Tq,\qquad e_p=p/\rho.

Thus :math:`q` and :math:`e_q` are **vectors**,
whereas :math:`p` and :math:`e_p` are **scalars**. The
first-order port-Hamiltonian dynamics are

.. math::

   \partial_t q=\nabla e_p,\qquad
   \partial_t p=\operatorname{div}e_q.

The physical Hamiltonian is

.. math::
   :label: wave_2d_conservative_energy

   \mathcal H(q,p)=\tfrac12\int_{\Omega}
   \left(q\cdot Tq+\frac{p^2}{\rho}\right)\,\mathrm dx.

In the **GetFEM Generic Weak Form Language** (GWFL),
``q.T.q`` denotes :math:`q\cdot Tq`: ``T`` is the
**tensor parameter's identifier**, not the transpose
operator from Python/NumPy. The actual code therefore
correctly uses ``"0.5*q.T.q"``.

With our boundary orientation, the collocated outputs are

.. math::

   y_N=e_p \quad\hbox{on }\Gamma_N,\qquad
   y_D=n\cdot e_q \quad\hbox{on }\Gamma_D.

The energy exchanged with the environment is

.. math::
   :label: wave_2d_conservative_power

   \frac{\mathrm d\mathcal H}{\mathrm dt}
   =\int_{\Gamma_N}u_N y_N\,\mathrm ds
    +\int_{\Gamma_D}u_D y_D\,\mathrm ds.

Our short demonstration uses **zero** Neumann tractions,
with a nontrivial prescribed left velocity. Thus the
left-hand boundary supplies the entire mechanical work;
the right, top and bottom boundary powers are zero.

From the strong equation to PFEM
--------------------------------

Integrating :math:`\partial_t p=\operatorname{div}e_q`
against a test function :math:`\varphi_p` gives

.. math::

   \int_\Omega\partial_t p\,\varphi_p
   =-\int_\Omega e_q\cdot\nabla\varphi_p
    +\int_{\Gamma_N}u_N\varphi_p
    +\int_{\Gamma_D}y_D\varphi_p.

The first equation is tested without integration by parts:

.. math::

   \int_\Omega\partial_t q\cdot\varphi_q
   =\int_\Omega\nabla e_p\cdot\varphi_q.

These two choices yield the partitioned finite-element
discretization (PFEM): ``DG1`` for vector strain
:math:`q`, ``CG2`` for scalar momentum :math:`p`,
and ``DG1`` boundary traces. The exact finite-element
spaces and material definitions are shown below.

In SCRIMP's assembly convention, the *flow-side*
algebraic variables for the Neumann velocity observations
are **``-Y_B``, ``-Y_R``, ``-Y_T``**, and the
velocity input on the left appears as **``-U_L``**.
The corresponding effort-side test terms are negative
velocity traces. The signs make the measured
powers consistent with :eq:`wave_2d_conservative_power`.
They have been checked by native GetFEM/PETSc simulations,
not merely assumed from the matrix layout.

Build the real SCRIMP model
---------------------------

To reproduce the snippets in a notebook, begin with these
imports and the *same tested defaults* (the command-line entry
point accepts other values):

.. code-block:: python

   import numpy as np
   import scrimp as S

   mesh_size = 0.2
   final_time = 0.08
   time_step = 0.008
   factor_solver = None  # Production default: MUMPS

The four boundary identifiers are taken directly from
the Python module, too:

.. literalinclude:: ../../../examples/wave.py
   :language: python
   :start-after: # TUTORIAL 2A-00:
   :end-before: def wave_eq(

The following seven sections are taken **literally** from
``examples/wave.py:wave_eq()``; no duplicate hand-written
Python excerpts are maintained in Sphinx. The boundary
index tuple ``BOUNDARIES`` is defined at module level.

**1. Instantiate the DPHS, domain, states and co-states.**

.. literalinclude:: ../../../examples/wave.py
   :language: python
   :start-after: # TUTORIAL 2A-01:
   :end-before: # TUTORIAL 2A-02:
   :dedent: 4

**2. Define the four mixed boundary ports.** Only the left port
uses ``position="flow"``, since it imposes the velocity
instead of the normal traction.

.. literalinclude:: ../../../examples/wave.py
   :language: python
   :start-after: # TUTORIAL 2A-02:
   :end-before: # TUTORIAL 2A-03:
   :dedent: 4

**3. Choose the PFEM spaces and material coefficients.**
GetFEM's nested GWFL matrix definition is symmetric and
represents the actual inhomogeneous tensor :math:`T`.

.. literalinclude:: ../../../examples/wave.py
   :language: python
   :start-after: # TUTORIAL 2A-03:
   :end-before: # TUTORIAL 2A-04:
   :dedent: 4

**4. Assemble the Dirac structure and constitutive laws.**
The rows with ``dt=True`` represent the differential
part of the port-Hamiltonian DAE. Boundary observations and
constraints are algebraic; the signs are visible here.

.. literalinclude:: ../../../examples/wave.py
   :language: python
   :start-after: # TUTORIAL 2A-04:
   :end-before: # TUTORIAL 2A-05:
   :dedent: 4

**5. Prescribe compatible initial data and controls.**
We set :math:`q_0=p_0=0`, so :math:`e_q=e_p=0`
at the initial instant, and apply

.. math::

   u_D(t,y)=0.2\sin(4\pi t)\sin(\pi y).

This input also vanishes at :math:`t=0` and at the
corners of the left boundary. It is already nonzero
during the ten-step installation smoke.

.. literalinclude:: ../../../examples/wave.py
   :language: python
   :start-after: # TUTORIAL 2A-05:
   :end-before: # TUTORIAL 2A-06:
   :dedent: 4

**6. Solve the actual DAE.** The ``algebraic`` initial
projection is strict and preserves the prescribed
differential data. For the quick installation check
we use fixed Crank--Nicolson,
:math:`t_f=0.08`, :math:`\Delta t=0.008`,
and ten accepted time steps, not a five-second
publication reproduction.

.. literalinclude:: ../../../examples/wave.py
   :language: python
   :start-after: # TUTORIAL 2A-06:
   :end-before: # TUTORIAL 2A-07:
   :dedent: 4

**7. Register the real anisotropic Hamiltonian.**
The two integrals are evaluated by GetFEM on the
actual finite-element fields.

.. literalinclude:: ../../../examples/wave.py
   :language: python
   :start-after: # TUTORIAL 2A-07:
   :end-before: # TUTORIAL 2A-08:
   :dedent: 4

To reproduce the complete program, use the canonical
entry point, which encapsulates the steps above
(and defines ``BOUNDARIES``):

.. code-block:: bash

   python -m examples.wave --solver mumps --check --export --output outputs/my_wave_2d

Installation and genuine scientific checks
------------------------------------------

The result must be a real **GetFEM/PETSc** solve, not an LU
factorization smoke. The common ``measure_wave_2d()``
function checks convergence, the strict initial DAE
projection, nonzero left power and the mechanical
energy/work equality **at every accepted CN step**:

.. math::

   \mathcal H(t_k)-\mathcal H(t_{k-1})
   -\Delta t_k\left(P_{B,k}+P_{R,k}+P_{T,k}+P_{L,k}\right)
   \simeq 0.

The port powers :math:`P_{\ast,k}` are the genuine
``cn_step`` GetFEM integrals sampled on the **mean FE
state** of the accepted step, not instantaneous endpoint
powers. Never trapezoid-integrate the ``P[0]`` sample.

The **initial MUMPS native pilot** passed in MPI 1
and MPI 2, with 2,230 DOFs, ten accepted steps and
numerical energy/power defects below
:math:`3\times10^{-18}`. Both ranks had identical
full trajectories; the MPI 1/MPI 2 state difference
was below :math:`4.2\times10^{-14}`. All four
fields (``q``, ``p``, ``e_q``, ``e_p``)
were exported as valid, root-only shared VTU/PVD,
with **11 time frames per field**. This is a
qualified **MUMPS MPI1/2** result only, not yet
a five-backend installation acceptance.

The complete installation acceptance requires:
MUMPS (MPI 1/2), sequential SuperLU (MPI 1),
and SuperLU_DIST (MPI 1/2). The three backends
must reproduce the same state, Hamiltonian,
physical powers, and ParaView outputs.
SuperLU standard is **serial-only**.

To save the four original GetFEM variable series,
open ``outputs/my_wave_2d/pv/q/q.pvd``,
``pv/p/p.pvd``, ``pv/e_q/e_q.pvd`` or
``pv/e_p/e_p.pvd`` in ParaView and play the
11 frames. The vector fields ``q`` and
``e_q`` may be displayed with arrows; the
velocity is ``e_p`` (not momentum ``p``).

Longer native showcase: 100-step energy and wave propagation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**PASS, 25 November 2026**, on the complete original
GetFEM/PETSc MUMPS MPI1 archive
``v130-wave2d-2a-showcase-bd6f5468-teEK1Ejz.zip``
(source ``bd6f5468e3a167f2b2a3b9e1b72bcc54e1948d45``,
native Docker image ``scrimp-mpi:v1.3.0-rc3``).
This longer visual run used ``h=0.1``, ``dt=0.01``
and ``tf=1``, with **8,138 DOFs**, **100/100
accepted fixed Crank–Nicolson steps**, zero
rejections, exact DAE initialization, and
**7/7 original physical acceptance gates passed**.

.. math::

   \mathcal H(0)=0,\qquad
   \mathcal H(1)=0.03707053141235693,\qquad
   W_{\partial}(1)=0.03707053141235698.

The **open conservative** boundary-work balance
``H(t)-H(0)=W_boundary(t)`` holds at each
accepted CN step (maximum absolute defect
``2.07e-16``) and cumulatively (maximum
``1.25e-16``). The *left imposed-velocity
port* supplies the energy; all three other
Neumann boundaries carry exactly zero power.

The original ``rank0.npz`` is accompanied by
**four PVD series and 404 authentic GetFEM VTU
files**, for ``q``, ``p``, ``e_q`` and ``e_p``,
each at **101 times** from 0 to 1.
The renderer generated from the actual FE data:
``visuals/field_snapshots.{png,svg}``,
``visuals/energy_balance.{png,svg}``,
``visuals/p_evolution.gif`` (101 frames),
``visuals/energy_and_ports.csv``, plus a
machine-readable provenance JSON. The animation
shows a boundary-generated disturbance
propagating through the anisotropic rectangle.

**Interpretation of fields matters:**
``q=Grad(w)`` is the *displacement gradient*
(not stress), while ``p=rho*dt(w)`` is *momentum*
(not velocity). The velocity is ``e_p=p/rho``,
and the stress-like energy co-variable is
``e_q=T*q``. The generic maps and GIF currently
display ``p`` and ``|q|``; they must not be
labeled velocity or stress. The original
``e_p`` and ``e_q`` VTK series permit these
physical quantities to be visualized separately.

To reproduce the independently qualified
long MUMPS MPI1 visual run:

.. code-block:: bash

   cd ~/Simulations/scrimp
   git pull --ff-only
   conda activate scrimp
   SCRIMP_DOCKER_NO_SUDO=1 \
   SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
   bash scripts/run-v130-wave2d-showcase.sh 2a

This long-run result is **not an MPI1/MPI2 parity
test or a mesh/time convergence study**. The short
MUMPS MPI1/2 and other solver installation
acceptance remains a separate experiment.
The displacement ``w`` is not an independent
state of this energy-based DAE; any animation
of it must identify its reconstruction
method and gauge.

Scientific references
---------------------

This introductory 2D example illustrates the
PFEM framework analyzed by **Haine, Matignon and
Serhani (2023)**, *Numerical Analysis of a
Structure-Preserving Space-Discretization for an
Anisotropic and Heterogeneous Boundary Controlled
N-Dimensional Wave Equation As a Port-Hamiltonian
System*, DOI `10.4208/ijnam2023-1005`, and
the mixed-boundary-control work of **Brugnoli,
Haine and Matignon (2022)**. The tutorial is
**not** a numerical reproduction of either paper.

The next section introduces a physically distinct
**resistive port** in ``examples/wave_dissipative.py`
while keeping the same material, domain and
conventions. The later :doc:`wave_coenergy`
chapter changes the state variables and adds
spatially localized damping.

2B — Adding a resistive port: heterogeneous damping
----------------------------------------------------

.. automodule:: examples.wave_dissipative
   :members:
   :undoc-members:
   :show-inheritance:

What changes relative to 2A?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This second example introduces **one new physical feature**:
distributed viscous damping of the 2D anisotropic membrane.
The independent, executable source is
`examples/wave_dissipative.py`. Its rectangle,
material tensor, density, normal directions, four boundary
control ports, finite-element spaces for `q` and `p`,
initial conditions and fixed Crank--Nicolson scheme
are *identical* to those of :ref:`Wave_2D`.

On :math:`\Omega=(0,2)\times(0,1)`, the equation now reads

.. math::

   \rho(x)\,\partial_{tt}w
   -\operatorname{div}(T(x,y)\nabla w)
   +\nu(x)\,\partial_t w=0,

with the same symmetric positive-definite tension and mass as
example 2A and the additional nonnegative, spatially varying
damping coefficient

.. math::

   \nu(x)=\tfrac12(2-x)\geq0.

The damping is nonzero throughout the interior except on the
right boundary. We may choose a nonnegative strength multiplier
in the executable; the default is `damping_strength=1.0`.
**Damping does not change the Hamiltonian:**

.. math::

   \mathcal H(q,p)=\tfrac12\int_\Omega
   \left(q\cdot Tq+\frac{p^2}{\rho}\right)\,\mathrm dx.

Why add a separate resistive port?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We retain the *skew-adjoint* interconnection of 2A and
expose the dissipation as a new scalar algebraic port
`(f_r,e_r)`. The state and co-energy variables are still

.. math::

   q=\nabla w,\quad p=\rho\,\partial_t w,\qquad
   e_q=Tq,\quad e_p=p/\rho.

The additional constitutive and coupling equations are

.. math::

   f_r=-e_p,\qquad e_r=\nu f_r,

.. math::

   \begin{pmatrix}\partial_t q\\ \partial_t p\\ f_r\end{pmatrix}
   =
   \begin{bmatrix}
   0&\nabla&0\\
   \operatorname{div}&0&I\\
   0&-I&0
   \end{bmatrix}
   \begin{pmatrix}e_q\\e_p\\e_r\end{pmatrix}.

Indeed, `e_r=-\nu e_p` and
`\partial_t p=\operatorname{div}e_q-\nu e_p`.
The discrete resistive port introduces **two algebraic
fields** `f_r` and `e_r` (both CG2), while
preserving the mechanical variables and boundary ports
from example 2A.

The dissipation power is *positive*:

.. math::
   :label: wave2b_dissipated_power

   D(t)=\int_\Omega f_r\,e_r\,\mathrm dx
       =\int_\Omega \nu\,f_r^2\,\mathrm dx\geq0.

Consequently, for the same boundary work
:math:`W_\partial(t)` as 2A,

.. math::
   :label: wave2b_energy_dissipation

   \mathcal H(t)-\mathcal H(0)
   = W_\partial(t)-\int_0^t D(s)\,\mathrm ds.

A dissipative system can **still gain mechanical energy**
if its boundary ports inject more power than it dissipates.
The correct identity is therefore
:math:`\Delta H=W_\partial-L_\nu`, *not*
the unconditional monotonicity of :math:`H`.

Build the actual GetFEM/PETSc model
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following code fragments are **read directly from**
`examples/wave_dissipative.py`; they are not second,
manually maintained versions of the solver.

**1. Construct the same 2D physical and co-energy states:**

.. literalinclude:: ../../../examples/wave_dissipative.py
   :language: python
   :start-after: # TUTORIAL 2B-01:
   :end-before: # TUTORIAL 2B-02:
   :dedent: 4

**2. Reuse the exact 2A mixed boundary control:**

.. literalinclude:: ../../../examples/wave_dissipative.py
   :language: python
   :start-after: # TUTORIAL 2B-02:
   :end-before: # TUTORIAL 2B-03:
   :dedent: 4

**3. Add precisely one scalar dissipative port:**

.. literalinclude:: ../../../examples/wave_dissipative.py
   :language: python
   :start-after: # TUTORIAL 2B-03:
   :end-before: # TUTORIAL 2B-04:
   :dedent: 4

**4. Reuse DG1/CG2 and add CG2 for damping,
with the same tensor and density:**

.. literalinclude:: ../../../examples/wave_dissipative.py
   :language: python
   :start-after: # TUTORIAL 2B-04:
   :end-before: # TUTORIAL 2B-05:
   :dedent: 4

**5. Keep all 20 original 2A weak-form bricks
and add exactly five resistive ones.**

.. literalinclude:: ../../../examples/wave_dissipative.py
   :language: python
   :start-after: # TUTORIAL 2B-05:
   :end-before: # TUTORIAL 2B-06:
   :dedent: 4

The extra bricks are `M_r`, `I_r`,
`-I_r^T`, `-M_e_r` and `CR_r`.
Their orientations implement
`f_r=-e_p` and `e_r=\nu f_r`.
In particular, adding a negative `e_r` term by hand
to the Hamiltonian would be wrong: the port itself
carries the positive dissipated power in
:eq:`wave2b_dissipated_power`.

**6. Start with the same compatible states and control pulse:**

.. literalinclude:: ../../../examples/wave_dissipative.py
   :language: python
   :start-after: # TUTORIAL 2B-06:
   :end-before: # TUTORIAL 2B-07:
   :dedent: 4

**7. Apply the same strict algebraic DAE initialization
and fixed CN step, without pseudo-time relaxation:**

.. literalinclude:: ../../../examples/wave_dissipative.py
   :language: python
   :start-after: # TUTORIAL 2B-07:
   :end-before: # TUTORIAL 2B-08:
   :dedent: 4

**8. Compute the unchanged physical Hamiltonian:**

.. literalinclude:: ../../../examples/wave_dissipative.py
   :language: python
   :start-after: # TUTORIAL 2B-08:
   :end-before: # TUTORIAL 2B-09:
   :dedent: 4

Physical qualification: original boundary work minus dissipation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For each *accepted* CN step, the program measures the original
GetFEM Hamiltonian and all four port powers on the **mean
finite-element state**. The resistive power is computed by
SCRIMP's native `Port.compute()`, and checked independently
by a second GetFEM volume quadrature of `nu*f_r*f_r`.
The intended stepwise residual is

.. math::

   \mathcal R_k
   =\mathcal H_k-\mathcal H_{k-1}
   -\Delta t_k\sum_{\Gamma}P_{\Gamma,k}
   +\Delta t_k D_k \simeq0.

The script also verifies the cumulative
:eq:`wave2b_energy_dissipation`, the nonnegativity
of the actual dissipated power, the resistive
flow equation, and the consistency of the
11-snapshot accepted FE history. Because
`\nu f_r` has higher polynomial degree than
the CG2 space of `e_r`, the constitutive
relationship is **discrete weak projection**;
a pointwise zero `L^2` error is **not**
an appropriate acceptance criterion.

Run a short, *real* GetFEM/PETSc check (not a
standalone LU smoke):

.. code-block:: console

   python -m examples.wave_dissipative --solver mumps \
       --mesh-size 0.2 --final-time 0.08 --time-step 0.008 \
       --damping-strength 1 --check --export \
       --output outputs/my_wave_2d_dissipative

The native qualification runner below checks **MUMPS MPI1/2**
and all six finite-element exports (`q`, `p`,
`e_q`, `e_p`, `f_r`, `e_r`).
Its physical results must be inspected **before**
declaring this example qualified:

.. code-block:: console

   SCRIMP_DOCKER_NO_SUDO=1 \
   SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
   bash scripts/run-v130-wave2d-dissipative-pilot.sh

Long native dissipative showcase — 100 CN steps
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**PASS, 25 November 2026**, on the authentic
GetFEM/PETSc archive
``v130-wave2d-2b-showcase-dbdc49cf-G1CKQAVG.zip``
(source ``dbdc49cfcac47d4daa1706b3fc2fa4b53646b0a2``,
pinned Docker ``scrimp-mpi:v1.3.0-rc3``).
The actual executable ``examples/wave_dissipative.py``
ran with ``h=0.1``, ``dt=0.01``, ``tf=1`` and
``damping_strength=1.0``, on **MUMPS MPI1**:
**10,204 DOFs**, **100/100 fixed CN steps**,
zero rejected steps or nonlinear failures,
and PASS on exact-algebraic DAE initialisation
(the algebraic residual is zero). All
**10 native physical gates passed**.

.. list-table:: Original long 2B power/energy audit
   :header-rows: 1
   :widths: 60 40

   * - Quantity
     - Measured value
   * - Hamiltonian at ``t=0``
     - ``0``
   * - Hamiltonian at ``t=1``
     - ``0.03205555671004583``
   * - Injected signed boundary work ``W(1)``
     - ``0.03710798721523879``
   * - Cumulative viscous loss ``L(1)``
     - ``0.005052430505192933``
   * - Maximum absolute defect per CN step
     - ``9.39e-17``
   * - Maximum absolute cumulative defect
     - ``9.28e-17``
   * - Native time series ``q,p,e_q,e_p,f_r,e_r``
     - 6 × 101 = **606 VTU**, six PVD

The original FE arrays satisfy the **open
dissipative balance** ``H-H0=W_boundary-L``,
not unconditional monotonicity of ``H``.
The three non-driven Neumann boundaries
carry exactly zero power. The integral
``D=integral f_r*e_r`` is positive and agrees
with an independently assembled GetFEM
volume integral ``integral nu*f_r*f_r``.
The flow equation ``f_r=-e_p`` is satisfied
to machine precision in its FE test.

The archive includes true 101-frame VTK/PVD
series for all six fields, plus energy CSV,
PNG/SVG figures and a 101-frame animation
``p_evolution.gif``. All 101 VTK timestamps
were verified against the original CN states.

2A versus 2B — same physical inputs, different dissipation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The independently inspected 2A and 2B long
archives have the **same rectangle, anisotropic
tensor, density, velocity boundary drive,
mesh size and CN time grid**. The VTK mesh
coordinates and cell connectivity for the
field ``p`` coincide exactly at ``t=0``,
``t=0.5``, ``t=1``. Adding two algebraic
resistive fields increases the model DOFs
from 8,138 (2A) to 10,204 (2B).

.. list-table:: Original FE energy comparison at ``tf=1``
   :header-rows: 1
   :widths: 40 30 30

   * - Quantity
     - Conservative 2A
     - Dissipative 2B
   * - Injected work ``W(1)``
     - ``0.03707053141235698``
     - ``0.03710798721523879``
   * - Dissipated energy ``L(1)``
     - ``0``
     - ``0.005052430505192933``
   * - Remaining mechanical energy ``H(1)``
     - ``0.03707053141235693``
     - ``0.03205555671004583``

**2B has 13.528% less mechanical energy**
than 2A at ``t=1``. The *applied control is
identical*, but the *work delivered* need not
be: it depends on the dynamical response.
Consequently the difference in stored energy
is not exactly the dissipated energy.
The exact comparison identity is

.. math::

   \mathcal H_{2A}-\mathcal H_{2B}
   = L_{2B} - (W_{2B}-W_{2A}),

which matches the real GetFEM values to
``2.08e-17``. This distinction matters when
illustrating dissipation in an **open**
port-Hamiltonian system.

A separate comparison figure of the native
fields ``p`` at ``t=0.5`` and ``t=1`` uses the
**same color scale** for both simulations.
It depicts *momentum* ``p``, not velocity
``e_p``; the original ``e_p`` VTK series
is available for a velocity-specific figure.
The resulting pointwise visual differences
are moderate on this time interval; the
energy deficit is the clearer diagnostic.

To reproduce the long 2B showcase:

.. code-block:: bash

   cd ~/Simulations/scrimp
   git pull --ff-only
   conda activate scrimp
   SCRIMP_DOCKER_NO_SUDO=1 \
   SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
   bash scripts/run-v130-wave2d-showcase.sh 2b

**Scope:** MUMPS MPI1 only at the long showcase
resolution. The short native MUMPS MPI1/MPI2
2B pilot passed separately; this archive
does **not** establish mesh/time convergence,
long-run MPI1/2 parity or SuperLU_DIST
qualification for 2B.

The progressive example **2C**, :doc:`wave_coenergy`,
will then introduce an explicitly different co-energy
formulation, instead of merely renaming these variables.
