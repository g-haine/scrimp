Maxwell equations in 3D — minimal Nédélec / Raviart–Thomas example
==================================================================

.. _Maxwell_3D:

Why a new Maxwell example?
--------------------------

The wave tutorials :doc:`wave` and :doc:`wave_coenergy` work
in two spatial dimensions. This is an intentionally **tiny
genuine 3D test**: a homogeneous sphere, just two vector
unknowns, no sources, no feedback controllers, and four
weak-form bricks. It exercises the corrected **intrinsic
vector finite elements** in the pinned SCRIMP GetFEM 5.5
image, not ordinary componentwise CG vector elements.

This example is **new**, independent of the published
``COMPEL_2025/Maxwell_2D.py`` research script (which is
not modified), and of the 2D Nédélec-independent
wave benchmark.

A normalized Maxwell system
---------------------------

Let :math:`\Omega` be the unit ball in
:math:`\mathbb R^3`, filled with a homogeneous,
isotropic, lossless material with **dimensionless**
:math:`\varepsilon=\mu=1`. The electric field
:math:`E` and the magnetic induction :math:`B` satisfy

.. math::

   \partial_t E = \operatorname{curl} B,\qquad
   \partial_t B = -\operatorname{curl} E,

without free currents or sources. Initially,

.. math::

   E(0,x,y,z)=(-y,x,0),\qquad B(0)=0.

The electric field has nonzero curl and
zero divergence initially. It is an exact
low-order edge-element-type polynomial and
excites a nontrivial magnetic field.

The true 3D domain and both co-energy/state
declarations come directly from the executable:

.. literalinclude:: ../../../examples/maxwell_3d.py
   :language: python
   :start-after: # TUTORIAL M3D-02:
   :end-before: # TUTORIAL M3D-03:
   :dedent: 4

.. literalinclude:: ../../../examples/maxwell_3d.py
   :language: python
   :start-after: # TUTORIAL M3D-03:
   :end-before: # TUTORIAL M3D-04:
   :dedent: 4

The coupling uses the **homogeneous natural
boundary closure** of the weak form, namely
zero magnetic tangential boundary work, not an
explicitly imposed electric PEC trace and not
a driven port. The discrete model closes the
boundary power balance by construction. This
small demonstration does not claim strong
boundary-trace or electromagnetic eigenmode
accuracy.

Compatible 3D spaces and their GetFEM patches
---------------------------------------------

We use the finite element de Rham sequence
on actual tetrahedra:

.. math::

   E_h \in \mathrm{N\acute ed\acute elec}_0
        \subset H(\mathrm{curl};\Omega),\qquad
   B_h \in \mathrm{RT}_0
        \subset H(\mathrm{div};\Omega).

**Nédélec** edge DOFs encode tangential continuity.
**Raviart–Thomas** face DOFs encode normal continuity.
Both are intrinsic vector shape functions in
GetFEM, not three unrelated scalar fields.

.. literalinclude:: ../../../examples/maxwell_3d.py
   :language: python
   :start-after: # TUTORIAL M3D-04:
   :end-before: # TUTORIAL M3D-05:
   :dedent: 4

A strict post-registration guard retains both
original ``S.FEM`` objects and verifies that
their **type and order are unchanged** after
``DPHS.add_FEM()``, and that the actual GetFEM
``MeshFem`` has vector dimension 3 and
nonzero degrees of freedom. This is important:
SCRIMP's existing FEM adapter *mutates the
original type/order to CG/1* on GetFEM's
``AssertionError`` fallback. A silent
replacement is thus explicitly detected,
without relying on the Python binding's
version-dependent ``MeshFem.fem()`` layout.
This runtime guard is **not yet a complete
proof of Maxwell numerical correctness**.

A skew-adjoint discrete Maxwell interconnection
------------------------------------------------

The weak formulation is deliberately minimal:

.. math::

   (\partial_t E_h,v_h)_{L^2(\Omega)}
      &= (B_h,\operatorname{curl}v_h)_{L^2(\Omega)},\\
   (\partial_t B_h,w_h)_{L^2(\Omega)}
      &= -(\operatorname{curl}E_h,w_h)_{L^2(\Omega)}.

The GetFEM generic weak-form language does
not supply a portable built-in 3D ``Curl`` in
this image. The canonical example defines
the operator explicitly from the 3×3 gradient
tensor, **without** imposing an ad hoc
stencil or swapping the curl orientation:

.. literalinclude:: ../../../examples/maxwell_3d.py
   :language: python
   :start-after: # TUTORIAL M3D-01:
   :end-before: def maxwell_3d(

There are only **four** variational bricks:
two mass blocks and the two mutually
opposite curl couplings. They are
skew-adjoint in the discrete GetFEM
mass inner products.

.. literalinclude:: ../../../examples/maxwell_3d.py
   :language: python
   :start-after: # TUTORIAL M3D-05:
   :end-before: # TUTORIAL M3D-06:
   :dedent: 4

The original electromagnetic Hamiltonian is

.. math::

   \mathcal H(t) = \frac12\int_{\Omega}
      (|E_h|^2+|B_h|^2)\,\mathrm dx.

The semidiscrete identity follows by
testing the first equation with ``E_h``
and the second with ``B_h``:

.. math::

   \frac{\mathrm d\mathcal H}{\mathrm dt}
   =(B_h,\operatorname{curl}E_h)
    -(\operatorname{curl}E_h,B_h)=0.

The initial data and CN integration
are explicitly defined in the runnable
source, including the absence of any
unnecessary algebraic DAE projection:

.. literalinclude:: ../../../examples/maxwell_3d.py
   :language: python
   :start-after: # TUTORIAL M3D-06:
   :end-before: # TUTORIAL M3D-07:
   :dedent: 4

The original two positive electromagnetic
energy terms are likewise the canonical
GetFEM Hamiltonian, not recomputed from
an invented displacement field:

.. literalinclude:: ../../../examples/maxwell_3d.py
   :language: python
   :start-after: # TUTORIAL M3D-07:
   :end-before: def measure_maxwell
   :dedent: 4

The **fixed Crank–Nicolson** scheme has
the corresponding *exact quadratic*
discrete identity:

.. math::

   \mathcal H^{n+1}-\mathcal H^n
    =\Delta t\left[
      (B_h^{n+1/2},\operatorname{curl}E_h^{n+1/2})
     -(\operatorname{curl}E_h^{n+1/2},B_h^{n+1/2})
     \right]=0.

Both variables are differential, so
there is **no algebraic DAE initial
projection** or pseudo-time spin-up.
Initial data are projected into their
exact native finite-element spaces.

How to run and qualify this 3D example
--------------------------------------

We separate a **source contract** from
real native numerical acceptance.
The qualification runner uses the same
pinned GetFEM/PETSc/MUMPS Docker image
as the 2D examples. It checks
10 fixed CN steps at ``h=0.7``,
``dt=0.02``, ``t_f=0.2``,
on **MPI 1 and MPI 2**:

.. code-block:: bash

   cd ~/Simulations/scrimp
   git pull --ff-only
   conda activate scrimp

   SCRIMP_DOCKER_NO_SUDO=1 \
   SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
   bash scripts/run-v130-maxwell3d-pilot.sh

The runner will only print
``PASS_MAXWELL3D_MUMPS_MPI1_MPI2`` if all
native checks succeed. It saves even
failed traces to
``sandbox/outputs/v130-maxwell3d-*.zip``.

The **physical gates** include the
presence of a genuine 3D mesh,
real Nédélec/RT0 spaces, exact ten-step
fixed CN integration, positive initial
electromagnetic energy, nontrivial magnetic
excitation, and agreement between GetFEM
Hamiltonian values and independent
electric/magnetic volume quadratures.
Every-step and cumulative energy defects
must be below ``1e-9`` times the
physical energy scale. All MPI rank traces
and MPI1/MPI2 results are compared.
The divergence of ``B`` is independently
measured by GetFEM's squared L2 norm
and **checked** against the compatible
discrete Gauss law. Its initial value is
zero and, for matching Nédélec/RT0 spaces,
the discrete curl maps into the
divergence-free RT subspace. This is an
actual additional acceptance gate,
**not** a claim of convergence of all
Maxwell constraints. The electric Gauss
law and electromagnetic boundary-mode
accuracy are outside this minimal pilot.

After a native PASS, the archive contains
true complete FE histories (``rank*.npz``),
original energy/diagnostic values, JSON
acceptance metadata, and two **genuine
3D ParaView PVD time series** for ``E``
and ``B`` (11 VTU snapshots each).
Open ``mumps_mpi1/pv/E/E.pvd`` and
``mumps_mpi1/pv/B/B.pvd`` in ParaView;
view the 3D vectors with a cut plane,
glyphs, or streamlines.

**Qualification status on preparation:
native pilot not yet executed.**
Passing static Python tests does not
establish correctness of patched GetFEM
3D assembly, physical energy balance,
MPI compatibility, VTK export,
mesh/time convergence or boundary
eigenmodes. Those claims require the
real archived native evidence.
