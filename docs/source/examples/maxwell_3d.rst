Maxwell equations in 3D — patched Nédélec-II / BDM1 example
============================================================

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

We use the **SCRIMP-patched** pair on genuine tetrahedra:

.. math::

   E_h \in \mathrm{N\acute ed\acute elec}_{\mathrm{II},1}
        \subset H(\mathrm{curl};\Omega),\qquad
   B_h \in \mathrm{BDM}_1
        \subset H(\mathrm{div};\Omega).

On each tetrahedron, both fields have local polynomial
shape functions in :math:`[P_1]^3`; the Nédélec-II
edge DOFs impose tangential continuity, while BDM1
face DOFs impose normal continuity. In particular
:math:`\operatorname{curl} E_h` is piecewise
constant and belongs to the divergence-free
subspace of BDM1, compatible with Maxwell's
magnetic Gauss law.

**Important distinction:** the original accepted
baseline at source commit ``daac18cb`` used
GetFEM's *built-in* lowest Nédélec
``FEM_NEDELEC(3)`` with built-in RT0. It
did **not** exercise SCRIMP's new element patches.
The current tutorial now explicitly uses
``FEM_NEDELEC2(3,1)`` (patch ``0004``) and
tetrahedral ``FEM_BDMK(3,1)`` (patch ``0001``).

These FE names use the new SCRIMP shortcuts
``S.FEM("E", 1, "NED2")`` and
``S.FEM("B", 1, "BDM")``. Additional
shortcuts ``NED1`` / ``NEDELEC1`` and
``NED2`` / ``NEDELEC2`` are mapped by
``scrimp/fem.py`` to the corresponding
GetFEM constructors. The number in
``NED1`` or ``NED2`` is the **family
(1st/2nd kind)**, while the numeric argument
is the **polynomial order k**. Patch 0004
supports only k=1 or 2.

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
real patched Nédélec-II/BDM1 spaces, exact ten-step
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
zero and, for the compatible Nédélec-II k1
and BDM1 spaces, the discrete curl maps into
the divergence-free BDM subspace. This is an
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

Native patched Nédélec-II / BDM1 qualification — PASS, 25 November 2026
-------------------------------------------------------------------------

**The patched pair, not just the historical native-element baseline,
has now passed the complete native Maxwell 3D pilot.**
The original archive is
``v130-maxwell3d-3687b11e-gSQ8jeC3.zip``,
source
``3687b11e473813c103cd4e9c02f82812677853d7``,
Docker ``scrimp-mpi:v1.3.0-rc3``.
The ``summary.json`` status is
``PASS_MAXWELL3D_MUMPS_MPI1_MPI2``,
with an empty error list.
Both MPI configurations used **genuine**
``FEM_NEDELEC2(3,1)`` for ``E``
and ``FEM_BDMK(3,1)`` for ``B``.

.. list-table:: Native patched Maxwell 3D qualification
   :header-rows: 1
   :widths: 54 23 23

   * - Check
     - MUMPS MPI1
     - MUMPS MPI2
   * - Tetrahedra / VTK points
     - **110 / 58**
     - **110 / 58**
   * - Nédélec-II k1 ``E`` DOFs
     - **444**
     - **444**
   * - BDM1 ``B`` DOFs
     - **825**
     - **825**
   * - Total FE DOFs
     - **1,269**
     - **1,269**
   * - CN steps / physical gates
     - **10/10, 9/9 PASS**
     - **10/10, 9/9 PASS**
   * - Maximum cumulative energy defect
     - ``4.1078e-15``
     - ``8.8818e-16``
   * - Maximum L2 norm of div(B)
     - ``2.0531e-14``
     - ``2.0387e-14``
   * - Original native VTU and PVD
     - **22 VTU, 2 PVD**
     - **22 VTU, 2 PVD**
   * - Automatic PNG/SVG/GIF
     - **PASS**
     - **PASS**

At final time, original GetFEM quadratures
give
``H_E(0.2)=0.436003804787049``
and
``H_B(0.2)=0.25824828616211315``
for MPI1. The total electromagnetic
Hamiltonian is conserved with a maximum
defect ``4.11e-15``, while about
**37.20%** of the initial electric energy
is stored magnetically at ``t=0.2``.
The independent GetFEM magnetic-divergence
integral :math:`\int_{\Omega}(\operatorname{div} B)^2\,\mathrm dx` remains
at most ``4.2153e-28`` on MPI1.

The maximum coefficient difference between
full MPI1 and MPI2 FE histories
(11 times, 1,269 unknowns) is
``2.1094237467877974e-15``.
Each MPI2 rank saved an identical
trajectory, with its SHA256 verified.
**All five graphics** are generated from
the *same run* by
``examples/maxwell_3d.py --export``,
including the **11-frame native magnetic
GIF**. The integrated renderer is
therefore **qualified in MPI1 and MPI2**,
not just source-tested. The true vector
VTU exports and original time samples
have been independently inspected.

The earlier source commit ``daac18cb``
was a valid **separate baseline** using
GetFEM native Nédélec0/RT0 (497 DOFs),
not a test of either patched family.
Neither the new or the old qualification
establishes mesh/time convergence,
electric Gauss-law validity, explicit PEC
boundary conditions, cavity eigenmodes,
or SuperLU_DIST acceptance.

Figures and animation generated by the example itself
---------------------------------------------------

The **same Python executable** now creates its
own scientific figures and animation whenever
``--export --output /artifacts`` is requested.
No separate notebook, rendering script or
``vtk`` Python package is needed: the existing
``matplotlib`` stack and the standard library
read the **original GetFEM 3D VTU/PVD**.
The energy curves are drawn from the
independently assembled native GetFEM
``H_E``, ``H_B`` and ``H`` arrays, not from
resampled nodal magnitudes.

For a complete run:

.. code-block:: bash

   python -m examples.maxwell_3d --check --export --output ./maxwell_results

The output folder contains:

- ``pv/E/E.pvd`` and ``pv/B/B.pvd``:
  original GetFEM vector time series, with
  one native VTU per accepted CN state.
- ``visuals/maxwell3d_echange_energie.png``
  and ``.svg``: electric, magnetic and
  total Hamiltonian over the original
  FE time samples.
- ``visuals/maxwell3d_vecteurs_E_tfinal.png``
  and ``maxwell3d_vecteurs_B_tfinal.png``:
  3D final-time arrows at the authentic
  exported VTK nodal coordinates.
- ``visuals/maxwell3d_B_3D_animation.gif``:
  all saved states of the *true* magnetic
  induction, with fixed camera and
  fixed vector scale (including the
  initially zero field).
- ``visuals/provenance.json`` and ``README.txt``:
  original FE trajectory SHA256,
  true tetrahedral/point counts,
  number of frames and scientific
  rendering caveats.

**Important numerical distinction:**
the VTK arrows are visualizations of
the native GetFEM vector *field values
at nodes*, not the edge or face degrees
of freedom of the **Nédélec-II/BDM1** spaces.
The mesh outline is extracted from the
true boundary faces of the tetrahedra;
no synthetic field, reconstructed
displacement or time interpolation is used.

The ParaView exports still involve all
MPI ranks, but **only rank 0** creates
PNG/SVG/GIF files. The native pilot
verifies the generated files, the entire
GIF frame count and provenance against
the original FE states; any rendering
failure is retained as an explicit
unqualified pilot rather than hidden.

**Qualification scope:** the original
Nédélec0/RT0 physics baseline passed
at ``daac18cb``. Both the **current
patched Nédélec-II/BDM1 model and its
integrated plotting** passed separately
at ``3687b11e`` on MPI1 and MPI2.
These are distinct archived proofs.

Historical *baseline* Maxwell 3D qualification — PASS, 25 November 2026
------------------------------------------------------------------------

**This historical PASS belongs exclusively to
the earlier built-in Nédélec0 + RT0 pair.**
The **current** patched Nédélec-II k1 + BDM1
variant and its integrated visual outputs
were **separately accepted** in the native
archive ``v130-maxwell3d-3687b11e-gSQ8jeC3.zip``.
Do not interchange their FE histories or DOF counts.

The earlier GetFEM/PETSc/MUMPS MPI1 and MPI2 pilot
**passed nine out of nine physical gates**,
using the original two-field Nédélec/RT0
model at ``h=0.7``, ``dt=0.02``, ``tf=0.2``.
The full evidence is the archive
``v130-maxwell3d-daac18cb-FCPngFrH.zip``,
source commit
``daac18cb36fbc404fcf5c4fe4f2f69121b2ee530``,
Docker image ``scrimp-mpi:v1.3.0-rc3``
(revision
``1d6d0c90115213298047b28237a7f38087563a35``).

.. list-table:: Original Maxwell 3D native evidence
   :header-rows: 1
   :widths: 52 24 24

   * - Check
     - MUMPS MPI1
     - MUMPS MPI2
   * - Real tetrahedra
     - **110**
     - **110**
   * - Baseline DOFs, native Nédélec0 E / RT0 B
     - **222 / 275**
     - **222 / 275**
   * - CN steps accepted
     - **10/10**
     - **10/10**
   * - Scientific gates
     - **9/9 PASS**
     - **9/9 PASS**
   * - Hamiltonian initial
     - ``0.6942520909491595``
     - ``0.6942520909491599``
   * - Hamiltonian final
     - ``0.6942520909491590``
     - ``0.6942520909491603``
   * - Final magnetic energy
     - ``0.2582694531373567``
     - ``0.2582694531373592``
   * - Maximum cumulative energy defect
     - ``2.33e-15``
     - ``2.00e-15``
   * - Maximum L2 norm of magnetic divergence
     - ``4.63e-15``
     - ``4.44e-15``

The two **full FE state histories** have
11 time samples × 497 coefficients.
The maximum absolute coefficient
difference MPI1/MPI2 is
``4.440892098500626e-16``. The MPI2 rank0
and rank1 binary histories are identical.
The saved state SHA256 values agree with
the original metadata.

**Electromagnetic energy transfer:** from
``H_E(0)=0.6942520909491595``, ``H_B(0)=0``,
we obtain ``H_E(0.2)=0.4359826378118022``
and ``H_B(0.2)=0.2582694531373567``,
so **37.20% of the original energy moves
from electric to magnetic storage**, with
no dissipative loss. The independent GetFEM
volume quadratures of electric and magnetic
energy add to ``H`` to within
``5.55e-17`` on MPI1; the full-step and
cumulative total Hamiltonian remains
constant at machine precision.

The magnetic Gauss law is checked by an
independent real GetFEM volume integral
of ``(div B)^2``, with
``||div B||_L2 ≤ 4.7e-15`` over all saved
times. This is a tested constraint for
the *discrete* field, not a claim that
the electric Gauss law was tested.

**VTK verification beyond XML:** all
**44 VTU files** (2 fields × 11 times ×
2 MPI configurations) were opened with
a genuine VTK unstructured-grid reader.
Every file contains 58 3D points,
110 native tetrahedral cells (VTK type 10)
and a finite vector field with three
components. Mesh geometry matches
across the eleven times and the two
MPI configurations. The 4 PVD catalogs
contain the corresponding 11 timestamps.
At final time the nodal magnetic-vector
norm ranges from about 0.367 to 0.375.

**Limitations:** these results qualify
only the *original native Nédélec0/RT0
3D Maxwell port-Hamiltonian baseline* at one
mesh size and one time step.
It establishes neither mesh/time convergence,
electric Gauss-law accuracy, explicit PEC
boundary traces, cavity-mode correctness
nor a SuperLU_DIST qualification. The
boundary remains the documented natural
zero-work closure, not an imposed PEC wall.

