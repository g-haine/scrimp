The wave equation — 2C: co-energy and localized damping
=========================================================

.. _Wave_2D_coenergy:

.. automodule:: examples.wave_coenergy
   :members:
   :undoc-members:
   :show-inheritance:

Why this third tutorial?
------------------------

Tutorial :doc:`wave` introduces the **anisotropic 2D membrane and
mixed boundary ports** (2A), then a **separate resistive flow/effort
pair** (2B). This tutorial introduces a genuinely different
construction: the **co-energy states** themselves, a **concentric
two-region domain**, inner damping and an outer boundary feedback.
It does not replace 2A or 2B or imply identical meshes.

The actual, executable code is
``examples/wave_coenergy.py``. Fragments below are included from
that file, rather than copied into the documentation. The short
default example (`h=0.2, dt=0.008, tf=0.08`) is a **native
acceptance pilot**, not a resolved wave-propagation demonstration
or a convergence study.

The true material inverse
-------------------------

The physical tension and density retain the definitions of 2A/2B:

.. math::

   T(x,y)=\begin{bmatrix}5+x&xy\\xy&2+y\end{bmatrix},
   \qquad \rho(x)=3-x.

However, the co-energy formulation needs the actual compliance,
**not** the same numeric tensor renamed ``Tinv``:

.. math::

   T^{-1}(x,y)=\frac{1}{(5+x)(2+y)-(xy)^2}
   \begin{bmatrix}2+y&-xy\\-xy&5+x\end{bmatrix}.

On the unit disk the determinant is strictly positive. The
GetFEM GWFL parameter ``Tinv`` below encodes this inverse;
``q.Tinv.q`` is a tensor contraction, not a Python transpose.

The co-energy unknowns are now **stress** `q=T grad(w)`
and **velocity** `p=dw/dt`. The conjugate co-states have
been substituted, so they do not introduce extra FE unknowns.
The mechanical energy is

.. math::

   \mathcal H(q,p)=\frac12\int_\Omega
   \left(q\cdot T^{-1}q+\rho p^2\right)\,\mathrm dx.

Two regions and substituted ports
---------------------------------

The geometry is a unit disk split into the inner disk
`\Omega_1` (radius 0.6, region 1) and an annulus
`\Omega_2` (region 2). Region 10 is their interface;
region 20 is the external boundary. The damping port
is present **only** on region 1 and the boundary port
**only** on region 20.

.. literalinclude:: ../../../examples/wave_coenergy.py
   :language: python
   :start-after: # TUTORIAL 2C-01:
   :end-before: # TUTORIAL 2C-02:

.. literalinclude:: ../../../examples/wave_coenergy.py
   :language: python
   :start-after: # TUTORIAL 2C-02:
   :end-before: # TUTORIAL 2C-03:

.. literalinclude:: ../../../examples/wave_coenergy.py
   :language: python
   :start-after: # TUTORIAL 2C-03:
   :end-before: # TUTORIAL 2C-04:

In the internal disk, the resistive port is substituted:
`nu^{-1} e_r = -p` in its finite-element weak form.
Unlike a non-substituted port, it has no independent
GetFEM flow variable. Consequently its power **cannot**
be obtained through `Port.get_power()`. The native audit
instead evaluates the original FE functions by independent
GetFEM quadrature.

The legacy coefficient `nu=10(0.36-r^2)` vanished at
the interface even though a brick integrates `1/nu`.
Here we deliberately use
`nu=0.2+10(0.36-r^2)` **on region 1 only**.
This is a strictly positive, bounded invertible
coefficient on the disk and still has exactly zero
damping on the outer annulus by region restriction.
It is a documented change to the legacy scenario,
not a numerical regularization hidden in SCRIMP.

.. literalinclude:: ../../../examples/wave_coenergy.py
   :language: python
   :start-after: # TUTORIAL 2C-04:
   :end-before: # TUTORIAL 2C-05:

The co-energy balance and feedback sign
--------------------------------------

The semidiscrete system is

.. math::

   T^{-1}\dot q=\nabla p,\qquad
   \rho\dot p=\operatorname{div}q+e_r
   \quad (e_r\text{ acts only on }\Omega_1),

with boundary observation `Y=p` and externally applied
traction `U`. The negative output feedback `U=-0.5Y`
is physically dissipative: its actual signed boundary
power is `P_b=integral UY <= 0`.
The original `U=+0.5Y` has the *opposite* sign and
must not be labeled feedback dissipation.

For the two dissipative mechanisms, define

.. math::

   D_i=\int_{\Omega_1}\frac{e_r^2}{\nu}\,\mathrm dx
      =-\int_{\Omega_1}p e_r\,\mathrm dx\ge0,
   \qquad
   D_b=-\int_{\partial\Omega} UY\,\mathrm ds\ge0.

The equality for `D_i` holds **in the FE weak
projection**. The pointwise integral `integral nu*p^2`
need not agree exactly, because `e_r` is DG1, and
`nu*p` generally is not. We record it as a diagnostic,
not a false exactness gate.

There is no independent boundary drive in this scenario;
therefore the physical balance to validate is

.. math::

   \boxed{\mathcal H(t)-\mathcal H(0)+
   \int_0^t (D_i(s)+D_b(s))\,\mathrm ds=0.}

.. literalinclude:: ../../../examples/wave_coenergy.py
   :language: python
   :start-after: # TUTORIAL 2C-05:
   :end-before: # TUTORIAL 2C-06:

A physically visible initial pulse
----------------------------------

Two smooth initial velocity bumps excite the inner
damping and outer feedback. The initial stress is zero.
The time integrator is **fixed Crank–Nicolson**, ten
accepted steps, with algebraic DAE initialization
and default MUMPS LU. The short time interval proves
the simulation and output pipeline works, but does
not by itself establish long-time spatial quality.

.. literalinclude:: ../../../examples/wave_coenergy.py
   :language: python
   :start-after: # TUTORIAL 2C-06:
   :end-before: # TUTORIAL 2C-07:

.. literalinclude:: ../../../examples/wave_coenergy.py
   :language: python
   :start-after: # TUTORIAL 2C-07:
   :end-before: # TUTORIAL 2C-08:

.. literalinclude:: ../../../examples/wave_coenergy.py
   :language: python
   :start-after: # TUTORIAL 2C-08:
   :end-before: # TUTORIAL 2C-09:

Native outputs and visualizations
---------------------------------

The pilot below records the original **full FE state
history** in NPZ, true GetFEM Hamiltonian and power
histories, detailed JSON physical gates, two independent
quadratures of inner resistance, and GetFEM/ParaView
time series `q`, `p`, `e_r` (one PVD and eleven VTU
per field and MPI case). It checks strict DAE, CN
work/loss balance at every step, nontrivial passive
feedback and damping, MPI replica/parity, and genuine
nonempty VTK fields.

Once a native case passes, the same script renders
**from the original FE NPZ/VTK**, not an invented
analytic image:

- ``visuals/energy_balance.png`` and ``.svg``:
  mechanical energy vs work/losses, powers, balance defect;
- ``visuals/field_snapshots.png`` and ``.svg``:
  spatial snapshots of `p`, `|q|`, `e_r`;
- ``visuals/p_evolution.gif``: animated accepted `p`
  snapshots (pilot: 11 frames);
- ``visuals/energy_and_ports.csv``: actual energetic
  histories with explicit columns;
- ``visuals/provenance.json``: commit, image revision,
  number of frames and Float32 VTK warning.

Open ``mumps_mpi1/pv/p/p.pvd`` (or q/e_r) in
ParaView to examine the whole native sequence. VTK
Float32 is for visualisation; the acceptance audit
uses full-precision original GetFEM FE values.

.. code-block:: bash

   cd ~/Simulations/scrimp
   git pull --ff-only
   conda activate scrimp

   SCRIMP_DOCKER_NO_SUDO=1 \
   SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
   bash scripts/run-v130-wave2d-coenergy-pilot.sh

The same native visualization tool can be run on
the already-qualified **2B** archive, without any
rerun of the solver:

.. code-block:: bash

   python scripts/render-v130-wave2d-native.py \
     --input sandbox/outputs/v130-wave2d-dissipative-e582874f-9ni2cPYZ.zip \
     --case mumps_mpi1 \
     --output sandbox/outputs/wave2d-2b-visuals

The renderer also supports compatible archived 2A
NPZ/PVD layouts. Use the case names in each archive.
A later long-time demonstration should be assessed
separately from this short physical acceptance,
and its validated images may subsequently be integrated
into the public Sphinx website.

**Qualification status:** 2C is a prepared native
candidate, **NOT PASS until real MUMPS MPI1/MPI2
logs, binary evidence and physical gates succeed**.
No generic solver regularization or hidden
SuperLU_DIST qualification is claimed.
