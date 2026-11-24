.. _saddle_point_solver_guide:

Mixed pH-DAEs: saddle points and sparse LU solvers
=================================================

Understanding the issue
-----------------------

Mixed finite-element formulations and algebraic ports can produce
**saddle-point** Jacobians containing an algebraic zero block:

.. math::

   \begin{bmatrix}A&B^\top\\ B&0\end{bmatrix}.

A zero diagonal does **not** imply a singular operator. Valid pivots
can occur off diagonal, but sparse direct solvers may respond differently
to the same matrix depending on permutation, pivoting and refinement.

**SCRIMP v1.3.0 does not automatically diagnose or repair these
factorizations.** MUMPS remains its default direct LU backend. Any
solver-specific options must be selected deliberately and verified
against the original numerical solution and the model's power balance.

Heat--Wave coupling: a related, but distinct case
-------------------------------------------------

The maintained coupled example is :doc:`examples/heat_wave`
(source: `examples/heat_wave.py`). It combines a heat flux variable,
wave variables, algebraic interface constraints, and a boundary port.
Its mixed formulation can exhibit the same *kind* of saddle-point
structure as the 2D wave model described below.

A previous **sequential SuperLU** experiment on Heat--Wave was
reported by the author as being resolved by **permutation/reordering**.
This is consistent with a pivot-ordering problem. However, the
**exact PETSc permutation setting and original matrix/factor evidence
for that historical run are not recorded in the current maintained
source**. It would be unjustified to present a particular option as
a newly qualified v1.3.0 remedy for this example.

Two other changes in the current Heat--Wave model are unrelated to
LU pivot ordering and should be kept separate:

* Switching `heat_region=1, wave_region=2` to the reverse assignment
  changes the **physical heat/wave subdomains**. It is not an
  algebraic row/column permutation performed by SuperLU.
* Correcting interface outward normals, boundary causality,
  variational signs and the mixed DG0/RT0/CG1 spaces is necessary
  for **physical consistency**, not a numerical pivot workaround.

The example supports optional `factor_solver` selection; without
it SCRIMP uses MUMPS. Heat--Wave's energy dissipation and interface
power checks need their own scientific qualification, regardless of
whether a factorization succeeds. See
`docs/V130_HEAT_WAVE_PHYSICAL_AUDIT.md` for the physical corrections.

Conservative 2D wave: exact matrix evidence
-------------------------------------------

See :doc:`examples/wave` and the canonical `examples/wave.py`.
Prescribed left-boundary velocity introduces the algebraic traction
variable `Y_L`. The first real CN factor matrix has size 2230 by
2230 and equals

.. math::

   P=K+\sigma M,\qquad \sigma=250,\quad t=0.008.

The native first-step matrix capture demonstrated:

* 25,801 stored entries, 25,791 numerically nonzero,
  no empty rows and no nonfinite coefficients.
* Exactly ten zero diagonals (indices 80--89), **all belonging
  to the algebraic variable Y_L**.
* Full **structural rank 2230/2230**, and **exactly the same**
  P, K and M matrices for MUMPS and SuperLU_DIST.
* Independent SciPy/SuperLU LU factorization **without
  regularization**, with residuals around :math:`10^{-15}`.

Thus the ten zero diagonals do not signify missing equations
or intrinsic singularity in this particular matrix.

With the default SuperLU_DIST factorization, the actual PETSc
failure at the first step was:

.. code-block:: text

   KSP_reason=-11               # KSP_DIVERGED_PC_FAILED
   PC_failed=2                  # PC_FACTOR_NUMERIC_ZEROPIVOT
   factor_solver=superlu_dist

**Permutation alone did not fix this case.**
The tested column permutation `MMD_AT_PLUS_A` retained the zero
pivot. The tested row permutation `LargeDiag_AWPM` entered a
CombBLAS-dependent path unavailable in the pinned Docker build
and segfaulted (exit 139). This failing AWPM setting must
not be recommended for that build.

Setting only `mat_superlu_dist_replacetinypivot=true` eliminated
the numeric zero-pivot failure in the isolated PETSc LU replay,
but the **real backward residual** remained about
:math:`5.08\cdot10^{-10}` (required below
:math:`10^{-11}`); right-preconditioned GMRES also failed
that residual threshold.

**The configuration actually qualified for this specific
conservative 2D example** is:

.. code-block:: text

   pc_type=lu
   pc_factor_mat_solver_type=superlu_dist
   mat_superlu_dist_replacetinypivot=true
   mat_superlu_dist_iterrefine=true

Both options are **explicit, opt-in PETSc choices**. Pivot replacement
modifies the LU factorization numerically; internal iterative refinement
recovers accurate solutions. On the *same archived first-CN matrix*,
the native PETSc replay on MPI 1/2 passed with backward residuals
:math:`9.54\cdot10^{-17}` and
:math:`1.27\cdot10^{-16}`, and known-solution errors
below :math:`2.3\cdot10^{-14}`.

The **complete unchanged model** was then qualified with these
options, not just the LU factor: ten accepted fixed CN steps,
zero rejected steps, strict algebraic DAE initialization, all
four physical boundary powers, anisotropic Hamiltonian,
eleven full finite-element snapshots, and eleven real ParaView
frames of each field `q`, `p`, `e_q`, `e_p` on
MPI1 and MPI2. The full FE trajectories, boundary powers and
VTK geometry agree with the previously qualified MUMPS MPI1
reference.

.. list-table:: Native full-model result with opt-in refined SuperLU_DIST
   :header-rows: 1
   :widths: 52 24 24

   * - Measurement
     - MPI1
     - MPI2
   * - Maximum per-step energy-work defect
     - :math:`3.9031\cdot10^{-18}`
     - :math:`1.8431\cdot10^{-18}`
   * - Maximum cumulative energy-work defect
     - :math:`6.0715\cdot10^{-18}`
     - :math:`3.4694\cdot10^{-18}`
   * - Maximum full-state difference vs. MUMPS MPI1
     - :math:`7.3830\cdot10^{-14}`
     - :math:`6.9500\cdot10^{-14}`
   * - Original VTU frames per field
     - 11
     - 11

Native evidence:
`v130-wave2d-refined-dist-physical-95fbd8b2-HfwaSFgK.zip`,
source commit `95fbd8b`. The reproducible opt-in
test runner invokes the original example without changing its
physics or SCRIMP:

.. code-block:: console

   SCRIMP_DOCKER_NO_SUDO=1 \
   SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
   bash scripts/run-v130-wave2d-refined-superlu-dist-physical.sh

The pilot injects the two options **after** SCRIMP constructs
the DPHS and **before** PETSc reads the solver options. That
order matters because creating a new DPHS clears the PETSc
options database in this version.

**Scope and status:** the explicitly **refined** SuperLU_DIST
configuration is validated here on MPI1/2. The five original
unmodified solver choices still count as **3/5 PASS**:
MUMPS MPI1/2 and sequential SuperLU MPI1 pass; default
SuperLU_DIST MPI1/2 fails on this case. A qualified **opt-in**
option pair does not imply general-purpose backend reliability
and does not change the library's default.

How to diagnose a similar failure
---------------------------------

#. First verify physical and variational correctness: FE spaces,
   boundary causality, normal orientation, constitutive laws and
   consistent initial algebraic constraints.
#. Read **TS/SNES/KSP/PC** convergence and failure reasons and confirm
   which LU backend actually ran. A PETSc "converged" status alone
   does not prove a sufficiently accurate solution.
#. Inspect the real factor matrix (sparsity, nonfinite values,
   empty rows, zero diagonals, structural rank). If possible,
   compare coefficients with a separately qualified backend.
   Distinguish a zero diagonal from a true singularity.
#. Consider a *documented, explicitly selected* solver-specific
   remedy: reorder first when there is evidence it helps; use pivot
   replacement and refinement only if the original-system
   residual and the solution accuracy justify them.
#. Finally, **requalify the original physical problem**:
   full FE trajectories, MPI consistency, Hamiltonian and
   signed boundary power balance, dissipation where relevant,
   and original ParaView exports.

These observations are specific to the actual matrix, solver
version, numerical precision and linked PETSc/SuperLU build.
They are **not** an automatic strategy-selection algorithm for
arbitrary pH systems, and no such feature is planned for v1.3.0.

Official PETSc references:

* `SuperLU, sequential <https://petsc.org/release/manualpages/Mat/MATSOLVERSUPERLU/>`_
* `SuperLU_DIST, distributed <https://petsc.org/release/manualpages/Mat/MATSOLVERSUPERLU_DIST/>`_
