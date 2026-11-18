**********************
Solver reuse and MPI
**********************

v1.3.0 development candidate
============================

For v1.3.0 SCRIMP focuses on MPI performance without changing the
port-Hamiltonian model, the PFEM discretization or the Crank--Nicolson
time integrator.

The default MUMPS reuse policy is ``on_change``. SCRIMP examines the
*actual assembled* PETSc preconditioner matrix on every rank:

* If the stored CSR graph (ownership range, row pointers or column
  indices) changes, MUMPS must recompute the symbolic factor.
* If the graph is unchanged but any stored numeric coefficient differs,
  only the numerical factor must be recomputed.
* If graph and numerical coefficients remain identical on every rank,
  the numerical LU can be reused.

This comparison is performed after construction of the actual Jacobian,
including the TS shift applied to the mass matrix and all nonlinear
tangent contributions. Linearity does not imply time invariance: an
adaptive time step can change the operator.

``SCRIMP_MUMPS_RESET_POLICY`` remains configurable through
``on_change`` (new default), ``always``, and ``never``. The older
``never`` merely avoids explicit PC resets; it does **not** guarantee
numerical factor reuse.

For MPI LU/MUMPS on more than one rank, SCRIMP also defaults to
``ICNTL(28)=2`` and ``ICNTL(29)=0`` (parallel symbolic analysis).
Explicit PETSc options always take precedence, and a one-rank run
retains the sequential analysis default.

**Qualification status:** a previous experimental (not production)
``static_cn`` benchmark demonstrated one instead of ten numerical
factorizations on ten fixed-CN wave steps, with matching state-block
norms. The new production ``on_change`` mechanism is different:
it uses exact coefficient comparisons and must pass its own numerical
and MPI regression before the version is released.

Reproducible performance qualification
=====================================

On the validated local Docker image ``scrimp-mpi:local``:

.. code-block:: bash

    cd ~/Simulations/scrimp
    conda activate scrimp
    bash sandbox/run_v130_wave_benchmark.sh 0.04 3

This executes **unmodified v1.2.0 sequentially** and the current code
with one, two, four and six MPI ranks, on the same unmodified wave
example, mesh and ten CN time steps (dt=0.002, final time=0.02).
The historical tag is checked out in a temporary Git worktree and
the script refuses to benchmark different wave.py sources.

The current version uses distributed COO, compact matrix memory, and
the *automatic*, unforced MUMPS policy. Both versions use the same
native Docker environment.

The report includes:

* model construction, solution and complete wall times, using medians
  across independent repeated runs;
* independently measured timing of the first CN step and later steps;
* current GetFEM, bridge, initial consistency and main TS phase timings;
* PETSc symbolic, numerical and triangular solve event counts/times,
  with initial consistency and main CN logged separately;
* MPI speedup and efficiency and exact full initial/final vector
  comparisons with the release reference.

Human-readable graphics are written under
``sandbox/outputs/<run>/png/``, alongside editable vector versions in
``svg/``. Current-version visualizations include complete solve/total
time, speedup, parallel efficiency, first versus steady CN steps,
SCRIMP phase costs and native PETSc/MUMPS symbolic/numeric/solve times.
The cross-release sequential bar chart is generated only when full
state parity with v1.2.0 passes. The short current-only numeric reuse
smoke (``sandbox/run_v130_numeric_reuse_smoke.sh 0.04 1,2``) additionally
creates ``png/01_smoke_numeric_factor_reuse.png`` and
``png/02_smoke_solve_time.png`` after all factor and state gates pass.

The runner fails when accepted steps are missing, current MPI
initial/final states disagree, or the production main CN stage
refactorizes numerically more than once. If current and native v1.2.0
states disagree, the current MPI scaling report is still valid,
but **no cross-version speedup claim is accepted** and the
release timing is excluded from the common scaling graph. Results
are stored under ignored ``sandbox/outputs/``: Markdown, CSV, JSON
and PNG.

**Caution:** SCRIMP phase timers are inclusive; nested phase totals
must not be added together. In particular, the TS time minus its
Python callbacks is *not* solely MUMPS computation.

Earlier repeated runs at h=0.01 (744,850 DOFs) showed 21--23 % total
benefit from parallel symbolic MUMPS ordering at 4--6 ranks, but still
performed ten numerical factorizations for ten CN steps. Those are
historical measurements and should not be advertised as v1.3.0
numeric-reuse benchmarks.

This performance qualification is separate from the physical
Hamiltonian, input-port power and dissipative power-balance regression;
those checks remain mandatory for scientific acceptance.
