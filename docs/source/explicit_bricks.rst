Fully explicit nonlinear bricks
================================

In SCRIMP v1.3.0, ``Brick(..., explicit=True)`` means a **fully
lagged weak-form forcing**, not merely an omitted Newton derivative.

Accepted-step contract
----------------------

For a generic nonlinear weak term `G(t,z)`, SCRIMP's implicit
PETSc TS residual on the interval starting from the last accepted
state `(t_n,z_n)` uses

.. math::

   F_{n+1}(z,\dot z) =
   M(z)\dot z + R_{\mathrm{implicit}}(t,z)
   + G(t_n,z_n),

with the sign of `G` determined by the brick's registered
port-Hamiltonian form. `G(t_n,z_n)` is assembled by GetFEM **once**
in PETSc's `PreStep` callback, copied into an owned PETSc RHS
vector, and held constant throughout Newton iterations and
line searches. A rejected time-step attempt must not update the
explicit forcing from an unaccepted Newton state. After a step
is accepted, the next `PreStep` refreshes the cache from the
new accepted state.

The corresponding implicit Jacobian is

.. math::

   J = \frac{\partial}{\partial z}
   (M(z)\dot z + R_{\mathrm{implicit}}(t,z))
   + \sigma M(z),

with the additional contracted derivative of `M(z)` in the
state-dependent-mass case. There is **no**
`D_z G(t,z)` contribution: the cached explicit term does not
depend on the current Newton state.

Usage and restrictions
----------------------

The original shallow-water gyroscopic term remains

.. code-block:: python

   S.Brick(
       "G", "(Gyro(p) * e_p) . Test_p", [1],
       linear=False, explicit=True, position="effort",
   )

SCRIMP's brick-registration convention negates the given form
for `explicit=True` (and for flow bricks); the actual GetFEM
registered form and physical sign are unchanged. When
`explicit=True` is requested, SCRIMP requires a **nonlinear,
non-`dt`, non-source** brick. Other combinations raise an
explicit error instead of silently producing an inconsistent
Jacobian.

The cache is distributed in PETSc, and its generic weak-form
GetFEM RHS follows the calibrated MPI generic-RHS scaling,
**not** the special scaling for prescribed source bricks.
`source_only` and `explicit_only` are distinct assembly
paths.

**Time accuracy:** the explicit term is held at
`(t_n,z_n)`. This is generally **first-order accurate**
for that contribution even if the main implicit scheme is
Crank–Nicolson or BDF2. In particular, the pH power identity
may acquire a splitting or lagging defect; a CN label does
not prove a second-order or structure-preserving *combined*
scheme. A second-order IMEX scheme or extrapolation is a
separate development and requires its own physics tests.

Native acceptance
-----------------

The shallow-water Jacobian regression uses the original
shipped weak forms and a genuinely distinct accepted-state
vector to freeze a **nonzero** gyroscopic forcing. It then
checks the production `IFunction` directional derivative
against `IJacobian` for independent `z` and `\dot z`
perturbations, and confirms zero drift in the cached RHS
through the Newton/FD evaluations. It compares the default
state-dependent-mass correction with the explicit
`legacy` negative control on MPI 1 and 2.

The native gate remains independent of source-only CI and
must pass before issue #27 can close:

.. code-block:: bash

   SCRIMP_DOCKER_NO_SUDO=1 \
   SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
   bash scripts/run-v130-shallow-mass-jacobian.sh

See `docs/V130_RELEASE_SPRINT_2026W48.md` for the
separate v1.3.0 release gates.
