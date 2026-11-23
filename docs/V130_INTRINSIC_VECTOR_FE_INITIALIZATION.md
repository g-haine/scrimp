# v1.3.0 — intrinsic vector finite-element initial states

**Status (23 November 2026): native small-mesh mass projection PASS; production
integration and general physical validation PENDING.**

## Immutable native evidence

Original user-supplied archive:
`v130-intrinsic-l2-1b46a613-reZMoWsJ.zip` from source
`1b46a613c2a7dbb27462df4349713de9f06a32c4`, qualified
Docker revision `1d6d0c90115213298047b28237a7f38087563a35`
(`scrimp-mpi:v1.3.0-rc3`). The pilot ran on one MPI rank and a
filtered region of a small 2-D simplicial mesh. It used the exact
constant vector `[1;2]` as an independent GetFEM GWFL source.

| FE space | Filtered DOFs | Relative mass-system residual | Old production |
| --- | ---: | ---: | --- |
| RT0 | 32 | 1.17850e-16 | REJECTED |
| RT1 | 88 | 3.49341e-16 | RETURNED |
| RT2 | 168 | 3.83978e-16 | REJECTED |
| BDM1 | 64 | 2.59574e-16 | RETURNED |
| BDM2 | 132 | 8.84826e-16 | REJECTED |
| lowest Nedelec | 32 | 1.09909e-16 | REJECTED |

The older `DPHS.set_initial_value()` relies on `MeshFem.eval()`.
For RT0/RT2/BDM2/lowest Nedelec, the output length is one more than
`MeshFem.nbdof()`, so production rejects. For RT1 and BDM1, the
sizes happen to match but there is **no physical-field evidence**
that these arrays represent the correct intrinsic coefficients.
Indeed, `MeshFem.is_lagrangian()` was false for all six spaces.

All six **separate GetFEM mass projections** solved
`M c = b` with a relative algebraic residual below 9e-16.
This proves native assembly and small dense linear-algebra
consistency, **not** physical reconstruction, general analytic
expression support, numerical stability across meshes or MPI, nor
PDE/DAE convergence. The pilot deliberately does NOT change
production `set_initial_value`.

## Next independent native gate

Source `9836f7de85a77cdbc6735b1a23be2b8208942e1a` extends
`sandbox/probe_v130_intrinsic_l2_projection.py` to check the
**physical** L2 error with GetFEM's own quadrature:

```text
||q_h - [1;2]||_{L2(region)} / ||[1;2]||_{L2(region)}
```

All six selected 2-D intrinsic spaces contain constant vectors.
The pilot requires relative physical error below `5e-8`, in
addition to the original matrix/linear-solver gates. It also
measures the old `set_initial_value`'s physical error for RT1
and BDM1, where that implementation previously returned without
error. These fields will be in the **new** native evidence archive;
they are not present in the original 1b46a613 evidence.

To reproduce with the pinned image:

```bash
git pull --ff-only
SCRIMP_DOCKER_NO_SUDO=1 \
SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
bash scripts/run-v130-intrinsic-l2-projection.sh
```

## Production requirements, still open

1. Preserve the legacy nodal interpolation semantics for supported
   Lagrangian spaces and avoid a silent behavioral change.
2. For intrinsic vector spaces, provide a documented initialization
   operator, preferably a PETSc-backed **sparse** L2 projection,
   not the diagnostic's small dense `numpy.linalg.solve`.
3. Test nonconstant analytically prescribed vector fields and the
   exact SCRIMP expression interface (Python/NumPy vs GetFEM GWFL),
   with clear fail-closed handling of unsupported expressions.
4. Check coefficient shape, mesh-region restriction, orientations,
   and collective ownership with MPI1/2; retain original prescribed
   differential initial coefficients once the DAE initializer
   is entered.
5. Preserve actual GetFEM `IFunction`-based DAE feasibility checks
   and reject incompatible Dirichlet/interface data explicitly.
6. Do not equate L2 projection with canonical commuting H(div) or
   H(curl) interpolation. Normal flux, tangent/curl, orientation and
   any PDE conservation claims require separate, physical tests.

No published article-era source, tagged version, or original native
output is changed by this work.
