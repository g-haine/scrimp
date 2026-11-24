# Changelog

Changes to SCRIMP are recorded here. See also the
[GitHub releases](https://github.com/g-haine/scrimp/releases).

## 1.3.0 — in development (not released)

### Promotion des interfaces conformes et des exports ParaView MPI (24 novembre 2026)

- **Branche `release/1.3.0-preparation` uniquement, sans fusion dans
  `master`, publication ni nouveau tag.** Les chemins de code par défaut
  utilisent désormais la compilation automatique prudente des interfaces
  internes conformes et l'export MPI partagé par le rang 0 si tous les
  rangs possèdent les mêmes états complets et la même géométrie.
  Les interfaces ambiguës et les données non répliquées **échouent
  explicitement** au lieu de produire un opérateur ou des VTU incomplets.
- Options de retour arrière :
  `SCRIMP_AUTO_INTERFACE=0` et
  `SCRIMP_PV_EXPORT_MODE=legacy` (ce dernier réintroduit l'ancien
  risque d'écritures concurrentes). Les anciens drapeaux
  `SCRIMP_V130_*_PILOT` restent des alias de compatibilité ;
  les nouvelles options ont priorité.
- **PASS native** dans l'image figée `scrimp-mpi:v1.3.0-rc3` :
  K manuel/automatique et MPI1/MPI2 ; initialisation DAE
  algébrique stricte ; 4/4 pas CN ; bilans énergétiques/échanges
  non triviaux ; géométrie et champs VTU/PVD ; export MPI2 partagé
  du rang 0. Deux affectations Heat-Wave sur le disque/anneau
  (`heat=1,wave=2` et `heat=2,wave=1`), 7689 / 8056 DOF,
  défauts énergétiques ≤ 2.69e-15 et 54/54 tests source précédents.
- **Échecs et correctifs consignés** : mesure initiale codant des
  facettes `30/31`, corrections des signes à normales extérieures,
  `scipy.sparse.vstack(generator)`, zéros CSR stockés en MPI2,
  Hamiltonien calculé deux fois, `host_error` non initialisé,
  sous-processus de tests tué par SIGKILL (cause externe inconnue),
  ancien libellé d'échec de comparateur, risque d'écritures VTK
  concurrentes. Voir
  `docs/V130_INTERFACE_PV_PRODUCTION_DECISION_20261124.md`.
- **PASS final des chemins par défaut** sans drapeaux expérimentaux
  au commit `241f58fb`, test natif
  `scripts/run-v130-production-defaults.sh` dans l'image RC3 figée :
  **71/71 tests source**, 8 instantanés publiés immuables,
  quatre simulations GetFEM/PETSc en MPI1/MPI2, strict DAE,
  quatre pas CN acceptés par simulation, zéro échec SNES,
  bilans énergétiques et puissances d'interface vérifiés,
  trajectoires MPI2 identiques entre rangs, fichiers
  ParaView partagés et valeurs physiques MPI1/MPI2 conformes.
  Archive `v130-auto-interface-MPI-CN-241f58fb-YEJXJhSK.zip`.
  Le rapport original conservait une ancienne métadonnée
  `production_auto_compiler_enabled_by_default=false`
  codée en dur ; sa correction est **report-only** (`942e23e`)
  et ne modifie aucune simulation.
- **Limites assumées** : les tests natifs ne couvrent pas la géométrie
  générique, les autres familles d'EF ni les trajectoires non
  répliquées. Le PASS final est spécifique aux deux orientations
  Heat-Wave en MPI1/2. Aucun tag, fusion ou publication.

### Performance and MPI
- Use exact `on_change` MUMPS factor reuse by default. Reuse the
  symbolic factor when the assembled PETSc CSR graph is unchanged and
  the numerical factor when the stored coefficients are identical on
  every MPI rank. Changes to the time-scheme shift or nonlinear
  tangents invalidate numerical reuse automatically.
- Prefer MUMPS parallel analysis `ICNTL(28)=2` for multi-rank LU/MUMPS
  runs, respecting user-provided PETSc options.
- Preserve the legacy `never`/`always` overrides and the physical
  model/discretization. No unconditional Jacobian freeze is introduced.
- Maintain the previously developed distributed GetFEM-to-PETSc COO
  bridge and compact-matrix-memory path.
- Add a release-to-current wave benchmark with native v1.2.0 serial
  reference, staged current PETSc profiling, complete state-vector
  comparisons and 1/2/4/6-rank strong scaling.
- Record the deferred MPI scalability, nonlinear-physics and
  experimental Schur results in `docs/V130_DEVELOPMENT_ROADMAP.md`.

### Truly explicit nonlinear weak-form bricks (directional FD and physical TS native PASS)

- Correct `Brick(explicit=True)` semantics: previously SCRIMP
  excluded such a nonlinear term from `IJacobian` while
  still reevaluating it at the *current* Newton state
  in `IFunction`. This produced an inconsistent
  Jacobian on `examples/shallow_water.py`
  (native zero-time-derivative directional FD
  discrepancy `2.02366589e-5`, MPI 1/2).
- Now the **entire nonlinear `dt=False` non-source brick**
  is assembled at `(t_n,z_n)`, the **last accepted**
  PETSc TS time/state, and held as a distributed,
  calibrated GetFEM/PETSc RHS through the following
  Newton and line-search iterations. A `TS.setPreStep`
  callback is installed for both the historical
  pseudo/BDF initialization TS and the main physical TS.
  The term has no derivative with respect to the
  current implicit trial state, and remains excluded
  from the implicit tangent. This is a true first-order
  lag of that contribution, not a fully implicit or
  midpoint evaluation.
- Preserve the **original gyroscopic G weak-form
  expression and sign** in the published shallow-water
  example; its `explicit=True` flag remains set.
  Requests for unsupported `explicit=True` brick
  classes (linear, `dt=True`, source) now fail
  explicitly instead of silently receiving
  half-implicit semantics.
- The underlying fully explicit split **does not
  automatically preserve the exact discrete
  Hamiltonian power balance** or BDF2 order for
  the lagged term. Validate its consequences
  separately on the short nonlinear SWE time test.
- **Full nonlinear mass-Jacobian directional-FD native gate PASSED**
  on original `examples/shallow_water.py` (1,246 DOFs),
  MPI 1 and 2, at source
  `fae066a64b649131351e686ac9383ea11c13c0a6`,
  immutable `scrimp-mpi:v1.3.0-rc3`:
  default maximal relative `IFunction`/`IJacobian`
  error **1.06105781e-8** versus
  **3.84378444e-6** for explicit `legacy`
  (~362× improvement). The gyroscopic explicit
  RHS was nonzero (`0.00536124096`) and
  exactly unchanged through all FD/Newton trials.
  The acceptance threshold was not relaxed.
  Host evidence:
  `sandbox/outputs/v130-shallow-mass-J-fae066a6-VvvbTRp4/`.
- **Actual native physical integration PASSED independently**
  at source `9ecd173983252795a3cba13b5879d29127b57638`
  on unchanged `scrimp-mpi:v1.3.0-rc3`. Four genuine
  BDF2 TS accepted steps per published/manufactured
  original-weak-form shallow-water case (MPI 1 and 2),
  zero rejected/SNES failures, nonzero accepted-state
  frozen gyroscopic `G` with zero in-Newton cache drift,
  maximum full-state MPI parity `6.39215556e-16`
  and Hamiltonian parity `1.26771841e-15`.
  Recorded signed final Hamiltonian changes:
  `-0.00206566661` and `-0.00206572620`
  (relative `~6.54e-8`); **not** an exact-power
  conservation or global order-two claim.
  Evidence:
  `sandbox/outputs/v130-shallow-explicit-TS-9ecd1739-LZjvda5C/`.
  Together with the independent native directional FD gate,
  the original scope of issue #27 is qualified.
  Longer-time energy/power, MPI rank 4 viscous SWE and
  general nonlinear DAE initialization remain separate.

### Consistent DAE initialization (new linear default; native two-model gate PASSED)

- Add a **distributed PETSc algebraic-only initialization** method for
  supported real, linear DAE models. It keeps the prescribed differential
  initial coefficients unchanged, solves only for algebraic efforts/ports
  with bounded LSQR refinements, and verifies the *actual*
  `IFunction` algebraic residual before starting the physical TS.
  Unsupported layouts or incompatible fixed initial data raise
  explicit errors; there is **no automatic pseudo-time fallback**.
- **Use `init_step_ts_type="algebraic"` as the implicit default**
  on the development branch for the supported linear class. The
  native **no-option** gate **PASSED** in the qualified
  `scrimp-mpi:v1.3.0-rc3` runtime at
  `3a17a65b7a17fbdade82bd3ea2bd3bb8cd8a53f8`
  on conservative and dissipative 2D acoustic pH-DAEs,
  MPI 1/2, ten CN steps, with exact prescribed dynamic
  states, independent full-trajectory and power references
  and negative Gaussian rejection.
  `init_step_ts_type="pseudo"` and `init_step_ts_type="bdf"`
  remain explicit alternatives; `init_step=False` remains valid.
  This is **not yet a tagged v1.3.0 release** or a claim of
  universal nonlinear DAE initialization.
- Document the distinction between **exact algebraic projection**,
  historical **pseudo-time relaxation** (which can move the
  prescribed differential initial data), and **small negative-start
  BDF2** (which also changes them, even with a very small step).
  See `docs/source/initialization.rst`.
- Preserve the historical pseudo behavior explicitly in the
  original top-level examples, legacy notebooks, published Maxwell
  example, GUI-generated legacy scripts and selected sandbox workloads.
  Their physical initial data and equations are not modified.
- **Native core opt-in qualification PASSED** on two independent
  port-Hamiltonian linear acoustic models, conservative and dissipative,
  on MPI 1/2 with ten fixed CN steps. Native algebraic residuals
  were about `1e-13`, initial dynamic drift zero, and full
  signed Hamiltonian + nonzero boundary work (+ positive internal
  damping on the dissipative model) defects about `1e-15`.
  The incompatible original Gaussian was explicitly rejected on
  both MPI ranks. The historical pseudo trajectory matched its
  frozen reference; tiny negative-time BDF2 remained usable with
  measured nonzero initial drift.
- **Remaining release gates:** the separate state-dependent-mass
  Jacobian on the original shipped `examples/shallow_water.py`
  (**#27**), MPI exceptional cleanup/memory, all shipped examples,
  3D Maxwell, GUI, Sphinx and packaging. The corrected nonlinear
  mass Jacobian **already defaults** to
  `SCRIMP_NONLINEAR_MASS_JACOBIAN=finite_difference`;
  the native directional-FD regression is pending.
  General nonlinear fixed-data DAE initialization by
  safeguarded Newton (**#28**) is a separate feature
  requiring independent native qualification.
  Do not tag or rebuild the release image before all
  release-blocking gates pass. See
  `docs/V130_RELEASE_SPRINT_2026W48.md`.

### Release qualification outstanding
- Confirm native numeric and symbolic LU call counts agree with the
  exact operator coefficient/graph changes detected collectively by
  `on_change`; qualify invalidation when coefficients change. The
  final fixed-CN MATCHSTEP may legitimately change the shift by
  floating-point roundoff (2 numeric LU for 10 steps observed).
- Execute the v1.2.0-vs-current serial and current MPI benchmark locally,
  validate complete states and discrete Hamiltonian/port-power balance,
  rerun the full test suite and build Sphinx HTML/LaTeX.
- Only then bump packaged/documentation versions and create a release tag.

## 1.2.0 — 2026-11-12

### Added

- **Global scalar (lumped) ports** through `DPHS.add_lumped_port()`:
  finite-dimensional state, co-state and algebraic variables with exactly
  one degree of freedom, independent of mesh resolution.
- Support for assembling coupled distributed PDE and global ODE weak-form
  contributions, with GetFEM variables synchronized to PETSc Newton
  trial iterates.
- Two new runnable, documented port-Hamiltonian ODE examples:
  `examples/mass_spring_damper.py` and
  `examples/nonlinear_pendulum.py`.
- `docs/source/lumped.rst` and worked Sphinx pages explaining the
  formulation, implementation and validation of the new examples.
- An **experimental** viscous shallow-water / moving-tank PDE–ODE
  prototype in `sandbox/sloshing_tank.py`, plus convergence diagnostics
  in `sandbox/tank_convergence.py`.

### Changed

- Corrected the nonlinear implicit residual assembly path for lumped
  systems having **linear time-derivative mass matrices**.
- Standardized the SCRIMP copyright/license header in tracked Python
  sources and the GUI code-generation template.
- Updated the test suite to use the currently supported public API,
  rather than removed symbolic/schema modules.
- Updated the package and documentation version to **1.2.0**.
- Adopted the standard setuptools PEP 517 build backend with
  `pyproject.toml`, retaining `setup.py` for package metadata.
- Used the root project README for the package's long description.

### Validation

- Python syntax check passed, excluding generated `outputs/` data.
- Full unit-test discovery: **57 tests passed**.
- Full default runs of both new finite-dimensional examples passed.
- Existing `examples/wave.py` completed: **501 steps**, no rejected
  steps and no nonlinear solver failures; Hamiltonian and boundary
  powers computed.
- Sphinx HTML and LaTeX built with `-W` and no warnings;
  the documentation PDF contains **100 pages**.
- The distribution wheel was built successfully. After the final
  README metadata correction, rebuild the wheel once more before
  tagging.

### Known limitations

- The sloshing tank remains **experimental**: wall impermeability and
  convergence of spatial fields have not yet been fully validated.
- The corrected nonlinear residual path does **not** address general
  state-dependent nonlinear mass matrices.
- PETSc may report an unused `-isset false` option; some real-valued
  runs emit complex-to-real conversion warnings.
- In the historical wave example, the last monitor entry was
  `t = 5.01` for `t_f = 5.0`; this endpoint discrepancy is not
  addressed in this release.
- Standalone examples have passed runtime checks; quantitative
  time-refinement and energy-balance studies for both new ODE
  examples remain future validation work.

## 1.1.1 — 2024-03-21

- Changed the output-directory behavior to use the current working
  directory.
- Corrected the sign of a source term and improved documentation.
- Added notebooks.

See the [v1.1.1 release notes](https://github.com/g-haine/scrimp/releases/tag/v1.1.1).
