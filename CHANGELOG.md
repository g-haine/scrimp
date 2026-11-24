# Changelog

Changes to SCRIMP are recorded here. See also the
[GitHub releases](https://github.com/g-haine/scrimp/releases).

## 1.3.0 — in development (not released)

### Onde 2D conservative : trois LU réussis, diagnostic SuperLU_DIST isolé (24 novembre 2026)

- **Campagne cinq configurations nativement exécutée**
  sur le code utilisateur `examples/wave.py`,
  source `99d99f2d`, archive
  `v130-wave2d-conservative-five-lu-99d99f2d-tEOjqgHy.zip`.
  MUMPS MPI1, MUMPS MPI2 et SuperLU MPI1 sont
  **PASS**, 2230 DOFs, 10/10 pas CN, trois
  bilans énergétiques par pas < `3.1e-18`,
  quatre séries VTU/PVD et puissance gauche
  `0.048996389`. Les différences de
  trajectoires à MUMPS MPI1 sont
  `9.25e-14` et `1.61e-13` respectivement.
- **SuperLU_DIST MPI1 et MPI2 échouent au
  tout premier pas** avec PETSc code 91
  `DIVERGED_STEP_REJECTED`, après
  assemblages GetFEM et initialisation
  DAE exacte (`||F_a||=0`). La cause
  profonde SNES/KSP/PC n'est pas visible
  dans cette archive ; ne pas en déduire
  une erreur de physique ou un pivot nul
  sans preuve. **Gate globale : FAIL 3/5.**
- Nouveau lanceur **strictement diagnostic**
  `scripts/run-v130-wave2d-superlu-dist-diagnostic.sh` :
  un seul pas `dt=0.008` avec SuperLU_DIST
  dans l'image figée, activation du diagnostic
  natif préexistant
  `SCRIMP_TS_FAILURE_DIAGNOSTICS=1`
  pour obtenir TS/SNES/KSP/PC en séquentiel,
  et archivage d'échec sans changer
  les équations, paramètres ou solveurs.
  **Nouvelle exécution native en attente.**
- Le rapport de l'acceptation à cinq backends
  est renforcé : après exécution des cinq,
  il consigne désormais **tous leurs
  codes de retour**, et non plus seulement
  le premier échec rencontré pendant
  la vérification scientifique.
  Aucun retour à MUMPS ni assouplissement
  des contrôles de conservation n'est
  autorisé silencieusement.

### Onde 2D conservative : modèle canonique et deux rangs MPI qualifiés (24 novembre 2026)

- **Second pilote natif MUMPS MPI1/MPI2 entièrement PASS**
  sur `cf42089d`, archive
  `v130-wave2d-mixed-conservative-cf42089d-pI8VXuBa.zip` :
  2230 DOFs, dix pas CN, initialisation DAE stricte,
  puissance gauche non nulle (`0.0489964`),
  défaut énergétique par pas < `3e-18`,
  tous les champs physiques et coénergies
  `q,p,e_q,e_p` exportés en **88 VTU**
  avec huit PVD. Les trajectoires des deux
  rangs MPI2 sont identiques ; écart
  MPI1/MPI2 `4.16e-14`.
- Le candidat temporaire de `sandbox/` a été
  **promu dans `examples/wave.py`**, en gardant
  l'API publique `wave_eq()`. Les briques,
  paramètres et données nativement qualifiés
  ont été conservés ; le double script de
  simulation et son lanceur exploratoire
  sont retirés (historiques Git conservés).
- **La page Sphinx `examples/wave.rst`
  est reconstruite** pour la conservative :
  vrai tenseur de tension anisotrope,
  traction scalaire `n·e_q`, commande
  gauche en vitesse, conservation à chaque
  pas. Les huit étapes `literalinclude`
  lisent le code Python réel et évitent
  la dérive d'expressions GWFL ;
  la section dissipative reste indépendante.
- **Acceptation cinq backends en attente** :
  `scripts/run-v130-wave2d-conservative-acceptance.sh`
  vérifie MUMPS MPI1/2, SuperLU MPI1,
  SuperLU_DIST MPI1/2, états, Hamiltonien,
  puissances des quatre frontières et
  géométries/fichiers VTK.
  Les simulations dissipatives et coénergie
  ne sont pas encore qualifiées dans ce lot.
  Aucun changement sur `master` ou dans
  `examples/published/`.

### Export ParaView MPI : coénergies et champs d'effort 2D (24 novembre 2026)

- Premier pilote natif onde 2D anisotrope **MUMPS MPI1 et MPI2** :
  véritable physique PASS (2230 DOFs, 10 pas CN sans
  rejet, projection DAE stricte, bilans énergétiques
  au niveau `1e-18`, puissance de bord gauche
  `≈0.049`) ; archive originale
  `v130-wave2d-mixed-conservative-abea0401-Jm7S0brY.zip`.
- **Première qualification complète FAIL** sur
  l'export MPI de `e_q` : le validateur de
  réplication recherchait `self.ports["e_q"]`,
  alors que le port propriétaire est `q` et
  que son coétat `e_q` est une variable
  GetFEM parfaitement exportable. MPI1 avait
  déjà écrit ses 11 VTU de `e_q` et `e_p`.
  Les **deux rangs MPI2 possèdent des
  trajectoires FE identiques bit-à-bit** ;
  écart maximal MPI1/MPI2 `4.08e-14`.
- Correction générique `DPHS._validate_root_only_pv_replication` :
  localiser **exactement un** port propriétaire
  d'une variable par ses noms `flow` ou `effort`.
  Les variables de ports d'état, leurs coénergies
  et les ports de contrôle reçoivent la même
  validation collective de géométrie/DOFs/FE.
  Champs inconnus et propriétaires ambigus
  restent explicitement rejetés.
- Régressions source du **vrai validateur**
  étendues dans
  `tests/test_v130_interface_pv_production_policy.py`.
  Le lanceur du 2D les exécute en préflight.
  **Nouveau run natif complet en attente** :
  ne pas présenter le premier échec de
  qualification d'export comme un PASS de bout en bout.
  Aucune promotion vers `examples/wave.py`
  ou modification de `master` pour cet essai.

### Synchronisation stricte des blocs Sphinx du tutoriel onde 1D (24 novembre 2026)

- Revue **bloc par bloc** de `docs/source/started.rst`
  face à la source scientifique qualifiée
  `examples/wave_1d.py`. Correction du maillage
  documenté (`h=0.05` au lieu de `0.01`),
  des noms Python, de la description de la tension
  physique, de `init_step=True`, des options
  explicites de solveur `gmres/lu`, de l'ordre des
  énergies et des paramètres de visualisation.
- **14/14 appels GWFL Brick** documentés sont
  désormais identiques à ceux de la source,
  y compris les deux ports signés. Les **11
  blocs Python** constituent une recette exécutable
  dans leur ordre d'apparition et conservent
  exactement les paramètres numériques du programme.
- Nouveau test `tests/test_v130_wave_1d_doc_parity.py`
  dans la CI : comparaison **AST** des constructeurs,
  briques, FEM, commandes, DAE/CN, solveurs et
  appels Hamiltonien ; échec explicite si
  la documentation dérive du code natif.
- Modifications **exclusivement documentaires
  et de tests** : les résultats physiques natifs
  des cinq solveurs et de la démonstration longue
  restent valides. Contrôle source supplémentaire
  à exécuter sur l'environnement de l'utilisateur.

### Tutoriel initial : corde vibrante 1D, cinq solveurs et figures natives (24 novembre 2026)

- **Exemple 01 qualifié** sur les chemins utilisateurs
  `examples/wave_1d.py`, avec les cinq vrais solveurs/tailles MPI
  admissibles : MUMPS (1/2), SuperLU (1), SuperLU_DIST (1/2).
  Conditions initiales P2 exactement compatibles,
  128 DOFs, dix pas CN fixes, ports de bord orientés
  et bilan énergétique à `1.791e-16`.
  Preuve native :
  `v130-first-tutorial-wave1d-2a170c82-cwkGUoXC.zip`.
- **Démonstration illustrative longue également PASS
  nativement** après le correctif d'archivage :
  `t_f=1.5`, 150 pas CN acceptés, 151 états enregistrés,
  MUMPS MPI1, aucun rejet/SNES, résidu d'initialisation
  DAE `2.481390789592951e-17`, défaut d'énergie
  par pas `3.434752482434078e-16`, cumulatif
  `7.216449660063518e-16`. Les **302 VTU**
  `q`/`p` ont été relus indépendamment.
  Preuve : `v130-wave1d-long-doc-7e2e95a4-75oCWG1T.zip`.
  Elle ne remplace pas la validation des cinq backends
  à horizon court.
- La documentation Sphinx `Getting started` comprend
  une **animation SVG légère de 1,5 s issue des vrais
  états acceptés** et un diagramme vectoriel du
  Hamiltonien, du travail et des puissances de bord.
  Le GIF complet de **151 images** et les cartes
  spatio-temporelles, profils, CSV et NPZ se
  régénèrent avec `scripts/render-v130-wave1d-doc.py`,
  sans versionner les nombreux fichiers binaires.
  Le GIF a fait l'objet d'une revue visuelle, et
  toute reconstruction de déplacement est identifiée
  explicitement, avec sa jauge.
- L'ancien fichier exploratoire
  `sandbox/wave_1D.py` a été supprimé ;
  `examples/wave_1d.py` reste l'unique tutoriel
  canonique, lié aux publications ISAE-SUPAERO
  pertinentes. Aucun changement dans
  `examples/published/` ni `master`.

### Moniteur PETSc TS : échéances d'archivage à l'ULP (24 novembre 2026)

- Correction d'un problème d'**arrondi flottant** dans
  `DPHS.monitor()` : le test strict
  `next_saved_at_t - t <= 0` pouvait ignorer un véritable
  pas CN accepté lorsque `last_saved_t+dt_save` dépassait
  `t` de **un ULP**.
- Reproduction native exacte : le tutoriel onde 1D prolongé
  jusqu'à `t=1.5` (`dt=dt_save=0.01`) a réussi
  **150/150 pas PETSc**, sans rejet ni échec SNES,
  avec un résidu DAE initial `2.48139079e-17`,
  mais n'a pas archivé l'état **i=149, t≈1.49**.
  La comparaison scientifique a **correctement refusé**
  l'archive de 150 états au lieu de 151, avant
  l'export VTK. Preuve de l'échec :
  `v130-wave1d-long-doc-f0b589a9-uOwc7fiY.zip`.
- La nouvelle fonction pure
  `scrimp/utils/trajectory_sampling.py:saved_time_is_due()`
  tolère au plus **16 ULP**, plafonnés à
  `1e-6 * dt_save`, sans fabriquer de solution
  ni déformer la grille physique ; le sous-échantillonnage
  explicite `dt_save>dt` est conservé.
  Régressions ciblées :
  `tests/test_v130_trajectory_sampling.py`,
  y compris l'oubli exact du pas 149.
- Le contrôle scientifique `last_ts_steps == len(saved_t)-1`
  et les bilans CN restent stricts ; l'exception
  explique désormais quelle quantité manque.
  En cas de nouvel échec, les états FE acceptés,
  l'Hamiltonien et les puissances sont conservés
  dans des fichiers `rank*.failed.npz/json`
  marqués **NON QUALIFIÉS**, pour diagnostic.
- **État final : correctif validé nativement.**
  La nouvelle simulation longue sauvegarde 151/151
  états, passe les portes énergétiques et produit
  302 fichiers VTU complets et le GIF vérifié.
  Archive `v130-wave1d-long-doc-7e2e95a4-75oCWG1T.zip`.
  La preuve de l'échec initial est conservée,
  sans assouplissement des contrôles physiques.

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
