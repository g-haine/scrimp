# Changelog

Changes to SCRIMP are recorded here. See also the
[GitHub releases](https://github.com/g-haine/scrimp/releases).

## 1.3.0 — in development (not released)

### Exemple 2C : identification de l'orientation énergétique réelle du port Grad–Grad (25 novembre 2026)

- Archive originale `v130-wave2d-coenergy-e1c5bc50-bRo68Qaf.zip` :
  MUMPS MPI1 et MPI2 terminent **10/10 pas CN**
  avec DAE algébrique exacte (2459 DOFs, zéro
  rejet), mais le **bilan énergétique FAIL** ;
  défaut maximal par pas `4.967881e-4`,
  cumulé `3.661328e-3`. Aucun résultat VTK
  validé ne doit être revendiqué pour cette exécution.
- Diagnostic sur la vraie trajectoire GetFEM :
  `dH + dt*(D_inner+P_boundary)=0` à
  `7.18e-16` près sur les dix pas.
  La convention originale `M_Y`, `-B^T` impose
  `Y=-p` au bord au sens faible ;
  la puissance `P_boundary=U*Y` est une **perte**
  (et non une entrée) dans cette formulation
  Grad–Grad. La commande testée `U=-0.5Y`
  était physiquement **active** ; l'écart précédent
  égalait exactement `2*|W_boundary|` cumulé.
- `examples/wave_coenergy.py` revient à
  **`U=+0.5Y`** avec dissipation
  `D_b=+integral UY=0.5*integral Y²`.
  Contrôles indépendants : `U-0.5Y=0`,
  `D_b=P_boundary>=0`,
  `H-H0+L_internal+L_boundary=0` ;
  les bilans par pas, la DAE, les tolérances
  et les vérifications MPI restent stricts.
- Le test de régression
  `tests/test_v130_coenergy_boundary_power_orientation.py`
  reprend **les vrais échantillons FE du pilote
  échoué** et valide le signe sans transformer
  cet échantillon en qualification du modèle
  réparé. CI et préflight 2C incluent ce test.
- **Statut : nouvelle simulation MUMPS MPI1/2 avec
  feedback corrigé attendue ; 2C reste NON QUALIFIÉ.**
  Pas de modification de `examples/published/`
  ni de libéralisation des critères énergétiques.


### Feedback à gain littéral : GetFEM linéaire et initialisation DAE 2C (25 novembre 2026)

- Deuxième tentative réelle du pilote 2C,
  archive `v130-wave2d-coenergy-0c44a1ae-rwh9fi2L.zip`,
  MUMPS MPI1 et MPI2 **FAIL avant le premier pas** :
  interfaces transparentes correctes, initialisation des champs
  et assemblage `M,K` effectués, mais
  `(-0.5*Y)*Test_U` classé à tort **non linéaire**
  par `DPHS.set_control`. L'initialiseur DAE
  linéaire refuse à juste titre cette classification.
- `scrimp/dphs.py` reconnaît désormais les
  feedbacks `-0.5*Y` et autres **gains littéraux
  numériques finis fois une seule inconnue déclarée**
  comme termes de raideur bilinéaires GetFEM, en
  conservant le Jacobien physique exact.
  Inconnues multipliées, paramètres symboliques,
  macros, fonctions et expressions ambiguës gardent
  leur traitement non linéaire fail-closed.
- Régressions positives/négatives
  `tests/test_v130_linear_scalar_feedback.py`,
  intégrées à la CI et au préflight 2C.
  Aucun changement de la physique, du gain
  `-0.5`, de l'initialiseur `algebraic`
  ou des seuils de passivité.
- **Statut : fix intégré, nouvelle exécution native
  et contrôle des sorties FE/VTK/visuels attendus**.


### Interfaces conformes : distinguer séparation matérielle et couplage par ports (25 novembre 2026)

- Première exécution native de `examples/wave_coenergy.py`,
  archive `v130-wave2d-coenergy-f110a3e3-2V4VQbbG.zip` :
  tests statiques PASS mais MUMPS MPI1/MPI2 **FAIL avant
  l'intégration**, car l'interface géométrique 10 du
  domaine concentrique ne comporte pas de ports déclarés.
  Aucun résultat physique 2C n'est donc revendiqué.
- Correction générique de `interface_deferred.py` :
  **interface matérielle transparente** si aucun port
  ni aucune brique n'utilise explicitement le label
  partagé ; FEM, paramètres et briques sont rejoués
  à l'identique, sans invention de ports, sans création
  de régions de faces orientées, avec statut traçable
  `SKIPPED_TRANSPARENT_MATERIAL_INTERFACE`.
- Les interfaces explicitement couplées gardent le
  **fail-closed** original ; déclarer seulement un port
  ou seulement une brique ne désactive aucun contrôle.
  Régressions ajoutées dans
  `tests/test_v130_interface_deferred.py` et exécutées
  avant le pilote 2C.
- **Qualification MUMPS MPI1/2 toujours en attente de
  nouvelle exécution réelle**, puis vérification des
  images GetFEM/ParaView. Pas de modification des
  équations, de `examples/published/` ou de `master`.


### Exemple 2C coénergie et visuels natifs 2D — préparés, non qualifiés (25 novembre 2026)

- `examples/wave_coenergy.py` utilise le **vrai**
  `T^{-1}` anisotrope et la densité `3-x` :
  coétats substitués, deux régions concentriques,
  damping sur le disque intérieur et feedback de signe
  **passif `U=+0.5Y`** au bord extérieur. Le signe
  avait été inversé à la préparation et a été corrigé
  après le diagnostic natif MPI1/2 du 25 novembre.
  `nu=0.2+10(0.36-r²)` reste strictement positive
  sur le seul disque intérieur pour `e_r/nu` ;
  c'est un changement **explicite** du scénario historique.
- Pilotage MUMPS MPI 1/2, DAE algébrique et dix pas
  CN avec bilan `H-H0+Li+Lb=0` :
  dissipation calculée indépendamment par
  `integral e_r²/nu` et `-integral p*e_r` ;
  perte de puissance au bord comparée à
  `+integral UY` (signe Grad–Grad `Y=-p`) ;
  défaut ponctuel de projection
  faible enregistré sans devenir un critère erroné.
  Traces originales `rank*.npz/json`, VTK/PVD
  `q,p,e_r`, parité MPI, archive d'échec
  et de succès : `run-v130-wave2d-coenergy-pilot.sh`.
- La visualisation commune
  `render-v130-wave2d-native.py` exploite
  les **vrais** VTU binaires et séries PVD : cartes
  `PNG/SVG` des champs, animation `GIF`
  de `p`, courbes énergétiques, CSV et provenance.
  Les pilotes 2A/2B produisent aussi ces sorties
  pour les cas qualifiés ; option de rendu après coup
  depuis leurs archives existantes.
- `run-v130-wave2d-showcase.sh {2a|2b|2c}`
  prépare une démonstration physique plus longue,
  toujours sur **MUMPS/FE natifs**, distincte
  de la qualification scientifique. Sphinx 2C
  intègre directement les blocs du code canonique.
- **Statut 2C : candidat uniquement, exécution native
  MPI1/2 et évaluation des visuels encore attendues.**
  Aucune assertion de convergence, aucune nouvelle
  qualification SuperLU_DIST et aucune modification
  de `examples/published` ou de `master`.


### Exemple 2B : qualification physique MUMPS MPI1/2 PASS (25 novembre 2026)

- Preuve native examinée :
  `v130-wave2d-dissipative-e582874f-9ni2cPYZ.zip`,
  source `e582874fe9d2fa296ab20a5963e4e4ce6437f092`,
  image `scrimp-mpi:v1.3.0-rc3` figée.
  Neuf tests source PASS, huit `examples/published/`
  inchangés. Les deux simulations MUMPS MPI1 et MPI2
  passent **10/10 pas CN**, 2804 DOFs, initialisation
  DAE algébrique exacte, sans rejet ni échec SNES.
- Bilan physique réel `H-H0=W_boundary-L_viscous` PASS :
  défaut par pas `2.70e-18 / 1.77e-18`,
  défaut cumulé `4.46e-18 / 2.06e-18`
  (MPI1/MPI2). Perte visqueuse finale
  `1.19096674652232e-5`, positive et non nulle ;
  puissance du port dissipatif et quadrature volumique
  GetFEM indépendante concordent à `1.36e-18 / 6.51e-19`.
- Parité sur les **vraies trajectoires FE complètes** :
  `max|z_MPI1-z_MPI2|=1.7142e-13`,
  rangs MPI2 identiques par SHA256, six séries
  `q,p,e_q,e_p,f_r,e_r` de onze VTK chacune
  vérifiées pour chaque cas. La loi `e_r=nu*f_r`
  reste une projection faible FEM, non une identité
  ponctuelle exigible.
- **Portée limitée : MUMPS MPI1/2, h=0.2, dt=0.008,
  tf=0.08.** Aucun résultat SuperLU_DIST 2B,
  étude de convergence ou temps long n'est revendiqué.
  Les équations, critères énergétiques, `master`
  et les exemples publiés restent inchangés.
  Prochaine étape : exemple 2C coénergie.


### Exemple 2B : onde 2D dissipative, extension stricte de 2A (24 novembre 2026)

- `examples/wave_dissipative.py` harmonisé avec le code canonique
  `examples/wave.py` : même rectangle, tenseur T, masse rho,
  FEM DG1/CG2, quatre ports mixtes, excitation gauche,
  initialisation DAE algébrique et 10 pas CN.
  Les **20 briques conservatives sont identiques**. L'extension
  ajoute uniquement le port résistif `f_r=-e_p`,
  `e_r=nu*f_r`, `nu=0.5*(2-x)>=0`,
  deux inconnues CG2 et **cinq briques** résistives.
- Audit énergétique basé sur **deux vraies intégrales GetFEM** :
  `D_port=∫f_r*e_r` via `Port.compute()` et
  `D_volume=∫nu*f_r*f_r` via assemblage indépendant,
  tous deux sur l'état CN moyen. Contrôle signé strict
  `H-H0-W_boundary+L_viscous=0` à chaque pas et en cumulé,
  avec dissipation non négative et non triviale.
  Les défauts constitutifs ponctuels ne sont pas
  confondus avec les équations faibles de projection FEM.
- Le nouveau pilote
  `scripts/run-v130-wave2d-dissipative-pilot.sh`
  demande **MUMPS MPI1/2** sur le véritable problème à 10 pas,
  compare les trajectoires complètes entre rangs,
  vérifie 6 sorties VTK/PVD × 11 images
  (`q,p,e_q,e_p,f_r,e_r`) et archive les erreurs
  **sans les requalifier**.
  Tests source `tests/test_v130_wave2d_dissipative_source.py`
  et Sphinx 2B reliés par **8 literalinclude** natifs.
- **Statut au 24 novembre** : qualification 2B native encore
  en attente à cette date ; **désormais PASS MUMPS MPI1/2 le
  25 novembre**, sur l'archive indiquée dans l'entrée ci-dessus.
  SCRIMP, `master` et les 8 exemples publiés restent inchangés.
  L'autodiagnostic des points-selles n'est pas ajouté en v1.3.0 ;
  la question historique Heat--Wave est différée.
- **2A précédemment validée** avec SuperLU_DIST
  `ReplaceTinyPivot+IterRefine` explicitement configuré,
  MPI1/2 : les 10 pas CN, les quatre puissances,
  les 2230 DOFs et les exports q/p/e_q/e_p sont PASS
  (archive `v130-wave2d-refined-dist-physical-95fbd8b2-HfwaSFgK.zip`).
  SuperLU_DIST *ordinaire* reste non qualifié sur 2A
  et le seuil des cinq configurations par défaut reste 3/5.


### Onde 2D : raffinement interne SuperLU_DIST à précision machine, prochain test physique (24 novembre 2026)

- Le dernier rejeu natif PETSc de la
  matrice CN originale, archive
  `v130-wave2d-superlu-dist-refinement-ad2d58a6-unyfltIn.zip`,
  confirme **ReplaceTinyPivot + IterRefine
  interne SuperLU_DIST PASS** sur 1 et
  2 rangs : `KSP_reason=4`, `PC_failed=0`,
  erreurs arrière sur `Px-b`
  `9.542e-17 / 1.272e-16`,
  erreurs sur solution connue
  `2.007e-14 / 2.256e-14`.
  Les seuils initiaux `1e-11 / 1e-7`
  sont strictement respectés.
- Le candidat alternatif GMRES
  `PC_RIGHT` reste **FAIL** :
  résidu arrière `1.97e-10 / 2.33e-10`
  même après convergence KSP.
  Ne pas confondre l'indication
  de convergence de PETSc avec un
  résidu réel suffisamment précis.
- Ajout du **pilotage opt-in uniquement**
  `sandbox/probe_v130_wave2d_refined_superlu_dist.py`
  et du lanceur
  `scripts/run-v130-wave2d-refined-superlu-dist-physical.sh`.
  Le candidat appelle **sans modification**
  `examples.wave.main()` ; il n'ajoute que
  les deux options de LU après la création
  de la DPHS, avant les solveurs PETSc,
  et exige `superlu_dist` authentique.
  La vraie validation physique
  **10 pas CN, 2230 DOFs, énergie
  anisotrope, 4 puissances signées,
  deux rangs MPI, 4×11 VTU**
  est comparée aux trajectoires
  complètes MUMPS MPI1 déjà
  archivées, sans re-simulation MUMPS.
- Sept tests de contrat source plus
  test de syntaxe AST du validateur
  scientifique, intégrés à la CI.
  **Qualification physique native
  de ce candidat toujours EN ATTENTE**.
  Aucun changement du solveur
  SuperLU_DIST par défaut, de
  `scrimp/dphs.py`, du modèle
  `examples/wave.py` ou de
  `master`. La qualification
  habituelle reste **3/5**.

### Onde 2D : remplacement de pivots opérationnel mais insuffisamment précis (24 novembre 2026)

- Archive native
  `v130-wave2d-superlu-dist-factor-replay-e87b20d1-H7IgKhB8.zip`
  examinée : **quatre vrais replays**
  SuperLU_DIST d'un même Jacobien
  original de 2230 DOFs, avec
  vérification stricte du SHA256,
  du graphe creux et des
  coefficients PETSc sur 1/2 rangs.
- SuperLU_DIST **standard** : échec
  `KSP=-11`, `PC_FAILED=2`
  (zéro pivot numérique) sur MPI1/2.
  Avec `mat_superlu_dist_replacetinypivot=true` :
  **factorisation opérationnelle**
  `KSP_reason=4`, `PC_FAILED=0`,
  sur les deux cas. MAIS
  **aucune réussite numérique**
  : résidu arrière `~5.084e-10`
  (tolérance < `1e-11`),
  erreur de solution `~7.3e-7`
  (tolérance < `1e-7`).
  Les sorties restent **FAIL**
  conformément aux exigences.
- Nouveau test **strictement isolé**
  `sandbox/replay_v130_wave2d_superlu_dist_refinement.py`
  et `scripts/run-v130-wave2d-superlu-dist-refinement.sh` :
  expérimente en MPI1/2 le
  `mat_superlu_dist_iterrefine=true`
  (raffinement interne SuperLU_DIST)
  et `KSPGMRES` avec `PC_RIGHT` et
  norme `UNPRECONDITIONED`.
  Les deux gardent le remplacement
  de pivots activé **uniquement**
  dans ce banc et contrôlent les
  vraies erreurs sur le Jacobien
  original, sans tolérance allégée.
  Six tests source et CI ajoutés.
  **Exécution native encore attendue.**
- Modèle canonique, noyau numérique,
  huit `examples/published/` et
  `master` inchangés.
  Installation onde 2D = **3/5
  configurations physiquement qualifiées**.

### Onde 2D : matrices CN égales bit à bit, rang plein et replay SuperLU_DIST de l'opérateur exact (24 novembre 2026)

- Archive native originale
  `v130-wave2d-first-cn-jacobian-574d9c71-NKqIZ0BZ.zip`
  examinée : **deux captures PETSc
  `P,K,M` complètement identiques**
  entre MUMPS et SuperLU_DIST (SHA256
  commun `24346856f358e4f81cf332ebc974e329dc084c60f5d6b98bfbcd788fdff37f54`),
  `P=K+250M` sans différence,
  2230 DOFs, 25801 coefficients stockés
  dont 25791 non nuls. Les **dix
  diagonales nulles sont les DOFs
  80–89 de `Y_L`** ; elles correspondent
  aux contraintes algébriques du
  port gauche, pas à des lignes vides.
- Le rang structurel réel est
  `2230/2230`. Un appariement biparti
  maximum couvre toutes les lignes.
  La vraie matrice est factorizable
  sans régularisation avec SciPy/SuperLU
  (`COLAMD` résidu relatif
  `2.28e-15`; `MMD_AT_PLUS_A`
  `5.54e-15`). Estimation
  `cond_1≈4.28e5`. L'échec
  `PC_FACTOR_NUMERIC_ZEROPIVOT`
  appartient donc au chemin numérique
  PETSc/SuperLU_DIST, pas à un Jacobien
  différent ou structurellement déficient.
- Nouveau **replay PETSc indépendant**
  `sandbox/replay_v130_wave2d_superlu_dist_factor.py`
  et lanceur
  `scripts/run-v130-wave2d-superlu-dist-factor-replay.sh` :
  charge uniquement le **CSR original
  SHA256 vérifié**, distribue les
  lignes locales MPI1/MPI2 et contrôle
  chaque coefficient. Compare
  SuperLU_DIST **par défaut** puis avec
  `mat_superlu_dist_replacetinypivot=true`,
  option volontairement **non activée
  dans SCRIMP**. Contrôles natifs
  du résidu arrière `Px-b` et
  d'une solution exacte, sans aucun
  changement de modèle, de CN ni
  recours caché à MUMPS.
- Ajout de six tests source et d'un
  préflight CI ; **replay natif
  EN ATTENTE**. Un PASS de
  factorisation isolée ne garantit
  ni les dix pas CN, ni les quatre
  puissances, ni la qualification
  MPI du programme canonique.
  **Statut d'installation toujours 3/5**.

### Onde 2D : permutations SuperLU_DIST rejetées et diagnostic exact du Jacobien (24 novembre 2026)

- **Étude native des permutations MPI1 FAIL**,
  preuve
  `v130-wave2d-superlu-dist-ordering-8cb2a51f-422NJxVt.zip`.
  Avec `LargeDiag_AWPM`, SuperLU_DIST
  imprime `CombBLAS is not available`
  puis subit un **segfault 139** avant tout
  pas accepté ; avec `MMD_AT_PLUS_A`,
  retour identique
  `PC_FACTOR_NUMERIC_ZEROPIVOT`.
  Aucun de ces essais ne constitue un
  contournement validé.
- Diagnostic PETSc des matrices réellement
  assemblées au premier pas CN :
  `J=K+250M`, **2230×2230**,
  **25801 coefficients stockés**,
  dix diagonales nulles (indices 80–89),
  **aucune ligne entièrement nulle**,
  pas de NaN/Inf. Ces diagonales
  ne prouvent pas la singularité d'une
  matrice de point-selle.
- Nouveau banc de diagnostic
  `scripts/run-v130-wave2d-superlu-dist-matrix-audit.sh`
  (code de capture
  `sandbox/probe_v130_wave2d_factor_matrix.py`,
  analyse
  `sandbox/analyze_v130_wave2d_factor_matrix.py`).
  Il conserve les **vraies matrices PETSc
  P, K, M** depuis le callback IJacobian
  avant la LU, comparées entre MUMPS et
  SuperLU_DIST, avec la liste des
  intervalles d'inconnues GetFEM,
  le rang structurel SciPy, les dix
  indices à diagonale nulle et un essai
  de factorisation SciPy **sans
  régularisation**. Les tests vérifient
  qu'un point-selle avec dix diagonales
  nulles peut rester de rang plein.
- **Ce nouveau diagnostic est EN
  ATTENTE de vraie exécution native.**
  L'acceptation de l'onde 2D reste
  **MUMPS MPI1/2 et SuperLU MPI1 PASS,
  SuperLU_DIST MPI1/2 FAIL**.
  `examples/wave.py`, `scrimp/dphs.py`,
  `master` et `examples/published/`
  restent inchangés par ce banc.

### Onde 2D : SuperLU_DIST PC numeric zero pivot et permutations sans régularisation (24 novembre 2026)

- Diagnostic **natif** d'un premier pas de `examples/wave.py`
  avec SuperLU_DIST MPI1, preuve
  `v130-wave2d-superlu-dist-first-step-62fb8765-QIaYkugY.zip` :
  `SNES_reason=-3`, `KSP_reason=-11`,
  `KSP_iter=0`, `PC_failed=2`,
  `factor_solver=superlu_dist`. Selon l'énumération
  officielle PETSc, **PC_failed=2 = PC_FACTOR_NUMERIC_ZEROPIVOT**.
  La TS a rejeté le premier pas jusqu'à
  `dt=1.9073486328125e-09`.
  La singularité mathématique de la matrice
  n'est **pas** démontrée par ce code de retour.
- Expérience de discrimination **non régularisante**
  prête : `sandbox/probe_v130_wave2d_superlu_dist_ordering.py`
  et `scripts/run-v130-wave2d-superlu-dist-ordering-study.sh`.
  Test de permutations SuperLU_DIST
  `LargeDiag_AWPM` et `MMD_AT_PLUS_A`
  sur le même modèle physique, sans remplacement
  de pivots, sans substitution de MUMPS ni
  changement d'énergie/DAE/CN.
  Vraie qualification de dix pas MPI1,
  quatre puissances, quatre exports VTK,
  comparaison aux **états FE complets MUMPS
  déjà archivés**, diagnostic de lignes nulles.
- **À exécuter nativement** : un succès de
  permutation n'est pas encore la qualification
  de SuperLU_DIST **par défaut** ni en MPI2.
  Le banc d'installation reste **3/5 PASS**.
  La CI vérifie les scripts, les contraintes
  physiques et l'absence de correctif
  de pivot appliqué aveuglément.
- `examples/wave.py`, `scrimp/dphs.py`,
  `master` et `examples/published/` inchangés
  dans ce diagnostic.

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
