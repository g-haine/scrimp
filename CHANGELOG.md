# Changelog

Changes to SCRIMP are recorded here. See also the
[GitHub releases](https://github.com/g-haine/scrimp/releases).

## 1.3.0 — in development (not released)

### Saint-Venant 2D non visqueux — qualification physique BDF2 MPI1/2 prolongée : pilote préparé (25 novembre 2026)

- **Poursuite des exemples après la qualification complète Heat–Wave.**
  Le vrai `examples/shallow_water.py` utilise une
  masse d'évolution `M_p(h)` non linéaire,
  une surface initiale discontinue,
  quatre ports de bord fermés, un
  terme gyroscopique explicitement
  retardé et l'intégrateur BDF2.
  Des campagnes **natives antérieures**
  ont déjà qualifié le jacobien de masse,
  l'initialisation algébrique
  `nonlinear_affine` et de très courtes
  simulations physiques avec raffinements
  en temps, en MPI1/2. Ces preuves
  sont conservées et ne sont pas
  confondues avec la qualification
  de l'horizon par défaut.
- Nouveau pilote de validation
  `scripts/run-v130-shallow-water-extended-physical.sh` :
  **quatre véritables exécutions MUMPS**
  sur le rectangle original,
  à temps final commun `T=0.002`
  et pas BDF2 `dt=0.0001`
  (**20 pas**) ou `0.00005`
  (**40 pas**), chacune en MPI1
  puis MPI2. Le probe d'origine
  `sandbox/probe_v130_issue28_shallow_short_physical.py`
  gagne un simple drapeau de sortie
  `--capture-full-fe` qui écrit les
  **vecteurs FE véritables par rang**,
  sans changement des équations, FEM,
  contrôles, données différentielles
  ni tolérances DAE.
- Vérificateur indépendant
  `scripts/verify-v130-shallow-water-extended.py` :
  parité des coefficients FE et des
  intégrales masse/H/puissance des
  deux rangs MPI2 et entre MPI1/MPI2 ;
  conservation de masse, travail aux
  frontières, défaut observé
  `H(t)-H(0)-int P_bord dt`,
  comparaison à temps final fixé
  pour les deux raffinements.
  **Aucun bilan Hamiltonien exact
  n'est revendiqué** pour BDF2 avec
  gyroscope explicitement retardé,
  ni aucun ordre de convergence
  présupposé pour la surface initiale
  discontinue. Toutes les sorties
  sont archivées, y compris les FAIL.
  Tests de source et contrôles négatifs
  intégrés à la CI ; les huit sources
  `examples/published/` restent immuables.
- **Statut : pilote préparé,
  exécution native encore en attente.**
  Le scénario \`t_f=0.5\` avec
  initialisation historique \`pseudo\`
  demeure distinct et non qualifié
  par ce nouveau pilote. Le Saint-Venant
  **visqueux tournant** avec contrainte
  algébrique d'indice supérieur
  est également un autre problème
  et fera l'objet de sa propre
  qualification. Référence :
  `docs/V130_SHALLOW_WATER_EXAMPLE_ACCEPTANCE.md`.


### Heat–Wave 2D — qualification native physique MPI1/2 et bilan énergétique CN complet sur 15 s : PASS (25 novembre 2026)

- **Qualification native finale obtenue** :
  archive réelle
  `v130-heat-wave-MPI2-longCN-97c7606b-7NBvjofW.zip`,
  source
  `97c7606bac14f7420140dd1fc2a1d88847301538`,
  image épinglée
  `scrimp-mpi:v1.3.0-rc3`, MUMPS.
  Statut de la preuve
  `PASS_HEAT_WAVE_MPI2_SHORT_AND_MPI1_LONG_CN_PHYSICAL`,
  `errors=[]`, contrôle des
  huit `examples/published/`
  inchangés PASS. SHA256 de
  l'archive :
  `9aff9f1f6cf3708f8736f48a77731ad4d814c08fd0d54c915a32f8aa1ac516fe`.
- **MPI1 et MPI2 (4 pas CN à
  `dt=.001`, 7 689 DOFs) : PASS.**
  Vecteurs FE complets identiques
  à `2.2297719226571644e-12`
  au plus entre MPI1 et MPI2,
  Hamiltoniens à
  `5.662137425588298e-15`
  au plus ; **égalité bit-à-bit
  des deux rangs MPI2** pour les
  états et toutes les observables.
  Les quadratures GetFEM des
  puissances d'interface orientées,
  dissipation et bilans chaleur/onde
  respectent les critères physiques
  sur chaque rang.
- **MPI1 temps long (3 000 pas CN,
  `dt=.005`, `t_f=15`) : PASS.**
  `preflight-rank0.json` vérifie
  **3 000 pas acceptés, 0 rejet,
  0 échec SNES et exactement
  3 001 états FE réellement
  enregistrés**, sans indice
  manquant. Le Hamiltonien
  gaussien diminue de
  `0.7836421211876778`
  à `0.008268389052222169`,
  sans augmentation sur aucun
  pas accepté. Défauts maximaux
  des bilans CN pas à pas :
  chaleur `1.3808398868775384e-15`,
  ondes `3.978506976989871e-15`,
  total indépendant
  `4.017619570362285e-15`.
  **Bilan énergétique cumulatif
  sur les 3 000 pas** :
  `5.697699100711716e-15` ;
  puissance nette d'interface
  orientée `3.1086244689504383e-15`.
  Les tableaux NPZ ont fait l'objet
  d'un contrôle indépendant
  des identités et de tous les
  instants CN réellement adjacents.
- **Régression générique de sortie
  corrigée et qualifiée :**
  le pas PETSc accepté `n=2999`
  avait été omis des sorties de
  l'archive précédente `42116b58`
  par un seuil de temps flottant.
  Le nouveau `DPHS.monitor`
  enregistre directement chaque
  pas réellement accepté pour
  un CN fixe lorsque
  `dt_save=dt`.
  Il ne modifie ni le calcul,
  ni l'initialisation DAE, ni
  les équations ou les bilans.
  Les autres politiques
  d'échantillonnage sont
  conservées. L'ancienne archive
  FAIL est conservée en tant que
  preuve historique ; la nouvelle
  archive `97c7606b` contient
  tous les 3 001 états et
  satisfait les portes physiques
  originales, sans relâchement.
- **Limites strictes :**
  qualification MPI1/MPI2 sur
  **4 pas**, temps long physique
  sur **MPI1 avec condition
  initiale gaussienne uniquement**.
  Les deux longues animations
  (`gaussian` et `gallery_packet`)
  sont des preuves distinctes
  de visualisation. MPI2 temps
  long, VTK parallèle,
  SuperLU_DIST, convergence
  h/dt et identification spectrale
  d'un mode *whispering gallery*
  ne sont pas revendiqués.
  Détails dans
  `docs/V130_HEAT_WAVE_PHYSICAL_AUDIT.md`
  et `docs/source/examples/heat_wave.rst`.


### Heat–Wave 2D — vrai PASS physique MPI1/2 et correctif générique de l’échantillonnage CN sur temps long (25 novembre 2026) — incident antérieur désormais résolu

- **Qualification native MPI1/MPI2 acquise :**
  archive originale
  `v130-heat-wave-MPI2-longCN-42116b58-FCxHdEia.zip`,
  source `42116b5803e5295225394582e6ddf930a4269afa`,
  MUMPS, GetFEM/PETSc, **7 689 DOFs**,
  quatre pas CN fixes `dt=.001`,
  `short-verification.json` :
  `PASS_HEAT_WAVE_MPI1_MPI2_SHORT_PHYSICAL`.
  Toutes les quadratures physiques
  chaleur/onde, puissance d'interface
  orientée, dissipation et identité
  énergétique CN ont passé sur
  **chacun des deux rangs MPI2** ;
  leurs trajectoires FE sont
  **strictement identiques**.
  Écart maximal MPI1/MPI2 sur
  tous les coefficients FE :
  `2.2297719226571644e-12` ;
  sur H : `5.662137425588298e-15`.
  Le plus grand défaut CN indépendant
  en MPI2 est `9.610368056911511e-16`.
- **Calcul CN long convergé, mais
  pas encore de PASS de bilan complet :**
  MUMPS a bien exécuté **3 000 pas CN**
  à `dt=.005` jusqu'à `t=15`
  (0 rejet, 0 échec non linéaire).
  Le calcul original GetFEM de H a
  lui aussi abouti, puis l'audit
  strict a rejeté une trajectoire
  **incomplète** : le moniteur avait
  conservé seulement 3 000 des
  **3 001 états acceptés**.
  Audit de tous les événements
  PETSc : **seul le pas n=2999,
  t≈14.995**, a été omis
  en raison d'un test flottant
  `t>=t_last_saved+dt_save`.
  Aucune valeur ne doit être
  interpolée ni inventée.
- **Correctif générique
  `scrimp/utils/trajectory_sampling.py`
  + `DPHS.monitor` :** lorsque
  `ts_type=cn`,
  `ts_adapt_type=none` et
  `dt_save=dt`, archiver directement
  **chaque pas TS accepté**, sans
  dépendance au seuil d'horloge
  arrondi ; conserver l'ancien
  échantillonnage contrôlé
  pour les méthodes adaptatives,
  les autres intégrateurs et
  les sorties sous-échantillonnées.
  La résolution PETSc, les
  équations, les données initiales
  et les tolérances physiques/DAE
  ne sont pas modifiées.
  Tests de non-régression contre
  le cas exact du pas 2999,
  avec vérification que les autres
  modes d'enregistrement restent
  inchangés. Diagnostics
  `preflight-rankN.json` conservés
  avant les assertions strictes.
- **Portée : MPI2 physique déjà PASS ;
  bilan complet t=15 encore PENDING**
  jusqu'à la nouvelle exécution
  native sur le même Docker épinglé.
  Sorties VTK concurrentes,
  SuperLU_DIST et modes spectraux
  « whispering gallery » restent
  hors du périmètre.


### Heat–Wave 2D — qualification native MPI1/2 et bilan complet à 15 s : campagne préparée (25 novembre 2026) — protocole exécuté et qualifié depuis

- Nouveau test **en deux phases**
  `scripts/run-v130-heat-wave-mpi-long-physics.sh`
  avec image GetFEM/PETSc/MUMPS
  `scrimp-mpi:v1.3.0-rc3`
  épinglée, échecs conservés et
  ZIP de preuves systématique :
  **4 pas CN dt=.001 en MPI1 et MPI2**,
  puis **3 000 pas CN dt=.005
  jusqu'à t=15 en MPI1** après
  vérification préalable de la parité.
  Aucun solveur alternatif ou
  changement de modèle n'est introduit.
- Option **`physics_only=True`**
  dans le vrai `examples/heat_wave.py` :
  seul le post-traitement graphique et
  VTK est évité pour ne pas lancer
  d'écrivains concurrents en MPI ;
  la PDE, les termes FE, la correction
  `clamp`, le solveur PETSc, les
  contraintes DAE, les deux côtés
  de l'interface et leurs intégrales
  GetFEM sont inchangés. Le comportement
  public normal (figures + ParaView)
  conserve sa valeur par défaut.
- Le collecteur physique orienté
  accepte maintenant **des nombres
  et pas de CN explicitement prescrits**,
  en conservant ses valeurs historiques
  de 4 pas à .001. L'audit long exige
  **dt_save=dt**, soit 3 001 états FE
  pour 3 000 pas effectivement acceptés,
  contrairement aux précédents
  visuels (301 états sous-échantillonnés).
  Les bilans CN chaleur/onde/total,
  la dissipation, les puissances
  d'interface et du bord extérieur,
  le cumul et la monotonie de H
  sont réévalués sur **chaque**
  vrai milieu de pas GetFEM.
- `sandbox/probe_v130_heat_wave_mpi_extended.py`
  mesure la physique sur **tous
  les rangs**, sauve les traces
  FE complètes pour la comparaison
  MPI et les séries énergétiques
  et de puissance à tous les pas
  pour l'horizon long. Vérificateur
  **indépendant** :
  `scripts/verify-v130-heat-wave-mpi-long-physics.py`,
  avec contrôle de parité des
  vecteurs FE, concordance exacte
  des rangs MPI2, et des bilans
  recalculés depuis les fichiers
  natifs. Des tests négatifs de
  falsification (rangs différents,
  bilans modifiés, horloge incomplète)
  protègent les conditions de PASS.
  Aucun calcul simulé n'est
  accepté comme preuve native.
- **État : campagne prête, PASS natif
  MPI2 et bilan complet t=15 ENCORE
  EN ATTENTE.** Le PASS MUMPS MPI1
  initial `b9d56fcc` et les deux
  simulations de visualisation
  longues demeurent acquis dans
  leur périmètre distinct.
  Les exports VTK parallèles,
  SuperLU_DIST et l'identification
  spectrale de modes de whispering
  gallery ne sont pas couverts.


### Heat–Wave 2D — bilan CN physique MUMPS MPI1 PASS et jeux de données 15 s complets, ré-audités sur archives (25 novembre 2026)

- **Premier PASS natif complet** du vrai
  `examples/heat_wave.py`, source
  `b9d56fccc18c1f9facacdd9b6c2645c4de25fc47`,
  Docker `scrimp-mpi:v1.3.0-rc3` :
  `v130-heat-wave-physical-b9d56fcc-3FQ8lFlC.zip`,
  `PASS_HEAT_WAVE_PHYSICAL_MUMPS_MPI1`,
  **10/10 contrôles physiques**,
  **7 689 DOFs**, quatre pas CN fixes
  `dt=0.001` à `t_f=0.004`, MUMPS MPI1,
  `errors=[]`. L'initialisation
  gaussienne `p` est strictement
  compatible après clamp des DOFs
  CG1 extérieurs, sans relâcher
  les contraintes DAE. Les deux
  normales physiques de l'interface
  conforme sont correctement utilisées.
- `H(0)=0.7836421211876778` ;
  `H(0.004)=0.6808984803105300`.
  Dissipation thermique positive et
  non triviale ; défaut maximal de
  puissance totale orientée à
  l'interface `8.881784197e-16`,
  travail au bord extérieur
  `7.733659e-24` ; résidus CN
  chaleur `4.91620633e-15`,
  onde `1.38517670e-15`,
  total `5.45744006e-15`,
  cumul `3.84370313e-15`.
  Conservation/dissipation de la
  puissance **vérifiée par des
  quadratures GetFEM indépendantes**,
  pas déduite du statut PETSc.
- **Deux simulations longues réelles**
  de `3 000 pas CN` chacune
  (`t_f=15`, `dt=0.005`,
  sauvegardes toutes `0.05`),
  MUMPS MPI1 et code d'exemple
  identique : gaussienne
  `H: 0.783642 → 0.00826839`
  (**−98.9449 %**), paquet
  angulaire annulaire
  `H: 0.462644 → 0.377942`
  (**−18.3084 %**).
  `H(t)` reste strictement
  décroissant aux **301 dates**
  réellement sauvegardées pour
  chaque excitation. Les archives
  ont chacune 301 VTU GetFEM pour
  chacun des champs `p,q,T`
  (**903 VTU par simulation**),
  quatre instantanés de l'onde,
  un GIF original de **48 images**,
  les courbes log-log PNG/SVG
  générées par
  `examples/heat_wave.py` et les
  tableaux FE énergétiques originaux.
- La commande native de ces
  simulations s'est bien terminée
  (`native_exit_code=0`), mais
  l'**ancien auditeur post-run**
  a marqué les deux archives
  `FAIL_OR_UNQUALIFIED_LONG_HEAT_WAVE_SHOWCASE`
  par erreur Python XML :
  `frame["file"]` au lieu de
  `frame.attrib["file"]`.
  Corrigé dans
  `scripts/run-v130-heat-wave-long-showcase.sh`
  avec test AST de non-régression.
  Nouveau
  `scripts/verify-v130-heat-wave-long-archive.py` :
  vérifie les **deux ZIP originaux
  sans recalcul ni modification**,
  écrit un sidecar `.reaudit.json`
  indépendant, limite le PASS à
  `PASS_ARCHIVED_LONG_NATIVE_VISUALS_NOT_WG_MODE_PROOF`.
  L'examen indépendant des ZIP
  fournis confirme la complétude
  des données, les temps VTK
  et les animations.
- La fraction cinétique nodale
  approximative dans la couronne
  `r>=0.82` est en moyenne
  `0.587` pour le paquet
  angulaire, contre `0.380`
  pour la gaussienne. Les images
  suggèrent une localisation
  préférentielle, **pas** un
  mode propre de
  *whispering gallery* qualifié :
  pas de mesure de propagation
  tangentielle, ni étude modale
  spectrale, ni preuve d'un
  exposant de décroissance.
- **Limites du PASS** :
  MPI1/MUMPS sur le cas
  CN court ; simulations longues
  et visuels requalifiés **à partir
  d'archives natives complètes**,
  sans nouveau test de bilan
  puissance/sous-systèmes sur
  tout l'horizon `t_f=15`.
  MPI2, SuperLU_DIST, ordre h/dt
  et modes propres restent
  des qualifications distinctes.
  Détails :
  `docs/V130_HEAT_WAVE_PHYSICAL_AUDIT.md`
  et
  `docs/source/examples/heat_wave.rst`.


### Heat–Wave 2D — initialisation au bord corrigée ; visuels log-log temps long et propagation annulaire (25 novembre 2026)

- Le premier pilote physique natif (`937718d2`) construit
  l'interface conforme automatique et assemble les **7 689 DOFs**
  mais **échoue dans la vraie initialisation DAE stricte** avant CN :
  `||F_a||=0.00762918411` pour tolérance
  `8.18210288e-09`, résidus dominants
  `Y_bnd=U_bnd≈0.005395`. Cause physique :
  projection `L²` de la gaussienne
  `p_0` incompatible avec `p|_{r=1}=0`.
  Archive exacte :
  `v130-heat-wave-physical-937718d2-9O1SrxEY.zip`.
  Aucun PASS CN n'en est déduit.
- Correction de l'**exemple canonique**
  `examples/heat_wave.py` : `initial_wave_boundary="clamp"`
  **par défaut** (DOFs CG1 extérieurs exactement nuls,
  DOFs intérieurs, `T_0` et `q_0` inchangés après
  projection FE `L²`), ou
  `"constrained_l2"` (minimisation
  dans la masse GetFEM, PETSc sparse/MUMPS,
  disponible en MPI1). Les deux méthodes
  avaient été validées séparément dans
  les audits natifs de compatibilité
  et DAE/CN/puissance. Les données
  `"unconstrained"` demeurent
  opt-in pour les audits historiques.
  **Aucune relaxation DAE, modification des
  équations, ni repli vers pseudo-TS.**
  Les diagnostics de l'ancien pilote
  ont été raccordés aux vraies deux
  faces physiques orientées plutôt
  qu'à la région logique 10 doublée.
- Le même exemple offre maintenant
  `--long-showcase` avec
  `--final-time 15 --time-step .005 --save-step .05`
  pour représenter **le Hamiltonien GetFEM
  réel en log-log** (PNG/SVG et NPZ de
  la trajectoire énergétique originale)
  et **les champs `p,q,T` originaux
  exportés en VTU/PVD**.
  Les images `p` sur l'anneau
  extérieur et un GIF d'au plus
  48 **vrais instants sauvegardés**,
  avec colorimétrie fixe, sont
  générés directement dans
  `examples/heat_wave.py`, sans
  script graphique externe.
  Un indicateur de concentration
  près de `r>=0.82` est explicitement
  limité à la **fraction approchée
  du terme cinétique scalaire `p²`**
  (pas à l'énergie d'onde complète).
- Nouvelle excitation **optionnelle**
  `initial_wave_profile="gallery_packet"`
  (enveloppe annulaire et `cos(14θ)`,
  chaleur initiale nulle) pour
  explorer les structures de type
  *whispering gallery* en temps long ;
  le cas gaussien historique reste
  **le défaut du modèle**. Il s'agit
  d'une excitation directionnelle
  **non qualifiée comme mode propre**.
  Nouveau pilote
  `scripts/run-v130-heat-wave-long-showcase.sh`
  (MUMPS MPI1, traces et ZIP,
  reprise sur erreur), source
  testée par
  `tests/test_v130_heat_wave_showcase_source.py`
  et `tests/test_v130_heat_wave_boundary_initial_data.py`.
- **Statut scientifique actuel** :
  les anciens PASS d'orientation,
  DAE/CN et parité MPI1/2 dans les
  recherches restent valables ;
  le pilote strict de l'exemple
  corrigé, et les nouveaux
  visuels de temps long, **attendent
  chacun une vraie exécution native**.
  Aucune affirmation de décroissance
  algébrique précise, de convergence
  h/dt ou de mode whispering gallery
  n'est présentée sans preuve.


### Heat–Wave 2D — reprise de la qualification physique native, MUMPS MPI1 strict (25 novembre 2026)

- **Prochain modèle du parcours Sphinx après Maxwell 3D.**
  Nouveau pilote
  `scripts/run-v130-heat-wave-physical-pilot.sh`,
  exécutant uniquement le véritable
  `examples/heat_wave.py` en **MUMPS MPI1**,
  Concentric heat=1/wave=2, **4 pas
  CN fixes** à `dt=0.001`,
  `t_f=0.004`. Il réutilise
  `sandbox/probe_v130_main_2d_examples.py`
  et `sandbox/measure_v130_heat_wave_physics.py`
  sans aucun changement physique.
- **10 portes de qualification, non une simple
  parité de solveurs** : intégration CN native,
  positivité des Hamiltoniens thermique
  et mécanique, parité de leur somme
  avec `H`, dissipation thermique
  positive et non triviale, annulation
  de la puissance signée d'interface
  et du travail extérieur, bilan CN
  pour chacun des deux sous-systèmes,
  bilan CN total et cumulatif,
  énergie globale non croissante.
  Vérification de véritables exports
  GetFEM `T`/`p` et de la figure `H`.
- **Fail closed et traçabilité** : le pilote
  enregistre `summary.json` et un ZIP
  même en cas de rejet PETSc/DAE,
  de problème de trace FE ou de bilan
  physique. Les tests
  `tests/test_v130_heat_wave_physical_pilot.py`
  analysent le programme Python embarqué
  et vérifient la présence des critères
  (source-only, intégrés à la CI).
  Les huit `examples/published/`
  restent strictement immuables.
- **Statut : exécution native encore à faire.**
  Aucun PASS Heat–Wave n'est revendiqué ;
  MPI2, SuperLU_DIST et le temps long
  seront traités séparément après
  acceptation du modèle physique.


### Maxwell 3D — deux patchs GetFEM Nédélec-II k1 / BDM1 et visuels intégrés : MUMPS MPI1/2 PASS (25 novembre 2026)

- **Première qualification native de la paire véritablement patchée**
  `S.FEM("E",1,"NED2")` → `FEM_NEDELEC2(3,1)`
  (patch `0004`) et `S.FEM("B",1,"BDM")`
  → `FEM_BDMK(3,1)` (patch tétraédrique
  `0001`) dans `examples/maxwell_3d.py`.
  Preuve immuable :
  `v130-maxwell3d-3687b11e-gSQ8jeC3.zip`,
  source `3687b11e473813c103cd4e9c02f82812677853d7`,
  image `scrimp-mpi:v1.3.0-rc3` figée.
  **6 tests aliases FEM PASS**, **11 tests
  source Maxwell PASS**, 8 fichiers
  `examples/published/` toujours immuables.
- Boule 3D tétraédrique : **58 points,
  110 tétraèdres**. **444 DOFs**
  Nédélec-II pour `E`, **825 DOFs**
  BDM1 pour `B` : **1 269 DOFs**,
  et non 497 (référence historique
  Nédélec0/RT0 natifs). **MUMPS MPI1
  et MPI2 : 10/10 pas Crank–Nicolson
  et 9/9 contrôles physiques PASS**,
  initialisation DAE absente.
- **Énergie électromagnétique FE conservée** :
  `H0=0.6942520909491591`,
  `Hf=0.6942520909491621` (MPI1).
  `H_E(tf)=0.436003804787049`,
  `H_B(tf)=0.25824828616211315`,
  soit **37.1981 %** de l'énergie
  initiale transférée vers l'induction
  magnétique ; défaut maximal
  `max|H-H0|=4.107825191113079e-15`
  en MPI1, `8.881784197001252e-16`
  en MPI2. Deux intégrales GetFEM
  indépendantes vérifient
  `max|H-H_E-H_B|=5.55e-17`.
- **Gauss magnétique (indépendante) PASS** :
  `max integral(div B)^2=4.2153e-28`
  en MPI1, soit
  `max||div B||_L2=2.0531e-14`.
  **Parité MPI1/2 complète** :
  `max|z_MPI1-z_MPI2|=2.1094237467877974e-15`
  sur les **11 × 1269** coefficients FE ;
  chaque rang MPI2 sauvegarde le même
  tableau et les empreintes SHA256
  concordent avec les métadonnées.
- **Génération de visuels directement
  dans `examples/maxwell_3d.py --export`
  maintenant qualifiée** : chacun
  des deux cas MPI produit les
  courbes `H_E/H_B/H` (PNG+SVG),
  les vecteurs `E` et `B` finaux
  en PNG et le GIF `B(t)` de **11 images**,
  avec `provenance.json` lié au SHA256
  FE original. Les 44 fichiers VTU
  et quatre PVD correspondent aux
  champs et temps originaux 3D ;
  les éléments exportés sont
  **110 véritables tétraèdres**.
- **Portée scientifique** : deux patchs
  GetFEM réellement employés par une
  simulation pH 3D, préservation de
  `H` et de la loi Gauss magnétique,
  reproductibilité MPI et sorties
  graphiques intégrées. Ne démontre
  **ni convergence h/dt, ni Gauss(E),
  ni condition PEC explicitement imposée,
  ni modes de cavité, ni SuperLU_DIST**.
  L'ancien PASS `daac18cb` pour
  Nédélec0/RT0 natifs reste
  une **référence historique distincte**.
  Détails :
  `docs/V130_MAXWELL_3D.md` et
  `docs/source/examples/maxwell_3d.rst`.


### Maxwell 3D : introduction de la paire patchée et des raccourcis FEM (historique avant qualification, 25 novembre 2026)

- **Correction scientifique importante :**
  le précédent `FEM_NEDELEC(3)` est
  le **Nédélec première espèce k=0 natif GetFEM**,
  et `RT0` est aussi natif. L'archive
  `v130-maxwell3d-daac18cb-FCPngFrH.zip`
  démontre authentiquement Maxwell 3D,
  la conservation du Hamiltonien,
  Gauss(B), la parité MPI1/2 et les
  VTK tétraédriques ; **elle ne teste pas
  les familles créées par nos patchs**.
  La qualification historique reste valide
  mais n'est plus attribuée à la paire courante.
- `scrimp/fem.py` accepte désormais
  **`S.FEM("E", k, "NED1")`** et
  **`S.FEM("E", k, "NED2")`**
  (synonymes `NEDELEC1` / `NEDELEC2`)
  qui construisent directement
  `FEM_NEDELEC1(N,k)` et
  `FEM_NEDELEC2(N,k)`.
  Numéro 1/2 = **espèce**, entier k =
  **indice polynomial**. NED1 k0 est
  le Nédélec natif, NED1 k>=1
  utilise le patch `0003` ;
  NED2 k=1,2 utilise le patch `0004`.
  Les ordres non implémentés pour
  NED2 sont rejetés explicitement.
  Les deux nouvelles familles
  **n'autorisent aucun repli silencieux CG1**
  en cas de patch GetFEM manquant.
  `RT` et `BDM` conservent leur API.
- Le code canonique `examples/maxwell_3d.py`
  utilise maintenant
  **`NED2`, k=1** (`FEM_NEDELEC2(3,1)`)
  pour `E`, et
  **`BDM`, k=1** (`FEM_BDMK(3,1)`,
  patch tétraédrique `0001`)
  pour `B`. La propriété locale
  `curl([P1]^3) ⊂ BDM1` maintient
  la compatibilité magnétique à tester.
  **Aucun changement des quatre briques pH,
  de la commande de bord nulle, du CN ou
  des neuf critères énergétiques/Gauss**.
- Ajout de
  `tests/test_v130_fem_nedelec_aliases.py`
  (descripteurs, synonymes, plages d'ordre,
  défaut de patch sans CG1, familles RT/BDM
  préexistantes), adaptation des
  tests source Maxwell et du pilote pour
  **exiger explicitement la paire patchée**,
  mise à jour de la CI et des docs Sphinx.
  Le pilote conserve l'archive d'échec
  et produit les PNG/SVG/GIF par
  `examples/maxwell_3d.py --export`.
- **Statut au moment de cette modification :
  qualification encore en attente.**
  Ce stade historique est désormais
  **dépassé par le PASS natif `3687b11e`**
  consigné ci-dessus ; les résultats
  de la référence Nédélec0/RT0
  (497 DOFs) ne sont jamais extrapolés.


### Maxwell 3D — visuels automatiques produits par l'exemple lui-même (25 novembre 2026)

- Le modèle qualifié `examples/maxwell_3d.py`
  crée désormais **sans script de post-traitement**
  ses propres figures et son animation dès
  `--export --output <répertoire>` :
  courbes `H_E/H_B/H` en **PNG/SVG** issues
  des intégrales FE originales ; vraies
  vues vectorielles 3D finales `E` et `B`
  en PNG ; **GIF B(t)** de toutes les
  images CN sauvegardées. Les flèches
  utilisent exclusivement les valeurs
  nodales des VTU/PVD tétraédriques
  **GetFEM authentiques** ; aucune
  interpolation en temps ou champ
  synthétique. Les échelles 3D et
  la caméra sont fixes dans le GIF.
- Le seul **rang MPI0** effectue le rendu
  après les exports GetFEM collectifs.
  Les fichiers sont déposés dans
  `<output>/visuals/` et enrichis d'un
  `provenance.json` référençant
  l'empreinte SHA256 de la trajectoire
  native, le nombre de nœuds/tétraèdres,
  les instants d'animation, et un
  avertissement sur les vecteurs nodaux
  (qui ne sont pas les DOFs Nédélec/RT0).
  Aucune nouvelle dépendance `vtk` :
  lecture VTK inline-base64 directement
  dans le même exemple avec la pile
  Matplotlib/Pillow déjà présente.
- Le pilote
  `scripts/run-v130-maxwell3d-pilot.sh`
  **échoue explicitement** si les figures,
  le GIF (11 images attendues) ou la
  provenance native manquent, si le
  lecteur rencontre une géométrie
  non tétraédrique, ou si les champs
  des VTK divergent de la série CN.
  Tests de protocole VTK natif et de
  générateur `--export` ajoutés à
  `tests/test_v130_maxwell3d_source.py`.
  **Aucune modification des équations,
  des seuils physiques ou du noyau.**
- **Qualification** : le calcul Maxwell
  3D MUMPS MPI1/2 du commit
  `daac18cb` est validé ; la
  **nouvelle sortie graphique intégrée**
  attend son propre pilote natif.


### Maxwell 3D Nédélec0/RT0 natifs — qualification historique de référence MUMPS MPI1/MPI2 PASS (25 novembre 2026)

- **PASS natif** sur `examples/maxwell_3d.py`,
  source `daac18cb36fbc404fcf5c4fe4f2f69121b2ee530`,
  image `scrimp-mpi:v1.3.0-rc3` figée,
  archive originale
  `v130-maxwell3d-daac18cb-FCPngFrH.zip`.
  La boule tridimensionnelle est réellement
  maillée en **110 tétraèdres**, 58 points.
  **222 DOFs Nédélec H(curl)** pour `E`,
  **275 DOFs RT0 H(div)** pour `B`,
  soit **497 DOFs**, et non des CG vectoriels.
- **MUMPS MPI1 et MPI2 : 10/10 pas CN**
  à `h=0.7, dt=0.02, tf=0.2`,
  aucun échec ni rejet, **9/9 critères
  physiques PASS** sur chaque configuration.
  Hamiltonien initial MPI1
  `H0=0.6942520909491595`, final
  `Hf=0.6942520909491590`.
  Défaut cumulé maximal : `2.33146835e-15`
  en MPI1, `1.99840144e-15` en MPI2.
- **Échange électromagnétique réel** :
  initialement `H_E=0.6942520909491595`,
  `H_B=0`; à `tf=0.2`,
  `H_E=0.4359826378118022` et
  `H_B=0.2582694531373567`
  (magnétique : **37.2011 %** de l'énergie
  initiale). Les deux intégrales énergétiques
  GetFEM indépendantes correspondent au
  Hamiltonien à `5.55e-17` près.
  `max ||div B||_L2=4.63e-15` en MPI1
  (Gauss magnétique discrète vérifiée).
- **Parité MPI** : les deux rangs MPI2
  donnent exactement la même trajectoire
  binaire sauvegardée ; écart maximal
  sur les trajectoires complètes MPI1/2
  `4.440892098500626e-16`,
  écart Hamiltonien `2.7755575616e-15`.
  Empreintes SHA256 des 3 trajectoires
  vérifiées sur les `rank*.npz` originaux.
- **VTK 3D authentiques** : **44 VTU**
  (deux vecteurs `E,B`, onze temps,
  MPI1 et MPI2) et quatre PVD,
  indépendamment ouverts par un vrai
  lecteur VTK. Chacun contient
  58 points, **110 tétraèdres VTK type 10**
  et trois composantes finies du
  champ correspondant. Le maillage
  des deux rangs est géométriquement
  identique et `B` final non nul.
- **Portée mesurée** : première démonstration
  port-Hamiltonienne 3D, conservatrice
  avec fermeture naturelle à puissance
  nulle, sans DAE artificielle, et
  éléments intrinsèques natifs GetFEM
  (Nédélec première espèce k=0 et RT0).
  **Ce PASS ne valide aucun élément ajouté
  par les patches SCRIMP** ; voir la nouvelle
  paire NED2/BDM1 en tête du changelog.
  Ni convergence h/dt, ni Gauss électrique,
  ni condition PEC imposée, ni modes de
  cavité, ni SuperLU_DIST qualifiés.
  Le pilote et les 8 exemples publiés
  restent inchangés. Détails :
  `docs/V130_MAXWELL_3D.md` et
  `docs/source/examples/maxwell_3d.rst`.


### Maxwell 3D — second pilote : introspection GetFEM retirée sans assouplissement du garde-fou (25 novembre 2026)

- Archive `v130-maxwell3d-ed9b9fd8-gnJLLDYD.zip`
  (source `ed9b9fd80f8ceb670085d7a58dc07f625dedbe99`) :
  **8/8 tests statiques PASS**, 8 exemples
  publiés intacts, boule de 58 points et
  110 tétraèdres construite ; Nédélec 3D
  **222 DOFs** et RT0 3D **275 DOFs**
  créés sur 1 et 2 rangs. Les exécutions
  n'ont **pas atteint la discrétisation
  temporelle**. Le vérificateur d'archive
  a correctement enregistré FAIL.
- Nouvelle erreur de l'**audit de l'exemple**
  et non preuve d'une faute du patch GetFEM :
  `MeshFem.fem()` n'exposait pas les objets
  `Fem` selon le format supposé après
  `np.asarray(...,dtype=object).ravel()`.
  Le modèle échouait avant assemblage avec
  `Invalid GetFEM MeshFem.fem() entries`.
- **Correctif robuste** : plus aucune
  dépendance à `MeshFem.fem()`.
  Conserver les objets originaux `S.FEM`
  passés à `DPHS.add_FEM()` ; leur
  **type et ordre** sont mutés en `CG/1`
  par l'ancienne politique de fallback
  de `scrimp/fem.py`. Rejeter tout
  écart à `FEM_NEDELEC(3)/-1` et
  `RT/0` ; contrôler les vraies
  dimensions de `MeshFem`.
  Tests de régression sur la mutation
  CG, la mauvaise famille, mauvais ordre
  et dimensions erronées. Aucune
  modification du noyau ni relaxation
  d'un critère physique.
- **Statut : Maxwell 3D MUMPS MPI1/2
  toujours NON QUALIFIÉ** ; une nouvelle
  exécution de `scripts/run-v130-maxwell3d-pilot.sh`
  est requise.


### Maxwell 3D : correctif de préqualification de l'audit FE/archives (25 novembre 2026)

- Première exécution native `acdb3c75` :
  **6 tests source PASS**, 8 exemples publiés intacts,
  58 sommets et 110 tétraèdres, **Nédélec 3D 222 DOFs**
  et **RT0 3D 275 DOFs** créés avec GetFEM 5.5
  sur MPI 1/2. **Aucun pas Maxwell CN exécuté**.
- Blocage dans le garde-fou Python :
  `MeshFem.fem()` rend un tableau NumPy de
  `Fem`; la boucle traitait une ligne comme
  un `Fem`, puis appelait `.char()` sur
  `numpy.ndarray`. Correctif : aplatissement
  du tableau, vérification de toutes les
  descriptions `Fem.char()` et rejet de
  toute substitution CG. Aucune modification
  du noyau, de la formulation ou des patchs.
- Correction distincte du script d'archives :
  `from pathlib import Path` au lieu du
  mauvais `from pathlib`. Nouveau test
  source : **AST du programme Python embarqué
  dans le heredoc Bash**, plus tests unitaires
  de la véritable fonction de contrôle FE
  sur des retours NumPy simulés.
- **Statut Maxwell 3D : toujours NON QUALIFIÉ**
  tant que le pilote MUMPS MPI1/2 n'a pas
  traversé les pas CN, les bilans de puissance,
  Gauss(B) et les vrais exports VTK.

### Maxwell 3D minimal avec éléments de Nédélec et RT0 — pilote préparé, non qualifié (25 novembre 2026)

- Nouvel exemple indépendant `examples/maxwell_3d.py` :
  **vraie géométrie tétraédrique 3D** de la boule,
  `E in H(curl)` (Nédélec 0),
  `B in H(div)` (Raviart–Thomas RT0).
  L'exemple exercera les éléments d'arête
  GetFEM corrigés pour SCRIMP v1.3.0,
  sans modifier `examples/published/`.
  L'adaptateur SCRIMP vérifie activement les
  familles FE effectivement assignées afin
  d'interdire un retour silencieux à CG1.
- Premier Maxwell pH **fermé** : matière
  homogène `eps=mu=1`, sans source ni
  contrôle, deux masses et deux briques
  de rotationnel transposées opposées,
  soit `M_E E_dot=C^T B`,
  `M_B B_dot=-C E`.
  Hamiltonien physique
  `H=1/2 integral(|E|²+|B|²)` :
  conservation visée `H(t)=H(0)`.
  Fermeture naturelle sans puissance de
  bord ; pas de fausse condition PEC ni
  de port fictif. Pas d'initialisation
  DAE : les deux états sont différentiels.
- Initialisation non triviale
  `E0=(-y,x,0)`, `B0=0`,
  rotationnel `curl(E0)=(0,0,2)`,
  projection L2 dans les vrais FE.
  Critères : conservation H par pas et
  en cumulé, deux quadratures GetFEM
  indépendantes, apparition du champ B,
  et Gauss magnétique `div(B)=0`.
  La loi Gauss électrique et la
  convergence ne sont pas qualifiées.
- Pilote natif
  `scripts/run-v130-maxwell3d-pilot.sh`,
  MUMPS MPI1/2, 10 pas CN fixes,
  `h=0.7`, `dt=0.02`, `tf=0.2`,
  9 portes scientifiques, parité complète
  FE entre rangs et MPI1/MPI2,
  vrais exports tétraédriques E/B
  (VTK types 10/24) et archive d'échec
  préservée. Contrat source/CI
  `tests/test_v130_maxwell3d_source.py`,
  documentation Sphinx
  `docs/source/examples/maxwell_3d.rst`
  et note `docs/V130_MAXWELL_3D.md`.
- **Statut : source prête, exécution native
  sur image `scrimp-mpi:v1.3.0-rc3`
  attendue. Aucun PASS de Maxwell 3D
  n'est encore revendiqué.**


### Exemple 2B — démonstration longue MUMPS sur 100 pas, dissipativité et comparaison native 2A/2B PASS (25 novembre 2026)

- Archive complète `v130-wave2d-2b-showcase-dbdc49cf-G1CKQAVG.zip`,
  source `dbdc49cfcac47d4daa1706b3fc2fa4b53646b0a2`,
  image `scrimp-mpi:v1.3.0-rc3` figée,
  MUMPS MPI1, `h=0.1`, `dt=0.01`,
  `tf=1`, `damping_strength=1` :
  **10204 DOFs**, **100/100 pas CN** sans rejet
  ni échec SNES, DAE algébrique PASS à résidu
  exactement nul, **10/10 portes physiques PASS**.
- **Bilan FE natif** `H-H0=W-L` sur 101 instants :
  `H(0)=0`, `H(1)=0.03205555671004583`,
  `W(1)=0.03710798721523879`,
  `L(1)=0.005052430505192933`.
  Défaut absolu maximal par pas
  `9.391901319e-17`, cumulé
  `9.280770596e-17`. Positivité et
  non-trivialité de `L`, puissance du port
  résistif identique à une intégrale GetFEM
  volumique indépendante, flux constitutif
  contraint faiblement. Trois ports de bord
  non forcés de puissance exactement nulle.
- **606 VTU** (6 champs `q,p,e_q,e_p,f_r,e_r`
  × 101 instants) et 6 PVD originaux,
  temps de 0 à 1 conformes aux tableaux FE ;
  SHA256 de trajectoire identique à
  `rank0.json`. Visuels GetFEM PNG/SVG,
  GIF `p_evolution.gif` (101 images),
  énergie CSV (101 lignes), provenance.
- **Comparaison indépendante 2A/2B** sur les
  deux archives longues (mêmes matériaux,
  maillage, excitation et pas de temps).
  Identité des VTK `p` en géométrie et
  connectivité vérifiée à `t=0,0.5,1`.
  À `tf=1`, `H_2A=0.03707053141235693`,
  `H_2B=0.03205555671004583` : **13.528 %**
  d'énergie en moins pour 2B.
  Les travaux injectés sont légèrement
  différents ; l'identité
  `H_2A-H_2B=L_2B-(W_2B-W_2A)` tient au
  défaut `2.08e-17` près. Deux figures de
  comparaison ont été produites depuis
  les VTK/NPZ natifs, dont une carte `p`
  **à échelle de couleurs identique**.
  `p` est la quantité de mouvement, non
  la vitesse `e_p`.
- **Portée distincte** : long showcase 2B
  validé en MUMPS MPI1 ; court pilote
  MUMPS MPI1/2 reste une autre preuve.
  Aucune convergence spatiale/temporelle,
  ni parité MPI longue, ni qualification
  SuperLU_DIST n'est déduite ici.
  Le modèle canonique `examples/wave_dissipative.py`,
  `examples/published/` et `master`
  n'ont pas été modifiés.


### Exemple 2A — démonstration longue native de l'onde conservative 2D : PASS (25 novembre 2026)

- Archive GetFEM authentique
  `v130-wave2d-2a-showcase-bd6f5468-teEK1Ejz.zip`,
  source `bd6f5468e3a167f2b2a3b9e1b72bcc54e1948d45`,
  Docker `scrimp-mpi:v1.3.0-rc3`, **MUMPS MPI1**,
  `h=0.1`, `dt=0.01`, `t_f=1`.
  **8138 DOFs, 100/100 pas CN** sans rejet,
  init DAE algébrique PASS avec résidu initial
  et final exactement nuls, **7/7 critères physiques PASS**.
- **Bilan énergétique signé original** :
  `H(0)=0`,
  `H(1)=0.03707053141235693`,
  `W_boundary(1)=0.03707053141235698` ;
  `H-H0=W_boundary` à tous les instants.
  Défaut maximal par pas `2.0654051e-16`
  et cumulé `1.2490009e-16`.
  Les trois bords Neumann non forcés ont une
  puissance exactement nulle ; l'énergie
  accumulée vient de la commande de vitesse
  au bord gauche (système conservatif **ouvert**).
- Trajectoire `z` complète vérifiée par SHA256.
  Les **404 vrais VTU** de `q,p,e_q,e_p`
  (101 instants chacun), leurs PVD, le CSV
  de 101 états, les cartes de champs PNG/SVG,
  le graphique de puissance et le GIF de
  101 images ont été inspectés. Le GIF montre
  une perturbation entrant depuis la gauche
  et se propageant dans le rectangle anisotrope.
- **Terminologie FE** : `q=∇w` = gradient
  du déplacement, `p=ρ∂_tw` = quantité de
  mouvement ; la vitesse réelle est
  `e_p=p/ρ` et la contrainte est `e_q=Tq`.
  Les cartes génériques `p,|q|` représentent
  momentum et gradient, non ces coénergies.
- **Portée limitée à cette simulation longue
  MUMPS MPI1**, distincte du pilote 2A court
  déjà qualifié MPI1/2 et autres solveurs.
  Pas de convergence h/dt ni de parité
  MPI longue revendiquée. Prochaine étape :
  long showcase **2B** avec mêmes paramètres
  géométriques, matériaux, excitation de bord,
  `h=0.1`, `dt=0.01`, `t_f=1`.


### Exemple 2C — démonstration longue native MUMPS validée sur 100 pas CN (25 novembre 2026)

- Archive complète `v130-wave2d-2c-showcase-bc880100-9P5gcxWf.zip`,
  source `bc880100f69f6f0dbbb5d31ee8fc229d0b32647d`,
  image `scrimp-mpi:v1.3.0-rc3`, MUMPS MPI1 :
  **7905 DOFs**, `h=0.1`, `dt=0.01`,
  `tf=1.0`, **100/100 pas CN** sans rejet
  ni échec non linéaire. Projection DAE algébrique
  PASS (résidu final `9.86e-14`).
- **10/10 critères physiques PASS** :
  `H(0)=0.2556143692106721`,
  `H(1)=0.10272364224930962`, baisse de
  `59.813 %` ; `Li=0.0882892889133343`,
  `Lb=0.06460143804802765` (respectivement
  `57.75 %` et `42.25 %` des pertes).
  Défauts absolus maximaux par pas et cumulés
  `1.32e-15` et `1.62e-15`; contrôle
  indépendant des 101 vrais états FE, énergie
  monotone, pertes positives, SHA256 identique
  à `rank0.json`.
- **303 VTU GetFEM** et 3 PVD originaux
  (`q,p,e_r` × 101 instants), confirmés en
  correspondance exacte avec les 101 temps
  sauvegardés ; les images PNG/SVG, le GIF
  `p_evolution.gif` de 101 frames et les
  101 lignes CSV ont été produits et examinés.
  Les visualisations démontrent la propagation
  et la montée progressive du feedback au bord.
- **Amélioration graphique à prévoir** :
  marquer ou masquer l'extérieur du disque
  `region=1` dans la carte `e_r`, qui
  apparaît aujourd'hui comme un zéro hors
  de son domaine de définition. Aucune erreur
  FE ou énergétique détectée.
- **Portée : MUMPS MPI1 uniquement pour cette
  démonstration longue** ; le court pilote 2C
  MUMPS MPI1/2 reste distinct et qualifié.
  Pas de convergence h/dt, de parité MPI longue,
  ni de qualification SuperLU_DIST déduite
  de cette archive. Les sources publiées
  et `master` restent inchangés.


### Exemple 2C : qualification physique complète MUMPS MPI1/MPI2 PASS (25 novembre 2026)

- Archive GetFEM/PETSc native `v130-wave2d-coenergy-eeae81f7-siWdW8ft.zip`,
  commit exact `eeae81f7099e61ed6edf62941c782e9d2ea56f2d`,
  image `scrimp-mpi:v1.3.0-rc3` :
  **10/10 pas CN sur 2459 DOFs**, MUMPS MPI1 et MPI2,
  aucun rejet ni échec non linéaire, initialisation
  DAE algébrique PASS (états différentiels préservés).
  Les dix critères scientifiques sont PASS des deux côtés.
- La vraie compliance anisotrope `T^{-1}`, les deux
  régions concentriques, la loi intérieure `e_r/nu=-p`
  (projection faible) et le feedback de bord passif
  `U=+0.5Y` sont simultanément validés :
  `H(tf)-H(0)+Li(tf)+Lb(tf)=0`.
  `Li=0.035508373983860454`,
  `Lb=0.0010393723378964645` ;
  `H(0)=0.25557110317334153`,
  `H(tf)=0.21902335685158436` (MPI1).
  Les pertes internes et externes sont positives et
  cohérentes avec les quadratures FE indépendantes.
- Erreurs énergétiques absolues maximales par pas :
  `9.81e-16` (MPI1), `4.39e-16` (MPI2).
  Défauts cumulés maximaux : `7.17e-16` /
  `2.60e-16`. Écart maximal sur toutes les
  trajectoires FE `z` MPI1/2 : `2.687e-11`.
  Les deux rangs MPI2 ont des tableaux FE
  numériquement identiques et leur SHA256 concorde.
- Les **66 VTU GetFEM authentiques** (3 champs
  `q,p,e_r` × 11 instants × 2 exécutions)
  et six séries PVD ont été vérifiés.
  Le pilote a produit et archivé les visuels
  `visuals/field_snapshots.{png,svg}`,
  `visuals/energy_balance.{png,svg}`,
  `visuals/p_evolution.gif`, le CSV énergétique,
  le README et la provenance. Le CSV et les images
  ont été examinés ; les figures utilisent les
  échantillons VTK Float32, tandis que les bilans
  reposent sur les intégrales GetFEM originales.
- **Portée : seulement MUMPS MPI1/2, h=0.2,
  dt=0.008, tf=0.08 et dix pas CN**.
  Reste à produire une démonstration longue
  (`scripts/run-v130-wave2d-showcase.sh 2c`),
  à examiner l'animation et à évaluer séparément
  convergence, temps long et autres solveurs.
  Pas de qualification SuperLU_DIST 2C,
  ni de changement de `master` ou
  `examples/published/`.


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
- **Statut historique de cette étape :** simulation
  réparée alors attendue ; **désormais PASS en
  MUMPS MPI1/2** (nouvelle entrée ci-dessus).
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
