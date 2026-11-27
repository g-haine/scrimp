# SCRIMP v1.3.0 — état des lieux et reste à faire (27 novembre 2026)

**Document de synthèse courant**, distinct de la chronologie des essais.
Références de portée : [plan de release](V130_RELEASE_PLAN.md),
[roadmap technique](V130_DEVELOPMENT_ROADMAP.md),
[sprint 2026-W48](V130_RELEASE_SPRINT_2026W48.md),
[politique des exemples](V130_EXAMPLES_SOLVER_POLICY.md).
Branche : `release/1.3.0-preparation`. **Aucun tag v1.3.0 publié.**
Les affirmations scientifiques ci-dessous concernent les configurations
explicitement testées, jamais une garantie universelle pour tout modèle.

## Dernière qualification de distribution — 27 novembre 2026

**PASS GHCR depuis une VM Ubuntu indépendante.** Le candidat
`ghcr.io/g-haine/scrimp:1.3.0-dist-rc3`, construit depuis la révision
`58dbe187de308979938153f9adfcf94dac4df85a`, a été publié puis
récupéré par digest OCI immuable
`sha256:0cbbce6ad42b23552aff9d12555342f3864b50522054fe039c161a7145b2732a`
sur un runner GitHub Ubuntu 24.04/AMD64 sans image préexistante.
Le ZIP utilisateur indépendant, la vraie simulation avec CSV persistant
sur l'hôte, le démarrage MPI2 et le modèle Heat--Wave natif MPI1/MPI2
ont tous PASSÉ. Parité relative L2 des 11 Hamiltoniens FE :
`4.275706381199414e-15`. Preuve :
[GitHub Actions #46025625273 (tentative 2)](
https://github.com/g-haine/scrimp/actions/runs/46025625273).

**PASS de la chaîne documentaire sur `0cdb2ad8`** :
tests sources, wheel, génération HTML et LaTeX strictes et compilation
PDF sur GitHub Actions (workflow #46026382133).
Le contrôle complémentaire `sdist` construit depuis une source
isolée a également passé son étape CI, avec reconstruction du wheel
et comparaison octet pour octet des fichiers Python.

Les preuves de performance, physique MPI et GUI Ubuntu sont détaillées
plus bas avec leur périmètre propre. **Restent à faire avant le tag
stable :** décider du sort des issues ouvertes, publier l'image
`1.3.0` et le ZIP officiel correspondant, vérifier leur cohérence,
et préparer la GitHub Release. La GUI SCRIMP et l'audit énergétique
du réservoir Saint--Venant sont reportés après la v1.3.0.

## Qualification Docker-first après installation en ZIP — 27 novembre

**PASS sur Ubuntu/amd64 pour l'image locale
`scrimp-mpi:v1.3.0-dist-rc2`**, construite depuis
`91d7eb1a12329d54854f13cd60d7860c36b856a4` :
contrôle de provenance OCI/versions, wheel v1.3.0, GetFEM
MPI2/PETSc/SLEPc/MUMPS, smoke MPI1/MPI2, ZIP autonome,
simulateur réel à 10 pas BDF2/11 états, CSV physiquement
cohérent et persistant sur l'hôte. Le gate de distribution a
rendu `PASS_STANDALONE_ZIP_IMAGE_SERIAL_MPI_STARTUP`.

**PASS du parcours Ubuntu VS Code Dev Containers sur le ZIP
autonome réel** : après fermeture puis réouverture de VS Code,
le projet a été exécuté ; l'utilisateur a confirmé **F5**
jusqu'à un point d'arrêt, **F10** instruction par instruction,
puis **F5** jusqu'à la fin de la simulation. Le fichier
`oscillator.csv` transmis permet la vérification de ses
11 états et de l'identité énergétique. Ce PASS concerne la
GUI locale du candidat testé, **pas** une installation depuis
GHCR ni le débogage MPI multiprocessus.

La CI Python source a obtenu **SUCCESS** sur
`7bdfb2685e9faae06a8d68a7197a6475a093a127`
après normalisation des seuls bandeaux de licence de
89 fichiers Python suivis (aucune modification des équations) ;
une nouvelle CI complète a également obtenu **SUCCESS**
sur `e81003b7ef32360c755a4db56b65d4093b744aea`.

**PASS du calcul Heat–Wave réellement embarqué dans
`scrimp-mpi:v1.3.0-dist-rc2`**, sans montage du dépôt :
maillage 3 154 éléments, 7 689 DoFs, CN de 0 à 0,02
par 10 pas `dt=0.002` acceptés, 0 échec de solveur,
11 Hamiltoniens FE en MPI1 et MPI2 ; écart relatif L2
`4.252240939683056e-15` sous la tolérance fixée
`1e-7`. Statut
`PASS_NATIVE_HEAT_WAVE_MPI1_MPI2_HAMILTONIAN_PARITY`.
La portée est la **parité de la trajectoire énergétique**,
et non une parité champ FE ou MPI générale. Preuves :
`sandbox/outputs/v130-distribution-heat-wave-mpi-wb9SrBRd/`.

**PASS documentaire 1.3.0 sur `e81003b7`** :
8/8 tests, Sphinx HTML et LaTeX avec
`-n -W --keep-going` et **zéro avertissement Sphinx**
de chaque côté, puis `latexmk` produisant le PDF réel
`SCRIMP, Release 1.3.0` de **168 pages**. Statut
`PASS_HTML_LATEX_PDF_BUILT`. L'archive de 258 entrées
a passé le CRC ; le PDF et le rapport sont présents.
Le journal LaTeX final ne contient plus de références
non résolues, mais conserve des avertissements typographiques
(`hyperref`, boîtes trop pleines). Preuves :
`sandbox/outputs/v130-docs-build-e81003b7-JoPsDG4H/`.

**Encore ouverts** : image définitive depuis la révision
de sources stabilisée, publication privée/prérelease GHCR
suivie d'un `docker pull` **par digest** sur hôte vierge
et du test ZIP externe, revue des issues et décision de
release. Aucun tag ni image publique à ce stade.

Voir [porte de qualification Docker-first](V130_DOCKER_FIRST_RELEASE_GATE.md).

## Décision Docker-first et compilation qualifiée — 27 novembre, après-midi

**Décision v1.3.0 : Docker précompilé devient l'installation recommandée.**
Le Docker construit localement reste l'option développeur et Conda
natif l'alternative avancée/existante. La branche de préparation
`feature/v130-docker-first-distribution` organise la distribution,
sans changer les équations ni les sources `examples/published/`.
La future image publique, son digest et le ZIP utilisateur **ne sont
pas publiés** ; le socle qualifié rc3 reste inchangé.

Le build documentaire **effectivement réalisé** sur la source
`cf75f4f761c2854788e5ef519375b31e7f040509` a passé
`test_v130_docs_release_source.py` (8/8), Sphinx HTML et LaTeX en
mode `-n -W --keep-going` (0 avertissement de chaque côté), et
`latexmk -pdf` : PDF réel **168 pages**, statut
`PASS_HTML_LATEX_PDF_BUILT`. Les quelques avertissements
de typographie pdfLaTeX ne sont pas des avertissements Sphinx et
n'empêchent pas le build. Le titre de la version dans le PDF
restait `1.2.0`, parce que Sphinx et `setup.py` étaient
encore à `1.2.0` ; tous deux sont maintenant passés à
`1.3.0` **sur la nouvelle branche** et requièrent un build
de confirmation avant la release.

Deux barrières de distribution sont préparées, sans PASS natif
revendiqué : `scripts/audit-v130-wheel.py` teste le wheel et
son installation isolée (métadonnées seulement, sans bibliothèques
natives), et `scripts/run-v130-docker-distribution-gate.sh` vérifie
la provenance exacte, le package installé dans l'image, le ZIP
autonome, le vrai oscillateur et son CSV sur l'hôte, et le
démarrage MPI à deux rangs. Voir
[V130_DOCKER_FIRST_RELEASE_GATE.md](V130_DOCKER_FIRST_RELEASE_GATE.md).
L'installation sur hôte vierge, l'acceptation GUI du nouvel artefact,
une qualification d'image finale et la publication GHCR
restent ouvertes. **Aucun tag v1.3.0 autorisé à ce stade.**

## Mise à jour documentation et sortie des exemples — 27 novembre, fin de matinée

Le guide public [`energy_and_power.rst`](source/energy_and_power.rst)
documente désormais les vraies conventions du code :
échantillonnage Hamiltonien sur les états GetFEM,
CN par pas *seulement* si chaque état TS accepté est sauvegardé,
CN sous-échantillonné, BDF/BDF2 et ARKIMEX avec **trapèzes
nodaux approximatifs** ; aucune conservation ou passivité
exacte n'est annoncée sans identité discrète propre au modèle.
Les pages `index.rst`, `examples.rst` et les tutoriels canoniques
ont été actualisés, ainsi que `install.rst`, `vscode.rst`,
`GUI.rst` et `README.md`. Le brouillon de GUI n'est plus
présenté comme une fonction qualifiée.

**Preuves natives distinctes et chronologiquement exactes :**

- Source `a471a637` : 9/9 **PASS** pour les vrais calculs
  d'exemples et leurs sorties, mais **dossiers parasites
  `outputs/{log,mesh,png,pv}` présents**.
- Source `ab6c6fd` : **Heat–Wave MPI2 UNIQUEMENT PASS**
  du correctif de structure ; exactement
  `runs/heat_wave_mpi2/outputs/heat_wave/`,
  aucun dossier générique voisin, exactement un
  `Concentric.msh` dans `mesh/`, deux PVD
  et 5 PNG / 5 SVG / 2 NPZ, ZIP CRC PASS.
- Source `20d79d0d` : **9/9 PASS FINAL** après
  correction de l'arborescence, y compris Heat–Wave
  MPI2 ; chaque dossier de travail possède
  uniquement `outputs/<system_name>/` et aucun
  `outputs/{log,mesh,png,pv}` générique ;
  maillage Gmsh unique, PVD/VTU valides, ZIP CRC PASS.
  **P0 exemples clos.**

**Compilation documentaire :** le workflow GitHub Actions
`.github/workflows/documentation.yml` est déclenché sur les
pushs et construit HTML, LaTeX puis PDF. Un test de
parité des textes, commandes et chemins est maintenant
exécuté avant le build. Un pilote local non destructif
`scripts/run-v130-docs-build.sh` reconstruit réellement
HTML, LaTeX et PDF dans `sandbox/outputs/`,
conserve les logs et produit un résumé avec le commit
source. Le **premier build effectif** du commit
`20d79d0d` a produit HTML et LaTeX (39
avertissements Sphinx), mais le PDF a **ÉCHOUÉ**
sur le caractère `U+2212` dans du code verbatim
(`pdflatex`). Les titres/docstrings/références
en cause ont été corrigés en source ; le **second build**
sur `cf75f4f7` a ensuite **PASSÉ HTML, LaTeX et PDF**,
avec 0 avertissement Sphinx (détails en tête du présent
document). Le pilote et la CI traitent les avertissements
Sphinx comme des échecs par défaut.

**Après ces deux gates :** vérifier le packaging/wheel,
l'installation propre et les versions Sphinx/`setup.py`
(qui sont encore `1.2.0`), puis décider des issues
et préparer la release. Les travaux #28 et
l'optimisation complète Docker restent hors chemin
critique par défaut.

## Historique — sortie des huit exemples (avant acceptation finale)

Un chantier transversal a été intégré à `release/1.3.0-preparation` :

- `DPHS.set_output_name(...)` et `get_output_path(...)` normalisent
  les fichiers des systèmes utilisateurs sous
  `working_dir/outputs/<system_name>/{png,svg,gif,mesh,log,npz,csv,json}`.
  Les anciennes API et les dossiers scientifiques de provenance
  restent lisibles et non déplacés.
- `plot_Hamiltonian` écrit PNG/SVG et traces NPZ ; le cache
  `compute_Hamiltonian` est recalculable sans duplication et
  respecte la véritable horloge du modèle GetFEM.
- `compute_powers` sélectionne `cn_step` **uniquement**
  quand tous les pas PETSc acceptés ont été sauvegardés,
  sinon une observation nodale et une intégration des
  puissances explicitement **approximative**.
- `examples/heat_wave.py` expose désormais directement
  les énergies chaleur/onde, les puissances sur les deux
  faces orientées, la dissipation, ainsi que le travail
  physique ; ses fichiers proviennent du modèle canonique
  et des vraies quadratures GetFEM, jamais d'un pilote.
- Les huit exemples canoniques ont été adaptés au nouvel
  agencement, sans modification de leurs formes faibles
  ni des huit sources immuables `examples/published/`.

La première campagne native 9/9 sur `a471a637`
a confirmé les sorties mais révélé des dossiers
parasites. Après le correctif de création des
dossiers et l'écriture directe des maillages,
la **seconde campagne complète 9/9** sur
`20d79d0d` a **PASSÉ** l'ensemble des
vérifications. Les anciennes preuves physiques
demeurent valides et distinctes.

Voir [contrat complet de sorties, quadratures et tests](V130_EXAMPLES_OUTPUT_CONTRACT.md).

### État des sorties utilisateur — 27 novembre, après la campagne native

La campagne originale de source `a471a637` a exécuté
**huit exemples canoniques en MPI1 et Heat–Wave en MPI2** :
`PASS_V130_ALL_NINE_REAL_USER_EXAMPLE_OUTPUT_CASES`, ZIP
CRC correct, réels fichiers PNG, SVG, GIF, NPZ, PVD
et VTU reproduits directement par les exemples.
Ce résultat prouve leurs sorties et leurs parcours
dans ces neuf configurations.

**La revue de leur arborescence a révélé un défaut
distinct** : des répertoires `outputs/log`,
`outputs/png`, `outputs/mesh`, `outputs/pv`
créés à l'import de SCRIMP et une écriture
transitoire supplémentaire des maillages Gmsh.
Un correctif source a depuis supprimé la création
automatique des répertoires génériques et fait écrire
et importer les maillages directement sous
`working_dir/outputs/<system_name>/mesh/`.
Les huit exemples canoniques sont raccordés
à ce nouveau chemin ; les chemins historiques
sans paramétrage explicite restent compatibles.
Le smoke-test est maintenant **fail-closed** sur
les anciens dossiers parasites et les doublons
de maillages. La campagne complète `20d79d0d`
a validé ces nouvelles barrières sur **9/9 cas**.

## 1. Décision générale

**Le socle numérique est très avancé ; la release n'est PAS encore prête
à être taguée.** Le prochain investissement de temps doit porter
prioritairement sur une distribution **utile à un utilisateur réel** :
exemples exécutables et compréhensibles, sorties reproductibles,
documentation, installation propre, et vérification finale de package.
Ne pas transformer une étude exhaustive de convergence d'ordre 2,
d'autres solveurs LU ou une nouvelle campagne de benchmarking en
blocage de sortie sans question physique identifiée.

### Acquis de référence — ne pas refaire sans changement pertinent

| Périmètre | État et réserve |
| --- | --- |
| GetFEM 5.5/MPI, PETSc/SLEPc et MUMPS | Image `scrimp-mpi:v1.3.0-rc3` validée sur Ubuntu ; COO distribué, réutilisation MUMPS `on_change`, physique acoustique MPI1/2/4/6 dans les campagnes bornées. |
| DAE | Initialisation algébrique **linéaire** qualifiée ; états prescrits et contraintes préservés, classes non supportées rejetées. Ne pas prétendre fournir un initialiseur non linéaire général. |
| Jacobien de masse non linéaire (#27) | Test réel sur `shallow_water.py`, MPI1/2 et comparaison `legacy` : **PASS** dans le périmètre documenté. |
| Checkpoint générique froid (#29) | API et preuve native sur Saint-Venant et une DAE linéaire, redémarrages MPI croisés : **PASS**. Ce n'est pas la sauvegarde de l'historique complet d'un intégrateur multistep. |
| Puissances discrètes (#30) | CN par pas corrigé ; validation wave CN/BDF2 : **PASS**, sans attribuer une conservation exacte à BDF2. |
| Interfaces conformes et sorties ParaView MPI | Chemins de production par défaut qualifiés sur les configurations Heat–Wave 2D testées ; n'implique pas l'ensemble des géométries, éléments et écritures distribuées. |
| Exemple 1D de corde | Physique/ports et figures natives, MUMPS/MPI et plusieurs LU : **PASS** dans le tutoriel validé. |
| Trois ondes 2D (2A, 2B, 2C) | Modèles conservatif, dissipatif, coénergie, physiques et démonstrations longues documentés comme **PASS** dans [leur rapport](V130_WAVE_2D_THREE_EXAMPLES.md). |
| `mass_spring_damper.py` | Horizon physique complet avec BDF2, MPI/parité et figures issues de l'exemple : **PASS**. |
| `nonlinear_pendulum.py` | Vrai sinus, CN complet 12 s, cas libre et forcé, MPI1/2 et animations/fichiers issus de l'exemple : **PASS**. |
| `heat_wave.py` (Gaussian) | 15 s, 3 000 CN, 3 001 états complets, MPI1/2, bilans chaleur/onde/interfaces, **PASS** en parité normique de coefficients. Ancien pointwise `allclose` historisé comme FAIL. |
| `heat_wave.py` : quatre pas de temps | Tous les calculs et bilans physiques **PASS** à maillage inchangé ; ratios auto-raffinement 3.0165/3.3382, ordres apparents 1.5929/1.7391. **Ce n'est PAS une preuve d'ordre 2** ; aucun besoin pédagogique de la réclamer maintenant. |
| Exemples historiques publiés | Les huit sources sous `examples/published/` sont **immutables** et protégées par test ; les sorties locales non suivies doivent rester intactes. |

## 2. La priorité scientifique immédiate : valeur des exemples pour l'utilisateur

Le critère pratique n'est plus « combien de campagnes scientifiques peut-on
lancer ? » mais « un nouvel utilisateur peut-il lancer le fichier d'exemple,
voir les phénomènes qu'il annonce et comprendre les ports/énergies ? ».

### Heat–Wave, prochain enrichissement recommandé

**Physique : suffisamment qualifiée pour le modèle Gaussian testé.**
L'exemple conserve la vraie géométrie concentrique, les flux RT0,
les traces orientées, l'initialisation DAE compatible et les
conditions extérieures homogènes ; il exporte déjà `T`, `p`
en PVD/VTU, le Hamiltonien total et une animation annulaire
optionnelle. Son manque pédagogique majeur est que l'échange
entre *deux sous-systèmes* n'apparaît pas côte à côte dans les
figures courantes.

Faire produire **par l'exemple canonique** (pas par un script
d'acceptation inaccessible à l'utilisateur), avec chemins de
sortie clairement indiqués :

1. Les courbes `H_heat(t)`, `H_wave(t)`, `H_total(t)`
   sur les vrais états sauvegardés, pour rendre visible le
   transfert et l'amortissement. Une énergie positive
   décroissante seule ne démontre pas la transmission.
2. Les puissances **des deux faces orientées**
   `P_heat^I(t)` et `P_wave^I(t)`, leur compensation
   et la dissipation `D_heat(t)>=0`, calculées avec
   GetFEM (sans remettre artificiellement le flux de somme à zéro).
3. Les champs de température et vitesse **aux vrais instants
   acceptés**, PVD/VTU et PNG/SVG reproductibles, plus
   animation GIF optionnelle. Ne pas interpoler/fabriquer un
   « état FEM ». Ajouter un résumé des paramètres
   `h, dt, t_f`, étiquettes des régions et signes des ports.
4. Un bilan CN **pas par pas** uniquement lorsque
   `dt_save=dt` ou qu'une véritable mesure des puissances
   de chaque pas est conservée. La configuration standard
   (`time_step=.001`, `save_step=.05`) est **sous-échantillonnée**
   et ne permet pas de reconstruire exactement le bilan à
   chaque étape à partir des seules sorties sauvées.
5. Expliquer les figures et l'API employée dans
   `docs/source/examples/heat_wave.rst`. L'utilisateur doit
   pouvoir reproduire ces fichiers via la ligne de commande
   documentée **sans importer un script de `sandbox/`**.

Les mesures physiques natives existent déjà dans
`sandbox/measure_v130_hw_dual_facet_physics.py`. Ne pas
copier sans revue toute cette instrumentation dans la
classe `DPHS` ni dans l'exemple. Préférer une petite
interface de présentation testable ou l'extraction d'un
helper générique seulement si une réutilisation réelle
est démontrée. **À qualifier nativement après
l'ajout effectif des sorties**, en vérifiant les fichiers,
leur provenance et leur lisibilité (PNG/SVG/VTU/GIF).

### Inventaire des exemples avant release

L'inventaire officiel comporte huit fichiers top-level :
`heat.py`, `heat_wave.py`, `mass_spring_damper.py`,
`nonlinear_pendulum.py`, `shallow_water.py`,
`wave.py`, `wave_coenergy.py`, `wave_dissipative.py`.
Les validations 2A/2B/2C, les deux exemples concentrés
et Heat–Wave couvrent désormais une grande part de la
liste. **Les trous de vérification les plus évidents**
restent la simulation utilisateur `heat.py` et la
démonstration complète et présentable de
`shallow_water.py` : le PASS du Jacobien (#27) n'est
pas automatiquement un PASS de toutes ses sorties
et de ses choix initiaux. Vérifier les scripts,
résultats, figures et notices avant de conclure ;
ne pas dupliquer des essais déjà archivés.

## 3. Issues GitHub : situation réelle au 27 novembre

**Quatre issues ouvertes :**

| Issue | Situation | Décision proposée v1.3.0 |
| --- | --- | --- |
| [#6 — MPI powers, GetFEM/PETSc, Saint-Venant](https://github.com/g-haine/scrimp/issues/6) | Les validations MPI et puissance en régime compatible actif, y compris MPI4, ont un sous-périmètre natif accepté ; le ticket initial garde un **champ plus large** que ces preuves. | **Vérifier les exigences résiduelles** : fermer uniquement avec une portée strictement reformulée et ses preuves, ou maintenir ouverte pour les cas généraux. Le développement tangentiel/LHMNC de longue durée sort de la release. |
| [#28 — initialisation DAE non linéaire](https://github.com/g-haine/scrimp/issues/28) | La structure haut indice/tractions et les incompatibilités de données rendent un Newton algébrique générique insuffisant. | **Reporter hors chemin critique v1.3.0** ; conserver le défaut linéaire qualifié et des rejets explicites. Ne pas promettre de résolution générale. |
| [#25 — Docker/context/cache/dev](https://github.com/g-haine/scrimp/issues/25) | Environnement rc3 et VS Code ont des PASS ciblés ; les ambitions complètes de cache/reconstruction depuis machine propre ne sont pas toutes démontrées. | **Séparer le minimum de reproduction/install nécessaire** du perfectionnement BuildKit/cache, reportable. Garder ouverte tant que ses critères complets ne sont pas couverts. |
| [#26 — documentation Sphinx/exemples](https://github.com/g-haine/scrimp/issues/26) | Tutoriels et documentation scientifique beaucoup améliorés, mais encore des références obsolètes et la validation HTML/LaTeX/PDF finale à faire. | **Priorité de release**, à fermer seulement après parcours utilisateur et builds finaux. |

**Issues récentes fermées, ne pas rouvrir sans régression :**
[#23](https://github.com/g-haine/scrimp/issues/23) — vieille anomalie
MPI GetFEM locale, solution de compatibilité exactement délimitée ;
[#27](https://github.com/g-haine/scrimp/issues/27) — masse/Jacobien
non linéaire ;
[#29](https://github.com/g-haine/scrimp/issues/29) — checkpoint froid ;
[#30](https://github.com/g-haine/scrimp/issues/30) — quadrature
des puissances. Les rapports de qualification sont dans
`docs/` et doivent garder leur portée.

## 4. Roadmap opérationnelle, ordonnée par valeur et risque

| Ordre | Chantier | Condition de sortie | Bloquant pour tag ? |
| --- | --- | --- | --- |
| **P0-A** | Exemplaires utilisateur | Les exemples principaux s'exécutent sans adaptation cachée ; chaque illustration proposée est créée par l'exemple et validée ; les modèles non qualifiés portent un avertissement explicite. Priorité `heat_wave.py`, `heat.py`, `shallow_water.py`. | **Oui**, pour une documentation qui revendique ces exemples comme disponibles. |
| **P0-B** | Installation propre Ubuntu | Image finale cohérente avec le code, projet externe VS Code et chemins de sortie persistants, `wheel`/`sdist` installables, vrai test d'import et un petit exemple physique. Sans toucher au Conda actuel. | **Oui**. |
| **P0-C** | Sphinx et README (#26) | Texte/figures conformes à l'API réelle, chemins et installation documentés, build HTML **et** LaTeX/PDF reproductibles, versions et liens justes. | **Oui**. |
| **P0-D** | Décision de portée MPI/#6 | Énoncer exactement les géométries/MPI/ports/solveurs qualifiés, maintenir ou scinder #6 sans prétendre avoir clos le Saint-Venant général. | **Oui**, pour la clarté du contrat, pas pour les travaux LHMNC. |
| **P0-E** | Tests de release et packaging | Unitaires, source, modèles physiques choisis, CI et smoke empaqueté ; `setup.py` et `docs/source/conf.py` portent encore `1.2.0` aujourd'hui : ne mettre `1.3.0` qu'au moment de figer le candidat. Notes, journal des limites, merge/tag/publish après verdict. | **Oui**. |
| **P1** | GUI, nettoyage et cohérence du dépôt | Smoke GUI→script si inclus comme fonction prête à l'emploi ; sinon limitation explicite. Reclasser les scripts `sandbox/` utilisés par les tests sans effacer les données locales. | **À arbitrer** : dégrader explicitement en expérimental si insuffisant. |
| **P2** | Études non nécessaires au contenu livré | Nouveau Maxwell 3D, initialiseur DAE non linéaire général (#28), longues études spatiales/temporelles, benchmarks supplémentaires, comparaison LU, optimisation complète Docker (#25), refonte interne de `DPHS`. | **Non**, sauf si une régression de l'existant l'exige. |

**Important :** l'objectif de simplifier `DPHS`
(et, notamment, déplacer dans `utils/linalg/` des fonctions
algébriques maintenant mêlées à la logique système)
est pertinent, mais c'est une refactorisation
risquée à engager **après gel de la v1.3.0**,
avec équivalence de comportement et tests
de non-régression. Elle ne doit pas repousser
indéfiniment une release scientifiquement qualifiée.

## 5. Contrôles concrets avant « release candidate »

- Le document d'exemples `V130_EXAMPLES_SOLVER_POLICY.md`
  contient encore des **états du 23 novembre** (« native
  pending ») devenus historiquement obsolètes. Mettre
  un **état courant en tête** et marquer les notes
  anciennes comme historiques, plutôt que supprimer
  les preuves et échecs datés.
- `setup.py` contient encore `version="1.2.0"`,
  `docs/source/conf.py` contient encore
  `release = version = '1.2.0'` : **pas prêt à taguer**.
- La documentation d'installation/README qualifie
  encore différentes parties de la chaîne de
  « preview/experimental ». Trier ce qui est
  désormais prouvé de ce qui reste expérimental.
- La roadmap de 2026-W48 et le plan de release
  comportent un long historique ; **ne pas lire
  un ancien « PENDING » comme état actuel**.
  Ce document résume les priorités mais ne remplace
  pas les rapports scientifiques d'origine.
- Vérifier `python -m build`, installation d'un
  wheel dans un environnement isolé, Sphinx HTML
  et PDF, imports, petit modèle physique,
  protection de `examples/published/`, et
  reproductibilité des visuels essentiels.
  Aucun de ces tests finaux n'est déclaré PASS
  tant qu'un journal d'exécution actuel ne l'a pas prouvé.

## 6. Décision / prochain travail

**Ne pas faire de nouvelles simulations juste pour
« prouver ordre 2 » sur Heat–Wave.** En priorité :

1. Construire la **vue pédagogique native**
   du transfert chaleur–onde, sans refaire la physique.
2. Finaliser `heat.py` et la démonstration utilisateur
   `shallow_water.py` dans leur périmètre réellement
   supporté, avec des sorties concrètes.
3. Consolider le tableau des exemples et les conditions
   d'installation/CI, puis livrer la Sphinx et le packaging.
4. Statuer explicitement sur #6/#25/#26/#28
   avant publication, sans demander leurs extensions
   non essentielles comme préalable systématique.

**Aucune issue n'est fermée et aucun tag n'est créé
par cette synthèse.** La release se décide sur
des preuves natives des fonctionnalités annoncées
et une expérience utilisateur vérifiée, pas sur
la quantité de diagnostics de recherche terminés.
