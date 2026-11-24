# SCRIMP v1.3.0 — Exemple 02 : onde 2D conservative, audit initial

**Branche unique : `release/1.3.0-preparation`.**  
**État : AUDIT DE SOURCE (24 novembre 2026), AUCUN NOUVEAU PASS NATIF.**  
**Sources examinées :** `docs/source/examples/wave.rst` (première partie,
avant « Adding Damping to the dphs »), `examples/wave.py`,
`examples/wave_dissipative.py`, `docs/source/biblio.rst` et
`docs/V130_PORT_POWER_QUADRATURE_ISSUE30.md`.

## Rectification explicite de l'audit

**Après réouverture des trois codes et de la syntaxe GetFEM,**
la conclusion antérieure sur `q.T.q` était incorrecte :
la chaîne GWFL contient bien le **nom du paramètre `T`**.
Le terme d'énergie potentielle est donc déjà la contraction
physique anisotrope. La mise en production ne doit
pas introduire une prétendue correction de ce terme.
La distinction avec `wave_coenergy.py` est essentielle :
ce dernier utilise le **paramètre de compliance**
`Tinv` en lieu et place de la tension `T`,
dans `q.Tinv.q` après substitution des coétats.

## Ordre pédagogique de Sphinx

Le premier chapitre « Getting started » avec la corde 1D est **terminé** :
4/4 tests de parité Sphinx/Python, 17/17 tests scientifiques du
tutoriel, cinq vraies configurations MUMPS/SuperLU et simulation
illustrative longue qualifiés, HTML/LaTeX compilés.

La **deuxième** étape est la membrane/onde **2D conservative**,
section « The wave equation », `examples/wave.py`.
La même page contient ensuite une **troisième étape physique**
« Adding Damping to the dphs » et le fichier
`examples/wave_dissipative.py`. Il faut travailler l'onde
conservative d'abord, puis le cas dissipatif ensuite, tout
en conservant le déroulé scientifique de la page.

## Ce qu'illustre réellement le modèle

Rectangle `Ω=(0,2)×(0,1)`, avec quatre régions de bord
(10 bas, 11 droite, 12 haut, 13 gauche).
`q=∇w` est un **champ vectoriel de déformation** ;
`p=ρ∂t w` une **quantité de mouvement scalaire**.
La tension anisotrope spatialement variable `T(x,y)`
est un tenseur symétrique, de chaîne GWFL
`[[5+x,x*y],[x*y,2+y]]`; la densité `ρ(x)=3−x`
reste positive sur le rectangle.
`q` : DG1 ; `p` : CG2 ; traces : DG1.

Les ports sur le bas, la droite et le haut imposent
la **traction normale scalaire** de Neumann
`u_N=n·e_q=n·(Tq)`, et observent la vitesse
`y_N=e_p`. Sur le bord gauche, le contrôle
de type « Dirichlet » est la **vitesse** `u_D=e_p`,
pas la déflexion `w` ; son observation colocalisée
est la **traction normale** `y_D=n·e_q`.
Le bilan continu de référence est

```text
H(q,p) = (1/2) ∫Ω [ q·Tq + p²/ρ ] dx
dH/dt = ∫ΓN u_N y_N ds + ∫ΓD u_D y_D ds.
```

Les signes effectifs des ports GetFEM, le traitement
de la normale et la relation algébrique DAE devront
être vérifiés par le calcul natif : **aucun PASS
n'est inféré de la seule écriture théorique**.

## Incohérences et questions scientifiques à traiter AVANT les solveurs

1. **ÉNERGIE TENSORIELLE : AUDIT RECTIFIÉ.**
   `S.Term("Potential energy", "0.5*q.T.q", [1])`
   est **correct** en syntaxe GetFEM GWFL :
   `q.T.q` est la contraction `q·T·q` avec le
   **paramètre tensoriel `T`**, et NON le produit
   `q^Tq` au sens d'une expression Python.
   L'audit initial identifiait à tort une omission
   de `T` dans ce terme. **Ne pas le remplacer
   par `q.T.T.q` ni par `q.q`.** La validation
   native du véritable Hamiltonien (positivité,
   anisotropie et bilan de puissance) reste à faire.
   Comparer ultérieurement le terme GetFEM aux
   intégrales alternatives `q.(T*q)` si le
   langage GWFL utilisé le permet, sans modifier
   la physique pour rendre un test artificiellement vert.
2. **ÉCRITURE DE BORD À PRÉCISER.**
   La documentation écrit notamment `e_q=u_N`
   et `y_D=e_q` sur les frontières ; ce sont des
   objets de nature vectorielle versus scalaire.
   Les équations correctes doivent faire apparaître
   `n·e_q`. La condition de Neumann dans la PDE
   exige `n·T∇w=u_N` ; une condition de vitesse
   `∂t w=u_D` n'est pas exactement `w=u_D`.
   Attention à la **normale extérieure** et
   aux signes induits par l'intégration par parties.
3. **INCOMPATIBILITÉ INITIALE POTENTIELLE.**
   Source actuelle : `q0=[0,0]`,
   `p0=3**(-20*((x-.5)^2+(y-.5)^2))`,
   `u_D(t,y)=0.1 sin(4t) sin(4πy) exp[-10(2.5−t)^2]`.
   Le `p0` gaussien n'est pas identiquement nul
   sur le bord contrôlé `x=0`, alors que
   `u_D(0,y)=0` ; le code force
   `init_step_ts_type="pseudo"` (non présent
   dans le bloc Sphinx). Ne pas présenter ce
   choix comme une initialisation DAE **strictement
   compatible**. Concevoir une initialisation
   simple et compatible avant le test physique.
4. **IMPULSION HORS HORIZON DU SMOKE.**
   L'ancienne sollicitation gaussienne est
   centrée à `t=2.5` et `t_f=5`. Une
   démonstration courte `t≤0.04` ne contient
   pratiquement aucun travail de bord :
   un PASS trivial sans puissance non nulle
   serait un mauvais test scientifique. Prévoir
   un court scénario de stimulation compatible
   et **une démonstration longue séparée** pour
   la propagation et les réflexions.
5. **ÉCART DOC/SOURCE ÉTABLI.**
   `docs/source/examples/wave.rst` appelle
   `set_time_scheme(ts_type="cn",t_f=t_f,dt_save=.01)`,
   tandis que `examples/wave.py` impose en plus
   `init_step_ts_type="pseudo"`. Le document
   omet cette modification de la projection DAE.
   Le bloc d'introduction source doit être copiable,
   piloté par une comparaison AST comme pour l'onde 1D.
6. **CONVENTIONS DE PUISSANCE DISTINCTES.**
   La correction de l'issue #30 distingue
   puissance `cn_step` et `instantaneous`.
   Ne pas interpréter `Port.get_power()`
   sans identifier les signes des ports, poids
   géométriques et schéma temporel. Le vrai test
   exige le bilan `ΔH−W_N−W_D≈0` pour les
   commandes non triviales, avec vérification
   indépendamment de `plot_Hamiltonian()`.
7. **CONTEXTE PHYSIQUE.** Parler ici d'une
   **membrane sous tension anisotrope**, pas
   simplement de module d'Young élastique scalaire.
   Les sorties `q`, `p`, `e_p=p/ρ` et
   `e_q=Tq` devront être clairement nommées
   dans figures et animations.

## Choix de livraison de l'exemple 02

**Étape A — validation de modèle conservative.**
Garder le rectangle, les quatre bords et les deux types de
conditions : leur contraste est la valeur pédagogique
principale par rapport à la corde 1D. Choisir
`q0,p0,u_D` de manière strictement compatibles
sans supprimer la physique anisotrope, et mesurer
`H` par la vraie énergie tensorielle. Prouver
que les données initiales sont préservées,
que les puissances ont les signes physiques et
que chaque pas CN satisfait la balance voulue.

**Étape B — matrice d'installation.** Sur la
**même configuration vérifiée**, MUMPS MPI1/2,
SuperLU MPI1, SuperLU_DIST MPI1/2 ; trajectoires
originales complètes, parité rangs/solveurs,
convergence PETSc, bilans physiques, export VTK.
Ne pas « corriger » les signes en adaptant
l'assertion à une formule numériquement fausse.

**Étape C — documentation et visualisations.**
Un programme canonique `examples/wave.py`,
blocs Sphinx générables/copiables, test de parité AST.
Générer à partir des **vrais** VTU/PVD :
instantanés 2D de `p` et `e_p`,
champ vectoriel `q` et contrainte `e_q`,
cartes de densité d'énergie potential/kinetic,
puissances intégrées sur les quatre frontières,
`H-H0-W` et comparaison des solveurs ;
GIF/SVG animé de la membrane uniquement si
le déplacement `w` est calculé ou reconstruit
avec une **jauge et un algorithme exposés**.
Conserver CSV/NPZ/provenance ; ne pas versionner
des centaines de sorties binaires.

**Étape D — dissipative.** Seulement après
qualification du cas conservative, passer à
`examples/wave_dissipative.py`, toujours dans
l'ordre de la page Sphinx. Pour ce dernier,
énergie `q·Tq` et perte par amortissement
port-Hamiltonienne non négative devront être
vérifiées **avec la bonne puissance de bord**.

## Publications de l'établissement

Déjà présentes dans `docs/source/biblio.rst` :
- **Haine, Matignon, Serhani (2023)**,
  *Numerical Analysis of a Structure-Preserving
  Space-Discretization for an Anisotropic and
  Heterogeneous Boundary Controlled N-Dimensional
  Wave Equation As a Port-Hamiltonian System*,
  DOI `10.4208/ijnam2023-1005`.
- **Brugnoli, Haine, Matignon (2022)**,
  *Explicit structure-preserving discretization of
  port-Hamiltonian systems with mixed boundary control*
  (MTNS 2022).
- **Serhani, Matignon, Haine (2019)**,
  *Partitioned Finite Element Method for
  port-Hamiltonian systems with Boundary Damping:
  Anisotropic Heterogeneous 2D wave equations*
  (pertinent surtout pour la seconde partie dissipative).

**Aucun de ces articles n'est encore revendiqué
comme exactement reproduit par l'exemple actuel.**

## Premier essai natif — 24 novembre 2026 — physique PASS, export de coénergie MPI à corriger

**Source du premier essai :** `abea0401`.  
**Archive exacte :**
`v130-wave2d-mixed-conservative-abea0401-Jm7S0brY.zip`
(fournie et relue ; 86 entrées).
Le statut global du lanceur est bien
`FAIL_OR_UNQUALIFIED_WAVE2D_MIXED_TENSOR_CN_MPI` :
**ne pas le rebaptiser PASS**.

### Mesures originales

| Contrôle réel GetFEM/PETSc | MUMPS MPI1 | MUMPS MPI2 |
|---|---:|---:|
| DOFs | 2230 | 2230 |
| Pas CN acceptés | 10/10 | 10/10 |
| Rejets TS / échecs SNES | 0 / 0 | 0 / 0 |
| Initialisation algébrique | PASS, résidu 0 | PASS, résidu 0 |
| Défaut énergétique maximal par pas | 2.9815559743351372e-18 | 1.2468324983583301e-18 |
| Défaut énergétique cumulé maximal | 5.421010862427522e-18 | 2.6020852139652106e-18 |
| Puissance absolue gauche maximale | 0.04899638915980502 | 0.048996389159805157 |
| Sept contrôles physiques | tous PASS | tous PASS |
| Export q et p | 11 VTU chacun | 11 VTU chacun |
| Export e_q et e_p | 11 VTU chacun | ÉCHEC avant e_q |
| Statut du processus | 0 | 1 |

L'archive contient les traces FE **complètes** des
deux rangs MPI 2. Une comparaison indépendante des
`rank0.npz` et `rank1.npz` montre une **égalité exacte**
de tous les champs enregistrés, y compris les
états de 2230 DOFs ; l'écart maximal de
**MPI1 vs MPI2** sur tous les états est
`4.08006961549745e-14`. Les puissances et
Hamiltoniens concordent également à l'arrondi.

### Origine de l'échec

Au moment de l'export MPI de `e_q`, le code
`DPHS._validate_root_only_pv_replication()`
utilisait `self.ports.get(name_variable)`.
Or le port d'état est indexé `q`, même si
son champ d'effort `e_q` est une **vraie
variable GetFEM** sur le même maillage.
Ce n'était **ni un champ manquant**, ni un
problème de réplication MPI. Le code séquentiel
exporte les 11 VTU de `e_q` et `e_p`
sans difficulté.

La correction **générique** résout le
port propriétaire par
`name_variable in (port.get_flow(), port.get_effort())`,
avec cardinalité obligatoirement égale à 1.
Les noms inconnus et les appartenances ambiguës
restent des erreurs collectives : **aucune
désactivation du contrôle MPI** ni retour
à la concurrence d'écriture n'est autorisé.
Les tests source du véritable validateur, dans
`tests/test_v130_interface_pv_production_policy.py`,
couvrent désormais `q, p, e_q, e_p, U_L, Y_L`,
ainsi que les contrôles négatifs de sécurité.

**Statut après correctif : SOURCE MODIFIÉE, NOUVEL
PASS NATIF MPI2 EN ATTENTE.** Relancer le même
lanceur avec l'image Docker figée, en conservant
les états de référence initialement fournis.
Si tout passe, seulement alors promouvoir
`examples/wave.py` et aligner Sphinx.

## Promotion de la source et préparation de l'acceptation LU (24 novembre)

**Le deuxième essai natif est PASS de bout en bout :**
`v130-wave2d-mixed-conservative-cf42089d-pI8VXuBa.zip`,
`PASS_WAVE2D_CONSERVATIVE_TENSOR_MIXED_MUMPS_MPI1_MPI2`.
Deux simulations physiques MUMPS (MPI 1/2), dix pas
CN acceptés, 2230 DOFs, données initiales DAE préservées,
puissance de bord non triviale, défaut énergétique
par pas **2.982e−18** (MPI1) et **1.789e−18**
(MPI2), avec **88 véritables VTU** (11 × 4 × 2)
et les huit PVD associés. Les états sauvegardés
par les deux rangs MPI2 sont identiques ; écart
maximal MPI1/MPI2 **4.16334e−14**.

À partir de cette preuve, l'ancien candidat
`sandbox/probe_v130_wave2d_conservative.py`
a été **promu dans `examples/wave.py`**, en
conservant le nom de fonction historique `wave_eq()`
et ses équations, paramètres et contrôles nativement
qualifiés. Seuls l'encapsulation en API canonique,
la documentation et les commentaires ont changé.
Le pilote temporaire et son lanceur ont été retirés :
leurs anciennes sources et leur preuve sont toujours
retrouvables dans l'historique Git/ZIP.

La partie conservative de
`docs/source/examples/wave.rst`
a été **reconstruite sans doublons des briques** :
huit `literalinclude` lient directement
les champs, ports, FEM, paramètres, briques,
contrôles, pas CN et Hamiltonien de la source
qualifiée. La partie **dissipative**, toujours
dans la même page, **n'est pas encore modifiée
ni qualifiée** ; la formulation coénergie est
également réservée à l'étape 2C.

Le nouvel outil
`scripts/run-v130-wave2d-conservative-acceptance.sh`
lancera directement `examples/wave.py` en
**MUMPS MPI 1/2, SuperLU MPI 1,
SuperLU_DIST MPI 1/2** et comparera
les trajectoires, Hamiltoniens, puissances
et quatre champs ParaView entre backends.

**Statut scientifique : MUMPS MPI 1/2 PASS
acquis. Cinq LU/MPI PENDING.**
Une nouvelle preuve est nécessaire pour
promouvoir le résultat *global* de l'exemple
2A à une qualification d'installation.

## Campagne cinq solveurs sur le code canonique : trois PASS, deux FAIL

**Date : 24 novembre 2026.**
Source `99d99f2d94c82fe95aac3b39a3aaf3a4c4c27eb8`,
archive `v130-wave2d-conservative-five-lu-99d99f2d-tEOjqgHy.zip`.
La campagne utilise bien `python -m examples.wave`,
sans exécuter le pilote abandonné.

**Trois configurations qualifiées :**
MUMPS MPI 1/2, SuperLU MPI 1. Chacune
effectue dix pas CN acceptés, 2230 DOFs,
projection DAE algébrique à résidu nul,
puissance de bord gauche non nulle et
44 VTU (quatre séries `q,p,e_q,e_p` de
11 instants). Les défauts énergétiques
par pas culminent respectivement à
`2.9815559743351372e-18`,
`1.6263032587282567e-18` et
`3.0899761915836876e-18`.

Comparaison indépendante depuis les
`rank0.npz` de l'archive : par rapport
à MUMPS MPI1, écarts maximaux FE sur
les 2230 DOFs
`9.248157795127554e-14` (MUMPS MPI2)
et `1.61259894326804e-13` (SuperLU MPI1).
Énergies finales
`0.0015763111255133718`,
`0.0015763111255133748` et
`0.0015763111255133716`.

**Deux configurations NON QUALIFIÉES :**
SuperLU_DIST MPI1 et MPI2.
Les deux réussissent l'assemblage GetFEM
et la projection DAE initiale à résidu
nul ; aucun pas CN n'est accepté. Le premier
`TSStep()` échoue avec
`DIVERGED_STEP_REJECTED` (PETSc code 91).
**Ce message ne détermine pas à lui seul
s'il s'agit de factorisation, KSP, SNES,
pivotage ou résidu non linéaire.**
Aucun reparamétrage ne sera appliqué
à l'aveugle.

**Correctif de diagnostic source uniquement :**
`scripts/run-v130-wave2d-superlu-dist-diagnostic.sh`
active le chemin **déjà existant**
`SCRIMP_TS_FAILURE_DIAGNOSTICS=1`,
sur SuperLU_DIST MPI1 pour un seul pas
`dt=0.008`. Il recueille les raisons
TS/SNES/KSP/PC et archive les erreurs,
sans modifier la formulation ou les
tolérances. L'extension MPI2
(`SCRIMP_V130_WAVE2D_DIAG_RANKS=2`)
est disponible, mais le diagnostic
interne détaillé n'est pas encore activé
en MPI : commencer par MPI1.

Le rapport d'acceptation à cinq cas a
également été amélioré pour **conserver
les deux échecs SuperLU_DIST** dans
`summary.json`, même si le premier
arrête l'analyse détaillée.

**État scientifique : 3/5 PASS,
2/5 FAIL, diagnostic natif approfondi PENDING.**
Les exemples 2B dissipatif et 2C coénergie
restent hors du périmètre de cette campagne.

## SuperLU_DIST : preuve native du pivot numérique nul

Le **diagnostic du premier pas** au commit `62fb8765`
est disponible dans
`v130-wave2d-superlu-dist-first-step-62fb8765-QIaYkugY.zip`.
Son archive native a été examinée intégralement,
y compris `native.host.log`, le manifest et
`summary.json` :

```text
MAIN TS FAILURE DIAGNOSTIC:
 TS_reason=-2, TS_step=0, TS_time=0
 TS_dt=1.9073486328125e-09
 SNES_reason=-3, SNES_iter=0
 SNES_norm=4.240793672320428e-09
 KSP_reason=-11, KSP_iter=0, KSP_norm=inf
 PC_type=lu, PC_failed=2, factor_solver=superlu_dist
```

D'après l'énumération officielle **PETSc PCFailedReason**,
`PC_failed=2` signifie
`PC_FACTOR_NUMERIC_ZEROPIVOT`. `KSP_reason=-11`
signifie `KSP_DIVERGED_PC_FAILED`.
La divergence TS n'est que la **conséquence** de
l'échec de factorisation. Ce constat ne démontre
**pas** la singularité mathématique de la
matrice DAE : MUMPS MPI1/2 et SuperLU MPI1
ont résolu exactement la même formulation,
avec bilans énergétiques vérifiés.

La bibliothèque PETSc permet plusieurs
**ordres de permutation** SuperLU_DIST
(`mat_superlu_dist_rowperm`,
`mat_superlu_dist_colperm`) ; le remplacement
des pivots minuscules est **désactivé par
défaut** dans son interface officielle.
Aucune pseudo-régularisation des pivots
n'est ici appliquée sans preuve.

### Expérience de discrimination, sans nouveau modèle

Le fichier
`sandbox/probe_v130_wave2d_superlu_dist_ordering.py`
injecte une option à la fois **après** l'initialisation
de `DPHS.set_time_scheme()`, car SCRIMP efface
la base PETSc au moment de construire une DPHS.
Il importe ensuite `examples.wave.main()`
**sans modifier** son code scientifique.

Le script
`scripts/run-v130-wave2d-superlu-dist-ordering-study.sh`
étudie en séquentiel :

1. `LargeDiag_AWPM` (permutation des lignes) ;
2. `MMD_AT_PLUS_A` (permutation des colonnes).

Chacun exécute les **dix pas CN** originaux
`dt=0.008`, `t_f=0.08`, `2230 DOFs`,
SuperLU_DIST réel. Le test continue de vérifier
la conservation de l'Hamiltonien anisotrope et
les **quatre puissances de bord**, enregistre
les quatre champs VTU, et confronte les vrais
états finis à l'archive MUMPS
`v130-wave2d-conservative-five-lu-99d99f2d-tEOjqgHy.zip`
**sans relancer MUMPS**. Il enregistre les
statistiques d'assemblage de la matrice
(`SCRIMP_MPI_TRACE=1`) et refuse une
configuration dont l'option de permutation
est signalée **non utilisée** par PETSc.

Commande native :

```bash
cd ~/Simulations/scrimp
git pull --ff-only
conda activate scrimp
SCRIMP_DOCKER_NO_SUDO=1 \
SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
bash scripts/run-v130-wave2d-superlu-dist-ordering-study.sh
```

**Statut : diagnostic ciblé prêt, exécution
native EN ATTENTE.** Même si une permutation
passe en MPI1, cela ne qualifie **pas**
SuperLU_DIST par défaut en MPI1/2 et ne
clôture pas la campagne cinq solveurs.
Le succès éventuel sera seulement un
**candidat non régularisant**, à confirmer
sur MPI2, avec tous les bilans physiques
et contrôles VTK, avant décision de
politique backend pour SCRIMP.

## Permutations non régularisantes : échec documenté et nouveau diagnostic matriciel

**Archive native analysée :**
`v130-wave2d-superlu-dist-ordering-8cb2a51f-422NJxVt.zip`,
source `8cb2a51f` ; les 13 tests source
associés ont réussi. Le résultat est
`FAIL_OR_UNQUALIFIED_SUPERLU_DIST_ORDERING_EXPERIMENT`.

Deux expériences **SuperLU_DIST MPI1** ont réellement échoué :

| Option | Échec | Indice objectif |
|---|---|---|
| `mat_superlu_dist_rowperm=LargeDiag_AWPM` | **139**, segmentation fault pendant la factorisation | Message natif `CombBLAS is not available` immédiatement avant le segfault ; le chemin AWPM est non opérationnel avec cette compilation SuperLU_DIST. |
| `mat_superlu_dist_colperm=MMD_AT_PLUS_A` | **1**, zéro pivot numérique | `KSP_reason=-11`, `PC_failed=2`, `SNES_reason=-3` et aucun pas accepté, comme sans permutation. |

Le log complet contient le **premier vrai
Jacobien CN** `J=K+250M`, de dimension
`2230×2230`, avec **25801 coefficients stockés**.
Aucune valeur infinie/NaN, **aucune ligne
numériquement nulle**, mais exactement **dix
diagonales nulles aux indices globaux 80–89**.
Les mêmes caractéristiques sont observées lors
des essais de pas raccourcis. Ces dix diagonales
nulles sont **compatibles avec une matrice de
point-selle pH-DAE** ; elles ne démontrent
en aucun cas que la matrice est singulière.

Le message `CombBLAS is not available` est
cohérent avec le code SuperLU_DIST qui tente
d'utiliser CombBLAS pour cette permutation ;
ce chemin a été vérifié dans les sources
SuperLU_DIST. **Ne plus lancer AWPM** sur
l'image Docker figée tant que cette capacité
n'a pas été rétablie et isolément qualifiée.
Il ne faut pas interpréter ce segfault
comme une défaillance de la formulation PFEM.

### Étape suivante : PREUVE sur l'opérateur exact

**Aucune nouvelle permutation ni remplacement
des petits pivots par défaut.** Les fichiers
suivants préparent un diagnostic reproductible :

- `sandbox/probe_v130_wave2d_factor_matrix.py`
  capture **après le vrai callback IJacobian
  et avant la factorisation** les CSR
  de `P` (préconditionneur réellement factorisé
  par PETSc), `K` et `M`, ainsi que
  l'instant `t=0.008`, le shift `250`
  et les intervalles des **variables GetFEM**.
  Cette capture, pour MPI1 seulement, ne
  change ni le solveur ni le code utilisateur.
- `sandbox/analyze_v130_wave2d_factor_matrix.py`
  vérifie indépendamment dans **SciPy**
  `P=K+250M`, compare numériquement les
  matrices effectivement construites avec
  **MUMPS et SuperLU_DIST**, calcule le
  **rang structurel** (sur les non-zéros
  numériques, sans compter les zéros stockés),
  identifie les variables propriétaires des
  dix indices, puis tente `splu` avec
  pivotage numérique actif,
  **sans régularisation**.
  Les résidus des factorisations réussies
  sont calculés sur un second membre
  déterministe compatible avec la matrice.
- `tests/test_v130_wave2d_factor_matrix.py`
  reproduit une matrice de point-selle
  de dimension 2230 avec **dix diagonales
  nulles mais un rang structurel plein**,
  ainsi qu'un cas réellement déficient.
  Il détecte également toute divergence
  entre matrices de référence.

Commande :

```bash
cd ~/Simulations/scrimp
git pull --ff-only
conda activate scrimp
SCRIMP_DOCKER_NO_SUDO=1 \
SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
bash scripts/run-v130-wave2d-superlu-dist-matrix-audit.sh
```

Le script effectue **un seul pas** avec
MUMPS MPI1 (contrôle) et SuperLU_DIST MPI1
(échec attendu), mais capture les matrices
**avant** l'éventuel échec de LU. Il
produit un JSON, les trois CSR dans chaque
répertoire et une archive ZIP même si la
factorisation SuperLU_DIST échoue.
**Ce n'est pas une acceptation cinq solveurs.**

**Statut : SOURCE prête ; diagnostic matriciel
natif EN ATTENTE.** Aucun changement du
modèle `examples/wave.py` ni de
`scrimp/dphs.py` n'est requis.
Les anciens trois PASS demeurent valables ;
**SuperLU_DIST MPI1/2 reste non qualifié**.

## Preuve matricielle : rang plein et LU sans régularisation (24 novembre 2026)

**Archive fournie et inspectée :**
`v130-wave2d-first-cn-jacobian-574d9c71-NKqIZ0BZ.zip`,
code `574d9c71`. Le script de capture a réellement
enregistré les matrices originales du premier
pas CN, **avant** la factorisation, pour MUMPS
MPI1 et SuperLU_DIST MPI1.

**Résultat du diagnostic (et non du solveur) :**
`PASS_EXACT_FIRST_CN_JACOBIAN_CAPTURE_ANALYSIS_NOT_SOLVER_ACCEPTANCE`.
Le solveur MUMPS accepte le pas ; le solveur
SuperLU_DIST conserve `PC_FAILED=2`.

Mesures vérifiées indépendamment depuis les
matrices originales (CSR) :

| Mesure | Valeur |
|---|---|
| Dimensions de `P=K+250M` | **2230 × 2230** |
| Non-zéros stockés / non-zéros numériques | **25801 / 25791** |
| Écart `P−(K+250M)` | **0 exactement** |
| Différence MUMPS/SuperLU_DIST sur `P`, `K` et `M` | **0 exactement** |
| Rang structurel des non-zéros numériques | **2230/2230** |
| Lignes vides | **0** |
| Diagonales nulles | **80–89**, exactement dix |
| Propriétaire des dix indices | **`Y_L`**, variable algébrique de traction gauche |

Ces dix lignes constituent les **contraintes du
port mixte**, elles ne sont pas des équations absentes.
L'inspection de leurs vrais coefficients révèle un
couplage non nul avec `U_L` et `e_p`.
Un appariement biparti maximal indépendant
sur la matrice réelle fournit **2230
correspondances sur 2230** et des colonnes
de pivot non diagonales pour les dix lignes
`Y_L`.

Une **véritable** factorisation SciPy/SuperLU
de ce même CSR (aucune modification des
coefficients) a réussi dans les deux cas :

| Réordonnancement | `min(abs(diag(U)))` | Résidu relatif |
|---|---:|---:|
| `COLAMD` | `5.6579925859605124e-05` | `2.2847109146099756e-15` |
| `MMD_AT_PLUS_A` | `0.00067651748410936` | `5.539095647629999e-15` |

Une estimation de la condition en norme 1,
à partir de la LU SciPy, donne
`cond_1(P)≈4.28e5` (estimation, pas un
calcul exact des valeurs singulières).

**Conclusion justifiée :** ce Jacobien est
structurellement de rang plein et possède
une factorisation LU numérique exploitable,
même avec ses dix diagonales nulles.
L'échec `PC_FACTOR_NUMERIC_ZEROPIVOT`
est spécifique au **chemin de factorisation**
SuperLU_DIST/PETSc employé dans l'image,
pas à une différence de données GetFEM.
On ne peut pas encore identifier la cause
exacte à l'intérieur de SuperLU_DIST
(scaling, permutations effectives, critère
de pivot ou autre chemin numérique).

### Expérience suivante, volontairement circonscrite

Le driver
`sandbox/replay_v130_wave2d_superlu_dist_factor.py`
et le lanceur
`scripts/run-v130-wave2d-superlu-dist-factor-replay.sh`
rejouent le **CSR original**, reconnu par le
SHA256 `24346856f358e4f81cf332ebc974e329dc084c60f5d6b98bfbcd788fdff37f54`.
Aucun assemblage GetFEM n'est relancé.
Les données sont distribuées uniquement
aux lignes locales PETSc, avec contrôle
de l'égalité de chaque coefficient.

Le banc compare en MPI1 **et** MPI2 la
factorisation SuperLU_DIST ordinaire et
l'option opt-in PETSc
`mat_superlu_dist_replacetinypivot=true`
(désactivée par défaut). Celle-ci peut
**modifier la factorisation numérique** ;
elle n'est donc pas présentée comme une
permutation non régularisante. Le banc
construit un vrai `b=P*x*` avec `x*`
connu, résout avec le **véritable
KSP/PC PETSc SuperLU_DIST** et impose
un résidu relatif arrière < `1e-11`
et une erreur relative de solution
< `1e-7`, sur **tous les rangs**.

```bash
cd ~/Simulations/scrimp
git pull --ff-only
conda activate scrimp
SCRIMP_DOCKER_NO_SUDO=1 \
SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
bash scripts/run-v130-wave2d-superlu-dist-factor-replay.sh
```

**Une réussite de ce replay isolé ne
qualifie pas la simulation des dix
pas CN, les puissances de bord ou
l'export MPI**. Elle permettrait
seulement d'envisager un essai complet
de la véritable onde 2D, avec une option
explicitement choisie et tracée.
La politique générale de SCRIMP
ne sera pas modifiée à ce stade.
**Statut : exécution native du replay PENDING ;
campagne originale toujours 3/5 PASS.**

## Rejeu PETSc SuperLU_DIST sur l'opérateur original : quatre calculs, aucun PASS

**Preuve native fournie :**
`v130-wave2d-superlu-dist-factor-replay-e87b20d1-H7IgKhB8.zip`,
source `e87b20d108d980ca9648cbf14d38aed96b37f4d3`.
L'archive a été ouverte et ses **six JSON de rang**
et quatre journaux natifs inspectés. Le test source
préalable était **6/6 PASS**. Le CSR original
`2230×2230` est vérifié avec
`SHA256=24346856f358e4f81cf332ebc974e329dc084c60f5d6b98bfbcd788fdff37f54`.
Tous les rangs déclarent le **même graphe CSR et
les mêmes coefficients** que l'archive initiale.

| PETSc réel (matrice inchangée) | MPI1 | MPI2 |
|---|---|---|
| SuperLU_DIST standard, `preonly` | `KSP=-11, PC_failed=2`, pas de solution | Même échec sur les 2 rangs |
| `ReplaceTinyPivot=true`, `preonly` | `KSP=4, PC_failed=0`, mais seuils de précision manqués | Même résultat sur les 2 rangs |

Pour le **remplacement des pivots** :

| Erreur réellement mesurée | MPI1 | MPI2 | Acceptation exigée |
|---|---:|---:|---:|
| Résidu arrière relatif sur `Px-b` | `5.083978939457479e-10` | `5.083978914609362e-10` | < `1e-11` |
| Erreur relative sur le vecteur connu `x*` | `7.280665330727204e-7` | `7.231143229884494e-7` | < `1e-7` |

L'option `ReplaceTinyPivot` permet donc **de
factoriser et de résoudre sans déclencher
`PC_FACTOR_NUMERIC_ZEROPIVOT`**. Elle ne
fournit toutefois **pas une solution suffisamment
précise** dans les tests numériques actuels.
PETSc annonce `KSP_reason=4` (PREONLY terminé),
ce qui ne constitue **pas** une validation des
résidus véritables ou de l'erreur de solution.
La qualification du replay est bien
`FAIL_OR_INCOMPLETE_OPTIONAL_PIVOT_PETSC_LU_REPLAY` ;
**aucun PASS n'a été réinterprété**.

### Expérience suivante : deux raffinement itératifs, seuils inchangés

La documentation PETSc expose
`mat_superlu_dist_iterrefine` (booléen,
raffinement interne en double précision), et
`KSPGMRES` accepte le préconditionnement
**à droite** avec contrôle du **résidu
non préconditionné**. Le nouveau driver
`sandbox/replay_v130_wave2d_superlu_dist_refinement.py`
réutilise strictement le CSR et le
second membre `b=P*x*` précédents,
avec `ReplaceTinyPivot=true` dans
les deux cas :

- `internal_refine` : SuperLU_DIST
  `IterRefine=true`, KSP `PREONLY` ;
- `right_gmres` : `KSPGMRES`,
  `PC_RIGHT`,
  `NORM_UNPRECONDITIONED`,
  `rtol=1e-12`, `atol=1e-14`,
  `max_it=30`.

Dans **chaque** cas, le résidu arrière
relatif original reste < `1e-11`
et l'erreur sur la solution exacte
< `1e-7`, sans possibilité de
PASS sur le seul `KSP_reason`.
Le test vérifie la configuration
**réellement active** de la méthode,
la présence des coefficients CSR
originaux (zéros structurels compris)
et la non-utilisation d'un solveur
de remplacement. Les versions MPI1
et MPI2 seront toutes deux exécutées.
Un cas n'est acceptable que si
les deux rangs réussissent strictement.

```bash
cd ~/Simulations/scrimp
git pull --ff-only
conda activate scrimp

SCRIMP_DOCKER_NO_SUDO=1 \
SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
bash scripts/run-v130-wave2d-superlu-dist-refinement.sh
```

**Statut : ce nouveau raffinement natif
n'a pas encore été exécuté.** Il ne
modifie ni l'onde conservative canonique,
ni le Jacobien, ni `scrimp/dphs.py`.
Une réussite serait seulement un
**candidat LU isolé** à qualifier ensuite
sur les dix pas CN, les quatre puissances
de bord, l'Hamiltonien et les VTU MPI.
**Acceptation actuelle de l'onde 2D : 3/5.**

## SuperLU_DIST : raffinement interne enfin PASS sur le Jacobien exact (24 novembre)

**Preuve native originale :**
`v130-wave2d-superlu-dist-refinement-ad2d58a6-unyfltIn.zip`,
source `ad2d58a60f2bc79f35c9bcff47968eaa5bd10461`.
Les tests source préalables ont réussi **6/6** ;
les quatre exécutions originales ont été réalisées
en PETSc/MPI dans `scrimp-mpi:v1.3.0-rc3`.
Tous les essais utilisent exactement le CSR
(P=K+250M) issu du premier pas CN, avec
le SHA256 précédemment qualifié et
`ReplaceTinyPivot=true`. Il s'agit
d'un **replay de facteur**, pas de la
simulation physique complète.

| Variante | MPI | KSP / PC | Résidu arrière relatif sur le vrai `Px-b` | Erreur relative `x-x*` | Résultat |
|---|---:|---|---:|---:|---|
| `IterRefine` interne SuperLU_DIST | 1 | `4 / 0` | `9.541996498414314e-17` | `2.007285126672514e-14` | **PASS** |
| `IterRefine` interne SuperLU_DIST | 2 | `4 / 0` sur les deux rangs | `1.272266199788575e-16` | `2.2559753193576045e-14` | **PASS** |
| GMRES `PC_RIGHT`, résidu non préconditionné | 1 | `2 / 0` | `1.9711896124742606e-10` | `1.8807418577168255e-08` | **FAIL résidu** |
| GMRES `PC_RIGHT`, résidu non préconditionné | 2 | `2 / 0` sur les deux rangs | `2.3300343615571625e-10` | `1.1213346996440907e-08` | **FAIL résidu** |

Les seuils d'acceptation restent strictement
`relative_backward_residual < 1e-11` et
`relative_known_solution_error < 1e-7`.
Le `KSP_reason>0` ne suffit donc pas
à lui seul : la branche GMRES est à
juste titre marquée **FAIL**.

**Conclusion expérimentale :** le mécanisme
`ReplaceTinyPivot=true` **plus** le raffinement
**interne** `mat_superlu_dist_iterrefine=true`
permet une solution d'une précision comparable
à celle d'une factorisation LU classique
sur la **matrice originale**, aussi bien
en MPI1 qu'en MPI2. C'est un résultat
solide sur le **rejeu isolé** et non
une validation du modèle 2D complet.

### Porte suivante : véritable simulation port-Hamiltonienne

Le candidat
`sandbox/probe_v130_wave2d_refined_superlu_dist.py`
ne modifie **aucune** brique, FEM,
paramètre, condition initiale ni commande
de `examples/wave.py`. Il intercepte
`DPHS.set_time_scheme()` uniquement pour
enregistrer ces **deux options PETSc**
explicitement et après la création
de la DPHS, avant `TS/SNES/KSP.setFromOptions`.
Il impose le solveur réel `superlu_dist`
et la demande de qualification par
`SCRIMP_V130_WAVE2D_ENABLE_REFINED_DIST=1`.
La modification des defaults de SCRIMP
est **hors périmètre**.

Le lanceur
`scripts/run-v130-wave2d-refined-superlu-dist-physical.sh`
exécute la **véritable onde anisotrope 2D**
en MPI 1 puis MPI 2, avec **dix pas CN
fixes**, initialisation DAE algébrique,
2230 inconnues, quatre puissances
physiques, énergie originale
`0.5*q.T.q+0.5*p*p/rho`,
rejets et échecs SNES nuls.
Il conserve les **11 états complets** et
les 11 sorties ParaView pour chacun
des champs `q,p,e_q,e_p`,
vérifie chaque rang, puis compare
les trajectoires, Hamiltoniens, puissances
et géométries au **MUMPS MPI1
précédemment qualifié**
(`v130-wave2d-conservative-five-lu-99d99f2d-tEOjqgHy.zip`)
**sans refaire ses calculs**.

```bash
cd ~/Simulations/scrimp
git pull --ff-only
conda activate scrimp
SCRIMP_DOCKER_NO_SUDO=1 \
SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
bash scripts/run-v130-wave2d-refined-superlu-dist-physical.sh
```

**Statut : SOURCE READY, validation physique
native MPI1/MPI2 EN ATTENTE.**
Un succès établirait un **candidat
SuperLU_DIST raffiné opt-in**
sur cette équation, **pas**
la qualification du solveur
SuperLU_DIST ordinaire. L'acceptation
cinq solveurs par défaut reste
**3/5**, tant que cette distinction
n'a pas été tranchée et documentée.

## Critère de passage à l'étape suivante

Conserver cette fiche comme **audit initial**, et enregistrer
séparément chaque PASS/FAIL GetFEM/PETSc natif.
Ne pas marquer l'exemple 02 PASS avant la
validation native de l'**énergie anisotrope déjà correctement écrite**,
de l'initialisation DAE, des commandes mixtes,
du bilan de puissance et des cinq backends/MPI.
La référence `examples/published/` ne doit pas être
modifiée. Aucun changement n'est à effectuer sur `master`.
