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

## Critère de passage à l'étape suivante

Conserver cette fiche comme **audit initial**, et enregistrer
séparément chaque PASS/FAIL GetFEM/PETSc natif.
Ne pas marquer l'exemple 02 PASS avant la
validation native de l'**énergie anisotrope déjà correctement écrite**,
de l'initialisation DAE, des commandes mixtes,
du bilan de puissance et des cinq backends/MPI.
La référence `examples/published/` ne doit pas être
modifiée. Aucun changement n'est à effectuer sur `master`.
