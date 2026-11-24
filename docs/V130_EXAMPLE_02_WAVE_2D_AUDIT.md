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

## Critère de passage à l'étape suivante

Conserver cette fiche comme **audit initial**, et enregistrer
séparément chaque PASS/FAIL GetFEM/PETSc natif.
Ne pas marquer l'exemple 02 PASS avant la
validation native de l'**énergie anisotrope déjà correctement écrite**,
de l'initialisation DAE, des commandes mixtes,
du bilan de puissance et des cinq backends/MPI.
La référence `examples/published/` ne doit pas être
modifiée. Aucun changement n'est à effectuer sur `master`.
