# SCRIMP v1.3.0 — Parcours progressif des trois ondes 2D

**Branche de travail :** `release/1.3.0-preparation`.  
**État du 24 novembre 2026 :** étape **2A CANDIDATE SOURCE**, aucune
nouvelle qualification GetFEM/PETSc.  
**Ordre dans Sphinx :** `examples/wave.rst` comporte la conservative
puis sa section « Adding Damping to the dphs » ; la page
`examples/wave_coenergy.rst` suit séparément. Cet ordre est conservé.

## Le fil conducteur scientifique

Les trois fichiers Python restent **exécutables et indépendants**,
pour ne pas cacher la construction de la DPHS derrière des helpers.
Les termes physiques, les signes et le choix des éléments finis
sont explicités de façon identique là où les formulations coïncident,
mais la nouvelle difficulté est isolée à chaque étape :

| Séquence | Programme canonique visé | Question nouvelle | Bilan physique visé |
|---|---|---|---|
| **2A** | `examples/wave.py` | Onde conservative anisotrope 2D, traction normale sur 3 côtés, vitesse imposée sur 1 côté, PFEM DG1/CG2 | `H(t)-H(0)=W_\partial(t)` |
| **2B** | `examples/wave_dissipative.py` | Même rectangle, même `T` et `rho`, mêmes 4 bords et excitation, plus un port résistif `(f_r,e_r)` | `H-H_0=W_\partial-D_\mathrm{volume}`, `D_\mathrm{volume}\ge0` |
| **2C** | `examples/wave_coenergy.py` | Même loi mécanique interprétée en coénergie : états `q=e_q` et `p=e_p`, substitutions de coétats, tenseur de compliance `T^{-1}`, domaine à 2 régions, amortissement localisé et feedback au bord | `H-H_0=W_\partial-D_\mathrm{interne}-D_\mathrm{feedback}` |

Les trois formulations ne doivent pas devenir trois variantes
d'une même visualisation sans explication. La 2A enseigne les
**ports et normales**, la 2B la **dissipation explicite**,
la 2C la **substitution physique et la localisation spatiale**.

### Loi matérielle

Conserver la loi testée dans les exemples de 2D,
sans la remplacer par un matériau homogène artificiel :

```text
T(x,y) = [[5+x, x*y], [x*y, 2+y]]
rho(x) = 3-x
```

Les expressions imbriquées GWFL du tenseur sont des
**colonnes** (documenté dans `docs/source/getfem_gwfl.rst`).
La matrice est ici symétrique, et SPD sur le rectangle
(et sur le disque unité envisagé en 2C).

**IMPORTANT :** `q.T.q` contracte le paramètre GetFEM
`T`, il représente bien `q·Tq` et n'est pas `|q|²`.
L'audit précédent accusant une omission de `T` était erroné
et a été rectifié dans
`docs/V130_EXAMPLE_02_WAVE_2D_AUDIT.md`.

En coénergie, **il ne suffit pas de renommer `T` en
`Tinv`**. Si l'on souhaite montrer la **même loi
matérielle**, il faut que `Tinv=T^{-1}`, c'est-à-dire
pour ce tenseur symétrique :

```text
detT = (5+x)*(2+y)-(x*y)^2
Tinv = [[(2+y)/detT, -(x*y)/detT],
        [-(x*y)/detT, (5+x)/detT]]
```

Ces expressions sont mathématiques ; leur représentation
**réelle en syntaxe GWFL** et l'intégration du modèle
doivent encore être qualifiées. La version actuelle
`examples/wave_coenergy.py` donne **le même
tenseur littéral** à `Tinv` que le cas 2A donne à
`T` ; elle ne représente donc pas le même matériau
tant que ce point n'est pas harmonisé. Son intérêt
pédagogique distinct (deux régions, dissipation locale)
reste pleinement pertinent.

### Normales et signes

Pour 2A et 2B, les Neumann portent des **tractions
normales scalaires** :

```text
u_N=n·e_q, y_N=e_p;  u_D=e_p, y_D=n·e_q
dH/dt = ∫Gamma_N u_N*y_N + ∫Gamma_D u_D*y_D
```

L'ancien `wave.rst` écrivait des égalités vecteur/scalaire
sans la contraction `n·e_q`. Les opérateurs GWFL doivent
aussi donner la **même orientation physique des ports
que celle décrite**. La correction proposée ne passe
en production qu'après mesure des quatre vraies puissances.

### Conditions initiales

Les anciennes données incluaient un pulse gaussien
de quantité de mouvement non exactement nul sur
le bord gauche contrôlé, avec projection `pseudo`.
Le scénario candidat utilise **q0=0, p0=0** et une
vitesse gauche `0.2*sin(4*pi*t)*sin(pi*y)`, les
autres tractions étant nulles. C'est **exactement
compatible à t=0**, y compris aux coins du rectangle,
et l'excitation devient non triviale pendant dix pas
CN. Cela dissocie proprement la validation des ports
et de l'initialiseur DAE de l'étude de propagation.

La simulation courte de validation
(`mesh_size=.2, tf=.08, dt=.008`) est séparée
d'une **illustration longue** destinée aux champs
2D, cartes de densité d'énergie et animations.

## État actuel : étape 2A — premier banc scientifique

Le candidat est enregistré dans
`sandbox/probe_v130_wave2d_conservative.py`,
**sans remplacer encore le code livré**
`examples/wave.py` ni ses blocs Sphinx.
La simulation effectue le **véritable** assemblage
GetFEM et la vraie intégration PETSc :

- q vectoriel DG1, p scalaire CG2, traces DG1,
  paramètres tensoriels et densité d'origine ;
- quatre véritables ports avec orientation
  Neumann/Dirichlet explicitée dans les briques ;
- initialisation DAE `algebraic` stricte à
  données différentielles conservées ;
- 10 pas CN fixes et enregistrement de tous
  les 11 états acceptés ;
- calcul original `H`, puissances séparées des
  quatre ports `cn_step`, intégrale du travail,
  défauts énergétiques **à chaque pas** ;
- états PETSc complets, diagnostics de rang,
  exports ParaView q, p, e_q et e_p.

**Premier gate : MUMPS MPI 1/2 seulement**, car
une erreur mathématique éventuelle doit être
comprise avant de multiplier les solveurs.

```bash
cd ~/Simulations/scrimp
git pull --ff-only
conda activate scrimp

SCRIMP_DOCKER_NO_SUDO=1 \
SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
bash scripts/run-v130-wave2d-conservative-pilot.sh
```

Résultat **recherché, pas encore obtenu** :

```text
PASS_WAVE2D_CONSERVATIVE_TENSOR_MIXED_MUMPS_MPI1_MPI2
```

Le script produit une archive indépendante comprenant
les journaux natifs même en cas d'échec ; aucun FAIL
n'est réinterprété comme un PASS, et aucun résultat
n'est inventé. Les tests source associés sont
`tests/test_v130_wave2d_conservative_pilot.py`.

Si la physique du candidat passe, **promouvoir 2A**
vers `examples/wave.py`, réécrire les seuls blocs
conservatifs de `wave.rst` en parfaite parité AST,
puis qualifier MUMPS/SuperLU/SuperLU_DIST, MPI1/2.
Produire SVG/PNG 2D, ParaView, GIF et CSV/NPZ depuis
les sorties originales. Une fois **2A achevé**,
reprendre les mêmes paramètres pour **2B** et
qualifier le terme `D`, puis adapter le domaine et
la compliance pour **2C**.

## Résultat du premier pilote 2A et correction de l'export MPI

**24 novembre 2026 : les deux vraies simulations MUMPS MPI1
et MPI2 réussissent la physique du modèle conservatif** :
10 pas CN acceptés et 11 états par rang,
2230 DOFs, DAE initialisée strictement
(résidu 0), puissance gauche non triviale
(`max|P_L|≈0.0489964`) et défaut d'énergie
par pas au plus `2.982e-18`. Le test de
réplication MPI donne **exactement les mêmes
états sur les deux rangs**, avec un écart
global MPI1/MPI2 de `4.08e-14`.

**Le premier lanceur reste FAILED** parce que
le validateur root-only de ParaView ignorait
les champs d'effort/coénergie `e_q` et
`e_p` en recherchant `self.ports.get("e_q")`
au lieu de reconnaître le propriétaire
`q↔e_q`. Le MPI1 avait pourtant exporté
tous les 44 VTU des quatre champs (11 chacun),
et le MPI2 avait déjà produit les 22 VTU
de `q` et `p`.

Le validateur est corrigé **génériquement**
pour résoudre l'unique propriétaire de
chaque variable par `get_flow()/get_effort()`,
sans neutraliser la vérification du maillage et
des états MPI ni autoriser des sorties concurrentes.
Régressions source élargies, dont champs
`q,p,e_q,e_p`, variables de bord, propriétaire
ambigu et inconnus. La **nouvelle validation
MPI1/MPI2 complète reste à exécuter** avant
de promouvoir le candidat.

Preuve initiale :
`v130-wave2d-mixed-conservative-abea0401-Jm7S0brY.zip`.
Détails chiffrés :
`docs/V130_EXAMPLE_02_WAVE_2D_AUDIT.md`.

## Publications liées (présentes dans la bibliographie Sphinx)

- **Haine, Matignon et Serhani (2023)**,
  *Numerical Analysis of a Structure-Preserving
  Space-Discretization for an Anisotropic and
  Heterogeneous Boundary Controlled N-Dimensional
  Wave Equation As a Port-Hamiltonian System*,
  DOI 10.4208/ijnam2023-1005.
- **Brugnoli, Haine et Matignon (2022)**,
  *Explicit structure-preserving discretization
  of port-Hamiltonian systems with mixed boundary
  control* — MTNS 2022.
- **Serhani, Matignon et Haine (2019)**,
  *Partitioned Finite Element Method for
  port-Hamiltonian systems with Boundary Damping:
  Anisotropic Heterogeneous 2D wave equations* —
  pour le cas dissipatif et ses déclinaisons.

Ces références justifient **l'approche** ;
les nouvelles simulations pédagogiques
ne sont pas des reproductions à l'identique
des courbes publiées.

**Ne modifier ni `master`, ni les huit
`examples/published/` ni des sorties de recherche
non suivies.**
