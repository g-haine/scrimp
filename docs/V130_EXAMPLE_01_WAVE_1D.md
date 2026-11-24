# v1.3.0 — Exemple 01 : corde vibrante 1D

**Branche : release/1.3.0-preparation uniquement.**  
**Statut scientifique : SOURCE STAGED / NATIVE PENDING.**  
Ne pas présenter le nouveau tutoriel comme validé avant l'exécution
du banc d'essai sur une installation GetFEM/PETSc/MPI réelle.

## Place dans le parcours utilisateur

L'ordre effectif de Sphinx est : présentation générale,
**Getting started — corde vibrante 1D**, puis la page d'exemples
2D (onde, chaleur, onde en coénergie, chaleur-onde, Saint-Venant,
oscillateurs lumped). Cet exemple doit être le tout premier
« hello world » scientifique, et non un benchmark de performance.

**Source canonique :** examples/wave_1d.py.  
**Documentation :** docs/source/started.rst.  
**Ancien témoin :** sandbox/wave_1D.py (inchangé, non recommandé
comme test v1.3.0). Les huit sources sous examples/published
restent absolument immuables.

## Ce qu'il enseigne

La corde 0<x<1 avec densité linéique hétérogène
rho(x)=1+x(1-x), tension T(x)=1, déformation
q=w_x et quantité de mouvement p=rho*w_t :

    H(q,p) = 1/2 integral_0^1 (T*q*q + p*p/rho) dx
    e_q=T*q, e_p=p/rho
    q_t = partial_x e_p
    p_t = partial_x e_q

Les forces de traction au bord, avec **normales extérieures**,
sont u_L=-e_q(0), u_R=e_q(1). Les sorties colocalisées
sont les vitesses y_L=e_p(0), y_R=e_p(1). Le
bilan de puissance attendu est :

    dH/dt = u_L*y_L + u_R*y_R.

Après intégration faible de p_t, les deux termes de bord
du code sont **+U_L*Test_p(0) et +U_R*Test_p(1)**.
Le signe moins de l'observation appartient à la
variable de **flux** (-Y) du système port-Hamiltonien :

    -Y_L = -e_p(0),    -Y_R = -e_p(1).

Ces signes explicites évitent de confondre :
(1) la traction de contrôle, (2) la normale extérieure,
(3) la convention d'assemblage interne de GetFEM et
(4) le signe de la puissance physique.

**Différence significative documentée :** l'ancien
sandbox/wave_1D.py utilisait -U_L côté dynamique,
+Y_L/+Y_R côté flow, et des signes -B^T non homogènes
sur les deux bords. Cette ancienne combinaison ne
correspond pas directement à la formule de traction
extérieure et d'observation positive de la documentation.
La nouvelle source aligne le système faible sur la
convention théorique ; il faut encore la valider
avec le calcul GetFEM réel et **ne pas supposer que
changer un signe vaut preuve physique**.

## Pourquoi ce scénario et ces données

On choisit q0=2*x*(1-x)*exp(-50*(x-.5)^2), p0=0.
Le facteur x*(1-x) annule la traction initiale aux
extrémités ; les commandes sont u_L(t)=sin(2*pi*t),
u_R=0. L'instant initial est ainsi conçu pour être
compatible avec l'initialiseur DAE strict du code.
Une courte simulation t=0..0.04, dt=0.004 avec
Crank–Nicolson fixe donne dix pas et onze snapshots.

L'énergie est une somme **des vrais termes d'intégration
GetFEM**, pas une norme de tableau arbitraire. En CN,
le terme Port.compute()[k>=1] est la puissance sur
l'état FE moyen de l'étape se terminant à t[k] :
on ne doit **pas** réintégrer ces valeurs par
trapèzes en réutilisant P[0]. Le test vérifie
pour **chaque** pas :

    H_k - H_(k-1) - (t_k-t_(k-1))*(P_L[k]+P_R[k]) ~= 0

et le défaut cumulatif, ainsi que le caractère
non trivial du travail de bord, la convergence TS
et la préservation exacte des données initiales
dynamiques lors de la projection algébrique.

Ce modèle est une **illustration pédagogique** de
Haine, Matignon et Serhani (2023),
*Numerical Analysis of a Structure-Preserving Space-Discretization
for an Anisotropic and Heterogeneous Boundary Controlled
N-Dimensional Wave Equation As a Port-Hamiltonian System*,
DOI 10.4208/ijnam2023-1005.
Il illustre aussi le cadre du logiciel de
Ferraro, Fournié et Haine (2024),
DOI 10.1016/j.ifacol.2024.08.267.
**Il ne reproduit pas** les figures ni les tests de convergence
de ces publications.

## Matrice de test d'installation

| PETSc LU | Série | MPI 2 | Statut |
|---|---|---|---|
| MUMPS | oui | oui | NON ENCORE QUALIFIÉ sur ce nouveau tutoriel |
| SuperLU | oui | **non** : solveur séquentiel | NON ENCORE QUALIFIÉ sur ce nouveau tutoriel |
| SuperLU_DIST | oui | oui | NON ENCORE QUALIFIÉ sur ce nouveau tutoriel |

Chaque installation doit satisfaire le même modèle et
les mêmes preuves indépendantes, pas seulement réussir
un factor solve PETSc. Le script natif crée et archive
les trajectoires FE, Hamiltonien, puissances gauche/droite,
les résidus énergétiques, rapports DAE et fichiers VTK.
Les résultats des cinq solveurs doivent être comparables
à la trajectoire MUMPS MPI1, avec une tolérance explicite.
Aucun benchmark de temps de calcul n'est revendiqué.

**Commande native, image qualifiée rc3 figée, sans reconstruire Docker :**

```bash
cd ~/Simulations/scrimp
git pull --ff-only
conda activate scrimp

SCRIMP_DOCKER_NO_SUDO=1 \
SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
bash scripts/run-v130-first-wave-1d-acceptance.sh
```

Statut recherché, et **non encore obtenu** :

    PASS_FIRST_WAVE_1D_ALL_AVAILABLE_LU_MPI_PHYSICS_VTK

En cas d'échec, le lanceur enregistre les cinq codes
de sortie, tous les journaux et les traces disponibles
dans une nouvelle archive ZIP. Ne modifier ni source
historique publiée ni équilibre énergétique pour faire
passer un comparateur.

## Conditions de clôture de cet exemple

1. Tests source sans GetFEM/PETSc réussis (signes,
   opérateurs, quadrature CN, contrôles négatifs).
2. Cinq exécutions natives réussies, même géométrie
   et paramètres ; tous les bilans physiques PASS.
3. VTK/PVD 1D réels lisibles, instantanés et
   géométries conformes entre backends et MPI.
4. Les fragments mathématiques et logiciels de
   docs/source/started.rst représentent exactement
   la source ; versions et lien de publication corrects.
5. Refaire si besoin l'illustration du bilan avec
   les vrais résultats de cette campagne et en
   consigner le commit, l'image et le paramétrage.

**Ne pas commencer l'exemple 2D avant la validation
ou le diagnostic explicite de ce premier jalon.**
