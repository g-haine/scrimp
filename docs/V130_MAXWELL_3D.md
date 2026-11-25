# SCRIMP v1.3.0 — Premier exemple Maxwell 3D Nédélec/RT0

**Branche :** `release/1.3.0-preparation`.  
**Statut initial (25 novembre 2026) :** implémentation et contrat statique
préparés ; **simulation GetFEM/PETSc MPI1/MPI2 non encore exécutée**.
Aucun résultat scientifique PASS n'est revendiqué avant les archives natives.

## Objectif et justification

Montrer en quelques lignes de physique et **quatre briques SCRIMP**
qu'une DPHS peut être construite sur un **vrai maillage tétraédrique 3D**
avec les éléments intrinsèques **Nédélec d'arête (Hcurl)** et
**Raviart–Thomas RT0 de face (Hdiv)**. L'exemple exerce notamment
la famille Nédélec disponible dans le GetFEM 5.5 corrigé de SCRIMP.
Ce n'est pas une simulation CG vectorielle, un modèle 2D extrudé
ou une reproduction du cas Maxwell 2D de
`examples/published/COMPEL_2025/` (répertoire intouchable).

**Choix minimal :** boule de rayon 1, `eps=mu=1` sans unités,
aucun courant ni commande, aucune dissipation. Champs états
`E` et `B` (induction magnétique). Interconnexion skew :

```text
Mass_E dE/dt =  C_transpose B
Mass_B dB/dt = -C E
H = 1/2 integral_Omega(E.E + B.B)
H(t) = H(0)
```

La forme faible exacte est
`(dE/dt,v)=(B,curl(v))` et
`(dB/dt,w)=-(curl(E),w)`.
L'identité pH est obtenue par les **deux termes
d'intégrale égaux et opposés**, sans utiliser une
source ou un port de bord artificiel. La première
équation correspond à une fermeture naturelle de
travail tangentiel magnétique nul dans cette
formulation faible ; ce n'est **pas** l'imposition
explicite de `n×E=0` (PEC). Une extension
multibords à puissance de Poynting est différée.

**Compatibilité FE :** `E_h` appartient aux Nédélec
d'ordre le plus bas, `B_h` à RT0. La donnée initiale
`E0=(-y,x,0)`, `B0=0` est solénoïdale en
continu et `curl(E0)=(0,0,2)`, donc excite
réellement le champ magnétique. Le curl discret
de Nédélec0 vit dans le sous-espace sans divergence
de RT0 ; le pilote vérifiera numériquement
`||div(B_h)||^2_L2` (loi de Gauss magnétique).

**Prudence scientifique :** l'énergie pH ne prouve
pas `div(E)=0` ni la compatibilité de toutes les
traces avec un problème Maxwell fort et régulier.
Les modes propres de cavité, la contrainte Gauss
électrique et la convergence de h/dt sont hors du
périmètre de ce premier exemple.

## Code, tutoriel, tests

- Modèle canonique : `examples/maxwell_3d.py` ;
  masse E/B, deux couplages rot explicites en GWFL,
  état initial projeté en L2 sur les vrais éléments.
  Pas d'initialiseur DAE : deux états à masse
  positive, sans port algébrique.
- Tutoriel Sphinx : `docs/source/examples/maxwell_3d.rst`
  (liens `literalinclude` sur le code canonique),
  indexé dans `docs/source/examples.rst`.
- Contrat source : `tests/test_v130_maxwell3d_source.py`,
  exécuté par la CI et le pilote ; aucun substitut
  à une exécution native 3D.
- Pilote : `scripts/run-v130-maxwell3d-pilot.sh`,
  image épinglée `scrimp-mpi:v1.3.0-rc3`,
  mémoire compact / transfert COO / solveur MUMPS,
  1 et 2 rangs MPI.

Le modèle SCRIMP emploie l'ancien adaptateur
`S.FEM(...)` qui peut revenir silencieusement
à CG1 après `AssertionError` GetFEM. Le
nouvel exemple **vérifie explicitement la famille
réellement installée pour chaque état** et
refuse l'exécution si un élément est remplacé.

## Qualification attendue (fail closed)

L'exécution native proposée dure **10 pas CN
fixes**, `h=0.7`, `dt=0.02`, `t_f=0.2`.
L'archive contient le vecteur de tous les
coefficients FE à chacun des 11 temps
(`rank*.npz`), les énergies originales
GetFEM et les quadratures indépendantes
`H_E`, `H_B` ; le champ `B` doit
devenir non trivial.

**Les neuf contrôles physiques** portent sur
l'authenticité 3D, CN fixe et sans rejet,
l'absence de DAE artificielle, l'énergie
initiale positive, les transferts E→B,
la parité énergie Hamiltonien/quadratures,
la conservation énergétique par pas,
la conservation cumulée et la loi
magnétique Gauss `div B=0`.

Le runner recalcule indépendamment les
bilans, contrôle les **deux rangs MPI2**,
compare toutes les trajectoires FE MPI1/2
et valide **11 instants VTU pour E et B**,
avec les vrais types VTK **tétraèdre**
(10 ou 24), pas seulement un titre « 3D ».

```bash
cd ~/Simulations/scrimp
git pull --ff-only
conda activate scrimp

SCRIMP_DOCKER_NO_SUDO=1 \
SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
bash scripts/run-v130-maxwell3d-pilot.sh
```

L'archive `sandbox/outputs/v130-maxwell3d-*.zip`
est **conservée même en cas d'échec** ; un
`PASS_MAXWELL3D_MUMPS_MPI1_MPI2` ne sera
accepté que sur les vrais résultats GetFEM.
Aucun test de convergence, SuperLU_DIST,
performance ou visuel pédagogique final
n'est revendiqué avant les qualifications.

**Prochaine étape :** exécuter sur l'hôte et
examiner la première archive. Si GetFEM
refuse le macro `Curl3`, l'initialisation
Nédélec/RT0 ou l'export 3D, diagnostiquer
le point précis dans les journaux de la
simulation, sans remplacer silencieusement
les FE ni détendre les critères.
