# SCRIMP v1.3.0 — Qualification de l'exemple Saint-Venant 2D non visqueux

**Statut au 25 novembre 2026 : PASS natif du pilote BDF2 borné,
MPI 1/2 et deux raffinements temporels au même temps final.**
Archive GetFEM/PETSc originale vérifiée indépendamment :
`v130-shallow-original-extended-9d7a4df7-c8ZnnOFN.zip`.
La qualification concerne la trajectoire `T=0.002` avec
initialisation `nonlinear_affine` choisie explicitement,
**pas** l'exécution publique historique `T=0.5` / `pseudo`. Ce document concerne
le vrai `examples/shallow_water.py` du tutoriel Sphinx, **pas** le
modèle tournant visqueux du répertoire `sandbox/`.

## Qualification native acquise — campagne MPI1/MPI2 à temps final fixe

**Provenance exacte :** source Git
`9d7a4df78f4fe48e9186e5d255cc9f95e5873291`,
image `scrimp-mpi:v1.3.0-rc3`, révision OCI
`1d6d0c90115213298047b28237a7f38087563a35`.
Archive utilisateur originale
`v130-shallow-original-extended-9d7a4df7-c8ZnnOFN.zip`,
**SHA256**
`47a2dd674d45372314b281e1624fd60e28d77f46b9db37b4b0f881ab4fc92d70`.
Le `summary.json` original indique
`PASS_SWE_BOUNDED_NATIVE_BDF2_MPI1_MPI2`,
`errors=[]` et
`published_examples_modified=false`.
Les quatre processus Docker se sont terminés avec
`exit_code=0`, sans aucun pas rejeté ni échec SNES.

**Relecture indépendante du ZIP :** les six fichiers
`shallow_fe_rank*.npz` contiennent les **vrais coefficients FE
complets de 4 270 DOFs**, les dates, le Hamiltonien, la
masse et les puissances des frontières, un fichier par rang.
La relecture vérifie le maillage temporel de 21 ou 41
états, la fin exacte à `t=0.002`,
la cohérence NPZ/rapports JSON,
`H(t)-H(0)-\int P_{bord}dt`,
la reproductibilité coefficient par coefficient
sur MPI1/MPI2, et l'identité
**bit-à-bit des deux rangs MPI2**.
Le schéma BDF2, les données prescrites et le
terme `G` explicitement retardé restent originaux ;
la sélection `nonlinear_affine` ne modifie
que l'initialisation DAE du test.

| Grandeur native (à `T=0.002`) | 20 pas / MPI1 | 20 pas / MPI2 | 40 pas / MPI1 | 40 pas / MPI2 |
|---|---:|---:|---:|---:|
| `dt` | `1e-4` | `1e-4` | `5e-5` | `5e-5` |
| Degrés de liberté | 4270 | 4270 | 4270 | 4270 |
| États FE réels enregistrés | 21 | 21 | 41 | 41 |
| Nombre de rejets / échecs SNES | 0 / 0 | 0 / 0 | 0 / 0 | 0 / 0 |
| Dérive relative de masse | `3.0199e-15` | `1.7764e-15` | `4.2634e-15` | `2.1317e-15` |
| `max abs(P_bord)` | 0 | 0 | 0 | 0 |
| `H(0)` | `31659.721892373385` | `31659.721892373353` | `31659.721892373385` | `31659.721892373353` |
| `H(T)` | `31659.714296554033` | `31659.714296554073` | `31659.718048003433` | `31659.718048003422` |
| `H(T)-H(0)-W_{bord}` | `-0.00759581935199094` | `-0.007595819279231364` | `-0.0038443699522758834` | `-0.0038443699304480106` |
| Dérive relative du Hamiltonien | `2.39921e-7` | `2.39921e-7` | `1.21428e-7` | `1.21428e-7` |

**Parité numérique (maximum absolu sur TOUS les états temporels) :**
pour 20 pas, `max |z_MPI1-z_MPI2|=9.094947017729282e-11`,
`max |H_MPI1-H_MPI2|=1.4551915228366852e-10` ;
pour 40 pas, respectivement
`1.2732925824820995e-10` et
`1.4915713109076023e-10`.
Les deux tableaux originaux des rangs MPI2
sont **strictement égaux** pour
`t,z,H,mass,boundary_power`.

**Tendance de la perte énergétique BDF2 :**
à temps final et maillage spatial identiques,
`|ε_H(1e-4)|/|ε_H(5e-5)|=1.9758294457312005`,
soit un exposant **observé** sur seulement deux points
`log2(ratio)=0.9824584183712797`.
Le défaut mesuré est **négatif** dans les quatre cas
(perte numérique d'énergie, sans travail au bord).
Avec le splitting gyroscopique retardé, ce
comportement est **compatible avec une erreur de bilan
au premier ordre**, sans établir l'ordre de convergence
de la solution, un exposant asymptotique, ni une propriété
exactement dissipative de tout le modèle.
Le champ initial `h_0` étant discontinu, une étude
d'erreur spatiale et plusieurs `dt`
sont nécessaires pour conclure davantage.
La sortie antérieure `relative_energy_defect_exchanged=1`
est **tautologique** dans ce problème fermé
(`W_{bord}=0`) et ne constitue pas un critère physique
de succès ; l'audit retient le **défaut signé**
et sa normalisation par `H(0)`.

**Classification finale :**
`PASS_SWE_BOUNDED_NATIVE_BDF2_MPI1_MPI2`.
Le PASS couvre exactement 4 calculs physiques
(MUMPS, 20/40 pas, MPI1/2, `T=0.002`),
le jacobien de masse et l'initialiseur ayant déjà
leurs preuves natives séparées.
Il **ne** qualifie **pas** `T=0.5` avec
initialisation `pseudo`, la convergence en
maillage, la conservation exacte BDF2, le Saint-Venant
visqueux tournant ni le solveur DAE d'indice supérieur.

## Périmètre scientifique exact

Le modèle original possède une surface libre `h`, un vecteur
`p = rho*v`, une masse d'évolution de quantité de mouvement
`M_p(h) = h*p.Test_p`, quatre commandes normales homogènes,
un terme gyroscopique `G` explicitement retardé, et une intégration
**BDF2**. Le maillage rectangle, les éléments finis, les conditions
initiales discontinu(es) `h0=3-(sign(x-.5)+1)/3` et `p0=0`,
et les équations des briques ne sont **pas modifiés** pour ce pilote.

**Point de rigueur :** la v1.3.0 corrige le jacobien de la masse
dépendant de `h`, mais `BDF2 + G explicit=True` n'est pas un
intégrateur à gradient discret garantissant `H(t)=H(0)`.
Même avec travail nul à la frontière, on doit **mesurer**
`H(T)-H(0)-∫P_bord dt` et son comportement en raffinant `dt`,
non exiger l'annulation à la précision machine comme pour une
formulation linéaire CN conservative.

## Preuves natives précédentes — à conserver distinctes

- Le jacobien de la matrice de masse, avec `G` complètement
  explicite au dernier état accepté, a été vérifié par différences
  finies directionnelles sur l'assemblage réel GetFEM/PETSc
  en MPI 1 et 2. Cette validation n'est pas une intégration en temps.
- L'initialisation DAE `nonlinear_affine` a passé les contrôles
  originaux sur le modèle non linéaire en `h`, dont les équations
  algébriques sont **affines dans les inconnues algébriques**.
  Cela ne qualifie pas un solveur DAE d'indice arbitrairement élevé.
- De courtes trajectoires physiques BDF2 MUMPS MPI 1/2 et une
  campagne à même temps final `T=0.0004`
  (`dt=0.0001/0.00005/0.000025`,
  4/8/16 pas) ont fait l'objet de vrais journaux GetFEM/PETSc.
  Leurs défauts énergétiques ont été **mesurés**, sans
  revendiquer un ordre temporel théorique ni un bilan
  exactement conservatif. Voir
  `docs/V130_RELEASE_SPRINT_2026W48.md`,
  la campagne `scripts/run-v130-issue28-shallow-energy-refinement.sh`
  et son comparateur.

## Procédure native de reproduction — campagne acquise

```bash
cd ~/Simulations/scrimp
git pull --ff-only
conda activate scrimp

SCRIMP_DOCKER_NO_SUDO=1 \
SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
bash scripts/run-v130-shallow-water-extended-physical.sh
```

Le pilote reprend le **véritable modèle canonique** :
`sandbox/probe_v130_issue28_shallow_short_physical.py`
ne neutralise que le lancement par défaut des 5 000 pas et
les figures pendant la construction. Il sélectionne ensuite
l'initialisation de production `nonlinear_affine`
avec préservation bit-à-bit des données différentielles
originales, et les mêmes BDF2, FEM, `G` retardé et
contrôles nuls.

**Quatre simulations indépendantes, toutes exécutées et qualifiées :**

| Temps final fixé | Pas de temps BDF2 | Pas acceptés visés | Nombre de rangs |
|---|---:|---:|---|
| `0.002` | `0.0001` | 20 | MPI 1 et MPI 2 |
| `0.002` | `0.00005` | 40 | MPI 1 et MPI 2 |

Le probe sauvegarde désormais **en option**
`--capture-full-fe`, qui écrit les véritables coefficients
FE à tous les instants originaux enregistrés
dans un NPZ **par rang**, sans changer la simulation.
Le comparateur externe
`scripts/verify-v130-shallow-water-extended.py`
réévalue depuis ces données :

1. Les étapes PETSc effectivement acceptées, rejets, SNES,
   initialisation DAE stricte et conservation des données initiales.
2. La parité des **coefficients FE complets** et des intégrales
   physiques entre MPI 1/2, et la concordance des deux rangs MPI 2.
3. La masse intégrale, le Hamiltonien `H_E+H_P`,
   le travail aux quatre frontières et le défaut indépendant
   `H(t)-H(0)-∫ P_bord dt`.
4. La comparaison **observée, non présupposée** du défaut
   entre 20 et 40 pas à temps final identique,
   sans annoncer un ordre BDF2 pour une donnée `h_0` discontinue.

Les bornes de contrôle physique du probe précédent
sont conservées : dérive de masse relative
`≤ 1e-5`, dérive de `H` relative `≤ 5e-3`,
puissance externe maximale `≤ 1e-7`.
La parité FE est contrôlée indépendamment avec
des tolérances explicites. Un échec de ces critères
reste un **FAIL** documenté, pas un prétexte pour
modifier la loi physique.

**Traçabilité :** le pilote utilise l'image épinglée
`scrimp-mpi:v1.3.0-rc3` avec sources actuelles en
lecture seule, crée un répertoire unique sous
`sandbox/outputs/v130-shallow-original-extended-*/`,
conserve rapports par rang, fichiers NPZ, journaux,
`summary.json` et une archive ZIP même en échec.
Le contrôle d'immutabilité des huit fichiers
`examples/published/` est exécuté au préalable.
Des tests négatifs hôtes valident la résistance
du comparateur à des coefficients FE incohérents ;
**ces seuls tests ne constituent pas un PASS natif**. Le PASS physique est désormais établi par l'archive native vérifiée ci-dessus.

## Qualifications encore ouvertes

Le vrai horizon historique par défaut de
`examples/shallow_water.py` (`t_f=0.5`,
`dt=0.0001` et `init_step_ts_type="pseudo"`)
n'est **pas qualifié par ce pilote**, qui utilise
l'initialiseur de production `nonlinear_affine`
**explicitement en option**. Un tableau de
décroissance du défaut BDF2 ne constitue pas
une preuve de convergence spatiale ou d'ordre
sans analyse et données supplémentaires.

Le **Saint-Venant visqueux tournant**
(`sandbox/viscous_swe_2d_mpi.py`, géométrie P2,
reprise de checkpoint, commande tangentielle et
retour normal) est un **autre exemple** :
son problème initial de contrainte d'indice
supérieur ne se déduit pas du cas rectangulaire.
Voir `docs/V130_ROTATING_CONTROLLED_INITIAL_STATE.md`
pour les validations et limites propres à ce modèle.
