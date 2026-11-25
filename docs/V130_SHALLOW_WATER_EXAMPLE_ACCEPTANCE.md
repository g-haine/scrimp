# SCRIMP v1.3.0 — Qualification de l'exemple Saint-Venant 2D non visqueux

**Statut au 25 novembre 2026 : nouveau pilote natif borné préparé ;
exécution et qualification encore EN ATTENTE.** Ce document concerne
le vrai `examples/shallow_water.py` du tutoriel Sphinx, **pas** le
modèle tournant visqueux du répertoire `sandbox/`.

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

## Nouvelle étape — simulation non visqueuse plus longue, MPI 1/2

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

**Quatre simulations indépendantes** sont proposées :

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
**ces tests ne sont pas un PASS natif**.

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
