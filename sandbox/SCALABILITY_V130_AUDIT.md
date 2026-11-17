# SCRIMP v1.3.0 — audit de scalabilité MPI, 17 novembre 2026

## Résumé décisionnel

**La cible de 50 % d'efficacité à 4 et 6 rangs n'est pas atteinte**
sur le solveur acoustique 2D linéaire, même à ~3 millions de DDL.
Les 96 runs analysés proviennent des quatre archives de l'utilisateur
(3 répétitions × 2 politiques × 4 nombres de rangs × 4 maillages).

**Priorité 1 : vérifier si MUMPS refactorise numériquement la
*Jacobienne constante* à chaque pas CN.** `never` empêche un
`PC.reset()` explicite, mais **n'implique pas** que PETSc/KSP
réutilise la factorisation LU numérique entre deux
`SNESSolve` / pas TS. Le code actuel reconstruit `P` à
chaque `IJacobian()`, même pour les opérateurs linéaires fixes.

**Priorité 2 : réduire le coût des DEUX factorisations coûteuses**
(initialisation de cohérence DAE, premier pas CN). Les
optimisations du régime établi ne suffiront pas à atteindre 50 %
sur les simulations de seulement dix pas.

**Priorité 3 : si nécessaire, passer de LU globale au Krylov
avec préconditionneur bloc port-hamiltonien** (`PCFIELDSPLIT`
+ Schur approximé + AMG/GAMG ou hypre sur les blocs appropriés).
C'est la piste algorithmique la plus crédible pour dépasser
les limitations de scalabilité d'une LU multifrontale globale.

## Provenance et paramètres (mesures utilisateur)

Archives :

- `mpi-coo-h0p04-159028b9-v3D0jn.zip`
- `mpi-coo-h0p02-159028b9-rKkPQi.zip`
- `mpi-coo-h0p01-159028b9-VYXNe1.zip`
- `h0p005.zip` : seulement JSON/logs/figures, sans les NPZ
  massifs. Le `performance_summary.json` inclus atteste que
  la validation des vecteurs complets a été exécutée chez
  l'utilisateur (relative L2 maximale indiquée :
  6.67e-13). **Impossible de revérifier indépendamment
  les états de h=0.005 depuis cette archive.**

Toutes les campagnes sont issues du commit `159028b9`,
avec la même équation acoustique
`examples/wave.py:wave_eq()`, `dt=0.002`,
`t_f=0.02`, **10 pas Crank–Nicolson**,
GetFEM distribué → PETSc **COO**,
MUMPS `never`/`on_change`, mode matrices `compact`,
`p=1,2,4,6`, 3 répétitions par cas. Temps en secondes,
médiane des répétitions pour une configuration.

### Résultats — politique MUMPS never

| h | DDL | p | T_solve médian (s) | T_pas établi médian (s) | E_global | E_pas |
|---|---:|---:|---:|---:|---:|---:|
| 0.04 | 47 942 | 1 | 1.986 | 0.1091 | 1.000 | 1.000 |
| 0.04 | 47 942 | 2 | 1.614 | 0.0790 | 0.615 | 0.690 |
| 0.04 | 47 942 | 4 | 1.396 | 0.0593 | 0.356 | 0.460 |
| 0.04 | 47 942 | 6 | 2.080 | 0.0844 | 0.159 | 0.215 |
| 0.02 | 187 978 | 1 | 8.989 | 0.4904 | 1.000 | 1.000 |
| 0.02 | 187 978 | 2 | 8.077 | 0.4213 | 0.556 | 0.582 |
| 0.02 | 187 978 | 4 | 8.159 | 0.3632 | 0.275 | 0.338 |
| 0.02 | 187 978 | 6 | 8.944 | 0.3548 | 0.168 | 0.230 |
| 0.01 | 744 850 | 1 | 42.387 | 2.4677 | 1.000 | 1.000 |
| 0.01 | 744 850 | 2 | 37.512 | 2.0276 | 0.565 | 0.609 |
| 0.01 | 744 850 | 4 | 36.303 | 1.6667 | 0.292 | 0.370 |
| 0.01 | 744 850 | 6 | 39.470 | 1.6717 | 0.179 | 0.246 |
| 0.005 | 2 966 786 | 1 | 215.666 | 13.3208 | 1.000 | 1.000 |
| 0.005 | 2 966 786 | 2 | 185.937 | 10.7332 | 0.580 | 0.621 |
| 0.005 | 2 966 786 | 4 | 172.187 | 8.3573 | 0.313 | 0.398 |
| 0.005 | 2 966 786 | 6 | 189.972 | 8.6999 | 0.189 | 0.255 |

`on_change` est très proche pour l'acoustique linéaire
(même ordre de coût aux quatre granularités). Ce n'est
donc pas le levier principal ici.

### Critère chiffré à h=0.005

```text
                actuel          cible E>=0.50     baisse nécessaire
p=4 T_solve     172.19 s        <=107.83 s         >=37.4 %
p=6 T_solve     189.97 s        <= 71.89 s         >=62.2 %
p=4 T_step        8.36 s        <=  6.66 s         >=20.3 %
p=6 T_step        8.70 s        <=  4.44 s         >=49.0 %
```

Les deux métriques (**global vs pas établi**) ne sont pas
interchangeables. Optimiser uniquement le pas établi
ne peut pas garantir la cible globale sur dix pas.

### Coûts qui plafonnent les performances

À h=0.005, politique `never` :

- **Premier pas** : ~37.90 s (1 rang), ~42.20 s
  (4 rangs), ~49.78 s (6 rangs).
- **Initialisation de cohérence DAE** (exemple rep1) :
  ~38.2 s (1 rang), ~43.41 s (4 rangs), ~50.27 s
  (6 rangs), selon les lignes `Initialisation done in ...`.
- **Pas établis (n=2..10)** : 13.32 / 8.36 / 8.70 s
  pour 1 / 4 / 6 rangs.
- **Assemblage GetFEM initial** (rep1) : masse ~3.68 s
  et raideur ~13.91 s à 1 rang ; masse ~2.02 s,
  raideur ~8.20 s à 6 rangs. Gains présents mais
  insuffisants au regard des ~190 s du solve global.
- La somme des **moniteurs PETSc mesurés** par rang est
  ~0.17 s sur dix pas et six rangs : la synchronisation
  dans `DPHS.monitor()` existe, mais ne constitue pas
  le goulet dominant de cette campagne.
- Une large part de `T_solve` hors pas établis est
  associée aux étapes initiales et non à une fuite
  mémoire ou au stockage des trajectoires.

Ces observations n'attribuent **pas encore** la durée
des pas établis à la factorisation numérique. Il faut
d'abord lire les compteurs `MatLUFactorNumeric`,
`MatLUFactorSymbolic`, `MatSolve`, `PCSetUp` et
`KSPSolve` dans `PETSc.Log.view()`.

## Revue du code : points concrets

### A. `scrimp/dphs.py` — Jacobienne linéaire reconstruite

`_solve_impl()` installe
`TS.setIJacobian(self._profile_ijacobian,self.J)`.
Dans `IJacobian()`, même lorsque
`linear_mass=True`, `linear_stiffness=True`,
`ts_adapt_type=none`, donc shift constant :

```python
self.tangent_stiffness.copy(P, structure=DIFFERENT_NONZERO_PATTERN)
P.axpy(sig, self.tangent_mass, structure=DIFFERENT_NONZERO_PATTERN)
```

Chaque modification de `P` peut être interprétée
par KSP/PC comme une nouvelle matrice **numérique**.
La logique `SCRIMP_MUMPS_RESET_POLICY=never` supprime
les `PC.reset()` explicites, mais ne demande pas
`KSPSetReusePreconditioner()`. Hypothèse forte :
**MUMPS effectue une LU numérique à chaque pas**.
Aucun compteur natif de ces quatre archives ne permet
encore de la confirmer. La cible expérimentale pour
une matrice réellement constante est **une LU
numérique initiale pour le TS principal**, puis
de simples applications de facteurs aux pas suivants.

PETSc documente :
- `KSPSetReusePreconditioner()` et option
  `-ksp_reuse_preconditioner` ;
- `SNESSetLagJacobian(-2)`, avec persistence
  entre appels `SNESSolve()` ; ne jamais utiliser
  `-1` *avant* la toute première construction ;
- `SNESSetLagPreconditioner()`.
Sources :
https://petsc.org/release/manualpages/KSP/KSPSetReusePreconditioner/
https://petsc.org/release/manualpages/SNES/SNESSetLagJacobian/
https://petsc.org/release/manualpages/SNES/SNESSetLagPreconditioner/

**Avertissement scientifique** : un « brick linéaire »
peut avoir des coefficients dépendant de `t`.
De même `sig` peut changer avec un pas adaptatif.
Il serait **incorrect** d'activer automatiquement
le cache LU pour toute EDP dite « linéaire ».
Démarrer par un test explicitement certifié
matrice+shift constants : modèle acoustique présent,
`dt` fixe, pas d'événement modifiant les coefficients,
et drapeau opt-in. Vérifier la parité de TOUS les états
et des bilans port-hamiltoniens avant adoption.

### B. `scrimp/dphs.py` — double LU liée à l'initialisation DAE

`_solve_impl()` exécute `init_step()` (un
`TS=pseudo` avec son **propre** KSP/PC/MUMPS),
puis crée un **deuxième** TS principal et
une autre LU MUMPS. La cohérence DAE doit être
préservée, mais il pourrait être plus efficace de
résoudre les contraintes algébriques par blocs
sans factoriser le Jacobien global de 3 millions
de DDL deux fois. Expérimenter une projection
exacte/cohérente des conditions initiales, et
démontrer qu'elle produit le même état projeté.

### C. `scrimp/utils/linalg.py` — COO n'est que le pont initial

`convert_gmm_to_petsc()` utilise la voie PETSc
`setPreallocationCOO()`/`setValuesCOO()` seulement
pour l'assemblage initial des masses/raideurs linéaires.
Il effectue sur chaque rang un
`gf.Spmat("copy",M)`, une transposition et
l'extraction CSR **globale** du modèle GetFEM,
puis un COO des contributions dont une grande
partie appartiennent à d'autres rangs.

Optimisation possible : transfert de **lignes
possédées** et partitionnement GetFEM plus étroit,
en validant les contributions et réductions MPI.
Mais les 8–14 s d'assemblage initial à h=0.005
ne peuvent expliquer à eux seuls une simulation
de 190 s. À traiter **après** la LU.

### D. `scrimp/dphs.py` — moniteur

Le monitor contient `Scatter.scatter()`,
`comm.barrier()`, `gf_model.set_time()`
et `gf_model.next_iter()`. En mode mémoire,
il réplique une trajectoire. Ici la somme mesurée
par rang est faible (~0.17 s à h=0.005, 6 rangs),
donc **pas** un levier pour améliorer de moitié
les temps de calcul de cette campagne.

### E. Solveur direct et exploitation des DDL

La formulation mixte `wave.py` contient les
inconnues dynamiques `q,p`, les coénergies
algébriques `e_q,e_p` et les variables de port.
Pour h=0.005, GetFEM annonce ~1 109 844 DDL
pour `q` (DG1) et autant pour `e_q`,
~371 149 pour `p` (CG2) et autant
pour `e_p`. La moitié de l'inconnue est
donc associée aux variables algébriques
de coénergie, hors ports. Une LU multifrontale
globale sur cette structure n'exploite pas
automatiquement la réduction de blocs.

Direction : conserver le système discrétisé
d'origine, mais utiliser `KSPGMRES`/`KSPFGMRES`
+ `PCFIELDSPLIT`, avec blocs cohérents
`(q,p)` et `(e_q,e_p,ports)`. Pour le sous-système
acoustique grad-div, après élimination des
contraintes appropriées, le Schur physique
ressemble à
`M_p + (dt^2/4) D^T M_q^-1 D`
(signes selon conventions). Ce bloc elliptique
peut se prêter à AMG/hypre/GAMG ; GetFEM DG
permet éventuellement de résoudre certains blocs
localement. Une condensation « exacte » de
`e_p` CG2 n'est **pas** locale par défaut :
l'inverse d'une matrice de masse consistante
couplée est global. Ne jamais la remplacer
silencieusement par une masse lumpée ni modifier
les bilans de puissance.

Références :
https://petsc.org/release/manualpages/PC/PCFIELDSPLIT/
https://petsc.org/release/manualpages/PC/PCGAMG/
https://petsc.org/release/manual/ksp/

## Matériel et topologie : vérification obligatoire

Le benchmark épingle Docker sur les CPU logiques
`0,2,4,6,8,10`, mais lance OpenMPI avec
`--bind-to none`. Rien ne garantit que ces
CPU correspondent à six **cœurs physiques
distincts**, ni que les rangs ne migrent pas
sur un ensemble autorisé. Vérifier impérativement :

```bash
lscpu -e=CPU,CORE,SOCKET,NODE,ONLINE
lscpu | grep -E 'Thread\\(s\\) per core|Socket\\(s\\)'
```

Puis établir une mesure A/B avec répartition
physique explicite (`mpiexec --map-by core
--bind-to core --report-bindings`, sous réserve
du hwloc et de l'image Docker). Les tests
actuels ne permettent pas de conclure que
la machine offre six cœurs physiques distincts
au processus Docker. Examiner également les
compteurs `MatLUFactorNumeric` et la bande
passante mémoire/NUMA.

## Roadmap v1.3.0 — ordonnée par rendement possible

1. **Gate 0 (diagnostic, pas de modification scientifique).**
   PETSc `Log.begin()` / `Log.view()` + profiler SCRIMP,
   h=0.02, 1/4/6 rangs, 5 pas CN.
   Nombre et durée de `MatLUFactorNumeric`,
   `MatLUFactorSymbolic`, `MatSolve`, `PCSetUp`,
   `MatAssemblyEnd`, collectives.
   Contrôle matériel des cœurs physiques.
2. **Gate 1 (quick win possible).** Prototype opt-in
   `constant_jacobian_factor_reuse` :
   `KSP.setReusePreconditioner(True)` et/ou
   `SNES` lagging persistant validés sur une
   Jacobienne certifiée constante, afin d'éviter
   les refactorisations MUMPS à chaque pas.
   Comparer séparément LU numérique et MatSolve,
   puis trajectoires complètes/Hamiltonien/puissances.
   Les matrices non linéaires et les pas adaptatifs
   restent sur le chemin sûr actuel.
3. **Gate 2.** Optimiser la projection initiale
   des contraintes DAE : éviter l'excès de
   factorisation globale sans modifier l'état
   initial physiquement admissible.
4. **Gate 3 (changement algorithmique).** Prototype
   Krylov + préconditionneur pHS bloc/Schur/AMG.
   Démontrer convergence, tolérances de résidu,
   précision des champs et bilans port-hamiltoniens.
   C'est la voie la plus plausible pour obtenir
   une **bonne scalabilité forte à 6 rangs** si LU
   direct plafonne, même avec réutilisation.
5. **Gate 4 (secondaire).** Accélérer GetFEM → PETSc
   initial, supprimer extraction CSR globale
   redondante et choisir ordonnancement MUMPS
   (`ICNTL(28)=2` analyse parallèle) si bénéfique.
   Évaluer également MPI/OpenMP hybride dans
   une construction PETSc/MUMPS adaptée.
6. **Gate 5.** Validation Saint-Venant non linéaire
   dissipatif, puis comparaisons acoustiques /
   Saint-Venant, durée plus longue si nécessaire.
   Éviter de déclarer un gain uniquement parce
   que la méthode directe était volontairement
   recalculée inutilement dans le baseline.

## Critères d'acceptation v1.3.0

- Préciser **quel** objectif `E_p>=0.5` :
  temps **global** (10 pas) ou temps **par pas**
  après amorçage, et pour quels `p`.
  Pour l'objectif le plus contraignant :
  `E_6>=0.5` sur les deux métriques, même maillage.
- 3 répétitions minimum, ordre randomisé/alterné,
  CPU réellement isolés, thread BLAS/OMP documentés.
- Au moins h=0.01 et h=0.005, p=1/2/4/6,
  CN fixe dt=0.002, même sortie et mêmes équations.
- Comparaison des **vecteurs complets** (l'utilisateur
  peut fournir les métriques/hachages pour h=0.005),
  résidu discrétisé, nombre d'itérations PETSc,
  contraintes algébriques, énergie et travail
  portuaire. Un simple état final proche n'est
  pas une preuve suffisante.
- La réutilisation d'une LU n'est légitime que si
  les valeurs de la Jacobienne, la structure et
  le shift `sig` n'ont pas changé.
- Les résultats sur une machine à six cœurs ne
  constituent pas une garantie d'efficacité
  universelle sur tout cluster : la version
  1.3.0 doit publier protocole, topologie,
  maillages et nombre de rangs.

**Décision : ne pas poursuivre en priorité
l'optimisation mémoire.** Déjà ~réglée pour
les besoins actuels ; les temps LU/Schur
constituent le vrai verrou de performance.
