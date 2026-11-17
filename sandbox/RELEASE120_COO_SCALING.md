# SCRIMP : v1.2.0 séquentielle vs COO/MPI (17 novembre 2026)

**Objectif : figures de performance reproductibles pour une
présentation**. Ce document décrit la campagne à exécuter et la
provenance des données ; il ne fournit **aucun gain inventé**.

## Deux questions différentes (à ne pas confondre)

1. **Gain de version, à un seul processus** : la release originale
   `v1.2.0` (exécutée avec `python` seul, **sans `mpiexec`**) face à
   la branche actuelle **sur un seul rang PETSc**, GetFEM distribué +
   transfert COO + `MUMPS=never` / `on_change`. Cette mesure inclut
   les changements logiciels entre versions ; elle ne saurait être
   attribuée à la seule technologie COO ou au seul MPI.
2. **Scalabilité forte** de la version actuelle : maillage **identique**
   sur `p=1,2,4,6` rangs MPI, séparément avec COO+`never` et
   COO+`on_change`. L'accélération interne à chaque politique vaut
   `S_p = median(T_1)/median(T_p)`, et l'efficacité `E_p=S_p/p`.
   Une accélération comparée à la release est une **autre** métrique :
   `median(T_release)/median(T_current,p)`.

**Le tag v1.2.0 ne possède pas les nouveaux réglages COO et
`SCRIMP_MUMPS_RESET_POLICY`** : on le mesure avec son code d'origine
et sa politique de factorisation d'origine. On ne prétend pas
qu'il tourne avec ces paramètres qui n'existaient pas.

## Reproductibilité et invariants

Le fichier `examples/wave.py` est **bit à bit identique** sur
`v1.2.0` et la branche courante (blob SHA
`aad91dd77734c3307b71da64821097014306a658`).
Le pilote `sandbox/bench_wave_release_scaling.py` est commun aux
deux variantes et charge la même EDP avant l'appel à `solve`.
La release originale est isolée dans un **worktree détaché** créé
temporairement (puis nettoyé).

Pour chaque cas : même h, même nombre de DDL, Crank–Nicolson,
pas fixe `dt=0.002`, `t_f=0.02` (**10 pas**, selon le lanceur fourni), même cadence
`dt_save=0.01`, même PETSc LU/MUMPS et **même image Docker
`scrimp-mpi:local`**, épinglage sur les cœurs physiques
`0,2,4,6,8,10`, un seul thread OpenBLAS/OpenMP.
Toutes les configurations restent sur le **même maillage**.

La validation exige que chaque run atteigne t_f, fasse 10 pas
acceptés sans rejet ni échec non linéaire. Le pilote conserve le
**vecteur initial et final complet** : l'analyse compare coefficient
par coefficient les deux états contre le calcul natif v1.2.0,
avec une erreur relative L2 et Linf < 1e-7 au final.
Si les calculs physiques ne sont pas équivalents, les graphes ne
sont pas publiés. La release v1.2.0 ne fournit pas une propriété
`last_ts_converged_reason` : le contrôle se fait via le moniteur
PETSc TS, en enregistrant steps/rejections/SNES failures et t_f.
La référence et les runs actuels sont identifiés par leurs SHAs.

La durée **`wall_solve_s`** mesure `set_time_scheme` puis
`DPHS.solve`, avec une barrière MPI aux limites de la mesure.
Le temps de construction du modèle est conservé séparément
(`wall_build_s`). La durée annoncée est donc une durée complète
de configuration/résolution, pas uniquement un temps noyau MUMPS.
Un échantillon unique ne donne aucune information statistique ;
avec ≥2 répétitions, le graphe utilise la médiane et représente
la plage min/max des mesures de scalabilité.

## Campagne `--mpi-only` : maillages plus fins sans release

Pour privilégier la **scalabilité forte** COO sur un maillage plus
fin, sans compiler, extraire ni lancer la release historique :

```bash
cd ~/Simulations/scrimp
git pull --ff-only

# 1, 2, 4 et 6 rangs ; COO + never et COO + on_change ;
# h = 0.025, 2 répétitions par configuration.
bash sandbox/run_release120_coo_scaling.sh --mpi-only 0.025 2

# La forme avec l'option à la fin fonctionne aussi :
bash sandbox/run_release120_coo_scaling.sh 0.025 2 --mpi-only
```

Le mode `--mpi-only` **n'accède pas au tag v1.2.0**,
ne crée aucun worktree temporaire et n'exécute aucun calcul de
la release. Il valide les états complets entre les différentes
configurations de la **branche actuelle**, toujours avec le même
maillage et les mêmes **10 pas CN**. Il produit **sept** figures :
trois globales (temps, accélération forte, efficacité) et quatre
sur les temps par pas (temps moyen en régime établi, accélération,
efficacité, coût du premier pas). Aucune accélération
**inter-versions** n'est calculée ni revendiquée.

Le répertoire produit est
`sandbox/outputs/mpi-coo-h0p025-<sha>-<identifiant>/`.
Les figures sont rangées ainsi :

```text
sandbox/outputs/<campagne>/
├── png/
│   ├── 02_strong_scaling_coo.png
│   ├── 03_speedup_coo.png
│   ├── 04_efficiency_coo.png
│   ├── 05_steady_step_time_coo.png
│   ├── 06_steady_step_speedup_coo.png
│   ├── 07_steady_step_efficiency_coo.png
│   └── 08_first_step_overhead_coo.png
├── svg/
│   ├── 02_strong_scaling_coo.svg
│   ├── 03_speedup_coo.svg
│   ├── 04_efficiency_coo.svg
│   ├── 05_steady_step_time_coo.svg
│   ├── 06_steady_step_speedup_coo.svg
│   ├── 07_steady_step_efficiency_coo.svg
│   └── 08_first_step_overhead_coo.svg
├── performance_medians.csv
├── performance_summary.json
└── current_*_rep*.json / current_*_rep*_state.npz / logs
```

Pour la campagne standard v1.2.0 + MPI, une **huitième**
figure `01_version120_vs_coo_serial` apparaît dans `png/`
et `svg/`, uniquement si la parité des états est vérifiée.

Pour **réanalyser sans aucun nouveau calcul** les résultats d'une
campagne déjà terminée (y compris `--mpi-only`) :

```bash
bash sandbox/analyze_existing_release120_scaling.sh \
  sandbox/outputs/mpi-coo-h0p025-<sha>-<identifiant>
```

Les données et rapports restent à la racine de chaque dossier
de campagne. Les figures ne sont plus déposées en vrac.

## Commande pour les figures de ce soir

Depuis le dépôt de développement :

```bash
cd ~/Simulations/scrimp
git pull --ff-only
bash -n sandbox/run_release120_coo_scaling.sh

# Courte campagne : environ 12 682 DDL et 10 pas CN par test ;
# 1 release + 8 COO/MPI cas, chaque cas répété 3 fois.
bash sandbox/run_release120_coo_scaling.sh 0.08 3
```

Un essai unique (sans médiane statistiquement robuste) :
`bash sandbox/run_release120_coo_scaling.sh 0.08 1`.
Pour un maillage plus grand, conservant exactement le même protocole :
`bash sandbox/run_release120_coo_scaling.sh 0.04 3`
(~47 942 DDL ; davantage de calcul).

Le lanceur créera `sandbox/outputs/release120-coo-h0p08-<sha>-<id>/`.
Il produit **huit figures** en PNG 240 dpi (sous `png/`)
**et SVG vectoriel** (sous `svg/`) lorsque la comparaison avec
v1.2.0 est admissible :

| Fichier | Question illustrée |
|---|---|
| `01_version120_vs_coo_serial` | gain algorithme actuel vs release, 1 processus / 1 rang |
| `02_strong_scaling_coo` | durée COO + never/on_change sur 1/2/4/6 rangs |
| `03_speedup_coo` | accélération forte, courbe idéale `S=p` |
| `04_efficiency_coo` | efficacité parallèle sur temps total, idéal 100 % |
| `05_steady_step_time_coo` | temps **réellement mesuré** d'un pas CN moyen (pas 2–10, sans sauvegarde) |
| `06_steady_step_speedup_coo` | accélération MPI **par pas établi**, référence 1 rang de la même politique |
| `07_steady_step_efficiency_coo` | efficacité MPI **par pas établi**, idéal 100 % |
| `08_first_step_overhead_coo` | pénalité du premier pas (dont facteur MUMPS) par rapport aux suivants |

Les données et contrôles sont conservés dans
`performance_medians.csv`, `performance_summary.json`, les
rapports JSON par essai, les états complets initiaux/finaux NPZ
et les logs par configuration.

**Précaution :** les résultats anciens CSR vs COO sur une même
révision ne sont pas des benchmarks de la release v1.2.0.
Ils ne seront pas mélangés à cette campagne.

## Méthodologie temporelle (mesure directe par pas CN)

**Ne jamais utiliser** `T_solve/N` comme approximation du coût
d'un pas en régime établi : cette division inclurait artificiellement
les assemblages, la mise en cohérence DAE et la première factorisation.

Le pilote `sandbox/bench_wave_release_scaling.py` place une horloge
monotone à l'**entrée** et à la **sortie** des moniteurs PETSc TS,
sans MPI collective supplémentaire et sans modifier les pas
numériques. Pour chaque rang et chaque pas accepté `n ≥ 1` :

```text
dt_step[n, rank] =
    perf_counter_at_monitor_entry(n)
    - perf_counter_at_monitor_exit(n-1)

T_step[n] = max_rank(dt_step[n, rank])
```

Les **durées des moniteurs et de leurs sauvegardes sont exclues** de
ces intervalles. Le premier pas `n=1` est mesuré séparément, car
il englobe l'amorçage du solveur et notamment la première
factorisation. En régime établi, avec **10 pas acceptés**, nous
calculons sur `n=2..10` :

```text
T_steady(p, policy) = mean(T_step[2..10])
S_step(p, policy) = median_repetitions(T_steady(1, policy))
                    / median_repetitions(T_steady(p, policy))
E_step(p, policy) = S_step(p, policy) / p
```

Nous affichons la **médiane des moyennes par run**, et enregistrons
aussi chaque pas/rang dans les JSON, le temps du premier pas,
et les moyennes séparées sur le début et la fin de l'intervalle
(pour détecter une dérive de coût avec le temps physique).
Le tableau CSV contient les valeurs et les accélérations en plus
des temps totaux. Pas de barres d'erreur inférées avec une
seule répétition : la zone min–max est tracée à partir d'au moins
deux runs.

Ces métriques rendent l'étude **moins sensible au temps final**
pour le modèle acoustique linéaire, mais ne constituent pas une
preuve générale d'indépendance vis-à-vis de l'état physique pour
les systèmes non linéaires. Une campagne plus longue est nécessaire
seulement pour étudier une éventuelle dérive, pas pour fabriquer
la scalabilité à partir d'une division par N.

**Compatibilité** : l'outil peut redessiner les trois anciennes
figures globales à partir des campagnes précédentes. Il ne trace
les nouvelles figures par pas **que si les chronométrages détaillés
sont réellement présents et validés** dans les JSON : aucune
extrapolation ni valeur artificielle ne sont admises.

## Note sur le chantier mémoire

Le chantier mémoire est mis en pause. Reprendre
`sandbox/MEMORY_NEXT_STEPS.md` : prochaine priorité, distinguer
GetFEM, PETSc, MUMPS et le coût de base des bibliothèques.
