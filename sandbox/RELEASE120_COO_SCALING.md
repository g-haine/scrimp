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
pas fixe `dt=0.002`, `t_f=0.04` (20 pas), même cadence
`dt_save=0.01`, même PETSc LU/MUMPS et **même image Docker
`scrimp-mpi:local`**, épinglage sur les cœurs physiques
`0,2,4,6,8,10`, un seul thread OpenBLAS/OpenMP.
Toutes les configurations restent sur le **même maillage**.

La validation exige que chaque run atteigne t_f, fasse 20 pas
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
maillage et les mêmes 20 pas CN. Il produit **trois** figures
(temps, accélération forte, efficacité). Aucune accélération
**inter-versions** n'est calculée ni revendiquée.

Le répertoire produit est
`sandbox/outputs/mpi-coo-h0p025-<sha>-<identifiant>/`.
Les figures sont rangées ainsi :

```text
sandbox/outputs/<campagne>/
├── png/
│   ├── 02_strong_scaling_coo.png
│   ├── 03_speedup_coo.png
│   └── 04_efficiency_coo.png
├── svg/
│   ├── 02_strong_scaling_coo.svg
│   ├── 03_speedup_coo.svg
│   └── 04_efficiency_coo.svg
├── performance_medians.csv
├── performance_summary.json
└── current_*_rep*.json / current_*_rep*_state.npz / logs
```

Pour la campagne standard v1.2.0 + MPI, une **quatrième**
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

# Courte campagne : environ 12 682 DDL et 20 pas CN par test ;
# 1 release + 8 COO/MPI cas, chaque cas répété 3 fois.
bash sandbox/run_release120_coo_scaling.sh 0.08 3
```

Un essai unique (sans médiane statistiquement robuste) :
`bash sandbox/run_release120_coo_scaling.sh 0.08 1`.
Pour un maillage plus grand, conservant exactement le même protocole :
`bash sandbox/run_release120_coo_scaling.sh 0.04 3`
(~47 942 DDL ; davantage de calcul).

Le lanceur créera `sandbox/outputs/release120-coo-h0p08-<sha>-<id>/`.
Il produit **quatre figures** en PNG 240 dpi (sous `png/`)
**et SVG vectoriel** (sous `svg/`) lorsque la comparaison avec
v1.2.0 est admissible :

| Fichier | Question illustrée |
|---|---|
| `01_version120_vs_coo_serial` | gain algorithme actuel vs release, 1 processus / 1 rang |
| `02_strong_scaling_coo` | durée COO + never/on_change sur 1/2/4/6 rangs |
| `03_speedup_coo` | accélération forte, courbe idéale `S=p` |
| `04_efficiency_coo` | efficacité parallèle, idéal 100 % |

Les données et contrôles sont conservés dans
`performance_medians.csv`, `performance_summary.json`, les
rapports JSON par essai, les états complets initiaux/finaux NPZ
et les logs par configuration.

**Précaution :** les résultats anciens CSR vs COO sur une même
révision ne sont pas des benchmarks de la release v1.2.0.
Ils ne seront pas mélangés à cette campagne.

## Note sur le chantier mémoire

Le chantier mémoire est mis en pause. Reprendre
`sandbox/MEMORY_NEXT_STEPS.md` : prochaine priorité, distinguer
GetFEM, PETSc, MUMPS et le coût de base des bibliothèques.
