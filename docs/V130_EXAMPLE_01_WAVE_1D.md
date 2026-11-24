# SCRIMP v1.3.0 — Exemple 01 : corde vibrante 1D

**Branche : `release/1.3.0-preparation` uniquement** — aucune fusion
vers `master`, ni publication.  
**Statut scientifique : PASS NATIF — 24 novembre 2026**  
**Code effectivement qualifié :**
`2a170c8247f0f8aff09331747a3bb4812a04324f`.  
**Environnement figé :** `scrimp-mpi:v1.3.0-rc3`,
révision `1d6d0c90115213298047b28237a7f38087563a35`.  
**Preuve complète :**
`v130-first-tutorial-wave1d-2a170c82-cwkGUoXC.zip`.

## 1. Pourquoi ce modèle est le premier exemple

L'ordre réel de Sphinx est **Getting started : corde vibrante 1D**,
puis les exemples spatiaux 2D. Le fichier maintenu est
`examples/wave_1d.py`, **unique source recommandée**.
L'ancien `sandbox/wave_1D.py` est retiré après qualification :
ne pas conserver deux sources concurrentes ni des extraits divergents.
Les huit `examples/published/` demeurent strictement inchangés.

Le modèle illustre simultanément :
variables d'énergie `q=w_x` (déformation) et `p=ρw_t`
(quantité de mouvement), coénergies `e_q=Tq` et `e_p=p/ρ`,
opérateur de Dirac antisymétrique, formulation PFEM,
contrôles de bord, bilan de puissance, initialisation
d'une vraie DAE et intégration Crank–Nicolson.

Sur `0<x<1`, `T(x)=1`, `ρ(x)=1+x(1-x)` :

```text
H(q,p) = 1/2 ∫ (T q² + p²/ρ) dx
q_t = ∂x e_p,   p_t = ∂x e_q
u_L = -e_q(0),  u_R = +e_q(1)
y_L = e_p(0),   y_R = e_p(1)
dH/dt = u_L y_L + u_R y_R.
```

Les briques faibles comportent `+U_L Test_p(0)` et
`+U_R Test_p(1)` : les commandes sont des **tractions
orientées vers l'extérieur**. Les observations entrent avec
la convention pH `-Y_L` / `-Y_R` côté flux ;
elles valent physiquement la **vitesse positive**
aux extrémités. Le fichier exploratoire antérieur
comportait des signes non homogènes avec cette
interprétation — il n'est pas reproduit. Les vrais tests
natif énergétiques ci-dessous valident la correction :
le seul changement de signe du code n'est jamais accepté
comme preuve.

Données initiales délibérément pédagogiques :

```text
q0(x) = 4x(1-x) ∈ P2,  p0(x)=0
u_L(t)=sin(2πt),     u_R(t)=0
```

`q0` s'annule aux deux bords, exactement dans le
même espace P2 que l'approximation. L'initialiseur DAE
**algébrique strict** conserve les données différentielles,
tout en ramenant les contraintes au résidu natif
`≈2.5e-17`. Le modèle comporte 128 DOFs dont 66
algébriques à `h=.05`.

## 2. Résultat de l'acceptation native (10 pas CN)

| Solveur PETSc LU | Rang(s) MPI | État | Défaut énergétique maximal par pas | Écart maximal de la trajectoire complète vs MUMPS MPI1 |
|---|---:|---|---:|---:|
| MUMPS | 1 | **PASS** | 1.790865e-16 | référence |
| MUMPS | 2 | **PASS** | 7.258734e-17 | 2.220446e-15 |
| SuperLU | 1 | **PASS** | 1.790899e-16 | 1.221245e-15 |
| SuperLU_DIST | 1 | **PASS** | 1.790865e-16 | 1.554312e-15 |
| SuperLU_DIST | 2 | **PASS** | 7.264155e-17 | 1.998401e-15 |

Ces cinq simulations utilisent **la même formulation faible
GetFEM et les mêmes 11 instants `0:0.004:0.04`**. Tous
les rangs MPI conservent les traces complètes, la
convergence TS et des bilans physiques vérifiés. Le
travail de bord est réellement non nul : maximum de
puissance `0.090934301708179`. Le bilan CN utilise
la **puissance calculée sur la moyenne des états FE**
de chaque pas, sans intégration trapézoïdale fictive
des échantillons de puissance `P[0]`.

Les quatre erreurs explicites de base sont contrôlées :
erreur de signe, quadrature CN inadéquate, perte de
la donnée différentielle pendant l'initialisation,
historique non fini/non convergé.

**Les VTK sont également natifs :** pour chacune des
cinq combinaisons, 11 VTU `q` (41 points, 20 cellules)
et 11 VTU `p` (21 points, 20 cellules), avec des PVD
correctement horodatés. Le rapport source indique
`PASS_FIRST_WAVE_1D_ALL_AVAILABLE_LU_MPI_PHYSICS_VTK`
et `errors=[]`, vérifié indépendamment dans l'archive.

**SuperLU standard est séquentiel et n'a volontairement
pas de case MPI2.** Aucun benchmark de rapidité ni de
convergence spatiale n'est revendiqué.

## 3. Documentation illustrée et exploitation des résultats

Source documentaire : `docs/source/started.rst`.
Le chapitre contient **deux assets SVG directement
issus des véritables données natives** :

- `docs/source/_static/examples/wave_1d/wave_1d_power_balance.svg`
  : superposition `H(t)-H(0)` et du travail signé des
  deux ports, plus les dix défauts par pas (ordre 1e-16).
- `docs/source/_static/examples/wave_1d/wave_1d_short_animation.svg`
  : **animation pour la documentation HTML** des 11
  trames effectivement calculées (déformation P2,
  puis incrément de déplacement reconstruit).

Le postprocesseur réutilisable
`scripts/render-v130-wave1d-doc.py` génère également,
**sans recalculer la simulation**, depuis l'archive ZIP
ou son dossier extrait :

| Sortie | Usage scientifique ou pédagogique |
|---|---|
| `wave_1d_field_profiles.png/.svg` | Profils de `q`, vitesse `p/ρ`, déplacement reconstruit |
| `wave_1d_spacetime.png/.svg` | Cartes espace-temps des trois champs |
| `wave_1d_energy_balance.png/.svg` | Hamiltonien, puissances des ports, bilan et défaut |
| `wave_1d_solver_parity.png/.svg` | Différences réelles MUMPS/SuperLU/SuperLU_DIST en MPI |
| `wave_1d_string.gif` | Animation native courte (GIF pour Sphinx en ligne) |
| `wave_1d_energy_and_ports.csv` | Séries temporelles et bilan quantitatif |
| `wave_1d_fields.csv` | Profils spatio-temporels exploitables |
| `wave_1d_postprocessed.npz` | Tableaux NumPy compactés |
| `provenance.json` | Commit, image et origine des valeurs |

L'énergie et les puissances proviennent **des
assemblages GetFEM originaux** et des états PETSc
acceptés ; les tracés `q,p` proviennent de véritables
fichiers `Float32` VTK, pas de données inventées.

**Attention au déplacement :** SCRIMP résout `q=w_x`
et `p=ρw_t` ; il ne résout pas directement `w`
dans cette formulation pH. Le postprocesseur
reconstruit `w` en intégrant *exactement* le P2
sur chaque élément spatial puis la vitesse gauche
en temps (trapèzes), avec la jauge `w(0,0)=0`.
Ces données reconstruites servent à la visualisation
et **ne constituent pas une nouvelle validation
de la conservation**.

La courte animation `t<=0.04` illustre le **début**
de l'excitation, et non la propagation à grande échelle.
Une simulation illustrée **optionnelle** sur 1.5 s
(150 pas CN, MUMPS MPI1 seulement), exécutée
séparément et portant ses propres preuves de bilan,
est préparée par :

```bash
SCRIMP_DOCKER_NO_SUDO=1 \
SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
bash scripts/run-v130-wave1d-long-doc-demo.sh
```

**Le scénario long est désormais qualifié nativement** à
`7e2e95a` (MUMPS MPI1, 150 pas CN, 151 états acceptés),
indépendamment des cinq simulations de validation
d'installation courtes. Le GIF de 151 trames et les
figures ont été examinés visuellement, et les **302 VTU**
ont été analysés avec leurs métadonnées PVD, leur contenu
`Float32` et les deux maillages d'éléments finis.
La documentation Sphinx intègre des **SVG animés et
statiques légers** issus de cette simulation, tandis
que le GIF complet reste reproductible depuis l'archive.
Voir la qualification détaillée ci-dessous.

## 3bis. Incident de la démonstration longue : t=1,49 non archivé

**24 novembre 2026 — PREMIER ESSAI LONG : ÉCHEC SCIENTIFIQUE
DU CONTRÔLE D'ARCHIVAGE, PAS DE LA RÉSOLUTION PETSc.**
Archive initiale exacte :
`v130-wave1d-long-doc-f0b589a9-uOwc7fiY.zip`,
source `f0b589a9c94185eef12f4048407b0a5456646288`.
Elle contient le journal natif et le manifeste, mais aucun
fichier de sortie `rank0.npz` ni VTK valide ; ne surtout
pas l'utiliser pour produire des illustrations.

Dans cette simulation `t_f=1,5`, `dt=dt_save=0,01`,
la vraie TS PETSc a réussi **150/150 pas**, sans rejet
ni échec non linéaire. L'initialisation algébrique
était compatible avec un résidu de
`2.48139079e-17`; GetFEM a pu calculer les termes
du Hamiltonien et les puissances des deux ports.

Le monitor original de `DPHS` évaluait
`next_saved_at_t - t <= 0` avec des flottants.
Au pas accepté **i=149, t≈1,49**, l'arrondi du
calcul `previous_saved_time + dt_save` donnait
une valeur supérieure de seulement **un ULP**
au temps PETSc accepté. Le moniteur l'a donc
**ignoré**, puis enregistré `i=150` ; les
états sauvegardés étaient au nombre de **150**
plutôt que **151**. Le contrôle strict du
tutoriel a alors levé une exception :
`Real PETSc CN convergence or energy positivity failed`.
Il n'est **pas** établi que l'énergie soit devenue
négative : cette erreur agrégée déclenche aussi
sur le nombre de snapshots incorrect.

**Correction du code commun** : ajout de
`scrimp/utils/trajectory_sampling.py` et appel
à `saved_time_is_due()` depuis
`DPHS.monitor()`. Le seuil s'autorise au
maximum **16 ULP** d'écart numérique, plafonnés
à `1e-6*dt_save`. Il **n'invente aucun état**
et ne doit pas avancer une échéance physique
d'un intervalle appréciable ; les sorties
avec `dt_save>dt` restent volontairement
sous-échantillonnées.

Le nouveau
`tests/test_v130_trajectory_sampling.py`
reproduit le flottant raté à `i=149`, prouve
que les 151 états des 150 pas sont de nouveau
retenus, et vérifie l'absence de déclenchement
anticipé/sous-échantillonnage modifié. Le
contrôle de l'exemple reste **strict** :
`last_ts_steps == len(times)-1`, défaut
énergétique pas-à-pas et bilan global.

**État final : RÉSOLU, PASS NATIF LONG OBTENU.**
Le second essai au commit `7e2e95a` conserve les
151 instants attendus et satisfait les bilans
énergétiques stricts. Le GIF et les figures ont
été vérifiés sur les vrais résultats de GetFEM.
L'archive d'échec initiale demeure conservée comme
preuve et test de non-régression ; le correctif de
l'archivage n'a pas nécessité d'assouplir les
critères de conservation.

## 3ter. Qualification native de la démonstration longue : PASS

**Source exacte :**
`7e2e95a46645762ac09027428b5bb2abe8cff7a0`.  
**Image native figée :** `scrimp-mpi:v1.3.0-rc3`,
révision `1d6d0c90115213298047b28237a7f38087563a35`.  
**Archive complète :**
`v130-wave1d-long-doc-7e2e95a4-75oCWG1T.zip`
(319 entrées, ~3,5 Mo compressés).

**Résultat effectif du script :
`PASS_REAL_NATIVE_LONG_WAVE_1D_DOC_FIGURES`.**
Les dix régressions d'archivage et les seize tests
du tutoriel passent dans le préflight.

| Grandeur contrôlée | Résultat natif |
|---|---|
| Modèle | Corde 1D, 128 DOFs, `h=0,05` |
| Résolution | PETSc CN fixe, MUMPS, MPI1 |
| Horizon | `t=[0,1,5]`, `dt=dt_save=0,01` |
| Pas de temps / états effectivement conservés | **150 / 151** |
| Rejets TS / échecs SNES | **0 / 0** |
| Initialisation algébrique finale | `2,481390789592951e-17` |
| Hamiltonien min / max | `0,2666666667 / 0,8271913855` |
| Puissance totale min / max | `-0,2030118 / 1,6987573` |
| Défaut maximal par pas CN | `3,434752482434078e-16` |
| Défaut maximal cumulé `H-H0-W` | `7,216449660063518e-16` |
| VTK `q` (P2) | **151 VTU, 41 points, 20 cellules** |
| VTK `p` (P1) | **151 VTU, 21 points, 20 cellules** |
| PVD | **2 séries, 151 instants corrects chacune** |

Les **302 fichiers VTU** ont été relus indépendamment
du message `PASS` : chaque XML contient un champ
`Float32` fini et un maillage de dimension correcte.
Les traces originales enregistrent les 151 vrais
temps de PETSc ; aucune interpolation de trajectoire
n'a été ajoutée par le correctif de moniteur.

**Revue visuelle :** le GIF contient 151 trames, de
dimensions 880×770, et montre clairement la
propagation, le changement de signe de la déformation
et les réflexions. Le déplacement présenté n'est
**pas** une variable supplémentaire de la DAE :
il est reconstruit de façon documentée, avec la
jauge initiale `w(0,0)=0`.

La documentation HTML intègre des fichiers
texte SVG peu volumineux :

- `wave_1d_long_native_animation.svg` : 16 vrais
  instantanés espacés de 0,1 s, extraits parmi les
  151 états acceptés (11 des 41 points P2
  affichés). La réduction est **uniquement visuelle**.
- `wave_1d_long_native_energy.svg` : bilan
  (H(t)-H(0)=W_{\partial\Omega}(t)) et puissance
  de bord signée, sous-échantillonnés pour la figure,
  avec seuil de conservation calculé sur toutes
  les 150 étapes originales.

Les ressources graphiques complètes
`wave_1d_string.gif`, PNG/SVG et CSV/NPZ
sont générables depuis la simulation validée avec
`scripts/render-v130-wave1d-doc.py`; elles n'ont
pas été copiées intégralement dans le dépôt.
Le GIF et les données riches restent dans
l'archive de preuve, afin de ne pas encombrer
les sources et la documentation Sphinx.

**Limite :** ceci qualifie la simulation longue
**MUMPS MPI1 seulement** ; les cinq combinaisons
LU/MPI demeurent qualifiées séparément sur le test
d'installation court. Aucun benchmark de temps
de calcul ou d'ordre de convergence n'est inféré.

## 4. Reproductibilité

```bash
cd ~/Simulations/scrimp
git pull --ff-only
conda activate scrimp

# Accès utilisateur direct à l'exemple
python -m examples.wave_1d --solver mumps --check \
  --export --output outputs/my_first_wave

# Retester les 5 combinaisons dans l'image MPI figée
SCRIMP_DOCKER_NO_SUDO=1 \
SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
bash scripts/run-v130-first-wave-1d-acceptance.sh

# Rendu séparé des sorties à partir de l'archive déjà qualifiée
wave1d_zip="$(ls -t sandbox/outputs/v130-first-tutorial-wave1d-*.zip | head -n 1)"
python scripts/render-v130-wave1d-doc.py \
  --input "$wave1d_zip" \
  --output outputs/wave1d-doc-figures
```

Les résultats sont conservés dans des dossiers
**possédés par l'utilisateur**, sans écrire dans
`examples/published/`. La documentation
Sphinx et la source `wave_1d.py` exposent
exactement la même configuration physique.

## 4bis. Revue de cohérence entre Sphinx et la source (24 novembre 2026)

**Audit demandé après qualification : des écarts étaient encore présents
dans les fragments pédagogiques**, même si les cinq exécutions natives
et les formulations faibles étaient correctes.

Différences corrigées, sans changement du programme qualifié :

| Élément | Ancien bloc `started.rst` | Source native `examples/wave_1d.py` |
|---|---|---|
| Taille du maillage | `h=0.01` | `mesh_size=0.05` |
| Identifiants d'objet | `wave`, `alpha_q`, `alpha_p` | `system`, `q`, `p` |
| Paramètre de tension | `Young's modulus` | `String tension` |
| Projection initiale | `init_step_ts_type="algebraic"` seulement | `init_step=True` également |
| Solveur linéaire | options implicites | `ksp_type="gmres"`, `pc_type="lu"`, backend optionnel |
| Terme et tracé d'énergie | nom générique, ordre des termes inversé, tracé minimal | potentielle puis cinétique, appel de tracé et nom de fichier explicites |

Les **14 briques GetFEM GWFL** ont été
resynchronisées **textuellement** avec les expressions
du code canonique, et les conditions initiales et commandes
restent exactement celles de la preuve native.
La documentation conserve une présentation progressive,
mais l'utilisateur peut désormais exécuter ses blocs
dans l'ordre pour construire **le même modèle**.
Les variables de configuration reprennent exactement
les valeurs par défaut de `wave_1d()`.

Un test source dédié,
`tests/test_v130_wave_1d_doc_parity.py`,
extrait les **11 blocs Python** de la section
`Coding within SCRIMP`, vérifie leur compilation
successive, compare les AST des appels réels SCRIMP
(états, ports, FEM, paramètres, briques, commandes,
DAE, schéma temporel, Hamiltonien, visualisation)
et confronte leurs paramètres aux défauts de la fonction.
La CI exécute désormais ce test. Les différences de
mise en page sont autorisées, mais pas les différences
d'arguments ni les chaînes GWFL.

**Périmètre :** cet audit concerne uniquement
les extraits du tutoriel d'onde 1D. Les exemples
indépendants de ports lumped en amont du chapitre
ne doivent pas être interprétés comme des fragments
du modèle de corde. Le test source supplémentaire
reste à exécuter dans l'environnement utilisateur
après ces changements documentaires.

## 5. Publications scientifiques associées

- **Haine, Matignon, Serhani (2023)**,
  *Numerical Analysis of a Structure-Preserving
  Space-Discretization for an Anisotropic and
  Heterogeneous Boundary Controlled N-Dimensional
  Wave Equation As a Port-Hamiltonian System*,
  DOI **10.4208/ijnam2023-1005**.
- **Ferraro, Fournié, Haine (2024)**,
  *Simulation and control of interactions in multi-physics:
  a Python package for port-Hamiltonian systems*,
  DOI **10.1016/j.ifacol.2024.08.267**.

C'est un **exemple d'enseignement illustrant le
cadre des publications**, pas une prétention de
reproduire numériquement leurs figures ou leurs
tests de convergence.

## Décision

**Premier exemple SCIENTIFIQUEMENT QUALIFIÉ
et correctement exposé dans le guide utilisateur.**
Les simulations courte et longue sont maintenant
qualifiées, avec des figures et animations effectivement
reproduites depuis les données natives. Le contrôle
Sphinx HTML/LaTeX de ces nouvelles références
graphiques reste à exécuter après cette intégration
documentaire, sans rouvrir la physique du modèle.
La suppression de l'ancien fichier
`sandbox/wave_1D.py` constitue un nettoyage
de dépôt **sans perte d'un code publié**.
