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

**Ce scénario long n'est pas encore qualifié nativement.**
Il ne remet pas en question les cinq PASS déjà obtenus.
La documentation en ligne peut employer l'animation
SVG courte dès maintenant ; la future animation GIF
longue ne sera publiée qu'après examen des preuves.

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
python scripts/render-v130-wave1d-doc.py \
  --input sandbox/outputs/<archive-native>.zip \
  --output outputs/wave1d-doc-figures
```

Les résultats sont conservés dans des dossiers
**possédés par l'utilisateur**, sans écrire dans
`examples/published/`. La documentation
Sphinx et la source `wave_1d.py` exposent
exactement la même configuration physique.

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
Restent à vérifier avant d'éventuelles figures de
propagation longue : l'exécution native séparée,
ses bilans énergétiques et les images produites.
La suppression de l'ancien fichier
`sandbox/wave_1D.py` constitue un nettoyage
de dépôt **sans perte d'un code publié**.
