# SCRIMP v1.3.0 — Premier exemple Maxwell 3D Nédélec/RT0

**Branche :** `release/1.3.0-preparation`.  
**Statut scientifique actuel (25 novembre 2026) :**
**PASS NATIF Maxwell 3D en MUMPS MPI 1 et MPI 2**,
sur dix pas CN, neuf contrôles physiques, parité
des trajectoires et véritables données VTK tétraédriques.
Cette qualification est distincte d'une étude
de convergence spatiotemporelle ou d'un cas PEC.

## Qualification native Maxwell 3D — PASS (25 novembre 2026)

**Preuve complète :** archive
`v130-maxwell3d-daac18cb-FCPngFrH.zip`,
source `daac18cb36fbc404fcf5c4fe4f2f69121b2ee530`,
image figée `scrimp-mpi:v1.3.0-rc3`
(révision `1d6d0c90115213298047b28237a7f38087563a35`).
`summary.json` porte
`PASS_MAXWELL3D_MUMPS_MPI1_MPI2`, sans erreur.
Les 9 tests source ont passé, les 8 exemples
publiés sont restés immuables.

| Mesure native | MUMPS MPI 1 | MUMPS MPI 2 |
|---|---:|---:|
| Domaine | Boule 3D, 58 points, **110 tétraèdres** | idem |
| DOFs `E` Nédélec 0 + `B` RT0 | **222 + 275 = 497** | **497** |
| Pas CN fixes acceptés | **10/10** | **10/10** |
| Critères physiques | **9/9 PASS** | **9/9 PASS** |
| Énergie initiale `H(0)` | `0.6942520909491595` | `0.6942520909491599` |
| Énergie finale `H(0.2)` | `0.6942520909491590` | `0.6942520909491603` |
| Énergie magnétique finale `H_B(0.2)` | `0.2582694531373567` | `0.2582694531373592` |
| Défaut énergétique cumulé maximal | `2.3314683517128287e-15` | `1.9984014443252818e-15` |
| Maximum `||div B||_L2` | `4.630808670018643e-15` | `4.4379315325408354e-15` |
| Export VTK | `E` + `B`, 11 instants chacun | idem |

**Contrôles indépendants sur `rank*.npz` :**
pour MPI1, l'énergie électrique finale est
`H_E(0.2)=0.4359826378118022` et l'énergie
magnétique finale est `H_B(0.2)=0.2582694531373567`,
soit **37.20110555 % de l'énergie initiale**
transférée vers le champ magnétique.
`max|H-H_E-H_B|=5.55e-17`, et
`max|H-H(0)|=2.33146835e-15`.
La loi de Gauss magnétique a été mesurée
par la quadrature volumique GetFEM de
`(div B)^2` ; sa racine est
inférieure à `4.7e-15` aux onze temps.
Aucun port de puissance, dissipation ou
initialisation DAE n'a été introduit.

Les trois archives `rank0.npz` MPI1,
`rank0.npz` MPI2 et `rank1.npz` MPI2
contiennent chacune **11 × 497 coefficients
FE réels** ; leurs SHA256 recalculés
correspondent aux empreintes JSON.
Les deux rangs MPI2 sauvegardent exactement
les mêmes coefficients. La différence
maximale MPI1/MPI2 sur les 11 états complets
est **`4.440892098500626e-16`** ;
les énergies varient au plus de
`2.7755575615628914e-15` entre ces deux exécutions.

**Validation supplémentaire avec un véritable
lecteur VTK :** les **44 fichiers VTU** du
pilote (2 champs × 11 temps × 2 configurations)
ont été ouverts, et non simplement analysés
en tant que XML. Chacun contient
**110 cellules VTK type 10 (tétraèdres)**,
**58 points en 3D** et les 3 composantes
numériques finies du champ attendu.
Leur géométrie est identique à tous
les instants et entre les configurations
MPI. À `t=0.2`, les normes nodales du
champ magnétique `B` vont de
`0.3673212` à `0.3750458`,
ce qui confirme un champ 3D non trivial.
Quatre fichiers `PVD` décrivent exactement
les onze dates `0,0.02,...,0.2`.

**Portée exacte du PASS :** la
conservation du Hamiltonien et la
contrainte `div B = 0` sont vérifiées
dans **cette** résolution 3D,
`h=0.7`, `dt=0.02`, `t_f=0.2`,
avec MUMPS MPI1/MPI2. Ce résultat ne
constitue ni une convergence h/dt,
ni une validation de la loi de Gauss
électrique, ni une résolution de modes
électromagnétiques sous condition PEC,
ni un test SuperLU_DIST.

## Visualisation intégrée au fichier de l'exemple (25 novembre 2026)

Le fichier canonique **`examples/maxwell_3d.py` génère désormais
lui-même tous les visuels** lors de son exécution avec
`--export --output <répertoire>`. Il n'est
plus nécessaire d'exécuter un post-traitement
externe après le pilote ; les résultats sont
produits et archivés par la même simulation.

Le calcul et les neuf portes scientifiques
sont **strictement inchangés**. Une fois les
exports collectifs `export_to_pv("E")` et
`export_to_pv("B")` terminés, **seul le
rang MPI 0** lit les données VTK et écrit
dans `<répertoire>/visuals/` :

- `maxwell3d_echange_energie.png` et
  `.svg` : énergies `H_E`, `H_B`
  et `H=H_E+H_B` directement issues des
  quadratures FE GetFEM, à chaque temps CN
  sauvegardé ;
- `maxwell3d_vecteurs_E_tfinal.png` et
  `maxwell3d_vecteurs_B_tfinal.png` :
  champ vectoriel 3D **réel**, au temps
  final, avec une trame constituée des
  arêtes extérieures des tétraèdres VTK ;
- `maxwell3d_B_3D_animation.gif` :
  évolution 3D de `B` sur **tous** les
  instants sauvegardés, caméra, domaine
  et longueur de flèches fixes ;
- `provenance.json` et `README.txt` :
  SHA256 du vecteur d'état FE original,
  effectifs de points et de tétraèdres,
  nombre d'images, avertissement
  pédagogique sur les flèches nodales.

**Interprétation rigoureuse des images :**
les flèches représentent des valeurs
vectorielles aux nœuds VTK calculées par
GetFEM, **pas** des degrés de liberté
intrinsèques d'arête Nédélec ou de face RT.
Aucun champ analytique ni état temporel
synthétique n'est introduit pour la
visualisation. Les bilans énergétiques
physiques restent construits sur les
quadratures FE originales, et non les
normes des flèches nodales. La lecture
VTU/PVD utilise `base64`, `ElementTree`,
`NumPy`, `matplotlib` et `Pillow` via
la pile graphique existante de l'image ;
elle ne dépend pas de `vtk` Python.

Le pilote MUMPS MPI1/2 exige maintenant la
présence de **tous ces fichiers**, vérifie
leur format, les **11 images du GIF** et
l'égalité du SHA256 de provenance avec
la trajectoire native de chaque
configuration. Un échec graphique
interrompt le PASS intégré et reste
enregistré dans l'archive d'échec.

**État au moment du changement :**
le modèle Maxwell 3D physique a déjà été
qualifié au commit `daac18cb`, mais la
nouvelle génération automatique des visuels
doit faire l'objet d'une **nouvelle
exécution native**. Ne pas annoncer le
rendu intégré comme validé avant le
PASS correspondant.

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

## Premier essai natif et correctifs préalables (25 novembre 2026)

L'exécution réelle sur l'image `scrimp-mpi:v1.3.0-rc3`
au commit `acdb3c75` a **échoué avant l'assemblage
Maxwell** en MPI 1 et MPI 2. Les 6 tests de source
ont passé et les 8 exemples publiés sont restés
immutables. GetFEM a effectivement construit une
boule 3D de **58 points et 110 éléments**, puis
assigné les éléments intrinsèques attendus :

- `FEM_NEDELEC(3)` à `E` : **222 DOFs** ;
- `FEM_RTK(3,0)` à `B` : **275 DOFs**.

Le blocage à ce stade était **une erreur Python
dans l'audit ajouté par le nouvel exemple**,
pas une erreur démontrée de l'implémentation
Nédélec/RT0 ou des équations :

```text
AttributeError: 'numpy.ndarray' object has no attribute 'char'
```

En effet, `MeshFem.fem()` renvoie un
**tableau NumPy de `Fem`** : l'itération
directe donnait des sous-tableaux, et
non les `Fem` individuels. Correction
dans `examples/maxwell_3d.py` :
`np.asarray(mesh_fem.fem(), dtype=object).ravel()`,
comparaison de l'ensemble des `Fem.char()`
avec le descripteur demandé, contrôles
`qdim=3` et DOFs non nuls. **Aucune
conversion vers CG n'est autorisée**.

Un second défaut concernait le **vérificateur
d'archive**, qui ne pouvait pas démarrer
à cause de `from pathlib` au lieu de
`from pathlib import Path`. Corrigé ;
les tests de source compilent maintenant
le programme Python *à l'intérieur du
heredoc Bash*, et exercent le garde-fou FE
sur des tableaux NumPy imbriqués, des cas
CG substitué et des cas mal formés.
Les erreurs initiales ont empêché la
publication d'une archive complète :
elles n'ont fourni **aucune preuve**
de conservation énergétique ou de
correction de la loi de Gauss.

**Statut à la date de ce premier essai :
non qualifié.** Cette étape est historique ;
voir le PASS natif `daac18cb` en tête
du présent document.

## Second essai natif et audit simplifié (25 novembre 2026)

Archive GetFEM du deuxième essai :
`v130-maxwell3d-ed9b9fd8-gnJLLDYD.zip`,
source exacte `ed9b9fd80f8ceb670085d7a58dc07f625dedbe99`.
**8/8 tests de source PASS**, huit exemples
publiés intacts, création confirmée des mêmes
**110 éléments, 222 DOFs Nédélec et 275 DOFs RT0**
en MPI1 et MPI2. Le pilote produit désormais son
`summary.json` et l'archive même quand les
exécutions natives échouent.

**Échec avant tout pas CN :**
`RuntimeError: Invalid GetFEM MeshFem.fem() entries
for FEM_NEDELEC(3)`. Le retour effectif de l'API
GetFEM ne satisfait pas l'hypothèse du code
sur des objets `Fem` pour chaque entrée, même
après `np.asarray(..., dtype=object).ravel()`.
Aucune donnée de matrice, énergie ou dynamique
Maxwell n'a été produite. **Il serait incorrect
de conclure à une défaillance du patch Nédélec.**

Le contrôle a été simplifié après inspection
de **`scrimp/fem.py`** : lors de l'unique
bascule CG1 de `FEM.set_fem()`, le même objet
SCRIMP `FEM` reçoit explicitement
`self.__type="CG"` et `self.__order=1`.
Le nouvel exemple conserve donc les deux
objets originaux `S.FEM("E",-1,"FEM_NEDELEC(3)")`
et `S.FEM("B",0,"RT")` après
`model.add_FEM()` : il compare **type et ordre**
avec la demande originale, puis vérifie
**`MeshFem.qdim()==3` et `nbdof()>0`**.
Le refus de l'ancien repli CG est conservé
sans introspection ambiguë du tableau
`MeshFem.fem()`. Des tests unitaires
simulent exactement cette mutation CG1,
les erreurs de famille et de dimensions.
Cette vérification porte sur **le comportement
de l'adaptateur de SCRIMP**, et non sur
l'identité spectrale ou la conformité complète
des matrices GetFEM (encore à qualifier).

La mise à jour est limitée au nouvel
exemple, au contrat source et à la
documentation. La formulation faible,
les patchs GetFEM, les paramètres du pilote
et ses seuils physiques sont inchangés.
**Statut à la date de ce second essai :
non qualifié.** Les erreurs ont depuis
été levées ; voir le PASS natif ci-dessus.

## Protocole de qualification exécuté (fail closed)

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

**Résultat du protocole :** les deux exécutions
MUMPS MPI1/2 ont passé tous les critères
détaillés ci-dessus sur les résultats natifs.
Les prochaines tâches sont la réalisation
d'une visualisation pédagogique des
véritables vecteurs 3D `E` et `B`, puis,
séparément, une éventuelle étude de
convergence et des conditions électromagnétiques
de bord plus explicites.
