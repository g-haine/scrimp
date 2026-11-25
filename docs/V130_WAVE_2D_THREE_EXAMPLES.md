# SCRIMP v1.3.0 — Parcours progressif des trois ondes 2D

**Branche de travail :** `release/1.3.0-preparation`.  
**État au 25 novembre 2026 :** étape **2A native PASS** en
MUMPS MPI1/2, SuperLU MPI1 et SuperLU_DIST explicitement raffiné
MPI1/2 ; l'option standard SuperLU_DIST reste FAIL (gate 3/5).
Étape **2B native PASS en MUMPS MPI1/2** sur le véritable modèle
dissipatif : 10 pas CN, bilan H-H0=W-D, quadrature indépendante,
trajectoires MPI et exports VTK validés. SuperLU_DIST non testé pour 2B.  
**Ordre dans Sphinx :** `examples/wave.rst` comporte la conservative
puis sa section « Adding Damping to the dphs » ; la page
`examples/wave_coenergy.rst` suit séparément. Cet ordre est conservé.

## Passage opérationnel de 2A à 2B

Le modèle **2A conservatif** a passé les 10 pas CN en MUMPS MPI1/2,
SuperLU MPI1 et également avec SuperLU_DIST **explicitement**
`ReplaceTinyPivot + IterRefine` en MPI1/2, sans changer
`examples/wave.py`. La dernière preuve est
`v130-wave2d-refined-dist-physical-95fbd8b2-HfwaSFgK.zip`.
Ce résultat qualifie l'option expérimentale, **pas**
SuperLU_DIST ordinaire. Les choix numériques du noyau SCRIMP
v1.3.0 restent inchangés.

L'étape suivante ne cherche **pas** à requalifier toutes les
politiques de factorisation : elle introduit une difficulté physique
nouvelle dans un **second fichier indépendant**,
`examples/wave_dissipative.py`.

**Invariant pédagogique :** les 20 briques, les quatre commandes
de bord, les données initiales exactes, le tenseur anisotrope,
la densité, DG1/CG2 et le pas CN de 2A sont identiques à 2B.
**Extension unique :** port dissipatif
`(f_r,e_r)`, deux champs CG2, paramètre
`nu(x)=0.5*(2-x)>=0` et cinq briques
`M_r`, `I_r`, `-I_r^T`, `-M_e_r`, `CR_r`.

La puissance volumique de dissipation est vérifiée
sur les **vraies valeurs de GetFEM** de deux manières
indépendantes :
`D_port=integral f_r*e_r` et
`D_volume=integral nu*f_r*f_r`.
Pour chaque pas CN, à l'état FE moyen,
le nouveau critère strict est

```text
deltaH_k - dt_k*(P_B+P_R+P_T+P_L)_k + dt_k*D_k = 0,
H(t)-H(0)-W_boundary(t)+L_viscous(t) = 0,
L_viscous(t)>=0.
```

Le diagnostic vérifie également le résidu du flux
`f_r+e_p`, la dissipation effectivement non nulle,
les **onze vrais états complets** et les **six**
séries GetFEM/ParaView
`q,p,e_q,e_p,f_r,e_r`. La relation faible
`e_r=nu*f_r` est une **projection FEM** :
`nu*f_r` n'est pas en général dans CG2 et il serait
incorrect d'exiger un résidu L2 pointwise nul.

La section 2B de `docs/source/examples/wave.rst`
utilise maintenant **huit `literalinclude`** depuis le
programme canonique, exactement comme 2A.
Le premier essai natif est volontairement limité
à **MUMPS MPI1/2** ; les variantes LU seront
étudiées ensuite, une fois la physique qualifiée.

```bash
cd ~/Simulations/scrimp
git pull --ff-only
conda activate scrimp
SCRIMP_DOCKER_NO_SUDO=1 \
SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
bash scripts/run-v130-wave2d-dissipative-pilot.sh
```

**Statut scientifique 2B (25 novembre 2026) : PASS natif
MUMPS MPI1/2**, après examen des journaux et tableaux natifs du ZIP
cité ci-dessous. L'exemple 2C coénergie reste indépendant
et constitue la prochaine étape pédagogique.

## Qualification native 2B — PASS MUMPS MPI1/2 (25 novembre 2026)

**Preuve authentique examinée :**
`sandbox/outputs/v130-wave2d-dissipative-e582874f-9ni2cPYZ.zip`
(archive utilisateur, **non** ajoutée au dépôt).
Source des simulations :
`e582874fe9d2fa296ab20a5963e4e4ce6437f092` ;
image `scrimp-mpi:v1.3.0-rc3`, révision interne figée
`1d6d0c90115213298047b28237a7f38087563a35`.
Les neuf tests statiques PASS, les huit exemples
`examples/published/` sont restés inchangés.

| Mesure native | MUMPS MPI1 | MUMPS MPI2 |
|---|---:|---:|
| DOFs, états sauvegardés | 2804, 11 | 2804, 11 |
| Pas CN acceptés, rejetés | 10, 0 | 10, 0 |
| DAE algébrique exacte, résidu initial | PASS, 0 | PASS, 0 |
| Défaut maximal énergétique par pas | 2.700341e-18 | 1.771569e-18 |
| Défaut énergétique cumulé maximal | 4.458781e-18 | 2.063372e-18 |
| Écart maximal `D_port-D_volume` | 1.355253e-18 | 6.505213e-19 |
| Énergie dissipée finale `L_viscous` | 1.1909667465223212e-5 | 1.1909667465223223e-5 |
| Exports GetFEM VTK | 6 champs × 11 instants | 6 champs × 11 instants |
| Statut | **PASS** | **PASS** |

À `t_f=0.08`, sur la trajectoire MPI1, l'énergie
tensorielle originale GetFEM vaut
`H=0.0015744391155611648`, le travail signé des quatre
bords `W=0.0015863487830263924` et la perte visqueuse
`L=0.000011909667465223212`.
La relation `H-H0=W-L` est vérifiée à chaque pas CN et
en cumulé, sans modifier les signes, le tenseur,
les vingt briques conservatrices ou les tolérances.

Parité des trajectoires FE originales :
`max|z_MPI1-z_MPI2|=1.7141843500212417e-13` ;
`max|H_MPI1-H_MPI2|=4.119968255444917e-18`.
Les rangs 0 et 1 de MPI2 ont des SHA256 identiques
sur les trajectoires complètes ; les 10 portes physiques
sont `true` sur tous les rangs. Le résidu quadratique
de flux `f_r+e_p` est voisin de `2.3e-35` ou moins.
Le résidu quadratique ponctuel `e_r-nu*f_r` n'est
pas strictement nul (maximum environ `1.33e-8`) :
c'est attendu pour une **loi projetée faiblement**
dans CG2, pas un défaut de passivité. On ne doit pas
le requalifier artificiellement en contrainte ponctuelle.

**Portée :** qualification de l'exemple 2B pour
`h=0.2`, `dt=0.008`, `t_f=0.08`,
`nu(x)=0.5*(2-x)`, et MUMPS sur 1 et 2 rangs.
Ce pilote ne démontre pas convergence en maillage
ou en temps, dissipation sur temps longs, performance,
ni qualification SuperLU_DIST. Le modèle 2A et
`examples/published/` restent intacts.
**Suite :** qualification indépendante de 2C coénergie,
notamment la compliance `T^{-1}`, la localisation des
ports dissipatifs et le bilan `H-H0=W-D_interne-D_feedback`.

## Deuxième pilote 2C : feedback linéaire à gain constant (25 novembre 2026)

Archive native `v130-wave2d-coenergy-0c44a1ae-rwh9fi2L.zip`
(source `0c44a1a`). Les **7 tests source** et les **7 tests
d'interfaces différées** ont passé, les huit exemples publiés
sont intacts. La nouvelle exécution MUMPS MPI1/MPI2 arrive
jusqu'à l'initalisation algébrique DAE, mais **aucun pas CN
n'est encore accepté** :

```text
Non-linear form '(-0.5*Y)*Test_U' has been added ...
ValueError: Algebraic initializer supports only linear mass and stiffness
```

Diagnostic établi directement par les traces de briques GetFEM :
la loi physique `U=-0.5*Y` est exactement linéaire,
mais `DPHS.set_control` ne reconnaissait que les
feedbacks à variable nue (`Y`, `-Y`) comme
`linear=True`. La constante numérique provoquait
une fausse classification `linear_stiffness=False` ;
l'initialiseur DAE algébrique a correctement refusé
cette classification, sans basculer sur une pseudo-TS.

**Correctif générique** dans `scrimp/dphs.py` :
les expressions constituées d'une unique variable
GetFEM déclarée, éventuellement précédée d'un
coefficient **littéral numérique fini** (`-0.5*Y`,
`2.5e-3*Y`, `-Y`) sont désormais certifiées comme
feedbacks linéaires, assemblées par la véritable
raideur GetFEM. Les autres expressions (fonction,
produit d'états, coefficient symbolique, macro ou
ambiguïté) conservent la branche non linéaire.
L'algèbre des opérateurs et le feedback physique
`-0.5` restent inchangés ; les critères CN,
passivité, DAE et MPI ne sont pas assouplis.

`tests/test_v130_linear_scalar_feedback.py` vérifie
cette frontière de classification sur la fonction
Python réelle sans dépendance GetFEM/PETSc ; les
nouveaux tests sont appelés dans la CI et le
préflight du pilote 2C.

**Statut : nouveau candidat commité, exécution native
MPI1/2 et sorties graphiques encore en attente**.
Le premier passage algébrique ne constitue pas encore
une validation du Hamiltonien ou du feedback.

## Premier pilote 2C : interface matérielle transparente (25 novembre 2026)

Première archive native :
`v130-wave2d-coenergy-f110a3e3-2V4VQbbG.zip`
(source `f110a3e3c6c661d865eb1eea3cfbbc5eba2afdac`,
image `scrimp-mpi:v1.3.0-rc3`). Les 7 tests source ont passé,
les 8 `examples/published/` sont intacts, mais les deux
simulations **MPI1 et MPI2 ont échoué avant tout pas CN** :

```text
InterfaceOwnershipError: No declared interface ports on region 10
```

Cause : le compilateur automatique détectait correctement
l'interface géométrique 10 entre les volumes 1 et 2, mais
confondait **frontière matérielle de champs globaux continus**
et **interconnexion explicite par des ports**. En 2C, les
variables `q`, `p` sont définies sur `[1,2]` et aucun
port ni aucune brique n'est déclaré sur 10 : inventer des
ports y modifierait le problème physique.

Correction générique localisée dans
`scrimp/interface_deferred.py` : lorsqu'une interface
conforme ne comporte **ni port ni brique sur son label**,
on rejoue les FEM, paramètres et briques originales
sans création de faces orientées ni de nouveaux ports.
La découverte géométrique reste effective et le rapport
marque `SKIPPED_TRANSPARENT_MATERIAL_INTERFACE`.
Si une interface contient **au moins un port ou une brique
explicite**, le contrôle original fail-closed de
propriété, d'orientation et de placement reste strict,
sans relaxation. Deux tests de régression couvrent ces
scénarios.

**État de 2C** : correctif source commité, **nouvelle
qualification GetFEM/PETSc MPI1/2 à exécuter**, sans
revendiquer de PASS prématuré. Le pilote inclut désormais
les tests d'interfaces en préflight, puis les sorties
numériques et graphiques obligatoires si la simulation
réussit. Aucune modification de la physique 2C, de
`examples/published/`, du choix MUMPS ou des bilans.

## 2C coénergie : candidat natif et visualisations réelles (25 novembre 2026)

**Statut : candidat 2C prêt à tester, PAS encore qualifié**. Aucune
simulation native MPI1/MPI2 2C n'a été constatée à ce stade. Cette
étape ne modifie ni la physique déjà qualifiée de 2A/2B, ni
`scrimp/dphs.py`, ni `examples/published/`.

Le code canonique `examples/wave_coenergy.py` conserve les
**états coénergétiques substitués** `q` (contrainte) et `p`
(vitesse), le disque intérieur `region=1`, la couronne
`region=2`, le bord `region=20` et les dix briques
coénergie. Les corrections scientifiques explicites sont :

- **Compliance authentique** : `Tinv=T^{-1}` pour
  `T=[[5+x,x*y],[x*y,2+y]]`, `rho=3-x`.
  Le code historique utilisait en réalité `T` sous le nom
  `Tinv` : les deux matériaux ne pouvaient pas être considérés
  équivalents.
- **Feedback physiquement dissipatif** :
  `U=-0.5*Y` au bord extérieur. L'ancien signe positif
  `U=+0.5*Y` injectait de la puissance `U*Y`.
- **Viscosité inversible dans le disque** :
  `nu=0.2+10*(0.36-(x*x+y*y))` sur `region=1`
  uniquement. Cela évite le `1/nu` potentiellement singulier
  de l'ancienne loi à l'interface. La modification est
  **physique et documentée**, non une régularisation cachée
  dans le solveur. Pas de port dissipatif dans la couronne.
- **Données initiales représentatives** : vitesse en deux
  bosses lisses, non nulles près du disque intérieur et du
  bord externe ; contrainte initiale nulle ; DAE algébrique
  et CN à pas fixe.

Le bilan exact visé est `H(t)-H(0)+Li(t)+Lb(t)=0`
(sans travail de forçage externe), avec
`Di=integral_region1(e_r²/nu)=-integral_region1(p*e_r)`
et `Db=-integral_boundary20(U*Y)=integral_boundary20(0.5*Y²)`.
Le modèle **substitué** ne possède pas de variable indépendante
`f_r` : appeler `Damping.get_power()` serait incorrect.
Le calcul `integral nu*p²` est consigné comme diagnostic
ponctuel, sans imposer abusivement l'égalité exacte à la
projection faible DG1. Ces mesures, les champs complets,
les dix pas et la parité MPI sont vérifiés sur les vraies
sorties GetFEM/PETSc.

Pilote scientifique natif (MUMPS 1 et 2 rangs, 10 pas,
preuve `NPZ/JSON/VTK/PVD`, CSV et figures après PASS) :

```bash
cd ~/Simulations/scrimp
git pull --ff-only
conda activate scrimp
SCRIMP_DOCKER_NO_SUDO=1 \
SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
bash scripts/run-v130-wave2d-coenergy-pilot.sh
```

**Livrables visuels réels pour les trois modèles 2A/2B/2C** :
`scripts/render-v130-wave2d-native.py` lit les `rank0.npz`
GetFEM et les véritables `pv/<champ>/<champ>.pvd` /
`*.vtu` GetFEM (y compris les données binaires).
Il écrit :

- `energy_balance.png`, `energy_balance.svg` :
  énergie, travail des ports, pertes et défaut énergétique ;
- `field_snapshots.png`, `field_snapshots.svg` :
  champs spatiaux `p`, `|q|`, éventuellement `e_r` ;
- `p_evolution.gif` : animation des états réellement sauvegardés ;
- `energy_and_ports.csv`, `provenance.json`, `README.txt`.
  Les VTU/PVD originaux restent présents et sont ouvrables
  sous ParaView, avec leurs valeurs Float32 explicitées.

La visualisation est intégrée dans les ZIP des campagnes
2A/2B lors des prochaines exécutions valides et du pilote
2C après PASS. Pour créer immédiatement les figures de
la preuve 2B déjà validée, **sans relancer la simulation** :

```bash
python scripts/render-v130-wave2d-native.py \
  --input sandbox/outputs/v130-wave2d-dissipative-e582874f-9ni2cPYZ.zip \
  --case mumps_mpi1 \
  --output sandbox/outputs/wave2d-2b-visuals
```

La qualification scientifique en **dix pas** est volontairement
distincte de la démonstration visuelle en temps plus long
(`h=0.1`, `t_f=1`, `dt=0.01`, par défaut), lancée à
la demande sur un vrai GetFEM/MUMPS MPI1 :

```bash
SCRIMP_DOCKER_NO_SUDO=1 \
SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
bash scripts/run-v130-wave2d-showcase.sh 2c
```

L'argument accepte aussi `2a` et `2b`.
**La démonstration longue ne qualifie ni parité MPI, ni
convergence spatiale/temporelle** ; ses figures ne seront
intégrées dans Sphinx qu'après examen de leurs vraies
données et vérification de leur lisibilité physique.

## Le fil conducteur scientifique

Les trois fichiers Python restent **exécutables et indépendants**,
pour ne pas cacher la construction de la DPHS derrière des helpers.
Les termes physiques, les signes et le choix des éléments finis
sont explicités de façon identique là où les formulations coïncident,
mais la nouvelle difficulté est isolée à chaque étape :

| Séquence | Programme canonique visé | Question nouvelle | Bilan physique visé |
|---|---|---|---|
| **2A** | `examples/wave.py` | Onde conservative anisotrope 2D, traction normale sur 3 côtés, vitesse imposée sur 1 côté, PFEM DG1/CG2 | `H(t)-H(0)=W_\partial(t)` |
| **2B** | `examples/wave_dissipative.py` | Même rectangle, même `T` et `rho`, mêmes 4 bords et excitation, plus un port résistif `(f_r,e_r)` | `H-H_0=W_\partial-D_\mathrm{volume}`, `D_\mathrm{volume}\ge0` |
| **2C** | `examples/wave_coenergy.py` | Même loi mécanique interprétée en coénergie : états `q=e_q` et `p=e_p`, substitutions de coétats, tenseur de compliance `T^{-1}`, domaine à 2 régions, amortissement localisé et feedback au bord | `H-H_0=W_\partial-D_\mathrm{interne}-D_\mathrm{feedback}` |

Les trois formulations ne doivent pas devenir trois variantes
d'une même visualisation sans explication. La 2A enseigne les
**ports et normales**, la 2B la **dissipation explicite**,
la 2C la **substitution physique et la localisation spatiale**.

### Loi matérielle

Conserver la loi testée dans les exemples de 2D,
sans la remplacer par un matériau homogène artificiel :

```text
T(x,y) = [[5+x, x*y], [x*y, 2+y]]
rho(x) = 3-x
```

Les expressions imbriquées GWFL du tenseur sont des
**colonnes** (documenté dans `docs/source/getfem_gwfl.rst`).
La matrice est ici symétrique, et SPD sur le rectangle
(et sur le disque unité envisagé en 2C).

**IMPORTANT :** `q.T.q` contracte le paramètre GetFEM
`T`, il représente bien `q·Tq` et n'est pas `|q|²`.
L'audit précédent accusant une omission de `T` était erroné
et a été rectifié dans
`docs/V130_EXAMPLE_02_WAVE_2D_AUDIT.md`.

En coénergie, **il ne suffit pas de renommer `T` en
`Tinv`**. Si l'on souhaite montrer la **même loi
matérielle**, il faut que `Tinv=T^{-1}`, c'est-à-dire
pour ce tenseur symétrique :

```text
detT = (5+x)*(2+y)-(x*y)^2
Tinv = [[(2+y)/detT, -(x*y)/detT],
        [-(x*y)/detT, (5+x)/detT]]
```

Ces expressions sont mathématiques ; leur représentation
**réelle en syntaxe GWFL** et l'intégration du modèle
doivent encore être qualifiées. La version actuelle
`examples/wave_coenergy.py` donne **le même
tenseur littéral** à `Tinv` que le cas 2A donne à
`T` ; elle ne représente donc pas le même matériau
tant que ce point n'est pas harmonisé. Son intérêt
pédagogique distinct (deux régions, dissipation locale)
reste pleinement pertinent.

### Normales et signes

Pour 2A et 2B, les Neumann portent des **tractions
normales scalaires** :

```text
u_N=n·e_q, y_N=e_p;  u_D=e_p, y_D=n·e_q
dH/dt = ∫Gamma_N u_N*y_N + ∫Gamma_D u_D*y_D
```

L'ancien `wave.rst` écrivait des égalités vecteur/scalaire
sans la contraction `n·e_q`. Les opérateurs GWFL doivent
aussi donner la **même orientation physique des ports
que celle décrite**. La correction proposée ne passe
en production qu'après mesure des quatre vraies puissances.

### Conditions initiales

Les anciennes données incluaient un pulse gaussien
de quantité de mouvement non exactement nul sur
le bord gauche contrôlé, avec projection `pseudo`.
Le scénario candidat utilise **q0=0, p0=0** et une
vitesse gauche `0.2*sin(4*pi*t)*sin(pi*y)`, les
autres tractions étant nulles. C'est **exactement
compatible à t=0**, y compris aux coins du rectangle,
et l'excitation devient non triviale pendant dix pas
CN. Cela dissocie proprement la validation des ports
et de l'initialiseur DAE de l'étude de propagation.

La simulation courte de validation
(`mesh_size=.2, tf=.08, dt=.008`) est séparée
d'une **illustration longue** destinée aux champs
2D, cartes de densité d'énergie et animations.

## État actuel : 2A conservative qualifiée avec MUMPS MPI1/2

Le **candidat scientifique a été promu** dans
`examples/wave.py`, avec sa fonction publique historique
`wave_eq()`, après un vrai pilote MUMPS sur 1 et 2 rangs.
Le double fichier exploratoire
`sandbox/probe_v130_wave2d_conservative.py` est retiré :
la source destinée aux utilisateurs est maintenant
**unique et autonome**, avec des options CLI.

**Preuve native confirmée du pilote :**
`v130-wave2d-mixed-conservative-cf42089d-pI8VXuBa.zip`.
Le pilote et le modèle promu partagent les **mêmes
briques, paramètres, conditions initiales et contrôles** ;
la source promue a seulement renommé la fonction
(`wave_2d_candidate` → `wave_eq`) et ajouté
des marqueurs/commentaires destinés à Sphinx.
Les fonctionnalités natives n'ont pas été modifiées.

| Mesure | MUMPS MPI1 | MUMPS MPI2 |
|---|---:|---:|
| États et DOFs | 11 × 2230 | 11 × 2230 sur chaque rang |
| Pas CN acceptés | 10 | 10 |
| Initialisation algébrique stricte | PASS, résidu 0 | PASS, résidu 0 |
| Défaut énergétique maximal par pas | 2.98156e−18 | 1.78893e−18 |
| Défaut énergétique cumulé | 5.42101e−18 | 3.03577e−18 |
| Puissance maximale au bord gauche | 0.0489963892 | 0.0489963892 |
| Champs VTU/PVD | q, p, e_q, e_p | q, p, e_q, e_p |
| Instantanés par champ | 11 | 11 |

**Les traces FE des deux rangs MPI2 sont
identiques exactement**, et l'écart maximal
des trajectoires globales MPI1/MPI2 reste
`4.16334e-14`. Il s'agit bien de la
**vraie énergie tensorielle GetFEM**
et des puissances `cn_step` des quatre ports.

Le tout premier pilote
`v130-wave2d-mixed-conservative-abea0401-Jm7S0brY.zip`
était **FAIL en export MPI du champ `e_q`** :
le validateur associait incorrectement les
variables de coénergie aux clés du dictionnaire
des ports. Le correctif générique
`get_flow()/get_effort()` conserve la protection
de réplication entre rangs ; le second pilote
`cf42089d` a **entièrement passé** l'export.

### Documentation harmonisée

La partie **2A** de `docs/source/examples/wave.rst`
présente désormais le problème de Cauchy
dynamique cohérent, le tenseur et la densité
hétérogènes, les tractions
`n·e_q`, la commande de vitesse gauche,
les signes et le bilan exact à chaque pas.
Les fragments de code utilisent directement
**huit inclusions `literalinclude`**
depuis `examples/wave.py` : la source Python
est la référence et aucune deuxième
copie des briques GetFEM n'est conservée
dans la documentation.

La partie **2B dissipative** de cette page
et `examples/wave_dissipative.py` sont
conservés pour leur propre qualification ;
l'exemple **2C coénergie**
(`examples/wave_coenergy.py`) suivra
après le port dissipatif.

### Prochaine étape native : cinq backends sur 2A

Le nouveau script
`scripts/run-v130-wave2d-conservative-acceptance.sh`
lance **le véritable programme publié**
`python -m examples.wave`
dans cinq configurations : MUMPS MPI 1/2,
SuperLU **séquentiel** MPI 1, et
SuperLU_DIST MPI 1/2. Il produit
les états FE complets, les vrais
Hamiltoniens/puissances, les sorties VTK/PVD
des quatre variables `q,p,e_q,e_p`,
les comparaisons inter-rangs et inter-solveurs,
et une archive en cas d'échec comme de succès.
La campagne native au commit `99d99f2d` a désormais
confirmé **SuperLU MPI1 PASS**, en plus des deux
cas MUMPS, mais **SuperLU_DIST MPI1/2 FAIL au
premier pas CN**. Les résultats ne doivent
pas être amalgamés.

```bash
cd ~/Simulations/scrimp
git pull --ff-only
conda activate scrimp

SCRIMP_DOCKER_NO_SUDO=1 \
SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
bash scripts/run-v130-wave2d-conservative-acceptance.sh
```

Résultat **recherché, pas encore obtenu** :

```text
PASS_WAVE2D_CONSERVATIVE_FIVE_LU_MPI_CN_PORTS_VTK
```

### Campagne cinq backends : résultats réels, novembre 2026

Preuve : `v130-wave2d-conservative-five-lu-99d99f2d-tEOjqgHy.zip`,
**source canonique** `99d99f2d94c82fe95aac3b39a3aaf3a4c4c27eb8`,
image native `scrimp-mpi:v1.3.0-rc3` inchangée.

| Solveur / rangs | Résultat natif 10 pas CN | Defaut maximal énergétique par pas |
|---|---|---:|
| MUMPS / 1 | **PASS** | 2.9815559743351372e-18 |
| MUMPS / 2 | **PASS** | 1.6263032587282567e-18 |
| SuperLU / 1 | **PASS** | 3.0899761915836876e-18 |
| SuperLU_DIST / 1 | **FAIL premier pas, PETSc error 91** | Indisponible |
| SuperLU_DIST / 2 | **FAIL premier pas, PETSc error 91** | Indisponible |

Les trois résultats PASS comprennent les 2230 DOFs,
la totalité des 11 états CN, une puissance gauche non
nulle (`max|P_L|=0.0489963892`), les quatre ports,
quatre séries de 11 VTU/PVD et l'initialisation
DAE algébrique exacte. Une comparaison directe des
états originaux à MUMPS MPI1 donne
`max|Δz|=9.248157795127554e-14`
(MUMPS MPI2) et
`1.61259894326804e-13`
(SuperLU MPI1). Les puissances concordent
à `1.3e-16` ou mieux.

Les deux scénarios SuperLU_DIST assemblent la masse et
la raideur, effectuent la projection initiale DAE
(`||F_a||=0`) et atteignent le moniteur
`i=0, t=0`, **sans accepter de pas physique**.
PETSc retourne `TSStep DIVERGED_STEP_REJECTED`
(code 91). Le journal ne renseigne **pas encore**
la raison profonde SNES/KSP/PC : aucune conclusion
définitive sur un pivot nul, l'ordonnancement
SuperLU_DIST ou le conditionnement de la
matrice n'est justifiée. Ne pas changer la
physique, les tolérances CN/DAE, ni substituer
silencieusement MUMPS à SuperLU_DIST.

Le diagnostic désormais disponible réutilise
`SCRIMP_TS_FAILURE_DIAGNOSTICS=1` dans
`scrimp/dphs.py` ; il inspecte le solveur
PETSc réel après l'échec **sans réassemblage,
relance du solveur ou changement de paramètres**.
Il exécute uniquement **le premier pas**
SuperLU_DIST MPI1 dans l'image figée :

```bash
cd ~/Simulations/scrimp
git pull --ff-only
conda activate scrimp

SCRIMP_DOCKER_NO_SUDO=1 \
SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
bash scripts/run-v130-wave2d-superlu-dist-diagnostic.sh
```

Le lanceur écrit
`v130-wave2d-superlu-dist-first-step-*.zip`,
conserve ses logs et retourne **échec réel**
en cas de divergence. La clé
`MAIN TS FAILURE DIAGNOSTIC:` devrait
renseigner `SNES_reason`, `KSP_reason`,
`PC_failed` et le type réel de
factorisation. La qualification
**cinq backends reste FAIL**, même si le
diagnostic identifie la cause : il faudra
un vrai run de dix pas avec toutes les preuves
pour l'accepter.

### Facteur de défaillance SuperLU_DIST identifié

Le second diagnostic nativement exécuté, preuve
`v130-wave2d-superlu-dist-first-step-62fb8765-QIaYkugY.zip`,
confirme **PC_FAILED=2 =
PC_FACTOR_NUMERIC_ZEROPIVOT** (énumération PETSc),
`KSP_DIVERGED_PC_FAILED=-11`, SNES `-3`,
TS `-2`. Le pas a été réessayé jusqu'à
`dt=1.9073486328125e-09` et n'a jamais été
accepté : il s'agit d'un blocage numérique
du facteur LU, **pas** d'une réfutation du
bilan port-Hamiltonien du problème.

Une expérience supplémentaire sur les **permutations
réelles de SuperLU_DIST** est préparée, sans
utiliser `ReplaceTinyPivot` ni modifier
la matrice. Elle compare le modèle complet
contre les **états FE MUMPS déjà archivés**,
avec `SCRIMP_MPI_TRACE=1` pour examiner les
lignes/diagonales nulles, et maintient les
vérifications énergétiques.

```bash
cd ~/Simulations/scrimp
git pull --ff-only
conda activate scrimp

SCRIMP_DOCKER_NO_SUDO=1 \
SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
bash scripts/run-v130-wave2d-superlu-dist-ordering-study.sh
```

Ce lanceur éprouve `LargeDiag_AWPM` sur
les lignes, puis `MMD_AT_PLUS_A` sur
les colonnes ; **un éventuel PASS ne
concernera que le candidat MPI1 testé,
et non l'acceptation cinq backends**.
Le résultat natif reste à obtenir. Voir
`docs/V130_EXAMPLE_02_WAVE_2D_AUDIT.md`
pour les valeurs exactes du diagnostic.

### Résultats des permutations et contrôle de rang

La campagne native `8cb2a51f`,
archive
`v130-wave2d-superlu-dist-ordering-8cb2a51f-422NJxVt.zip`,
confirme **deux échecs SuperLU_DIST MPI1** :

- `LargeDiag_AWPM` : **segfault 139**
  immédiatement après `CombBLAS is not available`
  (chemin de permutation indisponible
  dans la configuration binaire) ;
- `MMD_AT_PLUS_A` : **premier pas rejeté**
  avec `PC_FACTOR_NUMERIC_ZEROPIVOT`.

Le Jacobien `2230×2230` contient
**25801 coefficients CSR**, dix diagonales
numériquement nulles (80–89), **aucune ligne
entièrement nulle** et aucune entrée non finie.
La présence de ces diagonales nulles n'implique
**pas** une perte de rang dans un système
pH-DAE à structure de point-selle.

Il faut désormais capturer la **vraie matrice
de préconditionnement** que PETSc tente de
factoriser, comparer MUMPS et SuperLU_DIST,
identifier les dix variables GetFEM concernées,
et calculer le rang **structurel** séparément
de la factorisation **numérique**.

Un banc ciblé existe désormais :
`scripts/run-v130-wave2d-superlu-dist-matrix-audit.sh`,
avec capture native `P,K,M` à
`t=0.008`, `sigma=250`, étude indépendante
SciPy et tests synthétiques de matrices
à diagonales nulles néanmoins inversibles.
Il ne change ni `examples/wave.py` ni
le code de SCRIMP. **Exécution native
en attente.**

```bash
SCRIMP_DOCKER_NO_SUDO=1 \
SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
bash scripts/run-v130-wave2d-superlu-dist-matrix-audit.sh
```

### Analyse de matrice : échec SuperLU_DIST isolé sans singularité du Jacobien

La preuve native
`v130-wave2d-first-cn-jacobian-574d9c71-NKqIZ0BZ.zip`
montre que les opérateurs **MUMPS et
SuperLU_DIST sont exactement identiques**
(`P,K,M` : écart maximal 0,
`P=K+250M` exact). Rang structurel
`2230/2230`, dix diagonales nulles
**uniquement sur `Y_L`** (contrainte
de traction), aucune ligne vide.
Une factorisation SciPy/SuperLU du
vrai Jacobien réussit avec résidu
relatif < `6e-15` sans régularisation.
Condition en norme 1 estimée
`≈4.28e5`, non la valeur exacte.

Pour isoler définitivement le
chemin de factorisation, un replay
PETSc indépendant du système
`GetFEM→TS` a été préparé :
`scripts/run-v130-wave2d-superlu-dist-factor-replay.sh`.
Il charge **le vrai CSR conservé**,
vérifie son SHA256, le distribue
aux rangs MPI1/MPI2, teste
SuperLU_DIST standard puis avec
**`ReplaceTinyPivot` explicitement
activé**. Ce dernier est une
**modification numérique du LU**
et n'est pas introduit dans
`examples/wave.py`. La vérification
porte sur le vrai résidu `Px−b`
et l'erreur par rapport à une
solution connue, sur tous les rangs.
**Replay natif à exécuter.**

```bash
SCRIMP_DOCKER_NO_SUDO=1 \
SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
bash scripts/run-v130-wave2d-superlu-dist-factor-replay.sh
```

Un PASS du replay ne serait
pas un PASS de la simulation
ni une raison de modifier
immédiatement la politique LU
de SCRIMP. **3/5 backends
physiques seulement sont
qualifiés à ce jour.**

### SuperLU_DIST : pivot remplacé, précision encore insuffisante

Le rejeu natif de la matrice `P=K+250M`
réellement utilisée lors du premier pas
conservatif est consigné dans
`v130-wave2d-superlu-dist-factor-replay-e87b20d1-H7IgKhB8.zip`.
Quatre cas ont été réellement exécutés :
SuperLU_DIST standard MPI1/2 échoue
encore avec `PC_FAILED=2`, tandis que
`ReplaceTinyPivot=true` supprime
l'échec de LU et donne `KSP_reason=4`
sur les deux nombres de rangs.
**Le remplacement de pivot NE PASSE
pas les critères de précision** :
résidu arrière ≈ `5.084e-10`
(exigence < `1e-11`), erreur de
solution ≈ `7.3e-7`
(exigence < `1e-7`).
Ce résultat ne doit pas être désigné PASS.

L'étape suivante examine une **correction
de solution, pas de modèle** : le
raffinement interne de SuperLU_DIST,
et un GMRES à droite contrôlant le
véritable résidu non préconditionné.
Deux cas × MPI1/MPI2, avec
`ReplaceTinyPivot` activé uniquement
dans ces essais isolés, les mêmes
coefficients et les mêmes seuils.
Les contrôles sources sont dans
`tests/test_v130_wave2d_refinement.py`.
Commande native :

```bash
SCRIMP_DOCKER_NO_SUDO=1 \
SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
bash scripts/run-v130-wave2d-superlu-dist-refinement.sh
```

**Aucun changement du modèle pédagogique
conservatif, dissipatif ou coénergie**.
Les trois passes historiques restent
qualifiées ; SuperLU_DIST demeure
non qualifié jusqu'au passage de
l'intégralité des dix pas CN et des
contrôles physiques et MPI.

### SuperLU_DIST : récupération de précision sur le facteur CN

L'archive native
`v130-wave2d-superlu-dist-refinement-ad2d58a6-unyfltIn.zip`
confirme le succès **strict** de l'option
SuperLU_DIST `ReplaceTinyPivot` **associée
au raffinement interne `IterRefine`**,
sur le Jacobien original : résidus arrière
`9.54e-17` (MPI1) et `1.27e-16`
(MPI2), erreurs sur solution connue
`2.01e-14` et `2.26e-14`.
La variante GMRES préconditionnée à droite
échoue toujours au critère de résidu
(≈ `2e-10`), malgré un
`KSP_reason` convergé.

L'**essai physique intégral**,
désormais préparé, garde le
**modèle utilisateur original**
`examples/wave.py` complètement
inchangé. Seules les deux options
numériques, sélectionnées
**explicitement en mode recherche**,
sont injectées avant l'appel natif
`TS.setFromOptions`. Il vérifie
dix pas CN, les **2230 états FE**, le
bilan Hamiltonien/travail aux **quatre
ports**, les deux rangs MPI1/2 et
les quatre séries originales
VTK/PVD de 11 temps chacune.
Chaque grandeur sera comparée à
l'archive MUMPS MPI1 déjà validée,
sans relancer les autres solveurs.

```bash
SCRIMP_DOCKER_NO_SUDO=1 \
SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
bash scripts/run-v130-wave2d-refined-superlu-dist-physical.sh
```

**Exécution native encore requise.**
Même avec un succès complet, le
résultat concernera le **solveur
SuperLU_DIST explicitement raffiné**,
et non l'option ordinaire par défaut.
La qualification initiale de
l'installation demeure **3/5**.

Une fois ces cinq cas qualifiés et leurs
preuves inspectées, produire les
**figures scientifiques 2D** :
snapshots `q,p,e_q,e_p`, vitesse,
carte anisotrope d'énergie, puissance
de chaque frontière, bilan
`H-H0-W`, parité LU/MPI.
La simulation courte est un test
d'installation ; une simulation
plus longue apportera une **animation
physique plus lisible**. Ne pas
représenter `q=∇w` comme un
déplacement directement résolu.

## Publications liées (présentes dans la bibliographie Sphinx)

- **Haine, Matignon et Serhani (2023)**,
  *Numerical Analysis of a Structure-Preserving
  Space-Discretization for an Anisotropic and
  Heterogeneous Boundary Controlled N-Dimensional
  Wave Equation As a Port-Hamiltonian System*,
  DOI 10.4208/ijnam2023-1005.
- **Brugnoli, Haine et Matignon (2022)**,
  *Explicit structure-preserving discretization
  of port-Hamiltonian systems with mixed boundary
  control* — MTNS 2022.
- **Serhani, Matignon et Haine (2019)**,
  *Partitioned Finite Element Method for
  port-Hamiltonian systems with Boundary Damping:
  Anisotropic Heterogeneous 2D wave equations* —
  pour le cas dissipatif et ses déclinaisons.

Ces références justifient **l'approche** ;
les nouvelles simulations pédagogiques
ne sont pas des reproductions à l'identique
des courbes publiées.

**Ne modifier ni `master`, ni les huit
`examples/published/` ni des sorties de recherche
non suivies.**
