# SCRIMP v1.3.0 — Passage en chemins de production : interfaces et ParaView MPI

**Décision : 24 novembre 2026.**  
**Branche unique : `release/1.3.0-preparation`.**  
**Ni `master`, ni fusion, ni tag, ni publication PyPI/Docker.**  
**Image native de qualification : `scrimp-mpi:v1.3.0-rc3`**, révision immuable
`1d6d0c90115213298047b28237a7f38087563a35`.

## Décision de mise en production dans la branche de préparation

Sur demande explicite du responsable du projet et compte tenu du calendrier,
les deux développements sortent du mode pilote : **leurs chemins d'exécution
deviennent les chemins par défaut du CODE de la branche de préparation**.
Aucune revendication de validation universelle n'en découle.

1. **Compilation automatique des interfaces internes conformes.** L'étiquette
   physique de l'interface reste inchangée ; les deux faces effectives sont
   attachées aux convexes des domaines adjacents, avec leurs propres normales
   extérieures. Les FEM et briques sont rejoués dans l'ordre de déclaration.
   Les interfaces ambiguës, non conformes ou les formes impossibles à affecter
   à un domaine **échouent explicitement** plutôt que de produire un faux
   opérateur. Sans interface interne reconnue, les anciens modèles restent
   dans leur chemin classique.
2. **Export ParaView MPI avec un seul écrivain.** Au-delà d'un rang MPI,
   le rang 0 écrit les fichiers `.pvd` et `.vtu` communs, au lieu de
   laisser tous les rangs écrire simultanément les mêmes noms. Le nouveau
   chemin par défaut exige une **preuve collective d'états complets et
   identiques et de maillage géométrique identique** avant l'export.
   En cas d'absence de cette preuve, l'export échoue explicitement ;
   il ne prétend pas créer un `.pvtu` distribué. En MPI1, conserver le
   comportement d'export historique.

| Réglage | Valeur par défaut du code | Retour arrière |
|---|---|---|
| `SCRIMP_AUTO_INTERFACE` | `1` | `0` désactive la compilation automatique |
| `SCRIMP_PV_EXPORT_MODE` | `auto` (rang 0 pour MPI>1 si données répliquées) | `legacy` restaure les écritures historiques de tous les rangs |
| `SCRIMP_V130_AUTO_INTERFACE_PILOT` | alias historique (si aucune option nouvelle) | valeur explicite `0` ou `1` |
| `SCRIMP_V130_PV_ROOT_ONLY_PILOT` | alias historique (si aucune option nouvelle) | valeur explicite `0` ou `1` |

**Important :** le mode MPI `legacy` reste disponible pour compatibilité,
mais il **réintroduit le risque connu d'écritures concurrentes**.

## Réussites natives prouvées (étapes de qualification)

Les résultats ci-dessous proviennent de véritables exécutions
GetFEM/PETSc/MPI dans l'image figée, et non des tests synthétiques.

| Étape | Résultat et preuve |
|---|---|
| Heat=1, Wave=2, matrice de référence manuelle contre compilation automatique | **PASS**. Matrice 7689×7689, 46935 coefficients CSR stockés ; parité exacte des blocs et propriétaire de chacune des 78 facettes. |
| Heat=1, Wave=2, dynamique MPI1 | **PASS**. Deux initialisations compatibles (`clamp`, `constrained_l2`), 4 pas CN ; trajectoires et bilans de puissance identiques au cas manuel, défauts énergétiques ≈ 1e-15. Archive `v130-auto-interface-CN-b001f2eb-7wIwM8qq.zip`. |
| Heat=2, Wave=1, changement de convention de normales | **PASS**. K 8056×8056 ; exactement **390 coefficients changent de signe**, dans `e_Q←U_T` (156) et `U_w←Y_T` (234), **aucune autre différence**. Archive `v130-swapped-interface-sign-K-1c08c78b-0cTy2vOb.zip`. |
| Heat=2, Wave=1, DAE/CN et énergie MPI1 | **PASS**. Résidu DAE 9.758748179398196e-14, 4 pas acceptés, défaut énergétique total 1.6792123247455493e-15, puissance d'interface non nulle 4.390457971145539. Archive `v130-swapped-interface-CN-3cbc738a-djMlJFrk.zip`. |
| Assemblage MPI1 contre MPI2, les deux orientations | **PASS**. Matrices CSR canoniques **strictement identiques** : 41445 coefficients non nuls pour 7689×7689, 43189 pour 8056×8056. Répartitions PETSc correctes. Archive `v130-auto-interface-MPI-K-d00052f5-AyVWkC3j.zip`. |
| DAE stricte, intégration CN et bilan énergétique, MPI1 contre MPI2 | **PASS**. 4/4 pas acceptés dans les quatre simulations ; zéro échec SNES ; les états complets des deux rangs MPI2 sont identiques bit à bit et proches de MPI1 à 2.2e-12 / 1.45e-11. Archive `v130-auto-interface-MPI-CN-c22e2748-MTltxsjk.zip`. |
| Exports indépendants par rang MPI | **PASS**. PVD/VTU réels lisibles par VTK, cinq instants, chacun des deux rangs produit le **maillage complet**, contenu identique octet par octet à l'autre rang. `p` : 1641 points / 3154 cellules ; `T` : 9462 points / 3154 cellules. Archive `v130-auto-interface-MPI-CN-4bb09edc-yY0sPh1V.zip`. |
| Export partagé **rang 0 seul**, plus physique MPI1/MPI2 | **PASS**. Vraies séries VTU/PVD, géométrie/connectivité identiques à MPI1 et valeurs Float32 conformes, sans répertoire parasite de rang 1. DAE/CN/puissances inchangées ; **54/54 tests de qualification source**. Archive `v130-auto-interface-MPI-CN-03852372-55Ii0zEg.zip`. |

Détails des derniers bilans MPI2 : pour heat=1/wave=2, résidu DAE
`9.782179042671012e-14`, défaut énergétique total
`9.610368056911511e-16`, puissance d'interface maximale
`3.5690759218838686` ; pour heat=2/wave=1, respectivement
`9.853778732954403e-14`, `2.6853519408120974e-15` et
`4.3904579711455245`.

Les huit instantanés publiés de `examples/published/` sont demeurés
**strictement inchangés** pendant les campagnes de validation.

## Échecs rencontrés, diagnostic et résolution

**Les échecs ci-dessous sont consignés sans les convertir en PASS.
Les campagnes corrigées figurent dans le tableau précédent.**

| Incident | Impact exact | Correction / état |
|---|---|---|
| Mesure initiale des puissances sur les facettes codées `30/31` | Collecteur incorrect en mode automatique `21/22` ; **pas** un défaut d'assemblage GetFEM. | Résolution des régions effectives depuis les ports ; test de régression `test_v130_hw_facet_collector_regions.py`. **Corrigé, confirmé.** |
| Conventions de signes héritées de la normale commune en orientation inversée | Seuls deux blocs d'interface doivent changer de signe avec les normales extérieures mathématiques. | Mise à jour autorisée de `examples/heat_wave.py` ; preuve K native des 390 entrées et bilans DAE/CN physiques. **Corrigé, confirmé.** |
| `scipy.sparse.vstack(generator)` | `TypeError: iteration over a 0-d array` **dans le comparateur hôte**, après quatre assemblages MPI réels achevés. | Passage à une liste explicite. **Corrigé, confirmé.** |
| Nombre brut de coefficients CSR différent entre MPI1/MPI2 | Faux échec : MPI2 stocke 18/24 **zéros explicites supplémentaires** ; aucune différence numérique. | Comparaison canonique après `sum_duplicates`, `eliminate_zeros`, tri ; vrais défauts toujours rejetés. **Corrigé, confirmé.** |
| Second calcul du Hamiltonien dans le pilote MPI | Historiques de 10 valeurs au lieu de 5 ; `ValueError` à l'addition ; la vraie DAE/CN avait réussi. | Réutilisation de l'Hamiltonien déjà calculé par l'exemple ; longueur vérifiée. **Corrigé, confirmé.** |
| `host_error` non initialisé | `NameError` dans le comparateur et absence de résumé d'échec. | Initialisation avant le `try` ; archivage succès et erreurs testé. **Corrigé, confirmé.** |
| Sous-processus de test du comparateur terminé par `SIGKILL (-9)` | Échec **avant Docker** ; cause native du signal non établie. | Exécution du véritable comparateur dans le processus de test ; 7 puis 9 tests réussis. **Contournement confirmé, cause externe inconnue.** |
| Modification du libellé d'échec avec l'option ParaView | Ancien test attendait `FAIL_..._PHYSICS` au lieu de `FAIL_..._OR_PV_DIAGNOSTIC`. | Conserver l'ancien libellé hors mode PV ; nouveaux tests de la branche PV. **Corrigé, confirmé.** |
| Écriture partagée ancienne par tous les rangs | **Risque de concurrence** : mêmes noms de PVD/VTU ouverts par tous les processus ; pas un échec systématique prouvé dans les logs. | Un seul écrivain natif MPI2 testé, avec synchronisation et propagation d'erreur ; maintenant comportement par défaut **sur la branche de release**, sous condition de répliques complètes. |
| Régressions de compatibilité lors de l'activation par défaut | Les 54 tests source précédents validaient le **mode pilote**, pas les nouveaux chemins sans drapeaux. | Nouvelles régressions de politique + validateur réel via AST ; **71 tests et validation native sans drapeaux PASS** à `241f58fb`. |

## Validation finale — vrais chemins par défaut : PASS

**Exécution native du 24 novembre 2026, code source exactement
`241f58fb2097a16721ca1592211bd3248c88c88e`.**
Commande utilisée :
`bash scripts/run-v130-production-defaults.sh` avec image
`scrimp-mpi:v1.3.0-rc3` de révision verrouillée
`1d6d0c90115213298047b28237a7f38087563a35`.
Le lanceur a explicitement désactivé la transmission des anciens et
nouveaux drapeaux de fonctionnalités vers Docker : les comportements
ont donc été empruntés **par défaut** dans le code de la branche
`release/1.3.0-preparation`.

**Archive native vérifiée :**
`v130-auto-interface-MPI-CN-241f58fb-YEJXJhSK.zip`
(74 entrées, rapport `comparison.json`, 6 rapports de rang,
6 trajectoires `npz`, 4 journaux natifs et fichiers VTK/PVD).
Statut final **`PASS_MPI2_ROOT_ONLY_SHARED_VTU_PVD_COMPLETE_CN_ENERGY`** ;
`host_comparison_failure=null`.
**71/71 tests source PASS** (les 54 précédents + les 17 nouveaux
tests des politiques de production et du validateur de réplication).
Les **huit exemples publiés sont immuables**.

| Mesure | Heat=1 / Wave=2, MPI2 | Heat=2 / Wave=1, MPI2 |
|---|---:|---:|
| Dimension du modèle | 7689 | 8056 |
| Résidu algébrique DAE | 9.782179042671012e-14 | 9.853778732954403e-14 |
| Pas CN acceptés | 4/4 | 4/4 |
| Rejets / échecs SNES | 0 / 0 | 0 / 0 |
| Défaut du bilan énergétique total | 9.610368056911511e-16 | 2.6853519408120974e-15 |
| Défaut de fermeture de puissance à l'interface | 8.881784197001252e-16 | 4.440892098500626e-16 |
| Écart maximal des états MPI1/MPI2 | 2.2297719226571644e-12 | 1.4509282664221246e-11 |
| Traces MPI2 rang 0 / rang 1 | identiques bit à bit | identiques bit à bit |
| Export partagé `p` et `T`, cinq instants | PASS | PASS |

Les audits des VTU/PVD natifs confirment cinq instantanés par
champ et configuration, `p` avec **1641 points / 3154 cellules**
et `T` avec **9462 points / 3154 cellules**. La géométrie
et la connectivité binaires MPI1/MPI2 sont identiques.
L'écart absolu maximal des valeurs des champs exportés,
sur **tous les instants et les deux orientations**, est
`2.441726822362823e-17` pour `p` et
`6.462348535570529e-27` pour `T`.
Aucune sortie parasite de rang 1.

**Petite anomalie de compte rendu, sans impact sur la preuve native :**
le `comparison.json` de cette campagne contient encore
`"production_auto_compiler_enabled_by_default": false`,
valeur **codée en dur** et héritée de la période pilote.
Cela contredit le chemin réellement exécuté (aucun drapeau
transmis, compilateur actif et CN/VTK PASS) mais **ne signifie pas
qu'une compilation manuelle ou pilote a été utilisée**.
Le générateur de rapport a été corrigé au commit `942e23e`
pour les campagnes futures. Ce correctif modifie **uniquement
les métadonnées**, non la simulation native ; la présente preuve
reste attachée à son commit source d'origine `241f58fb`.

**Décision :** les deux développements restent **activés par défaut
dans le code de la branche de préparation de la v1.3.0**.
Les risques non qualifiés demeurent explicitement documentés.
**Aucune fusion vers `master`, aucun tag, aucune publication**.

## Limites connues et décisions de livraison

- **Périmètre physique établi :** Heat-Wave 2D concentrique, mesh conforme,
  RT0/DG0 côté chaleur, RT0/CG1 côté onde, deux affectations des régions,
  quatre pas CN, MPI1/2, métriques de puissance mesurées.
- **Non établis :** géométries arbitraires, interfaces partielles de
  rectangles en dynamique, intersections multiples, interfaces non conformes,
  faces courbes, import physique 1D/3D, H(div)/H(curl) générique,
  RT(k) tous ordres, BDM, Nédélec, données ParaView réellement
  partitionnées, séries temps très longues, MPI>2 pour ce couple.
- Le compilateur **échoue explicitement** lorsqu'il ne sait pas affecter
  une interface ; il ne doit jamais masquer un cas non supporté.
- L'export partagé n'est **pas** un export `.pvtu` multi-parties.
  Les rangs doivent posséder la même trajectoire complète et la même
  connectivité avant le premier fichier. Les cas réellement partitionnés
  sont rejetés au lieu d'écrire des résultats partiels.
- Le validateur collectif de répliques et l'activation **sans aucun
  drapeau expérimental** sont à présent qualifiés sur les deux
  orientations Heat-Wave en MPI 1/2 : archive native finale
  `v130-auto-interface-MPI-CN-241f58fb-YEJXJhSK.zip`, ci-dessous.
  Cette preuve ne généralise pas le comportement aux maillages non répliqués.
- Cette décision de livraison ne crée **ni tag v1.3.0, ni image Docker
  définitive, ni publication**, et ne fusionne PAS les nombreux
  autres développements de la branche avec `master`.

## Validation minimale des nouveaux chemins par défaut

```bash
cd ~/Simulations/scrimp
git pull --ff-only
conda activate scrimp
python tests/test_v130_interface_pv_production_policy.py

# Deux orientations × MPI1/MPI2 avec les chemins de code par défaut :
# aucun ancien SCRIMP_V130_*_PILOT n'est transmis au conteneur.
SCRIMP_DOCKER_NO_SUDO=1 \
SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
bash scripts/run-v130-production-defaults.sh
```

**Cette commande a été exécutée et a réussi au commit
`241f58fb` ; la présente section en conserve la procédure pour
reproductibilité.** Elle n'effectue ni publication ni fusion de
branches. Ce PASS valide les chemins par défaut **dans le périmètre
Heat-Wave décrit**, pas une généralité encore non testée.
