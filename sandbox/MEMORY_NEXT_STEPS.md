# SCRIMP — mémoire : état et pistes à reprendre

**État : chantier mis en pause le 17 novembre 2026.**
Référence de la dernière campagne native : `6ac3c007`,
`compact-matrix-6ac3c007-pSIPxG`, 2 rangs MPI, trois pas
Crank–Nicolson par mode et par modèle. Ne pas relancer de simulation
longue pour cette étude avant d'avoir choisi une cible mesurable.

## Acquis vérifiés

1. Les copies intégrales de trajectoire à chaque rang ont été
   remplacées, **en option**, par un processus externe d'E/S asynchrones
   écrivant les portions PETSc locales en binaire. Le mode `async_binary`
   conserve les historiques temporels, les états Newton/GetFEM
   transitoires et les interfaces de post-traitement par lecture
   différée. Validation acoustique linéaire et Saint-Venant 2D
   non linéaire **sur cas courts** ; format binaire v2 + CRC32 validé
   par test natif `1bcd9a49`. Le défaut reste `memory`.
2. Consolidation finale `_select_saved_solution_through_final_time` :
   réduction mesurée de **792,992 à 701,211 Mio** en mode sauvegardes
   fréquentes (47 942 DDL, 250 pas, 2 rangs), soit **−91,781 Mio
   par rang (−11,57 %) ;** aucune différence physique significative.
3. `SCRIMP_MATRIX_MEMORY_MODE=compact` (opt-in ; défaut `legacy`) :
   partage des tangentes linéaires avec leurs matrices assemblées
   et un seul `Mat` + un seul `Vec` persistants pour la correction
   Hessienne de masse non linéaire. Parité numérique complète
   acoustique/Saint-Venant validée sur trois pas. Sur Saint-Venant,
   **une** allocation de matrice Hessienne et **quatre** réutilisations
   par rang. Gain total en mémoire faible à ce stade.

## Dernier audit natif : phases mémoire (6ac3c007)

Mémoire maximale par rang (`VmHWM`, Mio) :

| Cas | legacy | compact | Écart compact − legacy |
|---|---:|---:|---:|
| Acoustique, h=0.04, 47 942 DDL | 508.73 | 506.30 | **−2.43** |
| Saint-Venant dissipatif, h=0.35, 1 170 DDL | 389.48 | 390.02 | **+0.55** |

Le +0,55 Mio de Saint-Venant sur une exécution n'établit
**ni fuite, ni régression significative**.

| Phase de la résolution (maximum RSS des rangs, Mio) | Acoustique legacy | Acoustique compact |
|---|---:|---:|
| Entrée `solve` | 376.94 | 378.36 |
| Après allocation des matrices PETSc | 388.21 | 387.46 |
| Après assemblage masse | 411.64 | 410.29 |
| Après assemblage raideur | 462.86 | 460.80 |
| Après RHS | 463.85 | 461.65 |
| Avant TS principal | 475.95 | 475.57 |
| Premier pas accepté | 502.09 | 496.77 |
| Après TS | 508.60 | 506.11 |
| Sortie `solve` | 508.73 | 506.30 |

Les **assemblages GetFEM** sont associés à une hausse RSS
d'environ **75,63 / 74,19 Mio** (legacy / compact).
Le segment début TS → premier pas accepté est associé à
**+26,14 / +21,21 Mio**, incluant les premières opérations
PETSc/SNES/KSP et la factorisation MUMPS. Ces mesures **ne
permettent pas d'attribuer chaque octet** à une bibliothèque.

Saint-Venant : assemblages **+2,79 / +2,49 Mio** et premier
pas TS **+1,19 / +1,98 Mio**. L'essentiel du coût (~370 Mio
dès l'entrée `solve`) est déjà résident avant ces phases.

## Pistes prioritaires pour reprise

1. **Distinguer mémoire du processus de base, GetFEM, PETSc et MUMPS.**
   Capturer RSS anonyme/fichier et, lorsque disponible, allocation
   PETSc/MUMPS réelle (facteurs symboliques/numériques).
   Mesurer également les phases **avant** `DPHS.solve()` :
   import Python/NumPy/gmsh/PETSc, construction maillage/FEM,
   construction du modèle GetFEM. Une RSS `RssFile` de ~190 Mio
   par rang n'est pas nécessairement de la mémoire privée
   additionnable entre rangs ; enregistrer PSS/USS via
   `/proc/self/smaps_rollup` lors d'un audit léger.
2. **Déterminer la taille et la durée de vie des copies temporaires
   GetFEM → PETSc** malgré le transfert COO, spécialement les
   matrices GMM/CSR intermédiaires de masse/raideur, leur conversion
   et la préallocation COO. Identifier si les matrices GetFEM
   restent utiles après l'assemblage initial.
3. **PETSc/MUMPS :** mesurer la mémoire des facteurs,
   l'analyse symbolique et l'allocation numérique à la première
   factorisation, puis comparer `never` (opérateurs linéaires,
   structure fixe) et `on_change` (Jacobiennes non linéaires).
   Ne pas optimiser par `PC.reset` systématique : son coût
   temporel a déjà été démontré. Toute mesure MUMPS doit être
   faite avec mêmes rangs, maillage, politique et matrice.
4. **Répliquer moins de données globales transitoires** si les
   bibliothèques le permettent : `GlobalVec` est encore nécessaire
   à la synchronisation GetFEM du **Newton non linéaire**, mais
   les rassemblements globaux faits uniquement pour le suivi
   (mode linéaire, `monitor`) sont candidats à suppression.
   Ne jamais supprimer les vecteurs de l'intégrateur PETSc,
   le trial-state Newton, les Jacobiennes ou la mémoire requise
   par les non-linéarités.
5. **Nettoyage et recyclage des matrices de travail** :
   valider propriétés PETSc/Mat et durée de vie des
   `J`, matrices tangent, `nl_mass`, `nl_stiffness`,
   et buffers ; éviter toute libération prématurée ou
   destruction d'objets référencés par TS/KSP/MUMPS.
6. **Évaluation longue durée seulement si nécessaire**, après
   identification d'un indice de croissance inexpliquée.
   Le précédent test acoustique à sauvegardes rares (250 pas)
   atteignait un plateau RSS : aucune fuite par pas n'y était
   démontrée. L'effet du disque/queue de `async_binary`
   sur l'empreinte mémoire et le débit reste à caractériser
   séparément.

## Exigences de validation

- **Une seule dimension modifiée à la fois** : maillage, temps,
  méthode COO/CSR, politique MUMPS, nombre de rangs et backend
  résultats doivent être documentés et constants par paire.
- Parité des **vecteurs complets**, étapes acceptées/rejetées,
  convergence PETSc, Hamiltoniens, puissances et bilans de
  volume/masse pour les systèmes concernés.
- Pas de conclusion sur le gain à partir de quelques
  dixièmes de Mio ou d'une seule mesure bruitée.
- Conserver les benchmarks de performance de **v1.2.0**
  séparés : cette release ne possède pas les mêmes options COO
  et politiques MUMPS ; les comparaisons doivent être datées
  et reproduites sur la même machine et le même exemple.

À reprendre après les figures de performance de la version actuelle.
