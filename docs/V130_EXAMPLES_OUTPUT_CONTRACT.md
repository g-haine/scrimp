# SCRIMP v1.3.0 — sorties des exemples et API Hamiltonien / puissances

**Statut au 27 novembre 2026 : modifications source intégrées sur
`release/1.3.0-preparation`, qualification native des nouveaux chemins
**EN ATTENTE**.** Les campagnes déjà archivées conservent leurs anciens
noms de fichiers et leur valeur scientifique. Ne jamais renommer ni
écraser leurs archives historiques.

## 1. Une structure prédictible pour chaque système

Lancer depuis le répertoire de travail de l'utilisateur :

```text
working_dir/
└── outputs/
    └── system_name/
        ├── png/               # figures matricielles réelles
        ├── svg/               # mêmes figures, vectorielles
        ├── gif/               # animations des vrais états sauvegardés
        ├── mesh/              # maillage Gmsh .msh et champs GetFEM
        │   ├── Rectangle.msh  # si le domaine Gmsh en génère un
        │   ├── q/q.pvd        # ParaView, avec VTU d'origine
        │   └── p/p.pvd
        ├── log/               # scrimp.log sur le rang MPI 0
        ├── npz/               # matrices/traces physiques sauvegardées
        ├── csv/               # séries tabulaires si disponibles
        └── json/              # provenance et vérification des exports
```

Les sous-répertoires sont créés à la demande. Le fichier historique
`outputs/mesh/<geometry>.msh` est **toujours utilisé pour l'import
GetFEM** ; une copie authentique est conservée dans
`outputs/<system_name>/mesh/` après construction du `Domain`
lorsque la géométrie l'a généré. Le maillage n'est ni regénéré
ni altéré pour la présentation. Certains modèles concentrés ne
génèrent naturellement aucun maillage Gmsh.

Les exemples de pendule forcé/libre ont deux racines distinctes :
`nonlinear_pendulum_free` et `nonlinear_pendulum_forced`.
Un utilisateur peut changer la racine avec
`DPHS.set_output_directory(path)` ou nommer son problème :

```python
system = S.DPHS("real")
system.set_output_name("my_system")
# ou, pour choisir explicitement le dossier racine :
system.set_output_name("my_system", output_dir="/tmp/mes_resultats")
```

Un appel `export_to_pv("p")` écrit alors
`<root>/mesh/p/p.pvd` et ses véritables VTU.
L'API historique `export_to_pv("p", path=...)`
et les `DPHS` sans nom de sortie conservent
leurs anciennes destinations par compatibilité.

Les **dossiers de preuve** `sandbox/outputs/<pilot>/work_<case>/`
restent séparés de cette convention : chaque pilote exécute
de vraies simulations dans des espaces isolés et doit conserver
leurs données rank-local/rang0 pour les auditeurs.
Ne pas confondre ces métadonnées de test avec les figures d'un
programme utilisateur.

## 2. Méthodes publiques : cas physique et honnêteté numérique

`DPHS.compute_Hamiltonian()` évalue chaque terme de
`Hamiltonian` sur les *vrais* états GetFEM sauvegardés
avec le modèle au **temps correspondant**.
Un recalcul remplace le cache antérieur
et ne double pas le nombre de valeurs.
Tout `DPHS.solve()` invalide son cache et
celui des puissances précédemment calculées.

`DPHS.plot_Hamiltonian(with_powers=True,save_figure=True,
filename="Hamiltonian.png",show=False)` :

- assemble collectivement chaque énergie/power originale GetFEM ;
- dessine uniquement sur le rang MPI 0 ;
- exporte l'image PNG **et** SVG et une copie numérique
  `npz/<stem>.npz` comprenant les vrais temps/termes
  et, si demandé, les puissances avec leurs conventions ;
- garde le `show=True` historique pour les usages interactifs,
  mais les exemples CLI utilisent `show=False` ;
- ne prétend pas que tous les systèmes ou tous les
  schémas PETSc satisfont une identité énergétique exacte.

`DPHS.compute_powers()` lit le **schéma PETSc déclaré**
et la **longueur réelle de l'historique de pas acceptés**.
En CN, `cn_step` n'est autorisé que si les
**deux extrémités de chaque pas accepté** ont été sauvegardées
(`dt_save=dt` en pas fixe, ou même nombre exact de
snapshots et de pas pour le cas adaptatif).
À défaut, les puissances sont évaluées aux
**véritables états nodaux sauvegardés** avec
`sampling="instantaneous"`.
Le `plot_powers` associé les intègre alors par
**trapèzes approximatifs**, avec l'étiquette adéquate.
Pour BDF2/ARKIMEX/autres, cette convention nodale
est également approximative : aucune quadrature de
stade caché n'est inventée. Le mode CN conserve
la convention GetFEM originale `t[k]`,
`0.5*(z[k-1]+z[k])` pour `k>=1`.
Le terme initial `P[0]` **n'est pas intégré**
comme une puissance de pas CN.

```python
system.compute_Hamiltonian()
H = system.get_Hamiltonian()
system.compute_powers()
system.plot_Hamiltonian(
    with_powers=True,
    save_figure=True,
    filename="Hamiltonian.png",
    show=False,
)
```

**Prudence pour un Hamiltonien non quadratique :**
PETSc TSCN aux extrémités ne donne pas en général
une identité de chaîne exacte sur `H(z)` sans
un calcul énergétique spécifique au modèle.
Les tracés génériques sont des observations,
pas une certification automatique de passivité.

## 3. Ce que créent les exemples eux-mêmes

| Exemple | Physique & sorties utilisateur |
| --- | --- |
| `mass_spring_damper.py` | BDF2 concentré original ; états, énergie, ports, PNG/SVG, GIF de schéma, NPZ. |
| `nonlinear_pendulum.py` | Vrai `sin(theta)`, CN ; variante libre/forcée, PNG/SVG, GIF et données avec défaut de chaîne non quadratique documenté. |
| `heat.py` | Température réelle GetFEM/VTP, Lyapunov, puissance/dissipation et GIF de diffusion (CN). |
| `heat_wave.py` | Température/onde en PVD/VTU, énergie chaleur/onde et totale, **deux puissances orientées et dissipation**, cumul `H-H0` contre `W-D`, PNG/SVG/NPZ. GIF annulaire en option ; n'est pas une preuve de mode propre. |
| `shallow_water.py` | Saint-Venant inviscide, hauteur et vitesse originales CG2/DG1, énergie, masse, PNG/SVG/GIF, avec limite BDF2 explicitée. |
| `wave.py` | Onde 2D anisotrope conservative 2A : Hamiltonien, ports, champs originaux `q,p,e_q,e_p` via PVD/VTU. |
| `wave_dissipative.py` | 2B : ajout de la perte résistive, Hamiltonien, champs `q,p,e_q,e_p,f_r,e_r`. |
| `wave_coenergy.py` | 2C : coénergies, amortissement interne/feedback, Hamiltonien, champs `q,p,e_r`. |

Le mode `physics_only` de Heat–Wave reste réservé
aux campagnes de qualification et ne produit
pas de visualisations, sans changer les équations.

## 4. Contrat de test, refus de « PASS » prématuré

Nouveaux tests purement Python :
`tests/test_v130_user_output_power_hamiltonian.py`,
`tests/test_v130_examples_output_layout_source.py`.

Le pilote `scripts/run-v130-example-user-output-smoke.sh`
exécute neuf cas **réels** sous la même image GetFEM/PETSc
que les qualifications antérieures : huit exemples
en MPI1 et un Heat–Wave MPI2 supplémentaire,
chaque exemple via son **point d'entrée Python public**.
Il vérifie le contenu PNG/SVG/NPZ, les séries PVD et
chaque VTU annoncée, les horodatages ParaView,
les logs et l'absence de fichiers fantômes.
Une archive ZIP conservant les résultats est produite
en cas de réussite complète ; les dossiers sont conservés
pour diagnostic si une exécution échoue.

```bash
cd ~/Simulations/scrimp
git pull --ff-only
SCRIMP_DOCKER_NO_SUDO=1 \
SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
bash scripts/run-v130-example-user-output-smoke.sh
```

**Statut actuel : SOURCE STAGED / NATIVE PENDING.**
L'ancienne physique numériquement acceptée n'est pas
remise en cause par les nouvelles sorties, mais
les nouveaux chemins, conventions de power sampling
et figures doivent être acceptés avec cette campagne
avant fermeture du volet « exemples v1.3.0 ».
