# Recherche de phases (search-match) pour GSAS-II

*[English version: README.md](README.md)*

Extension pour **GSAS-II** qui aide à identifier les **phases cristallines présentes dans un diagramme de diffraction sur poudre**, y compris dans un mélange. Elle compare votre diagramme à un index local de structures de référence, construit à partir de la CSD, de la COD et/ou de vos propres CIF. Les meilleurs candidats s'affichent superposés au diagramme. **C'est vous** qui décidez quelle phase retenir ; la recherche suivante porte alors sur ce qui reste à expliquer.

## Ce que fait l'extension

- Ajoute à GSAS-II le menu **Export → Entire project as → (user) Phase search (search-match)**.
- Ouvre une fenêtre interactive qui affiche :
  - le diagramme mesuré, fond soustrait ;
  - les phases retenues ;
  - le candidat sélectionné ;
  - les positions de raies ;
  - ce qui reste à expliquer.
- **Comment se déroule une recherche.** D'abord, un tri grossier de tout l'index vérifie si les raies fortes de chaque référence sont présentes dans le diagramme, en autorisant une dilatation isotrope de la maille. Ensuite, un tri fin ajuste la maille sur les pics observés. Enfin, les meilleurs candidats sont validés avec la structure complète : maille, échelle et largeur des pics sont ajustées sur le profil.
- Quand vous cliquez sur **Retenir cette phase**, toutes les phases retenues sont ajustées ensemble (facteurs d'échelle positifs), puis la recherche repart sur le résidu.
- Filtre facultatif par **éléments exigés**, par exemple `C H N O`, avec un mode strict (**aucun autre**).
- **Voir dans Mercury** ouvre un candidat dans Mercury (macOS).
- Les phases retenues sont exportées, avec un résumé enregistré dans `<projet>_phases.json`.

## Prérequis

| Quoi | Pourquoi | Où le trouver |
|---|---|---|
| **GSAS-II** | accueille l'extension ; fournit numpy, scipy, matplotlib et wxPython | [GSAS-II](https://advancedphotonsource.github.io/GSAS-II-tutorials/install.html) |
| **Un index de références** (construit une fois, voir plus bas) | les structures auxquelles comparer | CSD, COD et/ou vos dossiers de CIF |
| Logiciels du CCDC avec la CSD Python API *(facultatif)* | ajout de la CSD à l'index | [CCDC](https://www.ccdc.cam.ac.uk/) (licence CSD nécessaire) |
| gemmi *(facultatif)* | symétrie des fichiers COD sans opérateurs de symétrie | `pip install gemmi` |
| rsync *(facultatif)* | téléchargement et mise à jour rapides de la COD ; présent d'office sur macOS | sinon, utilisez `--archive` |

Seules les données de poudre à **longueur d'onde constante** (rayons X ou neutrons) sont prises en charge.

## Installation

### Étape 1 : copier les fichiers dans GSAS-II

**macOS (et Linux)**
1. Téléchargez ce dépôt : bouton vert **Code** → **Download ZIP**, puis double-cliquez sur le ZIP pour le décompresser.
2. Ouvrez le **Terminal**. Tapez `bash ` (avec un espace après), faites glisser `install_phasematch_plugin.sh` dans la fenêtre, puis appuyez sur Entrée. Les fichiers sont copiés dans `~/.GSASII/exports/` et gemmi est installé.

**Windows**
1. Téléchargez le ZIP comme ci-dessus, puis clic droit dessus → **Extraire tout**.
2. Double-cliquez sur **`install_windows.bat`**. Les fichiers sont copiés dans `C:\Users\<votre nom>\.GSASII\exports\`.
3. Facultatif : installez gemmi avec le Python de GSAS-II :
   `"<dossier de GSAS-II>\python.exe" -m pip install gemmi`

### Étape 2 : construire l'index de références (une fois)

L'index est un fichier unique, `~/.GSASII/pm_index.sqlite`. Il ne contient que les raies les plus fortes et une fiche descriptive de chaque structure. La structure complète est relue dans sa source quand c'est nécessaire.

Tapez ces commandes dans un terminal, avec le Python de GSAS-II. Ci-dessous, `PY` désigne ce Python, par exemple `~/g2main/bin/python` sur macOS, et `EXP` le dossier `~/.GSASII/exports` (`%USERPROFILE%\.GSASII\exports` sous Windows). Choisissez une ou plusieurs sources :

```
# Vos propres CIF (n'importe quel dossier, par exemple des exports ICSD)
PY EXP/pm_index.py build --cif-dir /chemin/vers/mes_CIF --source MesCIF

# La CSD (logiciels du CCDC nécessaires)
PY EXP/pm_index.py build --csd

# La COD : d'abord une copie locale des CIF (long la première fois), puis l'index
PY EXP/pm_index.py cod-sync              # macOS / Linux (rsync)
PY EXP/pm_index.py cod-sync --archive    # Windows, ou si rsync est absent
PY EXP/pm_index.py build --cod

# Voir ce que contient l'index
PY EXP/pm_index.py stats
```

`build` est **incrémental**. Vous pouvez l'interrompre avec Ctrl-C et le relancer plus tard. Après une mise à jour de la CSD ou de la COD, relancez les mêmes commandes : seules les nouvelles entrées sont ajoutées.

Si le Python du CCDC n'est pas trouvé, ajoutez `--csd-python` suivi de son chemin. Emplacements habituels :
- **Windows :** `"C:\Program Files\CCDC\ccdc-software\csd-python-api\miniconda\python.exe"` (pour le trouver : `dir /s /b "C:\Program Files\CCDC\python.exe"`)
- **macOS :** `~/CCDC/ccdc-software/csd-python-api/miniconda/bin/python`

### Étape 3 : redémarrer GSAS-II

## Utilisation

1. Ouvrez un projet GSAS-II contenant un diagramme de poudre.
2. Ouvrez **Export → Entire project as → (user) Phase search (search-match)**.
3. Examinez les candidats et leur superposition. Cliquez sur **Retenir cette phase** pour une phase qui correspond. La recherche continue alors sur le signal restant.
4. Recommencez jusqu'à ce que tout soit expliqué, puis cliquez sur **Exporter les phases retenues…**.

**En ligne de commande (sans la fenêtre) :**

```
PY EXP/pm_search.py projet.gpx --elements "C H N O" --max-phases 4
PY EXP/pm_search.py --help
```

## Limites

- Données à longueur d'onde constante uniquement (pas de temps de vol ni de dispersion en énergie).
- Le résultat est une **proposition**. La décision finale, et la confirmation par un affinement Rietveld, reviennent à l'utilisateur.
- La qualité de la recherche dépend de ce que contient l'index.
- Développé et utilisé sur macOS. Windows devrait fonctionner, mais n'a pas encore été testé complètement. Vos retours sont les bienvenus.

## Crédits et références

- **GSAS-II :** B. H. Toby, R. B. Von Dreele, *J. Appl. Cryst.* **2013**, 46, 544–549. Les coefficients de facteurs de diffusion X de `pm_f0.py` sont repris de GSAS-II.
- **COD :** Crystallography Open Database, <https://www.crystallography.net/>. **CSD :** Cambridge Crystallographic Data Centre. **gemmi :** <https://gemmi.readthedocs.io/>.
- Écrit par Yoann Rousselin avec l'aide de Claude (Anthropic) comme assistant de programmation.

## Licence

MIT. Voir [LICENSE](LICENSE).

## Auteur

**Yoann Rousselin**, ingénieur de recherche CNRS (cristallographie), ICMUB, Dijon.
Questions ou signalement d'un problème : ouvrez une *Issue* sur cette page.
