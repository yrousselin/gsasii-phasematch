#!/bin/bash
# Installe / met à jour la recherche de phases (search-match) dans GSAS-II -- macOS/Linux
set -e
HERE="$(cd "$(dirname "$0")" && pwd)"
DEST="$HOME/.GSASII/exports"
mkdir -p "$DEST"
cp "$HERE"/G2export_PhaseMatch.py "$HERE"/pm_*.py "$HERE"/csd_pm_worker.py "$DEST/"
echo "Fichiers copiés dans $DEST"
PY="$HOME/g2main/bin/python"
[ -x "$PY" ] || PY="$(command -v python3)"
echo "Installation du module gemmi (symétrie des fichiers COD sans opérateurs)…"
"$PY" -m pip install --quiet gemmi || echo "gemmi non installé : quelques pour cent des entrées COD seront laissées de côté"
cat <<MSG

Pour ajouter la COD à l'index (une fois) :
  $PY $DEST/pm_index.py cod-sync          # copie locale des CIF dans ~/COD/cif (long la première fois)
  $PY $DEST/pm_index.py build --cod       # calcul des raies, reprenable
Mises à jour ultérieures : relancer les deux mêmes commandes.
Redémarrez GSAS-II pour la fenêtre de recherche.
MSG
