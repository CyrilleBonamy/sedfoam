#!/usr/bin/env bash
# fix_swak4foam_openfoam2606.sh
# Corrige les constructeurs de copie de groovyBC devenus incompatibles
# avec le nouveau schéma fvPatchField (constructeur de copie "brut" supprimé)
# introduit dans les versions récentes d'OpenFOAM ESI (>= ~v2306, testé v2606).
#
# Usage: ./fix_swak4foam_openfoam2606.sh /chemin/vers/swak4Foam

set -euo pipefail

SWAK_SRC="${1:?Usage: $0 /chemin/vers/swak4Foam}"

if [ ! -d "$SWAK_SRC" ]; then
    echo "Répertoire introuvable : $SWAK_SRC" >&2
    exit 1
fi

# Motif : initialisation d'une classe de base par copie brute, du type
#   NomDeClasse<Type>(ptf),
# éventuellement suivi d'une virgule, sur sa propre ligne.
PATTERN='^([[:space:]]*)([A-Za-z_][A-Za-z0-9_]*)<Type>\(ptf\)([[:space:]]*,?[[:space:]]*)$'

FILES=$(grep -lrE "$PATTERN" --include='*.C' "$SWAK_SRC" || true)

if [ -z "$FILES" ]; then
    echo "Aucun fichier candidat trouvé dans $SWAK_SRC."
    exit 0
fi

echo "== Fichiers concernés =="
echo "$FILES"
echo

for f in $FILES; do
    cp "$f" "$f.bak"
    # On patche TOUTES les classes de base <Type>(ptf) SAUF groovyBCCommon,
    # qui a un vrai constructeur de copie (const groovyBCCommon<Type>&) et
    # doit rester appelé tel quel : groovyBCCommon<Type>(ptf)
    sed -i -E \
        "/groovyBCCommon/!s/${PATTERN}/\1\2<Type>(ptf, ptf.internalField())\3/" \
        "$f"
done

echo "== Diffs appliqués =="
for f in $FILES; do
    if ! diff -q "$f.bak" "$f" > /dev/null; then
        echo "--- $f ---"
        diff -u "$f.bak" "$f" || true
        echo
    fi
done

echo "== Vérification : aucune ligne groovyBCCommon ne doit avoir été modifiée =="
if grep -rn 'groovyBCCommon<Type>(ptf, ptf\.internalField())' "$SWAK_SRC"; then
    echo "ATTENTION : des lignes groovyBCCommon ont été modifiées par erreur, vérifiez ci-dessus." >&2
else
    echo "OK — groovyBCCommon<Type>(ptf) intact partout."
fi

echo
echo "Pour tout annuler :"
echo "  find \"$SWAK_SRC\" -name '*.bak' -exec bash -c 'mv \"\$0\" \"\${0%.bak}\"' {} \;"
