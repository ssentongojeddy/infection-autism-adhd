#!/usr/bin/env bash
# Reproduce every derived statistic, table and figure in the paper.
#
# The analysis scripts read their inputs by bare filename, so this script
# assembles a build directory with the code and data side by side, runs
# everything there, and copies the outputs back into figures/.
#
# Usage:  bash run_all.sh

set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BUILD="$ROOT/build"

echo "==> Preparing build directory"
rm -rf "$BUILD"
mkdir -p "$BUILD"
cp "$ROOT"/code/*.py "$BUILD"/
cp "$ROOT"/data/*.json "$ROOT"/data/*.txt "$ROOT"/data/*.csv "$BUILD"/
cd "$BUILD"

echo "==> Deriving statistics (P values, FDR, absolute risks, reference adjustment, pooling, E-values)"
python stats.py

echo "==> Building Table 1 and the supplementary tables"
python table1_cmv.py
python gen_tables.py

echo "==> Building figures"
for f in fig1.py fig2.py fig3.py fig4.py efigs.py; do
  echo "    $f"
  python "$f"
done

echo "==> Collecting outputs"
mkdir -p "$ROOT/figures"
find "$BUILD" -maxdepth 1 -name '*.pdf' -exec cp {} "$ROOT/figures/" \;
find "$BUILD" -maxdepth 1 -name '*.png' -exec cp {} "$ROOT/figures/" \; 2>/dev/null || true

echo
echo "Done. Figures are in figures/. Intermediate files remain in build/."
