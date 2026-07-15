#!/usr/bin/env bash
# paper/build.sh -- render both SJ and arXiv versions of the stpois paper.
#
# Prerequisites (Stata outputs already generated):
#   do paper/sj/sj_examples.do        -- regenerate paper/sj/ex_*.log.tex
#   do paper/benchmark.do             -- regenerate benchmark_results.csv
#   do paper/benchmark_hdfe.do        -- regenerate benchmark_hdfe_results.csv
#   do paper/benchmark_ppmlhdfe.do    -- regenerate benchmark_ppmlhdfe_results.csv
#
# Usage (from repo root):
#   bash paper/build.sh [--sj-only | --arxiv-only]

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
SJ_DIR="$REPO_ROOT/paper/sj"
ARXIV_DIR="$REPO_ROOT/paper/arxiv"
BUILD_SJ=1
BUILD_ARXIV=1

case "${1:-}" in
  --sj-only)    BUILD_ARXIV=0 ;;
  --arxiv-only) BUILD_SJ=0    ;;
esac

mkdir -p "$ARXIV_DIR"

# ── SJ version ────────────────────────────────────────────────────────────
if [ "$BUILD_SJ" -eq 1 ]; then
  echo "=== Rendering SJ version (stpois_sj.tex) ==="
  cd "$SJ_DIR"
  quarto render stpois_sj.qmd

  # pandoc turns markdown "~" (intended as a nonbreaking space, e.g. before
  # \ref or in "Stata~14") into a printed \textasciitilde{}; restore the
  # nonbreaking space.
  perl -i -pe 's/\\textasciitilde(\{\})?/~/g' stpois_sj.tex

  echo "=== Compiling SJ PDF (main.pdf) ==="
  latexmk -pdf -quiet -interaction=nonstopmode main.tex
fi

# ── arXiv version ─────────────────────────────────────────────────────────
if [ "$BUILD_ARXIV" -eq 1 ]; then
  echo "=== Rendering arXiv version (stpois_arxiv.tex) ==="
  cd "$SJ_DIR"
  quarto render stpois_sj.qmd \
    --to latex \
    --template "$ARXIV_DIR/arxiv-template.tex" \
    --output stpois_arxiv.tex

  # Same nonbreaking-space fix as the SJ build (see above).
  perl -i -pe 's/\\textasciitilde(\{\})?/~/g' stpois_arxiv.tex

  echo "=== Copying supporting files to paper/arxiv/ ==="
  cp stpois_arxiv.tex   "$ARXIV_DIR/"
  cp stpois.bib         "$ARXIV_DIR/"
  cp ex_*.log.tex       "$ARXIV_DIR/"
  cp stata.sty          "$ARXIV_DIR/"   # StataCorp log/table macros
  rm stpois_arxiv.tex   # remove from sj/ once copied

  echo "=== Compiling arXiv PDF (stpois_arxiv.pdf) ==="
  cd "$ARXIV_DIR"
  pdflatex -interaction=nonstopmode stpois_arxiv.tex || true
  bibtex stpois_arxiv
  pdflatex -interaction=nonstopmode stpois_arxiv.tex || true
  pdflatex -interaction=nonstopmode stpois_arxiv.tex
fi

echo ""
echo "=== build.sh complete ==="
[ "$BUILD_SJ"    -eq 1 ] && echo "  SJ PDF:    paper/sj/main.pdf"
[ "$BUILD_ARXIV" -eq 1 ] && echo "  arXiv PDF: paper/arxiv/stpois_arxiv.pdf"
echo ""
echo "To regenerate Stata outputs before rendering:"
echo "  do paper/sj/sj_examples.do       (SJ log output)"
echo "  do paper/benchmark.do            (benchmark_results.csv)"
echo "  do paper/benchmark_hdfe.do       (benchmark_hdfe_results.csv)"
echo "  do paper/benchmark_ppmlhdfe.do   (benchmark_ppmlhdfe_results.csv)"
