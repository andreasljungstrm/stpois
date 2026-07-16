# Theoretical / methodological paper (computational statistics)

This directory contains the manuscript

> **Computationally Efficient Estimation of High-Dimensional Exponential
> Family Models via Tilted Cell Moments and Alternating Projections**

targeted at a computational statistics journal (formatted in the style of
the *Journal of Computational and Graphical Statistics*).  It is the
theoretical companion to the `stpois` software article in `paper/sj/` and
`paper/arxiv/`: the Stata-specific material is stripped out, the
framework is generalized from Poisson event-history models to arbitrary
exponential dispersion families, and the paper's contributions are
theorems rather than commands:

1. **Theorem 1 (+ Corollary 1).** Exact representation of the score and
   information as exponentially tilted within-cell moments; per-iteration
   complexity O(N·p²+J·k²) instead of O(N·(p+k)²), with no approximation.
2. **Theorems 2–3 (+ Corollary 2).** Global convergence of the nonlinear
   block Gauss–Seidel fixed-effect absorption cycle (the gap noted in the
   applied HDFE-GLM literature), and its local linear rate: for two
   absorbed factors the contraction equals cos²(Friedrichs angle) = the
   squared second singular value of a weighted bipartite transition
   kernel.
3. **Theorems 4–6.** Asymptotics in three regimes: fixed dimension;
   diverging dimension with q³/N → 0 (self-contained proof in the
   intrinsic norm under a balanced-leverage condition); and the
   proportional regime J/N → c for Poisson, where profile = conditional
   likelihood implies no incidental-parameter bias, asymptotic normality,
   efficiency, and valid cluster-robust inference.
4. **Propositions 3–4.** Penalized (ridge/lasso, proximal Newton)
   estimation at cell-collapsed cost, and Firth's bias-reducing
   adjustment computed inside the tilted pass.

## Files

- `main.tex`, `references.bib` — the manuscript (compile with
  `make`, or `pdflatex → bibtex → pdflatex ×2`).
- `main.pdf` — compiled manuscript.
- `code/` — self-contained Python/NumPy reference implementation and all
  simulation/case-study scripts (see `code/README.md`).  Every number in
  Sections 6–7 of the paper regenerates from these scripts with fixed
  seeds.

## Build

```sh
make            # builds main.pdf (needs texlive incl. natbib/booktabs/algorithm)
make clean
```
