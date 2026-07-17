# Computational statistics paper (main + online supplement)

This directory contains the manuscript

> **Exact Estimation of Exponential Family Models on Massive Mixed
> Categorical–Continuous Data via Tilted Cell Moments**

targeted at a computational statistics journal (JCGS-style), together
with its **online supplement**.  It is the theoretical/methodological
companion to the `stpois` software article in `paper/sj/` and
`paper/arxiv/`.

## Structure

**Main paper (`main.tex` → `main.pdf`, ~19 pp.)** — one computational
idea developed end to end:

1. *Proposition 1.* Ordinary cell collapse fails with continuous
   covariates: evaluating the likelihood function requires the entire
   within-cell empirical distribution (Laplace-transform uniqueness).
2. *Theorem 1 + Corollary 1.* The exact score and Fisher information at
   any parameter value are three low-order moments of an exponentially
   tilted within-cell measure; each Newton/Fisher iteration costs one
   streaming O(N·p²) pass plus O(J·k²) cell algebra, and the iterates
   are identical to full-design Newton.  Proposition 2 covers
   non-canonical links (gamma–log: weights ≡ 1).
3. Algorithms and implementation: streaming pass, numerical safeguards,
   memory/parallelism, exact cell-accumulated robust/cluster VCEs, and
   the composition with absorbed fixed-effect factors.
4. Experiments: iterate-level exactness at 1e-13–1e-15 (Poisson,
   logistic, gamma); total-time benchmarks vs dense Newton,
   `statsmodels`, and a sparse-design IRLS; imbalance and rare-event
   sensitivity.
5. Applications: the public `nycflights13` data (327k flights, logistic
   on-time model, J=3,863 cells, verified against `statsmodels` at 10x
   less cost) and a register-scale synthetic study (5M episodes, 10,000
   absorbed municipality–year effects, seconds on one core).

**Online supplement (`supplement.tex` → `supplement.pdf`, ~26 pp.)** —
companion theory with proofs and simulations:

- *Theorems S1–S3, Corollary S1.* Global convergence of the nonlinear
  block Gauss–Seidel absorption cycle for canonical and log-concave
  links; local rate = cos²(Friedrichs angle) = bipartite spectral gap
  for two factors; for G ≥ 3, the exact rate as the spectral radius of
  a level-space sweep operator built from pairwise cross-tabulations,
  with a product-of-sines certificate.
- *Theorems S4–S6.* Asymptotics: fixed dimension; diverging dimension
  (q³/N → 0, balanced leverage); the proportional regime J/N → c for
  Poisson (profile = conditional likelihood ⇒ no incidental-parameter
  bias, valid cluster-robust inference), and the logit
  failure-and-repair (Andersen's 2γ limit; exact conditional logit).
- Firth's adjustment inside the tilted pass (Proposition S2), penalized
  estimation (Proposition S1), and simulation studies for all of the
  above.

## Files

- `main.tex`, `supplement.tex`, `references.bib` — sources (shared
  bibliography); `main.pdf`, `supplement.pdf` — compiled.
- `code/` — self-contained Python/NumPy reference implementation and
  all experiment scripts (see `code/README.md`).  Every number in both
  documents regenerates from these scripts with fixed seeds.

## Build

```sh
make            # builds main.pdf and supplement.pdf
make clean
```
