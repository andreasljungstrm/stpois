# Reproduction code

Reference implementation and numerical studies for the manuscript

> *Computationally Efficient Estimation of High-Dimensional Exponential
> Family Models via Tilted Cell Moments and Alternating Projections*

Everything is plain Python; the only dependencies are `numpy` and
`scipy`.  Every number in Sections 6 and 7 of the paper regenerates from
fixed seeds.

## Contents

| File | Role |
|---|---|
| `tilted_glm.py` | Library: exponential-dispersion families (Poisson, logistic, gamma/log), tilted-moment Newton (`fit_tilted`, Theorem 1), dense baseline (`fit_full`), nonlinear Gauss–Seidel absorption (Algorithm 1), profile/conditional Poisson (`profile_poisson`, Lemma 1 / Theorem 7), weighted alternating projections with the Friedrichs-angle rate predictor (Theorem 3 / Corollary 2) and the G ≥ 3 exact spectral rate and product certificate (`gs_spectral_rate`, `friedrichs_product_bound`, Theorem 4), Firth's adjusted score inside the tilted pass (Proposition 4), OIM / robust / cluster variance estimators. |
| `sim1_validation.py` | Study 1 (§6.1): tilted vs dense Newton — iterate-level exactness and per-iteration cost, three families. |
| `sim2_convergence.py` | Study 2 (§6.2): empirical alternating-projections contraction vs the spectral predictions, for G = 2 (bipartite σ₂²) and G = 3 (exact sweep-operator spectral radius + product bound). |
| `sim3_highdim.py` | Study 3 (§6.3): proportional regime $J/N \to c$ — bias, SE accuracy, coverage of profile-information and cell-clustered intervals; naive (effects-ignored) contrast. |
| `sim4_separation.py` | Study 4 (§6.4): quasi-separation frequency; MLE vs Firth-adjusted tilted iteration. |
| `sim5_logit_bias.py` | Study 5 (§6.5): logit incidental-parameter bias at J/N = 1/m (profile MLE vs exact conditional logit). |
| `app_register.py` | Section 7 case study: 5M person-year episodes, 10,000 absorbed municipality–year effects, 96 demographic cells, cluster-robust SEs; dense cross-check on a subpopulation. |
| `run_all.py` | Runs all of the above in order (~15–30 min single-core). |
| `output/` | CSV/text results written by the scripts. |

## Running

```sh
pip install numpy scipy
cd paper/stats/code
python3 run_all.py          # or any individual script
```

## Notes

- The library favours clarity over micro-optimisation: all cell
  accumulations are `np.bincount` segmented sums, so each pass is O(N)
  regardless of the number of cells or absorbed levels, and no
  N × (p+k) design matrix is ever formed.
- The Stata package `stpois` (repository root) implements the Poisson
  event-history case of the same algorithms in production form.
