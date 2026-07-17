# Reproduction code

Reference implementation and numerical studies for the manuscript

> *Exact Estimation of Exponential Family Models on Massive Mixed
> Categorical–Continuous Data via Tilted Cell Moments*

and its online supplement.  Everything is plain Python; the core
requires only `numpy` and `scipy`.  The competitor benchmarks additionally
need `statsmodels`, `glum`, and `pyfixest` (`pip install statsmodels glum
pyfixest`); each competitor is fitted in a fresh subprocess
(`_bench_worker.py`) so peak memory is isolated.  Every number in the
main paper (Sections 5–6) and in the supplement (Section S4) regenerates
from fixed seeds.

## Contents

| File | Role |
|---|---|
| `tilted_glm.py` | Library: exponential-dispersion families (Poisson, logistic, gamma/log), tilted-moment Newton (`fit_tilted`, Theorem 1), dense baseline (`fit_full`), nonlinear Gauss–Seidel absorption (Algorithm 1), profile/conditional Poisson (`profile_poisson`, Lemma 1 / Theorem 7), weighted alternating projections with the Friedrichs-angle rate predictor (Theorem 3 / Corollary 2) and the G ≥ 3 exact spectral rate and product certificate (`gs_spectral_rate`, `friedrichs_product_bound`, Theorem 4), Firth's adjusted score inside the tilted pass (Proposition 4), OIM / robust / cluster variance estimators. |
| `sim1_validation.py` | Study 1 (main §5.1): tilted vs dense Newton — iterate-level exactness and per-iteration cost, three families. |
| `sim1b_baselines.py` | Study 1f (main §5.6): sensitivity to severe cell imbalance and rare events. |
| `bench_baselines.py` | Studies 1b/1d (main §5.2–5.3): time + peak-memory benchmarks vs `statsmodels`, `glum`, a sparse-design IRLS, and (absorbed factors) `pyfixest`, with p-sweep; drives `_bench_worker.py`. |
| `scaling.py` | Study 1e (main §5.4): O(N) time and memory to N = 5×10⁷; dense-design size contrast. |
| `stability.py` | Study 1g (main §5.5): numerical-stability battery (collinearity, ill-scaling, extreme offsets, huge η, unbalanced cells, catastrophic cancellation). |
| `_bench_worker.py` | Subprocess worker: fits one method on one dataset, reports isolated peak RSS, time, iterations, coefficients. |
| `sim2_convergence.py` | Rate study (supplement §S4): empirical alternating-projections contraction vs the spectral predictions, for G = 2 (bipartite σ₂²) and G = 3 (exact sweep-operator spectral radius + product bound). |
| `sim3_highdim.py` | Proportional-regime study (supplement §S4): proportional regime $J/N \to c$ — bias, SE accuracy, coverage of profile-information and cell-clustered intervals; naive (effects-ignored) contrast. |
| `sim4_separation.py` | Separation study (supplement §S4): quasi-separation frequency; MLE vs Firth-adjusted tilted iteration. |
| `sim5_logit_bias.py` | Logit-bias study (supplement §S4): logit incidental-parameter bias at J/N = 1/m (profile MLE vs exact conditional logit). |
| `app_flights.py` | Main §6.1 application: on-time performance of `nycflights13` (327k flights, logistic, J=3,863 cells) verified against `statsmodels`; plus a high-dimensional block absorbing 3,692 aircraft, verified against `pyfixest` to 1e-9. Downloads `data/flights.csv` if absent. |
| `app_register.py` | Main §6.2 application: 5M person-year episodes, 10,000 absorbed municipality–year effects, 96 demographic cells, cluster-robust SEs; dense cross-check on a subpopulation. |
| `run_all.py` | Runs all of the above in order (~15–30 min single-core). |
| `output/` | CSV/text results written by the scripts. |

## Running

```sh
pip install numpy scipy statsmodels
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
