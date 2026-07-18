# Existing Simulations — Parsed Inventory

Mode: **AUDIT** (paper has theorems AND a full simulation section, main §5 + supplement §S4).
Reproduction code present under `code/` (per-study scripts, fixed seeds, saved CSVs in
`code/output/`) ⇒ reuse-legitimacy prerequisites are largely met (see ADEQUACY_AUDIT A2.6).

| Study | Loc | DGP | Grid | Methods | Metrics | Reps (B) | Verifies |
|-------|-----|-----|------|---------|---------|----------|----------|
| 1 | §5.1 `sim1_validation.py` | Poisson/logistic/gamma, cells | $J\in\{50,200,500\}$, $p\in\{2,5\}$, $N\le10^6$ | tilted vs dense Newton | max iterate dev, final dev, per-iter time | 1/config (deterministic exactness) | T1(iv) exactness, C1 cost |
| 1b | §5.2 `bench_baselines.py` | Poisson/logistic, $N{=}2{\times}10^5$, $J{=}200$, $p{=}2$ | fixed | tilted, sparse-IRLS, glum, dense, statsmodels | time, peak RSS, iters, max dev | median over repeats | C1 cost+memory, exactness vs production |
| 1p | §5.2 `bench_baselines.py` | Poisson, $p\in\{2,20,50\}$ | fixed $N,J$ | tilted, dense, sparse | time, peak | median | C1 $p$-scaling |
| 1d | §5.3 `bench_baselines.py` | Poisson, 1 absorbed factor 500–20k levels | fixed | tilted, pyfixest, dense, glum | time, peak, dev | median | absorbed-path cost, exactness |
| 1e | §5.4 `scaling.py` | Poisson, $N=10^6$–$5{\times}10^7$ | $J{=}500,p{=}2$ | tilted | time, peak, dense-design size | 1/N | C1 $O(N)$ time+memory |
| 1f | §5.6 `sim1b_baselines.py` | Zipf(1.6) cells, rare events ~0.2% | 8 family×scenario | tilted vs dense | max dev, iters | — | exactness under imbalance/rarity |
| 1g | §5.5 `stability.py` | 6 adversarial scenarios, $N{=}10^5$, $\kappa$ to $10^{17}$ | — | tilted vs dense vs standardized | $\kappa$, dev-vs-dense, dev-vs-std, converged | — | implementation conditioning |
| 2 (G=2) | §S4.1 `sim2_convergence.py` | 2 factors $J_1{=}100,J_2{=}80$, mixing $\rho$ | $\rho\in\{0..0.99\}$ | weighted AP demeaning | predicted $\sigma_2^2$ vs empirical rate, sweeps | 1/ρ (deterministic rate) | S2(ii), Cor S1 |
| 2 (G=3) | §S4.1 `sim2_convergence.py` | 3 factors, chained mixing $\rho$ | $\rho\in\{0..0.95\}$ | 3-factor demeaning | exact $\rho_{GS}$ vs empirical, norm bound, sweeps | 1/ρ | S3(ii),(iii) |
| 3 | §S4.2 `sim3_highdim.py` | Poisson, $m\in\{2,5,20\}$, 1 FE/cell, effect-correlated covar | $N\in\{20k,80k\}$ | profile Poisson; naive contrast | bias, SD, mean SE, cover_P, cover_C | **1000** | S6 (no bias, coverage), Lem S2 |
| 4 | §S4.3 `sim4_separation.py` | logistic, $J{=}10$, rare exposure 3% | $N\in\{300,600,1200,2400\}$ | MLE vs Firth-tilted | separation freq, med est, MAE | **500** | Prop S2 (Firth), separation |
| 5 | §S4.4 `sim5_logit_bias.py` | logit, $m\in\{2,4,8\}$, 1 FE/cell | $N{=}8000$ | profile MLE vs conditional MLE | mean est, bias | **500** | logit incidental-param bias (contrast to S6) |
| App-flights | §6.1 `app_flights.py` | nycflights13 327k, logistic | real | tilted vs statsmodels; +absorb 3692 aircraft vs pyfixest | coef dev, time | real | external validation |
| App-register | §6.2 `app_register.py` | 5M episodes, 10k absorbed FE | synthetic-truth | tilted; dense subpop cross-check | recovery vs truth, cluster SE, dev | real-scale | full-composition demo |

## Stated conclusions
Exactness is algebraic ⇒ simulations *confirm the implementation inherits it* (main §5.1,5.6
say this explicitly). Rate studies confirm the spectral predictions to 3 decimals. Study 3
confirms the no-incidental-bias theorem; Study 5 exhibits the contrasting logit bias.
