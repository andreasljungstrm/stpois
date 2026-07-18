# Per-Experiment Adequacy Audit

General note: many headline studies verify an **algebraic exactness** claim (T1(iv)), not a
stochastic rate. For those, MCSE / replication counts are not the right yardstick — a single
deterministic fit per configuration with a reported max deviation is the correct and
sufficient design. MC-precision criteria apply only to the genuinely stochastic studies
(3, 4, 5).

## Study 1 / 1b / 1p / 1d / 1e / 1f (exactness + cost + scaling)
| Criterion | Status | Note |
|-----------|--------|------|
| Loss object matches claim | ✅ | max abs iterate/coef deviation is exactly the exactness target |
| Grid breadth | ✅ | 3 families × several $J,p,N$; $p$-sweep; scaling over 2 orders of $N$ |
| MCSE | n/a | deterministic exactness (not a MC estimand) |
| Comparators | ✅ | dense (ground truth), sparse-IRLS, glum, statsmodels, pyfixest |
| Failure modes reported | ✅ | OOM / singular / timeout stated explicitly (Table 1d) |
| Computational diagnostics (A2.10) | ✅ | runtime, peak RSS, iters all reported — **COMP adequacy fully met** |
| Captions content-bearing | ✅ | DGP, $N$, $J$, $p$, tolerance stated |
| Selection risk (A2.8) | LOW | exactness is algebraic; hard to cherry-pick; imbalance/rare regimes explicitly added (1f) |
**Verdict: adequate (strong).**

## Study 1g (stability battery)
| Criterion | Status | Note |
|-----------|--------|------|
| Adversarial coverage | ✅ | collinearity, ill-scaling, extreme offsets, huge η, Zipf, pedestal |
| Honesty about failure | ✅✅ | the $10^8$-pedestal divergence is reported, not hidden; dense diverges identically; remedy (standardize) shown |
| Diagnostic separation | ✅ | vs-dense and vs-standardized columns separate conditioning from reproduction |
**Verdict: adequate (exemplary honesty).**

## Study 2 (G=2 and G=3 rates)
| Criterion | Status | Note |
|-----------|--------|------|
| Metric matches theorem | ✅ | empirical per-sweep contraction vs predicted spectral rate |
| Grid | ✅ | $\rho\in\{0,.3,.6,.8,.9,.95,.99\}$ spans expander→near-disconnected |
| Precision | ✅ | rate is deterministic given optimum weights; sweeps-to-$10^{-13}$ reported |
| Both halves of S3 | ✅ | exact rate (ii) AND product bound (iii) shown, bound's conservativeness demonstrated |
**Verdict: adequate (strong).** Minor: the empirical contraction is one realization per $\rho$;
acceptable since the object is a deterministic operator norm, but a sentence noting the weights
are fixed at the Poisson optimum would preempt a referee question.

## Study 3 (proportional regime) — stochastic
| Criterion | Status | Note |
|-----------|--------|------|
| Reps | ✅ | B = 1000 |
| Asymptotic path | ✅ | $J/N=1/m$ held; $c=1/2$ extreme case included |
| Metric matches theorem | ✅ | bias, SD, mean SE, coverage of both estimators |
| MCSE on coverage | ⚠ | coverage reported to 3 decimals but **no Wilson/binomial CI** on the coverage estimate (with B=1000, ±~1.4% at 95%) `[precision-reporting]` |
| Truth source (A2.7) | ✅ | analytic $\gamma_{01}=0.3$ |
| Negative control | ✅✅ | effects-ignored contrast → zero coverage, showing FE are non-ignorable |
**Verdict: adequate; single minor reporting fix (add MCSE/Wilson band on coverage).**

## Study 4 (Firth / separation) — stochastic
| Criterion | Status | Note |
|-----------|--------|------|
| Reps | ✅ | B = 500 |
| Conditioning honesty | ✅ | MLE summaries condition on non-separated reps, stated |
| Grid | ✅ | $N\in\{300,600,1200,2400\}$ shows separation vanishing |
| Cost reported | ✅ | ~2× per-iteration overhead, matches Prop S2 |
**Verdict: adequate.**

## Study 5 (logit bias) — stochastic
| Criterion | Status | Note |
|-----------|--------|------|
| Reps | ✅ | B = 500 |
| Metric matches theorem | ✅ | mean/bias vs Andersen $2\gamma_0$ at $m=2$; $O(1/m)$ decay |
| Comparator | ✅ | exact conditional MLE |
**Verdict: adequate.** Minor: add SD/MCSE columns so the bias claim carries uncertainty.

## Applications (flights, register)
| Criterion | Status | Note |
|-----------|--------|------|
| External validation | ✅✅ | statsmodels ($10^{-13}$), pyfixest ($10^{-9}$) — real data |
| Reuse legitimacy (A2.6) | ✅ | scripts + fixed seeds + saved outputs present; deterministic |
| Truth source | ✅ | register truth is the known generator; flights validated cross-implementation |
**Verdict: adequate (strong external checks).**

## Cross-cutting
- **Tuning audit (A2.9):** no oracle tuning issue — the MLE has no tuning knob; tolerance
  $10^{-9}$ stated. Firth is a defined estimand, not a tuned fallback (stated). No `TUNING_GAP`.
- **Reuse (A2.6):** per-study scripts, fixed seeds, saved CSVs ⇒ extension rather than rerun is
  legitimate for all studies.
