# Gap Analysis

## A3.1 — Claims with NO experimental evidence (most serious)
| Claim | Priority | Why it matters | Required experiment |
|-------|----------|----------------|---------------------|
| **S5 diverging dimension** ($q_N^3/N\to0$: $\|\bF_N^{1/2}(\hat\btheta-\btheta_0)\|=O_p(\sqrt{q_N})$ + contrast CLT + relative-consistent plug-in variance) | SECONDARY (asymptotics chain) | It is the only theorem in the paper with **zero** simulation support; the register-scale motivation ($k_N$ in the thousands) lives in exactly this regime | **E-new1**: Poisson/logistic with $k_N$ growing with $N$ along $q_N^3/N=$ const (and a violating path $q_N^3/N\to\infty$ as negative control); verify (a) fixed-contrast coverage → nominal, (b) intrinsic-norm error scales as $\sqrt{q_N}$, (c) plug-in variance relative error → 0 |
| Prop S1 penalized (ridge/lasso path on cell moments) | PERIPHERAL | Stated as an extension; no numbers | **E-new2** (optional): one lasso-path timing+correctness check vs a dense penalized solver on a small cell design; OR relabel Prop S1 explicitly as an un-simulated methodological extension |

## A3.2 — Experiments with adequacy problems (medium)
| Experiment | Flaw | Fix |
|-----------|------|-----|
| S4 classical asymptotics | No dedicated fixed-$q$ normality/coverage experiment; only exercised indirectly | Add a small fixed-$q$ coverage/QQ panel, or state explicitly that Study 3 (its $m=20$, small-$J/N$ rows) + the application SEs stand in for it |
| S1 "converges from any starting value" | Only default starts tested | Add 2–3 dispersed/adversarial initializations to one absorbed-factor DGP and show identical optimum |

## A3.3 — Reporting / discipline (revision quality)
| Issue | Where | Fix |
|-------|-------|-----|
| Coverage estimates lack MC uncertainty | Study 3 (and implicitly 4) | Add Wilson/binomial CI (±~1.4% at B=1000) or an MCSE column |
| Bias claims without SD/MCSE columns | Study 5 | Add SD or MCSE so the bias figure carries uncertainty |
| Deterministic-rate realizations | Study 2 | One sentence noting weights are fixed at the Poisson optimum (pre-empts "is this one draw?") |

## A3.4 — Selection-bias risks
None material. Exactness is algebraic (hard to cherry-pick); imbalance/rare-event and
adversarial-stability regimes are *added* rather than omitted; failure modes (OOM, singular,
separation, pedestal divergence) are reported, not suppressed. No `SELECTION_RISK`.

## A3.5 — Tuning / procedure gaps
None. No oracle-tuning dependence; MLE has no tuning knob; Firth is a declared estimand. No
`TUNING_GAP`.

## A3.6 — Computational adequacy
**Fully met** (rare). The paper markets speed/scalability and backs it with runtime, peak RSS,
iteration counts, empirical $O(N)$ scaling, and honest failure reporting. No `COMP_GAP`.

## Summary
The simulation section is unusually complete and honest. **One real gap** (S5 diverging
dimension, a theorem with no experiment) and **one peripheral gap** (Prop S1 penalized). Two
minor reporting fixes (MCSE on coverage/bias). No contradictions, no selection or tuning
problems, exemplary computational reporting.
