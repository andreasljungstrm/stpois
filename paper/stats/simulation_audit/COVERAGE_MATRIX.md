# Coverage Matrix: Claims × Existing Evidence

Two axes: **Coverage** (is there an experiment aimed at the claim?) and **Evidentiary
strength** (does it identify the claim to top-journal standards?).

| Claim | Priority | Coverage | Strength | Final tag |
|-------|----------|----------|----------|-----------|
| T1(iv) iterate-level exactness (= dense Newton) | PRIMARY | Study 1 (+1b,1d,1f, App) | $10^{-13}$–$10^{-15}$ across 3 families, $N$, $J$, $p$; real-data $10^{-9}$–$10^{-13}$ | **YES[strong]** |
| C1 per-iteration cost / speedup | PRIMARY | Study 1, 1b, 1p | measured 12×–117×; FLOP vs measured explained honestly | **YES[strong]** |
| C1 $O(N)$ time & memory | PRIMARY | Study 1e (to $5{\times}10^7$) | linear fit; dense-design contrast | **YES[strong]** |
| C1 memory footprint (no $Nk$) | PRIMARY | Study 1b, 1e | peak RSS 12–53× lighter | **YES[strong]** |
| Exactness under imbalance/rare events | SECONDARY | Study 1f | 8 scenarios, dev $10^{-15}$–$10^{-14}$ | **YES[strong]** |
| Implementation conditioning (no amplification) | SECONDARY | Study 1g | $\kappa$ to $10^{17}$; honest pedestal failure | **YES[strong]** |
| S1 global convergence (any start) | PRIMARY | every absorbed fit converges; App-register | convergence observed, but *"from any starting value"* not adversarially stress-tested | **YES[weak]** `[stress-coverage]` |
| S2(ii)/Cor S1 two-factor rate $c_F^2$ | PRIMARY | Study 2 (G=2) | 3-decimal match over $\rho\in[0,0.99]$ | **YES[strong]** |
| S3(ii) exact multi-$G$ spectral rate | PRIMARY | Study 2 (G=3) | 3-decimal match; identification diagnostic | **YES[strong]** |
| S3(iii) product-of-sines certificate | SECONDARY | Study 2 (G=3) | shown valid + conservative, as theory says | **YES[strong]** |
| S4 classical $\sqrt N$ normality + efficient + sandwich | SECONDARY | indirectly via Study 3 coverage; App SEs | no dedicated fixed-$q$ normality/coverage plot | **PARTIAL** `[grid, reporting]` |
| S5 diverging dimension ($q_N^3/N\to0$; $\sqrt{q_N}$ rate; contrast CLT) | SECONDARY | **none** | Study 3 is *proportional* ($J/N\to c$), a different regime | **NO** |
| Lem S2 profile = conditional (Poisson) | SECONDARY | underlies Study 3, 5 | algebraic; exercised numerically | **YES[strong]** |
| S6 proportional: no asymptotic bias | PRIMARY | Study 3 | bias ≪ SD at $c=1/2$; 1000 reps | **YES[strong]** |
| S6 coverage of profile-info & cluster CIs | PRIMARY | Study 3 | 94–96% both estimators | **YES[strong]** |
| S6 vs logit failure mode | SECONDARY | Study 5 | profile bias → Andersen $2\gamma_0$; conditional unbiased | **YES[strong]** |
| Prop S1 penalized (ridge/lasso path) | PERIPHERAL | **none** | stated only | **NO** |
| Prop S2 Firth inside tilted pass | SECONDARY | Study 4 | separation freq + MAE, 500 reps | **YES[strong]** |
| P1 impossibility (fixed collapse fails) | SECONDARY | n/a (structural negative result) | not simulable | **N/A** (correctly untested) |

## Roll-up
- PRIMARY claims: 8 covered YES[strong], 1 YES[weak] (S1 "any start"). **0 PRIMARY at NO or CONTRADICTED.**
- SECONDARY: mostly strong; **S5 diverging-dimension = NO** (the one real coverage gap on a theorem); S4 PARTIAL.
- PERIPHERAL: Prop S1 penalized = NO (acceptable for an extension proposition, flag as illustrative).
- **No CONTRADICTED tags** — every covered prediction matches (exactness algebraic; rates to 3 decimals; coverage nominal). The CONTRADICTED protocol (A2.5) is not triggered.
