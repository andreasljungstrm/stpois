# Targeted Simulation Improvement Plan

Minimal work to close the gaps in GAP_ANALYSIS.md — not a redesign. The existing section is
strong; these are additive.

## Priority 1 — close the one real claim-gap (NEW experiment)
- **E-new1: diverging-dimension study (verifies S5).**
  DGP: Poisson (and logistic) with a cell design whose parameter count $q_N=p+k_N$ grows with
  $N$ along the path $q_N^3/N = \text{const}$ (e.g. $k_N \propto N^{1/3}$), plus a **violating**
  path $k_N\propto N^{1/2}$ ($q_N^3/N\to\infty$) as negative control. Balanced cells so
  Assumption S-diverging(ii) holds.
  Report, over $\ge 6$ points on each path, $B\ge500$ reps:
  1. empirical coverage of a fixed low-dimensional contrast $\bC_N'\btheta$ → nominal on the
     valid path, degrading on the violating path;
  2. intrinsic-norm error $\|\bF_N^{1/2}(\hat\btheta-\btheta_0)\|/\sqrt{q_N}$ → constant;
  3. relative error of the plug-in variance $\bC_N'\bF_N(\hat\btheta)^{-1}\bC_N$ → 0.
  Reuse: build on `tilted_glm.py`; new script `sim6_diverging.py`. ~30 lines.

## Priority 2 — strengthen existing / peripheral
- **E-new2 (optional): penalized check (Prop S1).** One lasso-path timing+correctness check
  vs a dense penalized solver (`glmnet`/`sklearn`) on a small cell design; or, if descoped,
  relabel Prop S1 in the text as an un-simulated methodological extension.
- **S1 robustness:** add 2–3 dispersed/adversarial starts to one absorbed-factor DGP; show
  identical optimum (supports "from any starting value").
- **S4 stand-in:** add a small fixed-$q$ coverage/QQ panel, or state that Study 3's small-$J/N$
  rows + application SEs cover it.

## Priority 3 — reporting fixes (NO new runs; recompute from saved outputs)
- Add Wilson/binomial CIs (or MCSE) to the coverage columns of Study 3 (and Study 4).
- Add SD/MCSE columns to Study 5's bias table.
- Add one sentence to Study 2 noting the weights are fixed at the Poisson optimum.

## Reuse vs rerun
- **Reuse (extend from saved `code/output/*.csv`, seeds present):** Studies 1–5 reporting fixes;
  Study 3/4/5 MCSE columns are pure post-processing of existing replicate outputs.
- **New runs required:** E-new1 (diverging), E-new2 (penalized), S1 dispersed-start check —
  these are genuinely new cells, but all reuse the existing library.

## Bottom line
The simulation section already meets top-journal computational and honesty standards. The
single substantive addition a referee is likely to demand is **E-new1** (a diverging-dimension
experiment for S5); everything else is polish.
