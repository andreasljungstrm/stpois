---
artifact: local_check
scope: local
theorem_ids: [prop:sufficiency]
assumption_ids: [MOD]
generated: 2026-07-18
---

## Proof Unit: P1 / Proposition 1 (Failure of fixed collapse)
- Location: main.tex L406 (statement), L1411 (proof, Appendix A)
- Type: proposition (impossibility / negative result)
- **Sketch class**: COMPLETE
- Proof Strategy: (i) direct; (ii) reduction to Laplace-transform uniqueness via
  real-analyticity + boundary expansion; (iii) Jensen
- Provability: (i),(iii) PROVABLE AS STATED; (ii) **PROVABLE AFTER WEAKENING** (see issue)
- Status: **Gap found** (part ii)
- Confidence: HIGH

### Claim Normalization
Part (ii) claims: *any* statistic from which the likelihood function can be evaluated
for $(\gamma,\delta)$ in an open set determines each within-cell cumulant functional
$G_j(\cdot,\gamma)$ on an open set, hence (by analyticity + boundary expansion)
$\mathcal L_{\nu_j}$, hence $\nu_j$; therefore the minimal statistic has dimension growing
with $\sum_j |\mathrm{supp}(\nu_j)|$ — "no fixed-dimensional cell summary can support the
likelihood function."

### Step-by-Step Verification
| Step | Location | Claim | Verdict | Notes |
|------|----------|-------|---------|-------|
| ii.1 | L1425 | ℓ computable from $\{S_j,D_j,G_j\}$ | Valid | direct grouping |
| ii.2 | L1434–1442 | evaluating ℓ on open set recovers each $G_j(\cdot,\gamma)$ on an open $s$-interval | **Gap** | requires $\delta\mapsto(s_1,\dots,s_J)$ to have **open range in $\R^J$**; see below |
| ii.3 | L1443–1453 | $G_j$ real-analytic ⇒ analytic continuation to $\R$; $e^{-s}G_j\to\mathcal L_{\nu_j}$ as $s\to-\infty$ | Valid | boundary hypothesis $b(u)=e^u(1+o(1))$ holds for Poisson (exact) and Bernoulli ($\log(1+\varepsilon)=\varepsilon(1+o(1))$) |
| ii.4 | L1456–1460 | $\mathcal L_{\nu_j}$ entire (finite discrete $\nu_j$) determined everywhere; Laplace transform determines the measure | Valid | standard |
| iii | L1464–1468 | Jensen on strictly convex $e^z$ under $\nu_j/\nu_j(\R^p)$ | Valid | equality iff $\gamma'x$ $\nu_j$-a.s. constant |

### The gap (Step ii.2)
The proof writes $s_j=\delta'w_j$ and asserts one recovers each $G_j(\cdot,\gamma)$ on an
open $s$-interval "provided $\delta\mapsto(s_1,\dots,s_J)$ has open range … as holds
whenever the $w_j$ are affinely independent, and in general after passing to the identified
subfamily."

- $\delta\in\R^k$ and $(s_1,\dots,s_J)=W\delta$ where $W$ is $J\times k$. Its range is a
  subspace of dimension $\le k$. For that range to be **open in $\R^J$** one needs $J\le k$
  and $W$ of full row rank. Affine independence of $\{w_j\}\subset\R^k$ likewise caps $J\le k+1$.
- The **entire motivating regime of the paper is $J\gg k$** (many cells, few cell-level
  parameters). There the $s_j$ are confined to a $k$-dimensional subspace and the individual
  $G_j$ cannot be disentangled from the accessible values of $\sum_j G_j(\delta'w_j,\gamma)$
  by varying $\delta$. Varying $\gamma$ does not help: all cells share $\gamma$.
- The escape clause "in general after passing to the identified subfamily" is not justified;
  passing to the identified subfamily does not restore independent variation of the $s_j$.

Consequently, as written, the recovery-of-each-$\nu_j$ argument establishes the impossibility
only for configurations with **enough free cell-offset directions** ($J\lesssim k$, e.g. a
single cell / pure-continuous design), not "in general."

### Why this is S2 and not S0/S1
- The **conclusion** (continuous covariates break fixed collapse) is almost certainly true
  and is fixable: the object $\{S_j,D_j,G_j(\cdot,\cdot)\}$ genuinely IS the minimal
  likelihood-evaluating statistic, and $G_j$ is an infinite-dimensional (function) object;
  a degrees-of-freedom / dimension-counting argument (or restricting the disentangling claim
  to per-cell sub-experiments) closes it without touching any other result.
- **Nothing downstream depends on P1** (dependency graph: P1 is standalone motivation). No
  theorem, corollary, or numerical claim inherits the gap.
- But it is the paper's advertised impossibility theorem ("a Laplace-transform uniqueness
  argument", abstract + intro), so a top-journal referee would require the $J>k$ case to be
  argued honestly. → **S2**, could be escalated to S1 by a strict referee.

### Issues
| Severity | Confidence | Description | Evidence | Proposed repair |
|----------|-----------|-------------|----------|-----------------|
| S2 | HIGH | Part (ii) recovers individual $G_j$/$\nu_j$ only when the cell-offset map $\delta\mapsto(s_j)$ has open range ($J\lesssim k$); the "in general" claim is unjustified in the paper's $J\gg k$ regime | L1436–1442 | Replace individual-recovery route with a dimension lower-bound argument, or restrict the disentangling step and prove the impossibility via total sufficient-statistic dimension |
| S3 | HIGH | Part (iii): $\bar\bx_j$ is the $\nu_j$-weighted mean here but denotes the plain/cell mean elsewhere (e.g. part i) | L1466 vs L413 | disambiguate notation |

### Final Verdict
**Conditionally verified** — parts (i),(iii) verified; part (ii) has a real proof gap in the
regime of interest but a repairable one with no downstream impact.
