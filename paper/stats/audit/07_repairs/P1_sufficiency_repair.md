---
artifact: repair_spec
scope: local
theorem_ids: [prop:sufficiency]
issue_ids: [I-01]
generated: 2026-07-18
---

# Repair: I-01 — Proposition 1(ii) disentangling gap

## Selected Strategy: Replace-Technique on the disentangling step (run the recovery in the cell-saturated parametrization)

### Reason for selection
The conclusion (each within-cell measure $\nu_j$ is determined, hence no fixed-dimensional
cell summary supports the likelihood function) is correct; only the *route* to it is flawed.
The flaw is that the current route needs the cell-offset map $\delta\mapsto(s_1,\dots,s_J)$ to
have open range in $\R^J$, which fails when $J>k$. The fix reorganizes the argument so the
required independent variation of the $s_j$ is supplied *within the existing model* by the
cell-saturated parametrization — the categorical apparatus that defines the cells already
carries a free per-cell level in the collapse setting. No new external assumption and no new
literature are required; the mechanism (real-analyticity + Laplace-transform uniqueness) is
unchanged.

## Repair Ladder Defense
- Chosen ladder level: **L3** (Alternative technique under existing assumptions)
- Chosen repair class: **Replace-Technique** (of the disentangling sub-argument only)
- Claim preserved: **yes**
- Assumptions preserved: **yes** (the cell-saturated parametrization is the natural collapse
  setting — each cell carries a free effect — not a new external hypothesis; see note)

### Phase A Exhaustion Record
| Branch | Tried? | Concrete attempt | Specific obstacle | Why the obstacle is genuine | Verdict |
|--------|--------|------------------|-------------------|-----------------------------|---------|
| L1 Internal correction | yes | Keep the "vary $\delta$" route, just tighten wording | The route needs open range of $\delta\mapsto(s_j)$ in $\R^J$; false for $J>k$ | Linear-algebra fact: $\mathrm{rank}(W)\le k<J$ | ruled out |
| L2 Supporting lemma from existing assumptions | not relevant | — | — | The gap is in the recovery technique, not a missing helper fact | not relevant |
| L3 Alternative technique | yes | Run the per-cell Laplace recovery using a free per-cell offset (cell-saturated design), then note a coarser design is a submodel | none — the saturated design gives $s_j$ independently variable and the boundary expansion recovers each $\nu_j$ | succeeds | **succeeded** |

### Note on "assumptions preserved"
Two readings are possible and I state both honestly:
- **Reading used here (claim + assumptions preserved, L3):** the collapse device operates
  cell-by-cell with a free cell effect; the cell-saturated parametrization ($\bw_j$ containing
  the cell indicators, so $s_j=\delta_j$ free) is the *native* setting of the impossibility
  question, not an added restriction. The proof simply makes this explicit.
- **Strict reading (would be L4):** if one insists Proposition 1(ii) claimed the impossibility
  for an *arbitrary, possibly non-saturated* cell design, then restricting to the saturated
  design narrows the hypothesis. Under that reading the change is an Add-Assumption (L4). I
  record the Assumption-Extension Change Log entry below so the plan is valid under *either*
  reading; the re-audit can accept the L3 defense and treat the L4 log as belt-and-suspenders.

### Assumption-Extension Change Log (recorded for the strict reading)
| Issue ID | Original assumption set | Added assumption (verbatim) | Natural weaker variant considered | Why the weaker variant fails | Scientific-scope impact | Propagation |
|----------|------------------------|-----------------------------|-----------------------------------|------------------------------|-------------------------|-------------|
| I-01 | model (edf)+(linpred); $b$ real-analytic; boundary $b(u)=e^u(1+o(1))$ | "the cell-level design contains a free effect for each cell (equivalently, $\{\bw_j\}$ spans enough directions that $\delta\mapsto(s_1,\dots,s_J)$ is locally onto)" | Keep an arbitrary non-saturated $\bw_j$ and recover each $\nu_j$ | With $\mathrm{rank}(W)<J$ the $s_j$ move in a proper subspace; two distinct $\{\nu_j\}$ can give the same likelihood on the accessible directions, so individual $\nu_j$ are not identified by the likelihood alone | Restricts the *proof* to the saturated collapse setting; the paper's message ("continuous covariates break the fixed collapse") is a per-cell statement and is unaffected — cells carry free effects in every collapse application in the paper | Prop 1 statement + Appendix A proof; nothing downstream (P1 has no dependents) |

### Semantic-Edit Log Pointer
- L4 (strict-reading) Assumption-Extension Change Log row: I-01 (above).
- L3 (primary reading): NA.

## Complete Repaired Proof (part ii disentangling step)

**Repaired claim (normalized).** Consider the model (edf)–(linpred) with $b$ real-analytic and
$b(u)=e^u(1+o(1))$ as $u\to-\infty$, in the cell-saturated parametrization: each cell $j$
carries a free scalar level $s_j$ (so $(s_1,\dots,s_J)$ ranges over an open subset of $\R^J$ as
the cell-level parameters vary). Then any statistic from which the log-likelihood can be
evaluated for parameters in an open set determines each within-cell weighted covariate measure
$\nu_j$, and hence the minimal such statistic has dimension $\gtrsim\sum_j|\mathrm{supp}\,\nu_j|$,
which is unbounded as the number of distinct within-cell covariate values grows.

**Proof strategy:** direct recovery per cell + Laplace-transform uniqueness (unchanged
mechanism), now with the independent-variation premise supplied explicitly.

**Dependency map.**
1. Uses: $\ell(\gamma,\delta)=\sum_j[\gamma'S_j+s_jD_j+\sum_{i\in j}o_iy_i-G_j(s_j,\gamma)]$
   with $G_j(s,\gamma)=\sum_{i\in j}b(s+\gamma'\bx_i+o_i)$ (already in the paper).
2. Cell-saturated parametrization ⇒ the map $(\text{cell params})\mapsto(s_1,\dots,s_J)$ is
   locally onto an open set (the point previously left implicit).
3. Real-analyticity of $b$ ⇒ each $G_j(\cdot,\gamma)$ real-analytic in $s$.
4. Boundary expansion $e^{-s}G_j(s,\gamma)\to\mathcal L_{\nu_j}(\gamma)$ (already in the paper).
5. Laplace transform of a finite discrete measure determines the measure.

**Proof.**
- *Step 1 (isolate each cell).* In the saturated parametrization the cell levels $s_1,\dots,s_J$
  vary independently over an open box $\prod_j I_j$. Holding $\gamma$ fixed and differencing
  $\ell$ in the single coordinate $s_j$ (all other $s_{j'}$ fixed) isolates
  $\partial\ell/\partial s_j=D_j-\partial_s G_j(s_j,\gamma)$; integrating in $s_j$ over $I_j$
  recovers $G_j(\cdot,\gamma)$ up to an additive constant on the open interval $I_j$, for each
  $\gamma$ in the (open) $\gamma$-slice. (The additive constant is fixed by the known $D_j$
  term; only $G_j$ up to a constant is needed for Steps 2–3.)
- *Step 2 (analytic continuation).* $G_j(\cdot,\gamma)$ is a finite sum of translates of the
  real-analytic $b$, hence real-analytic on $\R$; knowledge on the open interval $I_j$
  determines it on all of $\R$ by analytic continuation.
- *Step 3 (Laplace recovery).* By the boundary hypothesis,
  $e^{-s}G_j(s,\gamma)\to\sum_{i\in j}e^{\gamma'\bx_i+o_i}=\mathcal L_{\nu_j}(\gamma)$ as
  $s\to-\infty$ (exact for Poisson; via $\log(1+\varepsilon)=\varepsilon(1+o(1))$ for Bernoulli).
  Thus the statistic determines $\mathcal L_{\nu_j}$ on the open $\gamma$-slice;
  $\mathcal L_{\nu_j}$ is entire (finite discrete $\nu_j$), so it is determined everywhere by
  analytic continuation, and the Laplace transform of a finite positive measure determines the
  measure. Hence $\nu_j$ is determined.
- *Step 4 (dimension).* $\nu_j$ is a finite atomic measure on $\R^p$ with one atom per distinct
  within-cell value of $(\bx_i,o_i)$; its minimal description has dimension of order
  $(1+p)\cdot|\mathrm{supp}\,\nu_j|$. Summing over cells, any statistic permitting likelihood
  evaluation has dimension $\gtrsim\sum_j|\mathrm{supp}\,\nu_j|$, which grows with the number of
  distinct within-cell covariate configurations and is therefore not fixed-dimensional. ∎

**Coarser (non-saturated) designs.** A non-saturated cell-level design is a constrained
submodel of the saturated one; the practical claim the paper makes — that *the fixed cell
collapse* (which assigns each cell its own effect) cannot support the likelihood once a
continuous covariate is present — is exactly the saturated statement, since the collapse device
gives each cell a free level by construction. No generality relevant to the paper is lost.

### Verification checklist
- [x] Statement matches what is proved (not stronger) — now explicitly the saturated-cell case
- [x] Every assumption used is listed (added: independent variation of $s_j$, supplied by saturation)
- [x] Every nontrivial implication justified (differencing/integrating in $s_j$; analytic continuation)
- [x] Boundary expansion direction correct (Poisson exact; Bernoulli via $\log(1+\varepsilon)$)
- [x] Edge case $J\le k$ (already covered) and $J>k$ (now covered) both handled
- [x] No downstream dependence (P1 is standalone)

### Repair provability status
**PROVABLE AS STATED** (in the saturated-cell parametrization, which is the native collapse
setting). No new external reference required — the cited machinery (real-analyticity, Laplace
uniqueness) already appears in the paper.
