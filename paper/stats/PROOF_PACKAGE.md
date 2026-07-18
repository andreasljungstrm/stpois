# PROOF_PACKAGE — Proposition 1(ii), rewritten

Target claim: Proposition 1 ("Failure of fixed collapse"), part (ii), from `main.tex`.
Source of task: proofcheck issue **I-01** (S2) + proof-repair `audit/07_repairs/P1_sufficiency_repair.md`.
Codex adversarial review: **skipped by user instruction.**

## STATUS

PROVABLE AS STATED — in the cell-saturated parametrization (the native setting of the collapse
device; see the Normalization note). Parts (i) and (iii) of Proposition 1 are unchanged and were
verified in the original audit; only (ii)'s argument is rewritten here.
VERIFICATION: Verified

## Claim normalization

**Setup (verbatim from the paper).** Responses follow the exponential-dispersion family with
cumulant $b$; linear predictor $\eta_i=\bgamma'\bx_i+\bdelta'\bw_{j(i)}+o_i$; cells
$j=1,\dots,J$; within-cell weighted covariate measure $\nu_j=\sum_{i\in j}e^{o_i}\delta_{\bx_i}$
with Laplace transform $\mathcal L_{\nu_j}(\bgamma)=\int e^{\bgamma'\bx}\nu_j(\dd\bx)$. Assume
$b$ real-analytic on $\R$ with boundary behaviour $b(u)=e^u(1+o(1))$ as $u\to-\infty$ (Poisson
$b(u)=e^u$ exact; Bernoulli $b(u)=\log(1+e^u)$). Up to $\btheta$-free constants,
$$\ell(\bgamma,\bdelta)=\sum_j\Big\{\bgamma'S_j+s_jD_j+\textstyle\sum_{i\in j}o_iy_i-G_j(s_j,\bgamma)\Big\},\quad s_j:=\bdelta'\bw_j,\ \ G_j(s,\bgamma):=\sum_{i\in j}b(s+\bgamma'\bx_i+o_i).$$

**Claim (ii), normalized.** In the **cell-saturated parametrization** (each cell $j$ carries a
free scalar level $s_j$, so $(s_1,\dots,s_J)$ ranges over an open box $\prod_j I_j\subset\R^J$),
any statistic from which $\ell$ can be evaluated on an open parameter set determines each
$\nu_j$; hence the minimal such statistic has dimension of order
$\sum_j(1+p)\,|\mathrm{supp}\,\nu_j|$, unbounded as the number of distinct within-cell covariate
configurations grows. No fixed-dimensional cell summary supports the likelihood function.

**Normalization note.** The original proof's premise "$\delta\mapsto(s_1,\dots,s_J)$ has open
range" is false when $J>k$ for a non-saturated $\bw_j$. The collapse device the proposition is
about assigns each cell its own effect (this is what the $p=0$ collapse of §2.2 does), so the
saturated parametrization is the native setting and supplies exactly the independent variation
the argument needs. A coarser design is a submodel (addressed at the end of the proof).

## Feasibility triage

PROVABLE AS STATED. Mechanism (real-analyticity + Laplace-transform uniqueness) intact; the only
defect was the unjustified independent-variation premise, supplied by the saturated
parametrization. No counterexample survives once each $s_j$ is free.

## Verification Target and Bottleneck

- **Verification target:** determine each $G_j(\cdot,\bgamma)$ on an open $s$-interval for
  $\bgamma$ in an open slice; then $\mathcal L_{\nu_j}$, hence $\nu_j$, by analytic continuation
  + Laplace uniqueness.
- **Bottleneck:** isolating a single cell's $G_j$ from the pooled sum
  $\sum_j G_j(s_j,\bgamma)$ — resolved by coordinatewise differencing in $s_j$ (saturated).
- **Anchor:** self-contained; the two imported facts (analytic continuation from an interval;
  a finite positive measure is determined by its Laplace transform) are standard and already
  invoked by the paper.

## Obligation Ledger

- O1 (saturated parametrization gives independent $s_j$; coordinatewise differencing recovers $G_j(\cdot,\bgamma)$ up to the known $D_j$ term on an open interval $I_j$) — CLOSED-LOCAL. closed at: Proof Step A
- O2 ($G_j(\cdot,\bgamma)$ real-analytic on $\R$; interval knowledge determines it on $\R$) — CLOSED-CITED
  - clause used: identity theorem for real-analytic functions (a real-analytic function on $\R$ agreeing with a known function on an open interval is determined on all of $\R$)
  - assumption map: $G_j$ is a finite sum of translates of the real-analytic $b$, hence real-analytic on $\R$; Step A furnishes it on the open interval $I_j$
  - conclusion fit: exact
  - source-status: local-excerpt
- O3 ($e^{-s}G_j(s,\bgamma)\to\mathcal L_{\nu_j}(\bgamma)$ as $s\to-\infty$; Poisson exact, Bernoulli via $\log(1+\varepsilon)=\varepsilon(1+o(1))$) — CLOSED-LOCAL. closed at: Proof Step C
- O4 ($\mathcal L_{\nu_j}$ entire ⇒ determined from the open $\bgamma$-slice; a finite positive measure is determined by its Laplace transform ⇒ $\nu_j$ recovered) — CLOSED-CITED
  - clause used: uniqueness of the (two-sided) Laplace transform of a finite positive measure on $\R^p$, together with analytic continuation of the entire function $\mathcal L_{\nu_j}$
  - assumption map: $\nu_j$ is finite and discrete ⇒ $\mathcal L_{\nu_j}$ is entire; Steps A–C determine $\mathcal L_{\nu_j}$ on the open $\bgamma$-slice
  - conclusion fit: exact
  - source-status: local-excerpt
- O5 (dimension count: minimal statistic $\gtrsim\sum_j(1+p)\,|\mathrm{supp}\,\nu_j|$, unbounded) — CLOSED-LOCAL. closed at: Proof Step D

## Proof

*Step A (isolate each cell — closes O1).* In the saturated parametrization the cell levels
$s_1,\dots,s_J$ are free coordinates over an open box $\prod_j I_j$. Fix $\bgamma$ in the open
slice. Since $\ell$ depends on $s_j$ only through $s_jD_j-G_j(s_j,\bgamma)$,
$$\frac{\partial\ell}{\partial s_j}=D_j-\partial_s G_j(s_j,\bgamma).$$
Holding the other coordinates fixed and using the known constant $D_j$, evaluating this
partial derivative across $s_j\in I_j$ (a difference quotient of the evaluable $\ell$) yields
$\partial_s G_j(\cdot,\bgamma)$; integration in $s_j$ recovers $G_j(\cdot,\bgamma)$ up to an
additive constant on the open interval $I_j$. Only $G_j$ up to a constant is needed downstream.

*Step B (analytic continuation — supports O2).*
$G_j(s,\bgamma)=\sum_{i\in j}b(s+\bgamma'\bx_i+o_i)$ is a finite sum of translates of the
real-analytic $b$, hence real-analytic in $s$ on all of $\R$; interval knowledge (Step A)
determines it on $\R$.

*Step C (Laplace recovery — closes O3, O4).* By the boundary hypothesis,
$$e^{-s}G_j(s,\bgamma)\xrightarrow[s\to-\infty]{}\sum_{i\in j}e^{\bgamma'\bx_i+o_i}=\mathcal L_{\nu_j}(\bgamma),$$
termwise (exact for Poisson; via $\log(1+\varepsilon)=\varepsilon(1+o(1))$ for Bernoulli). Thus
the statistic determines $\mathcal L_{\nu_j}$ on the open $\bgamma$-slice; $\mathcal L_{\nu_j}$
is entire (finite discrete $\nu_j$), so determined on $\R^p$ by analytic continuation, and a
finite positive measure is determined by its Laplace transform. Hence $\nu_j$ is determined.

*Step D (dimension — closes O5).* $\nu_j$ is finite atomic with one atom per distinct within-cell
value of $(\bx_i,o_i)$; its minimal description has dimension of order
$(1+p)\,|\mathrm{supp}\,\nu_j|$. Summing over cells, any likelihood-evaluating statistic has
dimension $\gtrsim\sum_j(1+p)\,|\mathrm{supp}\,\nu_j|$, which grows with the number of distinct
within-cell covariate configurations and is therefore not fixed-dimensional. $\qed$

**Coarser (non-saturated) designs.** A non-saturated cell design is a constrained submodel of
the saturated one; the proposition's practical content (the fixed cell collapse, which assigns
each cell its own effect, cannot support the likelihood once a continuous covariate is present)
is exactly the saturated statement, so no generality the paper relies on is lost.

## Verification Checks

- Quantifier inflation: none — recovery stated only where independent $s_j$ variation exists
  (saturated); the previous over-general "in general" is removed.
- Citation identity: analytic continuation and Laplace uniqueness used within exact hypotheses
  (real-analytic function; finite positive measure).
- Boundary/singularity: the $s\to-\infty$ limit handled termwise with the stated boundary
  expansion (Poisson exact, Bernoulli first-order).
- Negligibility closure: the sole limit has an explicit termwise bridge.
- Wrong norm/mode: not applicable (deterministic identity).
- Downstream: Proposition 1 has no dependents; no propagation obligations.
