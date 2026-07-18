---
artifact: assumption_ledger
scope: global
source_files: [main.tex, supplement.tex]
generated: 2026-07-18
---

# Assumption Ledger

| ID | Assumption | Loc | Scope | Used by | Strength needed | Status |
|----|-----------|-----|-------|---------|-----------------|--------|
| MOD | EDF likelihood, $b$ strictly convex & $C^3$, canonical link until §2.4 | main L338–347 | global | T1,C1,P1 | as stated | OK |
| A-can `ass:canonical` | Canonical link, $b$ strictly convex & twice $C$-diff on $\Theta=\R$ | supp L219 | S1,S2,S4,S5 | as stated; S2/S5 additionally use $b\in C^3$ | OK; S1(iv) relaxes to log-concave link |
| A-mle `ass:mle-exists` | MLE exists (no separation), unique mod $\mathcal K$ | supp L224 | S1,S2 | data condition | OK; enforced by drop rule for absorbed blocks |
| A-blk `ass:blocks` | (i) each absorbed group has an interior-support obs; (ii) $F_{\beta\beta}$ nonsingular on initial superlevel set | supp L229 | S1,S2 | (i) ⇒ finite block maximizers; (ii) ⇒ bounded preconditioner | OK; (i) enforced constructively |
| A-cl `ass:classical` | $q$ fixed, bounded covariates; $\lambda_{\min}(F_N)\to\infty$; $N^{-1}F_N\to F_\infty\succ0$; interior $\theta_0$, $b\in C^3$ | supp L511 | S4 | matches Fahrmeir–Kaufmann | OK (delegation) |
| A-div `ass:diverging` | bounded covariates, compact true predictors; balanced leverage $\max_i h_i\le Cq_N/N$; $q_N^3/N\to0$; $b\in C^3$, $b''$ bounded away from $0,\infty$, $b'''$ bounded | supp L547 | S5 | all four load-bearing in proof | OK — each used explicitly (Steps 1–5) |
| A-prop `ass:proportional` | $2\le m_j\le\bar m$, $J/N\to c\in(0,1)$; bounded covariates & $\sup_j\alpha_j\le C$; $\bar F_J\to\bar F\succ0$; interior $\gamma_0$ in compact $\Gamma$ | supp L647 | S6 | all load-bearing | OK |

## Assumption-propagation check (Pass 3)
- Every cited assumption is in scope where used. No assumption is used beyond its stated
  scope. No contradictions among assumptions.
- **A-blocks(ii)** (nonsingular $F_{\beta\beta}$ on the *initial* superlevel set) is the
  one assumption doing quiet work in S1: it supplies the bounded-preconditioner constant
  $\bar\lambda$ for the sufficient-ascent step. Verified it is invoked correctly (S1 proof,
  "limit points stationary", β-block).
- No unused assumptions detected. A-div's four parts each map to a proof step (leverage→ε_N,
  $q_N^3/N$→remainder, $b'''$ bounded→Lipschitz $b''$, $b''$ bounded→$\underline\omega,\bar\omega$).
