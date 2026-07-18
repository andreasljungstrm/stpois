---
artifact: adversarial_review
scope: global
generated: 2026-07-18
note: Codex cross-review intentionally skipped (per user instruction).
---

# Pass 4 — Adversarial review (counterexample + hidden-assumption search)

## Counterexample / edge-case sweep
- **$d=1$, single cell, $J=1$** (P1): the impossibility argument *works* here (one free
  offset direction) — which is precisely why the general-$J$ gap (issue I-01) is easy to miss.
- **Zero-variance / singular $\bF_{\beta\beta}$** (S1, S2): excluded by A-blocks(ii). The
  drop rule enforces A-blocks(i). Handled honestly.
- **Disconnected factor graph** (CS1, S3): $\sigma_2=1$ / $\rho_{GS}=1$ — the theorems
  correctly report non-identification rather than a false rate; power iteration returns 1.
- **$m_j=1$ cells** (S6): excluded by A-proportional(i) ($m_j\ge2$) — necessary, since a
  size-1 cell carries no within-cell contrast; correctly assumed.
- **$D_j=0$ cells** (LS2, S6): $\hat\alpha_j\to-\infty$, cell dropped, summation over
  $D_j>0$ — correctly accounted for in both the algorithm and the theory.
- **Heavy tails / unbounded covariates** (S4,S5,S6): all three regimes assume bounded
  covariates and (S5) bounded $b'',b'''$. No overreach; the boundedness is used, not decorative.

## Hidden-assumption search
- **T1**: implicit positivity $b'>0$ for the "$b'$-tilt" to be a probability measure
  (S3, I-02). No effect on the score/info identity itself.
- **S1(iv)**: relies on $\ell$ concave on $\R^D$ for "stationary ⇒ global." For log-concave
  (non-canonical) links this still holds (sum of concave per-obs contributions). Not hidden —
  correctly flagged in the text as the boundary where the global claim stops.
- **S5 Step 4**: the tightness $O_p(q_N^{3/2}/\sqrt N)=o_p(1)$ *requires* $q_N^3/N\to0$, i.e.
  the assumption is exactly the necessary rate, not a convenience. No slack, no hidden
  stronger condition.
- **S6(iii)** cluster-robust: consistency under within-cell dependence misspecification uses
  correct conditional mean only; the $G/(G-1)$ correction called asymptotically negligible —
  verified consistent.

## External-theorem-misuse check (citation failure mode)
- Fahrmeir–Kaufmann (1985), Gouriéroux–Monfort–Trognon (1984): used within their canonical /
  linear-exponential-family scope (S4). Matching of A-classical to their conditions asserted
  (I-04, S3).
- Smith–Solmon–Wagner (1977) product bound, Aronszajn / Kayalar–Weinert angle results,
  Rockafellar Thm 8.4, Bertsekas Prop 1.2.1, Tseng (2001), van der Vaart Ch. 25: all invoked
  within their stated hypotheses. No misuse found.

## Quantifier / uniformity check (Pass 3 item 5)
- S5 Step 2 upgrades a pointwise curvature bound to a **uniform** one over the ball
  $\|\bs\|\le r\sqrt{q_N}$ via $\max_i|\eta_i-\eta_{0i}|\le\varepsilon_N$ — the uniformity is
  earned through balanced leverage, not assumed. ✓
- S6(ii) Hessian convergence stated $\sup_{\|\gamma-\gamma_0\|\le\epsilon_J}$ — uniform on a
  shrinking neighborhood, which is what the Taylor step needs. ✓
- No pointwise-used-as-uniform defects found.

## Negligibility / dropped-term ledger (Pass 3 item 8)
Only S5 and S6 drop asymptotic terms.
| Unit | Term dropped | Needed scale | Support | Bridge | Verdict |
|------|-------------|--------------|---------|--------|---------|
| S5 | Taylor remainder $\hat\bs-\nabla g(0)$ | $o_p(1)$ | Step-2 Hessian sandwich | $\frac{\delta_N}{1-\delta_N}\|\nabla g(0)\|=O_p(q_N^{3/2}/\sqrt N)$ under $q_N^3/N\to0$ | OK (explicit bridge) |
| S5 | Lindeberg negligible cells | $\max|c_{N,i}|\to0$ | balanced leverage | $O(\sqrt{q_N/N})\to0$ | OK |
| S6 | $L_J-\E L_J$ | $o_p(1)$ uniform | $\Var\le C/J$ + convexity lemma | pointwise→uniform for concave (Andersen 1982) | OK |
| S6 | eval-at-$\hat\gamma$ minus at-$\gamma_0$ (variance est.) | $o_p(1)$ | Lipschitz envelope × $\|\hat\gamma-\gamma_0\|$ | consistency (i) | OK |

No missing bridges, no wrong comparison scale, no unjustified pointwise→uniform upgrade.
