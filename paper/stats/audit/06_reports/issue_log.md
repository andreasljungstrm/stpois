---
artifact: issue_log
scope: global
source_files: [main.tex, supplement.tex]
theorem_ids: [prop:sufficiency, thm:tilted, cor:complexity, prop:noncanonical, thm:global-main, thm:global, thm:rate, cor:bipartite, thm:multiG, thm:classical, thm:diverging, lem:conditional, thm:proportional, prop:penalized, prop:firth, lem:levelsets]
issue_ids: [I-01, I-02, I-03, I-04, I-05]
commit: pending
generated: 2026-07-18
generator: proofcheck v1.7 (manual pipeline, Codex disabled)
---

# Issue Log

| ID | Severity | Confidence | Unit | Summary |
|----|----------|-----------|------|---------|
| I-01 | **S2** | HIGH | P1 `prop:sufficiency` (ii) | Laplace-transform recovery of each $\nu_j$ needs the cell-offset map $\delta\mapsto(s_1,\dots,s_J)$ to have open range in $\R^J$; this holds only for $J\lesssim k$, but the paper's regime is $J\gg k$. The "in general after passing to the identified subfamily" clause is unjustified. Conclusion likely true but the general-case argument is a gap. No downstream dependents. |
| I-02 | S3 | HIGH | T1 `thm:tilted` (i) | The "$b'$-tilt" $P^{b'}_{j,\theta}$ is a probability measure only if $b'>0$; true for in-scope families but positivity is left implicit. |
| I-03 | S3 | HIGH | P1 `prop:sufficiency` (iii) | $\bar\bx_j$ denotes the $\nu_j$-weighted mean in (iii) but the plain cell mean in (i); notation drift. |
| I-04 | S3 | MED | S4 `thm:classical` | Correspondence between Assumption S-classical and the exact regularity conditions of Fahrmeir–Kaufmann (1985) is asserted, not re-derived. |
| I-05 | S3 | MED | S1 `thm:global` | Algorithm S1 runs the inner α-loop to tolerance then a single damped-Newton β-step, whereas the proof analyzes essentially-cyclic single block updates; the equivalence is correct but stated only implicitly. |

## Severity roll-up
- S0 (fatal): **0**
- S1 (major, breaks a theorem): **0**
- S2 (moderate): **1** (I-01)
- S3 (minor): **4** (I-02, I-03, I-04, I-05)

No issue lies on a load-bearing dependency edge: the single S2 (I-01) is in a standalone
motivating proposition with no downstream consumers.
