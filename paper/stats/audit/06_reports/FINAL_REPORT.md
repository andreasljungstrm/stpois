---
artifact: audit_final_report
scope: global
source_files: [main.tex, supplement.tex]
theorem_ids: [prop:sufficiency, thm:tilted, cor:complexity, prop:noncanonical, thm:global-main, thm:global, thm:rate, cor:bipartite, thm:multiG, thm:classical, thm:diverging, lem:conditional, thm:proportional, prop:penalized, prop:firth, lem:levelsets]
assumption_ids: [MOD, ass:canonical, ass:mle-exists, ass:blocks, ass:classical, ass:diverging, ass:proportional]
issue_ids: [I-01, I-02, I-03, I-04, I-05]
commit: pending
generated: 2026-07-18
generator: proofcheck v1.7 Pass 5 (manual pipeline run; Codex cross-review disabled by request)
---

# Final Proof-Check Report

**Paper:** *Tilted Cell Moments: A One-Pass Exact Newton Method for Generalized Linear
Models at Scale* (main + online supplement).

## Executive Summary
- **Overall verdict: Correct modulo one minor revision.** The central engine result (T1)
  and every convergence/rate/asymptotic theorem in the supplement (S1–S6, Lemmas S1–S2,
  Props S1–S2) are mathematically sound as written.
- **Main theorem support:** T1 (tilted-moment representation) — Verified. T2/S1 (global
  convergence) — Verified. S3 (multi-factor exact rate) — Verified, including the delicate
  semisimplicity argument. S5 (diverging dimension) and S6 (proportional regime, "no
  incidental-parameter bias") — Verified with tight, honest use of every assumption.
- **Highest-severity issue: S2** (one), in the standalone impossibility Proposition 1(ii).
- **Checked units: 16 / 16.**
- **Open issues: 5** — S0: 0, S1: 0, S2: 1, S3: 4.
- **Sketches detected: 0** (all proofs COMPLETE; S4 is complete-by-reduction).

```
Sketches detected: 0
├── Expanded to complete proof: 0
├── Determined unprovable (blockage): 0
└── Outstanding: 0
```

## Checked Scope
- **Checked:** all theorem-like environments in both files (main P1, T1, C1, P2, T2;
  supplement S1, S2, Cor S1, S3, S4, S5, LS2, S6, PS1, PS2, LS1) with full per-step
  verification, step-completeness/skip audit, edge-case sweep, and a negligibility ledger
  for the asymptotic proofs.
- **Not deeply checked:** the numerical-experiment tables (§5 main, §S4 supplement) as
  *proof* objects — they are evidence, audited in the Stage-5 theory-simulation pass, not
  here. Bibliographic accuracy of the ~40 external citations was spot-checked for
  proof-load-bearing citations only.

## Main Dependency Chain
Def-tilt → **T1** → C1 (engine). Lemma S1 → **S1**(=T2) → S2 → {Cor S1, S3} (absorption).
Lemma S2 → **S6**; **S5** self-contained; S4 delegated (asymptotics). No circularity; the one
S2 issue (P1) is off every chain.

## Verified Results
| Result | Status | Confidence | Notes |
|--------|--------|-----------|-------|
| T1 tilted representation | Verified | HIGH | exact score/info identity + Newton-step exactness |
| C1 complexity | Verified | HIGH | operation counting exact |
| P2 non-canonical | Verified | HIGH | gamma–log $\omega_i\equiv1$ ✓ |
| S1 / T2 global convergence | Verified | HIGH | BCD on quotient; log-concave extension holds |
| S2 local rate | Verified | HIGH | $c_F^2$ two-factor rate |
| Cor S1 bipartite gap | Verified | HIGH | $\sigma_2$ of normalized biadjacency |
| S3 multi-$G$ exact rate | Verified | HIGH | semisimplicity of eigenvalue 1 fully re-derived |
| S4 classical asymptotics | Conditionally verified | MED | delegation (I-04) |
| S5 diverging dimension | Verified | HIGH | every step tight; $q_N^3/N\to0$ is exactly necessary |
| LS2 profile=conditional | Verified | HIGH | Poisson–multinomial + Schur form |
| S6 proportional regime | Verified | HIGH | no-bias mechanism + efficiency projection re-derived |
| PS1 penalized / PS2 Firth | Verified | HIGH | cell-decomposability confirmed |
| P1 impossibility (i),(iii) | Verified | HIGH | — |
| P1 impossibility (ii) | **Gap found** | HIGH | issue I-01 |

## Open Issues (ranked)
See `issue_log.md`. Headline: **I-01 (S2)** — Proposition 1(ii)'s recovery of each within-cell
measure $\nu_j$ from likelihood evaluability requires independent variation of the cell
offsets $s_j=\delta'\bw_j$, available only when $J\lesssim k$; the paper's regime is $J\gg k$,
where the individual-recovery route breaks and the "in general" clause is not justified.

## Conditional Results
| Result | Condition | Evidence |
|--------|-----------|----------|
| S4 classical asymptotics | A-classical matches Fahrmeir–Kaufmann conditions | assumptions align but not re-derived (I-04) |

## Recommended Repairs
1. **I-01 (S2):** replace the individual-$\nu_j$ recovery in P1(ii) with a dimension /
   degrees-of-freedom lower bound on any likelihood-evaluating statistic, OR restrict the
   disentangling claim to configurations with $\mathrm{rank}(W)=J$ and argue the general case
   separately. (Handled in Stage 2 repair + Stage 6 proof-writer.)
2. **I-02, I-03 (S3):** state $b'>0$ for the $b'$-tilt; disambiguate $\bar\bx_j$.
3. **I-04 (S3):** name the specific Fahrmeir–Kaufmann conditions being invoked.
4. **I-05 (S3):** state the inner-loop/one-sweep equivalence explicitly in the S1 proof.

## Final Judgment
**Correct modulo repairs.** No S0/S1. The paper's mathematical core is solid and, in several
places (S3 semisimplicity, S5 rate accounting, S6 efficiency), unusually careful. The single
substantive finding is a repairable gap in a non-load-bearing motivating proposition. Because
the original audit contains an S2 (not S0/S1), the post-repair convergence test (Stage 3) is
**recommended but not a hard gate**.
