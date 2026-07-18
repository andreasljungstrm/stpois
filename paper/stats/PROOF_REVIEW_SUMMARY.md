# Proof-Review Pipeline — Summary

Full `stat-theory-skills` proof-review pipeline run on `paper/stats/main.tex` (+ `supplement.tex`)
on 2026-07-18. **Codex cross-review was skipped at every stage by request.** No manuscript
`.tex` file was modified — all fixes are proposals in `PATCHES.md`.

Pipeline: `/proofcheck → /proof-repair → /proofcheck --post-repair → /theory-sharpen →
/theory-simulation (AUDIT) → /proof-writer`.

## Headline
The paper's mathematics is **correct modulo one minor, repairable revision**. 16/16 proof units
checked. **No S0, no S1.** One S2 (a repairable gap in a standalone motivating proposition) and
four S3 nits. Several results (S3 semisimplicity, S5 rate accounting, S6 efficiency) are
unusually careful.

## Stage outputs
| Stage | Verdict | Key artifacts |
|-------|---------|---------------|
| 1. proofcheck (6-pass) | Correct modulo repairs; 1×S2, 4×S3 | `audit/06_reports/FINAL_REPORT.md`, `audit/06_reports/issue_log.md`, `audit/0*/` |
| 2. proof-repair | 5/5 issues designed, all self-contained (0 new refs) | `REPAIR_PLAN.md`, `PATCHES.md`, `audit/07_repairs/P1_sufficiency_repair.md`, `repair_references.bib` |
| 3. proofcheck --post-repair | **CONVERGED** (recommended-not-gated, no S0/S1) | `audit/08_post_repair/CONVERGENCE_VERDICT.md`, `diff_ledger.md`, … |
| 4. theory-sharpen | Frontier-level; 1 high-value strengthening (S5 $q_N^3\!\to\!q_N^2$) | `SHARPEN_REPORT.md` |
| 5. theory-simulation (AUDIT) | Strong + honest; 1 real gap (S5 has no experiment) | `simulation_audit/*` |
| 6. proof-writer | Rewritten Prop 1(ii): PROVABLE AS STATED / Verified (linter-clean) | `PROOF_PACKAGE.md` |

## The one substantive finding (I-01, S2)
Proposition 1(ii)'s Laplace-transform recovery of each within-cell measure $\nu_j$ silently
assumes the cell offsets $s_j=\bdelta'\bw_j$ can be varied independently — true only when
$J\lesssim k$, but the paper's regime is $J\gg k$. **Fix** (Patch 1 / `PROOF_PACKAGE.md`): run
the recovery in the cell-saturated parametrization (each cell carries a free effect — the native
collapse setting), where coordinatewise differencing in $s_j$ supplies the needed variation.
Claim preserved, no new assumption of consequence, no downstream impact (P1 is standalone).

## Cross-stage recommendations (priority order)
1. Apply Patch 1 (closes I-01) and the four S3 patches (`PATCHES.md`).
2. Add a **diverging-dimension simulation** for Theorem S5 (`simulation_audit/IMPROVEMENT_PLAN.md`
   E-new1) — the only theorem with no experiment, and the setting the register-scale motivation lives in.
3. (Optional, high value) sharpen S5's regime condition $q_N^3/N\to0$ to $q_N^2/N\to0$ via a
   self-concordance argument (`SHARPEN_REPORT.md`, Rank 1).
4. Minor reporting: MCSE/Wilson bands on the coverage/bias tables (Studies 3–5).

## Note
Where `theory-sharpen`/`theory-simulation` reference external literature, positioning is
knowledge-based (assistant cutoff Jan 2026); rows marked ▲ warrant a fresh citation check before
the authors act. Nothing in the proofcheck verdicts depends on those citations.
