---
artifact: per_issue_closure
scope: dependency_expanded
theorem_ids: [prop:sufficiency, thm:tilted, thm:classical, thm:global]
issue_ids: [I-01, I-02, I-03, I-04, I-05]
generated: 2026-07-18
generator: proofcheck --post-repair (manual; Codex disabled)
---

# Per-issue closure verification

| Issue | Orig sev | Patch | Change type | Closure status | Notes |
|-------|----------|-------|-------------|----------------|-------|
| I-01 | S2 | Patch 1 | SEMANTIC (claim-preserving; scope of proof premise clarified) | **CLOSED-VERIFIED** | Re-verified: the reorganized proof isolates each cell via free per-cell offset $s_j$ (saturated parametrization), then recovers $\nu_j$ by the unchanged analytic-continuation + Laplace-uniqueness chain. The premise now holds in the $J\gg k$ regime. Scope narrowing to the saturated collapse setting is documented in the Assumption-Extension Change Log (per-issue file). |
| I-02 | S3 | Patch 2 | STRUCTURAL (parenthetical) | **CLOSED-VERIFIED** | Positivity $b'=\mu>0$ stated; the $b'$-tilt is now explicitly a probability measure. |
| I-03 | S3 | Patch 3 | STRUCTURAL (notation) | **CLOSED-VERIFIED** | $\bar\bx_j^{\nu}$ defined as the $\nu_j$-weighted mean; coincidence with the plain cell mean scoped to constant-offset case. |
| I-04 | S3 | Patch 4 | STRUCTURAL (citation detail) | **CLOSED-VERIFIED (minor residual)** | The invoked Fahrmeir–Kaufmann conditions are now named and tied to Assumption S-classical. Residual: the exact condition *labels* ("(N1)–(N3)") are proposed and should be checked against the source's own numbering before submission — see new_issues N-01 (S3, non-blocking). |
| I-05 | S3 | Patch 5 | STRUCTURAL (step fill) | **CLOSED-VERIFIED** | Inner-loop-to-tolerance is now explicitly framed as an essentially-cyclic ordering with finite $B$; uses Thm S2(i) R-linear inner convergence (same file, valid `\ref`). |

All five original issues reach a terminal closure status. No STILL-OPEN issue.
