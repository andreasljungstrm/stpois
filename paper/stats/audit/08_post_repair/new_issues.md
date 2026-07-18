---
artifact: new_issues
scope: dependency_expanded
issue_ids: [N-01]
generated: 2026-07-18
---

# New-issue scan on touched units

Focused scan (not a full adversarial re-run) of the five patched passages, asking only:
did a patch introduce a defect absent from the original?

| ID | Severity | Unit | Description | Source patch |
|----|----------|------|-------------|--------------|
| N-01 | **S3** (non-blocking) | Thm S4 proof | Patch 4 proposes the labels "(N1)–(N3)" for the Fahrmeir–Kaufmann (1985) conditions; the labels are illustrative and should be reconciled with that paper's own condition numbering before submission. Does not affect correctness — the conditions themselves (divergence of $\lambda_{\min}$, normalization) are the right ones. | Patch 4 (I-04) |

## NEW-S0 / NEW-S1 scan
- New hidden assumption? Patch 1 introduces the "free per-cell offset" premise, but it is
  **documented** (Assumption-Extension Change Log, I-01) and consistent with the model's
  categorical structure ⇒ not a NEW-S0. No other patch adds an assumption.
- New quantifier / rate / constant defect? None — no patch touches a rate, constant,
  probability level, or norm.
- New circular dependency? None — no new lemma inserted; no new cross-unit edge.
- New notation drift? Patch 3 *removes* drift (I-03); introduces none.
- New cross-file `\ref` breakage (Mode B)? None — every patch stays within one file.

**NEW-S0: 0. NEW-S1: 0. NEW-S2: 0. NEW-S3: 1 (N-01, non-blocking documentation nit).**
