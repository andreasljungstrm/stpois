---
artifact: diff_ledger
scope: dependency_expanded
generated: 2026-07-18
generator: proofcheck --post-repair
---

# Assumption / Rate Diff Ledger

Compares the patched paper against the pre-patch audit baseline.

## Assumption diff
| Assumption | Pre-patch scope | Post-patch scope | Change | Justified? |
|---|---|---|---|---|
| Prop 1(ii) proof premise: independent variation of cell offsets $s_j$ | implicit / claimed "in general" (false for $J>k$) | explicit: cell-saturated parametrization (free per-cell effect) | CLARIFIED + NARROWED to the native collapse setting | **yes** — Assumption-Extension Change Log row (I-01); no downstream consumer |
| all other assumptions (A-canonical … A-proportional) | — | unchanged | none | n/a |

## Rate / constant diff
| Object | Pre-patch | Post-patch | Change |
|---|---|---|---|
| (all rates, constants, complexities) | — | identical | **none** |

## Probability level diff
None (no high-probability statement touched).

## Norm / metric diff
None.

## Sample-size regime diff
None (no asymptotic regime touched; S4/S5/S6 statements unchanged).

## Dependency requirement diff
None — no new edges; no lemma inserted; no reordering.

## Summary
Total diffs: **1** (the I-01 proof-premise clarification).
Justified: **1**.
**Unjustified or unpropagated: 0** ⇒ the diff ledger is clean.
