# Repair Plan: *Tilted Cell Moments* (main + supplement)

Generated from the /proofcheck audit (`paper/stats/audit/`) on 2026-07-18.
Codex per-repair stress-test: **skipped by user instruction** (recorded, not gated —
see Consistency Verification).

## Reference Mode
Mode: **B — two-file** (`main.tex` + `supplement.tex` compile to separate PDFs).
Numbering: main text 1,2,3…; supplement S1,S2,S3… (counters renewed with `S\arabic`).
Cross-file citation style: hard-coded numbers ("Theorem 1 of the main paper", "Lemma S1")
— never `\ref{}` across files. All patches below stay **within a single file**, so no
cross-file `\ref` is introduced.

## Executive Summary
- Total issues found: **5** (S2: 1, S3: 4). No S0, no S1.
- Repairable issues: **5 / 5**. All fixes are **self-contained** (no new literature).
- New references needed: **0**.
- Main theorem status after repair: **Preserved** (no theorem statement changes; the only
  content change clarifies the parametrization of a standalone motivating proposition).
- **Convergence status: NOT YET RE-AUDITED (S2/S3-only — re-audit recommended but not
  required).** Post-repair delta audit is run as Stage 3 for completeness.

## Repair Priority Order

### Phase 1 — the one substantive repair
| Issue | Unit | Repair Class | Ladder | Strategy | New refs |
|-------|------|--------------|--------|----------|----------|
| I-01 | Prop 1(ii) `prop:sufficiency` | Replace-Technique | L3 | Run per-cell Laplace recovery in the cell-saturated parametrization | none (self-contained) |

### Phase 2 — minor clarifications (S3)
| Issue | Unit | Repair Class | Ladder | Strategy | New refs |
|-------|------|--------------|--------|----------|----------|
| I-02 | Thm 1(i) `thm:tilted` | Notation-Fix | L1 | State $b'>0$ on the fitted range so the $b'$-tilt is a probability measure | none |
| I-03 | Prop 1(iii) `prop:sufficiency` | Notation-Fix | L1 | Mark $\bar\bx_j$ in (iii) as the $\nu_j$-weighted mean | none |
| I-04 | Thm S4 `thm:classical` | Citation-Fix | L1 | Name the Fahrmeir–Kaufmann (1985) divergence/normalization conditions invoked | none |
| I-05 | Thm S1 `thm:global` | Fill-Skipped-Steps | L1 | State the inner-loop-to-tolerance / essentially-cyclic equivalence explicitly | none |

## Repair Ladder Summary
| Issue ID | Unit | Chosen repair class | Ladder level | Claim preserved? | Assumptions preserved? | Escalation justified? | Pointer |
|----------|------|--------------------|--------------|------------------|------------------------|-----------------------|---------|
| I-01 | Prop 1(ii) | Replace-Technique | L3 | yes | yes (see repair note; L4 log also recorded) | NA (Phase A) | `audit/07_repairs/P1_sufficiency_repair.md` |
| I-02 | Thm 1(i) | Notation-Fix | L1 | yes | yes | NA | PATCHES Patch 2 |
| I-03 | Prop 1(iii) | Notation-Fix | L1 | yes | yes | NA | PATCHES Patch 3 |
| I-04 | Thm S4 | Citation-Fix | L1 | yes | yes | NA | PATCHES Patch 4 |
| I-05 | Thm S1 | Fill-Skipped-Steps | L1 | yes | yes | NA | PATCHES Patch 5 |

## Repair Closure Matrix
| Issue ID | Orig sev | Unit | Repair class | Patch ID | Touched units | Closure status | Post-repair status | Downstream affected |
|----------|----------|------|--------------|----------|---------------|----------------|--------------------|--------------------|
| I-01 | S2 | Prop 1(ii) | Replace-Technique | Patch 1 | Prop 1 stmt + App A proof | DESIGNED | (set by re-audit) | none |
| I-02 | S3 | Thm 1(i) | Notation-Fix | Patch 2 | Thm 1 stmt | DESIGNED | (set by re-audit) | none |
| I-03 | S3 | Prop 1(iii) | Notation-Fix | Patch 3 | Prop 1 stmt + App A proof | DESIGNED | (set by re-audit) | none |
| I-04 | S3 | Thm S4 | Citation-Fix | Patch 4 | App S-B proof | DESIGNED | (set by re-audit) | none |
| I-05 | S3 | Thm S1 | Fill-Skipped-Steps | Patch 5 | App S-A proof | DESIGNED | (set by re-audit) | none |

## Weaken-Claim Change Log
None — no repair weakens any claim. (I-01 preserves the claim; it repairs the argument.)

## Assumption-Extension Change Log
The only (belt-and-suspenders, strict-reading) Add-Assumption entry lives in the per-issue
file `audit/07_repairs/P1_sufficiency_repair.md`. The primary L3 reading adds no assumption.

## New References Summary
None. All repairs are self-contained; every technique invoked (real-analyticity, Laplace
uniqueness, Fahrmeir–Kaufmann conditions, block-coordinate essential cyclicity) is already
cited in the paper.

## Reference Quality Summary
- T1 references: n/a (no new references)
- Self-proved / clarified in place: 5 / 5

## Consistency Verification
- [x] Assumption matrix: no contradictions (only I-01 touches a hypothesis, and only for the
      standalone Prop 1; no downstream unit consumes it).
- [x] Rate propagation: no rate or constant changes anywhere; main theorems untouched.
- [x] New references: none, so trivially compatible.
- [x] Downstream units: none affected by any repair.
- [x] Repair Closure Matrix complete (all 5 issues have a row).
- [x] Weaken-Claim Change Log complete (empty — no Weaken-Claim repairs).
- [x] Assumption-Extension Change Log present in the I-01 per-issue file.
- [x] Every per-issue/patch has a ladder record; the one Phase-B-eligible item (I-01) carries a
      Phase A Exhaustion Record.
- [ ] Codex per-repair stress-test — **skipped by user instruction** (not gated for S2/S3-only).

## Residual Issues
None. All five issues are designed for closure.

## Hard-Gate Completion Rule status
Original audit had **no S0/S1**, so `/proofcheck --post-repair` is *strongly recommended but
not a hard gate*. It is run anyway (Stage 3). Outstanding sketches = 0. Consistency checklist
complete except the deliberately-skipped Codex step.
