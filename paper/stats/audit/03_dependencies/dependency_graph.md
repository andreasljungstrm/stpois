---
artifact: dependency_graph
scope: global
generated: 2026-07-18
---

# Dependency Graph

(The mechanical indexer reported 0 edges because proofs cross-reference via `\eqref`
and hard-coded "Theorem 1 of the main paper" text rather than `\ref`; edges below are
reconstructed by reading each proof body.)

```
Def tilt ──► T1 (tilted representation) ──► C1 (complexity)
                     │                        └─► PS1, PS2 (penalized, Firth)
                     └─► P2 (non-canonical Fisher scoring)

P1 (impossibility)  : standalone motivation; nothing depends on it.

Lemma S1 (compact superlevel) ──► S1 (global convergence) ═ T2 (main)
S1 ──► S2 (local rate) ──► Cor S1 (bipartite) 
                    └────► S3 (multi-G rate)   [also uses Cor S1 for G=2 reduction]

Lemma S2 (profile=conditional) ──► S6 (proportional regime)

S4 (classical)   : reduction to Fahrmeir–Kaufmann 1985 + Gouriéroux et al 1984
S5 (diverging)   : self-contained (uses only A-canonical + A-diverging)
```

## Critical proof chains
1. **Engine chain (main contribution):** Def tilt → **T1** → C1. Load-bearing = T1.
2. **Absorption convergence chain:** Lemma S1 → **S1** (= main T2) → S2 → {Cor S1, S3}.
   Load-bearing = S1 and S2.
3. **Asymptotics:** S4 (delegated), **S5** (self-contained), Lemma S2 → **S6**.
   Headline of this chain = S6 (no incidental-parameter bias for Poisson).

## Circularity check (Pass 3)
No circular dependencies. Every chain terminates at base assumptions / external cited
results (Fahrmeir–Kaufmann; Smith 1977; Aronszajn–Kayalar–Weinert; Bertsekas/Tseng BCD;
Rockafellar 8.4; van der Vaart Ch. 25). Topological order is well-founded.
