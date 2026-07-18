---
artifact: theorem_inventory
scope: global
source_files: [main.tex, supplement.tex]
generated: 2026-07-18
generator: proofcheck (manual pipeline run, Codex disabled)
reference_mode: B (two-file, separate self-contained; supplement uses S-prefix counters)
---

# Theorem Inventory

Paper: *Tilted Cell Moments: A One-Pass Exact Newton Method for GLMs at Scale*
(main.tex) + online supplement (supplement.tex).

## Main text (proofs in Appendix A of main.tex)

| ID | Type | Loc (main.tex) | Short name | Claim summary | Depends on | Used by | Proof loc |
|----|------|------|------------|---------------|------------|---------|-----------|
| P1 `prop:sufficiency` | Prop | L406 | Failure of fixed collapse | (i) constant-within-cell ⇒ finite collapse exact; (ii) with continuous covariates, any statistic that evaluates the likelihood *function* determines the within-cell Laplace transform hence $\nu_j$ ⇒ no fixed-dim summary; (iii) Jensen strictness | b real-analytic + boundary $b(u)=e^u(1+o(1))$ | Motivation only (nothing downstream) | App A L1411 |
| T1 `thm:tilted` | Thm | L499 | Exact tilted-moment representation | Score & Fisher info at any θ are linear combinations of 3 tilted cell moments; Newton step identical to dense | eq (edf),(linpred); Def tilt | Cor1, Prop2, all supplement algs | App A L1470 |
| C1 `cor:complexity` | Cor | L573 | Scaling | Per-iter cost $O(N(p^2{+}p{+}1))+O(J(k^2{+}kp))+O((p{+}k)^3)$ | T1 | complexity claims throughout | App A L1495 |
| P2 `prop:noncanonical` | Prop | L625 | Tilted Fisher scoring | Non-canonical links: same representation with IRLS weight $\omega_i$, residual $r_i$; gamma–log ⇒ $\omega_i\equiv1$ | T1, eq (fisher-weight) | gamma–log claims | App A L1510 |
| T2 `thm:global-main` | Thm | L864 | Global convergence (absorbed factors) | Restatement of Supplement Thm S1 | = S1 | applications | Supp App S-A |

## Supplement (proofs in Appendices S-A … S-C)

| ID | Type | Loc (supp.tex) | Short name | Claim summary | Depends on | Used by | Proof loc |
|----|------|------|------------|---------------|------------|---------|-----------|
| S1 `thm:global` | Thm | L247 | Global convergence of nonlinear alternating projections | Algorithm S1 iterates increase ℓ monotonically & converge to MLE from any start, any cyclic order; (iv) extends to log-concave links + Fisher-scoring β-step | A-canonical, A-mle-exists, A-blocks; Lemma S1 | T2 (main) | L1112 |
| S2 `thm:rate` | Thm | L332 | Local linear convergence | Inner GS Jacobian = block-GS operator; $G{=}2$ rate $=c_F^2$; outer R-linear | S1 assumptions, $b\in C^3$ | Cor S1, S3 | L1191 |
| CS1 `cor:bipartite` | Cor | L360 | Bipartite spectral gap | $c_F=\sigma_2(D_1^{-1/2}BD_2^{-1/2})$; $c_F<1$ ⇔ graph connected | S2(ii) | rate intuition | L1238 |
| S3 `thm:multiG` | Thm | L406 | Exact rate & certificates, $G\ge2$ | (i) pairwise sufficiency; (ii) exact rate = max non-unit eig of sweep op, eig 1 semisimple; (iii) product-of-sines certificate; (iv) $G{=}2$ reduction | S2, Cor S1 | rate practice | L1258 |
| S4 `thm:classical` | Thm | L519 | Classical asymptotics | Fixed $q$: $\sqrt N(\hat\theta-\theta_0)\to N(0,\phi F_\infty^{-1})$, efficient; sandwich under mean-only correct | A-canonical, A-classical | asymptotics | L1355 (reduction to Fahrmeir–Kaufmann 1985 + Gouriéroux et al 1984) |
| S5 `thm:diverging` | Thm | L562 | Diverging dimension | $q_N^3/N\to0$, balanced leverage ⇒ rate $O_p(\sqrt{q_N})$ + contrast CLT | A-canonical, A-diverging | asymptotics | L1373 (self-contained) |
| LS2 `lem:conditional` | Lemma | L610 | Profile = conditional (Poisson) | Poisson profile ℓ = conditional multinomial ℓ; Schur-complement score/Hessian | Poisson model | S6 | L1455 |
| S6 `thm:proportional` | Thm | L656 | Proportional regime $J/N\to c$ | Poisson profile MLE consistent, $\sqrt J$-normal, **no asymptotic bias**, efficient; variance estimators consistent | A-proportional; LS2 | headline incidental-param result | L1481 |
| PS1 `prop:penalized` | Prop | L766 | Penalized proximal Newton | Ridge/lasso/elastic-net updates run on cell moments only | T1 | penalized use | L1569 |
| PS2 `prop:firth` | Prop | L800 | Tilted Firth | Firth adjusted score & leverages cell-decomposable at same order | T1 | separation remedy | L1588 |
| LS1 `lem:levelsets` | Lemma | L1090 | Compact superlevel sets | strictly concave + attains max ⇒ superlevel sets compact | — | S1 | L1096 |

## Sketch-vs-Complete classification (Pass 1 mandatory)

Every unit provides a **COMPLETE** proof (rigorous step-by-step derivation) except:
- **S4 (classical asymptotics)**: reduction-to-known-result. Proof delegates to
  Fahrmeir–Kaufmann (1985, Thms 1–2) and Gouriéroux–Monfort–Trognon (1984, Thm 3).
  Classified **COMPLETE-BY-REDUCTION** — acceptable, but the matching of the paper's
  Assumption S-classical to those papers' exact regularity conditions is asserted, not
  re-derived (see issue I-04, S3).
- **T2 (main-paper Thm 2)** is a restatement pointer to S1; the actual proof is S1.

No SKETCH-ONLY or PARTIAL-SKETCH units. No outstanding sketches. Sketches detected: 0.
