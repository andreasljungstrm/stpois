# Theory Sharpening Report: *Tilted Cell Moments*

Codex independent assessment: **skipped by user instruction**. Literature positioning
below is knowledge-based (assistant cutoff Jan 2026), drawn from established T1 work;
citation counts are not asserted. A fresh Semantic-Scholar pass is recommended before
the authors act on any benchmarking row (flagged ▲).

## Executive Summary
- Assumptions analyzed: 7 blocks. Relaxable with real payoff: **2** (A-diverging leverage;
  A-diverging rate condition $q_N^3/N\to0$).
- Rates analyzed: the diverging-dimension rate (S5) and the proportional-regime rate (S6).
  Sharpenable: **1** (S5's dimension condition $q_N^3/N\to0 \Rightarrow q_N^2/N\to0$).
- Theory–model gaps: 0 substantive (theory matches the model it defines).
- Theory–experiment gaps: 0 contradictions; experiments track theory to 3 decimals.
- Top priorities: (1) sharpen S5's $q_N^3/N$ to $q_N^2/N$; (2) relax balanced-leverage to an
  effective-dimension condition; (3) add observable diagnostics for A-diverging/A-proportional.

## Framework Classification (presented for confirmation; proceeding at HIGH confidence)

*This paper is a multi-regime methodological+theory paper; the classification is unambiguous
from its own statements, so the pipeline proceeds without blocking. Correct any axis if desired.*

- **Axis 1 — Data structure:** `PANEL` (grouped / cell-structured data with absorbed
  high-dimensional fixed effects), with conditionally-`IID` GLM responses within the design.
  Evidence: cells = categorical cross-classification; absorbed factors §S1; cluster-robust VCE.
- **Axis 2 — Modeling framework:** `PAR` (parametric exponential-dispersion GLM) with a
  `SEMI` component in the proportional regime (Thm S6(iv) treats the $\alpha_j$ as unrestricted
  nuisance and attains the semiparametric efficiency bound).
- **Axis 3 — Asymptotic regime:** **MULTI** — `CLA` (S4, fixed $q$), a diverging-dimension
  regime `q_N^3/N→0` (S5), and the proportional/incidental-parameter regime `PROP` `J/N→c`
  (S6). The explicit `PROP` result is the theoretically distinctive one.
- **Cross-axis:** (PANEL, PAR/SEMI, PROP) is the incidental-parameters combination of the
  panel-econometrics line (Neyman–Scott; Andersen; Hahn–Newey; Fernández-Val–Weidner); the
  paper sits squarely in that lineage while adding the exact computational identity.

## Theory–Model–Experiment Alignment Matrix
| Assumption / Claim | Theory requires | Model provides | Experiments test | Alignment |
|---|---|---|---|---|
| Canonical / log-concave link | S1 | canonical (Poisson, logit) + gamma-log | all three tested | ALIGNED |
| Bounded covariates | S4,S5,S6 | synthetic bounded; flights bounded | yes | ALIGNED |
| Balanced leverage $\max h_i\le Cq_N/N$ | S5 | balanced cells give $C\approx1$ | Zipf-imbalance study stresses it | ALIGNED (imbalance tested, §5.6) |
| $J/N\to c$, bounded cells | S6 | synthetic $m\in\{2,5,20\}$ | Study 3 at $c=1/2$ | ALIGNED |
| Profile = conditional (Poisson) | S6 | exact for Poisson-log | Study 3 (Poisson unbiased) + Study 5 (logit biased) | ALIGNED — both sides shown |
| Exactness of Newton step | T1 | algebraic identity | $10^{-13}$ agreement | ALIGNED |
| Spectral rate $=c_F^2$ / multi-G | S2,S3 | weighted biadjacency | Study 2 (3-decimal match) | ALIGNED |

No `THEORY-MODEL` or `THEORY-EXP` gaps; one `EXPLOITABLE` item (below).

## Assumption Relaxation Opportunities

### Feasible (T1 technique exists)
| Assumption | Current form | Relaxation | Which step breaks | Rate change | Priority |
|---|---|---|---|---|---|
| A-diverging(iii) | $q_N^3/N\to0$ | $q_N^2/N\to0$ (up to logs) | S5 Step 4 remainder $O_p(q_N^{3/2}/\sqrt N)$ — the cube comes from a global third-order Taylor bound | none (same $\sqrt{q_N}$ rate) | **HIGH** |
| A-diverging(ii) | uniform $\max_i h_i\le Cq_N/N$ | average/effective-dimension leverage (e.g. $\sum_i h_i^2$ control) | S5 Step 2 uniform curvature | possibly extra log | MEDIUM |
| A-proportional(i) | $m_j\le\bar m$ (bounded cells) | slowly growing $\bar m_J\to\infty$ | S6(i) $\Var(L_J)\le C/J$ and Lyapunov moments | mild; needs $m$-uniform moment control | LOW–MEDIUM |

**Detail — sharpening $q_N^3/N\to0$ (the one genuinely worth doing).**
The cubic condition is the classical Portnoy (1988, *AoS*) threshold for M-estimation with
diverging dimension. The subsequent literature reduced it: He & Shao (2000, *AoS*) and,
via self-concordance / finite-sample parametric theory, Spokoiny (2012, *AoS*) obtain
asymptotic normality of GLM-type MLEs under roughly $q_N^2/N\to0$ (sometimes
$q_N\log q_N/N\to0$) by replacing the global third-order Taylor bound with a
self-concordant or leave-one-out expansion. Because the paper's canonical-link objective is
*exactly* self-concordant-like (the log-partition $b$ controls all higher derivatives), S5
Step 4 is the natural place to import this: the remainder can be bounded by a
self-concordance argument that costs $q_N^2/N$ rather than $q_N^3/N$. This is a clean
reviewer-friendly upgrade that widens the register-scale regime the paper targets (where
$k_N$ is in the thousands). ▲ verify current best threshold.

### Infeasible / correctly bounded
| Item | Why not relaxable |
|---|---|
| S1(iv) log-concave-link generality | The paper already stops exactly at the honest boundary (log-concave links); non-log-concave links lose the global claim — this is a *feature*, not a gap. |
| S6 Poisson-only no-bias | The mechanism (profile = conditional) is algebraic to Poisson; logit genuinely has $O(1/m)$ bias (Study 5). Not relaxable — correctly scoped. |

## Rate Sharpening Opportunities
| Result | Claimed | Best known / minimax | Optimal? | Sharpening |
|---|---|---|---|---|
| S4 classical | $\sqrt N$, efficient, Cramér–Rao | parametric efficiency | **YES — minimax optimal** (attains CR bound) | none (strength) |
| S6 proportional | $\sqrt J$, variance $\bar F^{-1}$, semiparametric-efficient | conditional-MLE efficiency bound | **YES — attains the semiparametric bound** (Thm S6(iv)) | none (strength) |
| S5 diverging | $O_p(\sqrt{q_N})$ in intrinsic norm under $q_N^3/N\to0$ | same rate under weaker $q_N^2/N$ | rate optimal; *condition* loose | tighten the regime, not the rate (above) |
| S2/S3 GS rate | $c_F^2$ / exact spectral | this IS the exact asymptotic rate | **YES — exact, not a bound** | (iii) product-of-sines bound is admittedly loose for $G\ge3$; (ii) exact rate is the right object — already stated |

**Minimax status:** the paper is unusually strong here — S4 and S6 attain classical and
semiparametric efficiency bounds respectively, and S2/S3 give the *exact* linear rate (not an
upper bound). There is essentially no rate slack to close; the only opportunity is the S5
*regime condition*.

## Reviewer-Critical Dimensions Audit
| Dimension | Addressed? | How | Gap | Suggestion |
|---|---|---|---|---|
| Lower bounds / optimality | **Yes** | CR bound (S4), semiparametric bound (S6(iv)), exact GS rate (S3(ii)) | none | — |
| Necessity of assumptions | Partial | P1 impossibility; A-blocks enforced by drop rule | I-01 gap (now repaired) | keep repaired P1 |
| Inference / UQ | **Yes** | sandwich + cluster-robust, cell-accumulated; valid CIs proven (S4,S6) | none | — |
| Identification | **Yes** | connectedness ⟺ $\sigma_2<1$ (Cor S1); spectral condition (S3(ii)) | none | — |
| Adaptivity / tuning-free | **Yes** | a-priori spectral rate `gs_spectral_rate`, no oracle tuning | none | — |
| Structural guarantees | **Yes** | exact recovery = MLE; separation drop rule | none | — |
| Computational attainability | **Yes (the paper's core)** | one-pass $O(N)$; benchmarks | none | — |
| Robustness to misspecification | **Yes** | PML/sandwich under mean-only correctness (S4); cluster-robust (S6(iii)) | none | — |
| Uniformity / honesty | Partial | S5 contrasts are uniform; stability battery honest about the $10^8$ pedestal | minor | — |
| **Assumption verifiability** | Partial | balanced-leverage & $\bar F_J\to\bar F$ are somewhat abstract | small | give observable diagnostics: report $\max_i h_i / (q_N/N)$ and the empirical $\bar F_J$ conditioning as fit output |

This paper scores strikingly well on the reviewer-critical axes — the usual "where is the
lower bound / inference / identification" demands are already met in-text.

## Exploitable Model Property
| Property | Theory uses | Could exploit | Improvement |
|---|---|---|---|
| Canonical link = exact self-concordance (log-partition $b$ bounds all derivatives) | generic $C^3$ + Lipschitz $b''$ (S5) | self-concordance | powers the S5 $q_N^3\to q_N^2$ sharpening above — the single highest-value item |

## Competitive Positioning (knowledge-based; ▲ = verify counts/venue before use)
| Work | Venue | Relation to this paper |
|---|---|---|
| Abowd–Kramarz–Margolis (1999) | *Econometrica* | two-way FE identification via connectedness — this paper generalizes to the exact $G\ge2$ spectral condition (S3) |
| Gaure (2013) | *Comput. Stat. Data Anal.* | alternating projections for linear FE — this paper adds the **nonlinear** combined-cycle global proof (S1) and exact rate (S2–S3) |
| Correia–Guimarães–Zylkin `ppmlhdfe` (2020) | *Stata J.* | Poisson-HDFE by IRLS demeaning, validated by simulation — this paper supplies the missing unified convergence theorem + exact tilted step |
| Bergé `fixest`/`FENmlm` (2018) | *CREA WP* ▲ | fast HDFE-GLM; benchmarked against (`pyfixest`) — this paper matches to $10^{-9}$ at lower memory |
| Andersen (1970); Hahn–Newey (2004); Fernández-Val–Weidner (2016) | *JRSS-B / Ectrica / J. Econometrics* | incidental-parameters — this paper's S6 is a clean profile=conditional (Poisson, no bias) with a $J/N\to c$ limit + variance theory; logit failure (Study 5) matches Andersen's $2\gamma_0$ |
| Sur–Candès (2019, *PNAS*); Candès–Sur (2020, *AoS*) ▲ | proportional-regime **logistic** MLE phase transition | complementary: their regime is exactly where the paper's Poisson coincidence *fails*; the paper correctly scopes S6 to Poisson and points to conditional logit as the repair |
| Portnoy (1988); He–Shao (2000); Spokoiny (2012) | *AoS* | diverging-dimension M-estimation — the frontier for relaxing S5's $q_N^3/N$ condition |

**Positioning verdict:** the paper is *ahead of or on par with* the frontier on every axis it
claims — exact identity (novel), unified nonlinear convergence (fills a real gap), exact
multi-factor rate (novel), and efficiency/inference (matches best known). Its one condition
that lags the current frontier is S5's $q_N^3/N$, which the diverging-dimension literature has
since improved.

## Improvement Roadmap (prioritized)
| Rank | Improvement | Type | Impact | Feas. | Lit | Align | Reviewer | Score | Ref |
|---|---|---|---|---|---|---|---|---|---|
| 1 | S5 regime $q_N^3/N\to q_N^2/N$ via self-concordance | Regime-extend | 4 | 4 | 4 | 3 | 4 | 3·4+2·4+2·4+1·3+1.5·4 =37 | Spokoiny (2012), He–Shao (2000) |
| 2 | Relax balanced-leverage to effective-dimension | Assumption-relax | 3 | 3 | 3 | 3 | 3 | 3·3+2·3+2·3+3+1.5·3=28.5 | Koltchinskii–Lounici (2017) |
| 3 | Observable diagnostics for A-div/A-prop | Verifiability | 2 | 5 | 3 | 4 | 3 | 3·2+2·5+2·3+4+1.5·3=30.5 | Crump et al. (2009) |
| 4 | Growing cell size $\bar m_J\to\infty$ in S6 | Regime-extend | 2 | 3 | 2 | 2 | 2 | 3·2+2·3+2·2+2+1.5·2=23 | Fernández-Val–Weidner (2016) |

## Detailed Improvement Spec — Rank 1 (S5 dimension condition)
- **Current:** S5 needs $q_N^3/N\to0$; remainder in Step 4 is $O_p(q_N^{3/2}/\sqrt N)$ from a
  global cubic Taylor term $\sum_i b'''(\eta_i)(\ba_i'\bs)^3$.
- **Target:** $q_N^2/N\to0$ (up to logs), doubling the admissible $k_N$ growth exponent.
- **Technique:** replace the global cubic bound by a self-concordance inequality for the
  canonical exponential family (the log-partition $b$ satisfies $|b'''|\le c\,b''$ on compacts),
  or a leave-one-out (Sur–Candès-style) expansion of $\hat\bs$.
- **Steps changed:** only S5 Step 2 (curvature) and Step 4 (remainder); Steps 1,3,5 unchanged.
- **Downstream:** none — S5 is a leaf (no theorem consumes it); the register-scale claim in the
  intro benefits (wider valid $k_N$).
- **Effort:** MEDIUM. Feed to `/proof-repair` as a *voluntary improvement*, then `/proof-writer`.

## Recommended Actions for Authors
- **Quick win:** report the two observable diagnostics (max-leverage ratio; $\bar F_J$
  conditioning) as fit output — closes the only verifiability gap cheaply (Rank 3).
- **Medium effort:** the self-concordant S5 upgrade (Rank 1) — the single highest-value
  theoretical strengthening, and reviewer-likely for the diverging-dimension claim.
- **Future work:** growing-cell-size proportional regime; effective-dimension leverage.

## New References (for the sharpening directions; not yet in the paper)
| Key | Citation | Venue (Tier) | Supports |
|---|---|---|---|
| spokoiny2012 | Spokoiny, "Parametric estimation. Finite sample theory" (2012) | *AoS* (T1) | S5 self-concordant sharpening |
| heshao2000 | He & Shao, "On parameters of increasing dimensions" (2000) | *AoS* (T1) | S5 $q_N^2/N$ threshold |
| portnoy1988 | Portnoy, "Asymptotic behavior of likelihood methods... dimension" (1988) | *AoS* (T1) | S5 baseline (already cited) |
| koltchinskii2017 | Koltchinskii & Lounici, "Concentration inequalities... effective rank" (2017) | *Bernoulli* (T1) | leverage relaxation |

(sharpen_references.bib mirrors these; heshao2000/portnoy1988 already appear in the paper's
own bibliography.)
