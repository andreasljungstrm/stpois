# LaTeX Patches — Apply in Order

## Reference Mode
- Mode: **B — two-file** (`main.tex` + `supplement.tex`).
- Supplement numbering: S1, S2, S3, … (S-prefixed counters).
- Every patch below is **within a single file** and introduces **no cross-file `\ref{}`**.
- These patches are a *proposal*; they are recorded here, not applied to the manuscript,
  so the author can review before editing. (No manuscript `.tex` file is modified by the
  pipeline.)

All changes are **STRUCTURAL / DOCUMENTATION EDITS** except Patch 1, which is a
claim-preserving Replace-Technique on a proof sub-argument (semantic status: claim
unchanged; see REPAIR_PLAN Repair Ladder Summary).

---

## Patch 1 — Prop 1(ii): repair the disentangling step (Issue I-01)
- File: `main.tex`
- Statement location: L429–430 (end of part ii). Proof location: L1436–1442 (Appendix A).

**1a. Statement — replace (L429–430):**
> The support size of $\nu_j$, and with it the dimension of any such statistic, grows with
> the number of distinct covariate values in the cell: no fixed-dimensional cell summary can
> support the likelihood function.

with:

> In the cell-saturated parametrization---each cell carrying its own free effect, as the
> collapse device supplies by construction---the support size of $\nu_j$, and with it the
> dimension of any such statistic, grows with the number of distinct covariate values in the
> cell: no fixed-dimensional cell summary can support the likelihood function.

**1b. Proof — replace (L1436–1442), the passage beginning "Provided the map $j\mapsto\bw_j$
is injective … one recovers each $G_j(\cdot,\gamma)$ on an open interval of $s$ for each
$\gamma$ in the slice." with:**

> Work in the cell-saturated parametrization, in which each cell $j$ carries a free scalar
> level $s_j$ (the categorical cross-classification that defines the cells supplies one effect
> per cell; a coarser cell-level design is a submodel of this one). Then $(s_1,\dots,s_J)$
> ranges over an open box $\prod_j I_j\subset\R^J$ as the cell-level parameters vary, and
> $\gamma$ ranges over an open $\gamma$-slice. Holding $\gamma$ fixed and differencing $\ell$
> in the single coordinate $s_j$ isolates $\partial\ell/\partial s_j = D_j -
> \partial_s G_j(s_j,\gamma)$; integrating in $s_j$ over $I_j$ recovers $G_j(\cdot,\gamma)$,
> up to the known additive $D_j$ term, on the open interval $I_j$, for each $\gamma$ in the
> slice.

(The remainder of the proof — analytic continuation of $G_j$ to $\R$, the boundary expansion
$e^{-s}G_j\to\mathcal L_{\nu_j}$, and Laplace-transform uniqueness — is unchanged and now
rests on a premise that holds in the paper's $J\gg k$ regime.)

- Repair class: Replace-Technique (L3). Claim preserved. No cross-file ref. No new `\cite`.
- Full corrected proof: `audit/07_repairs/P1_sufficiency_repair.md`.

---

## Patch 2 — Thm 1(i): positivity of the $b'$-tilt (Issue I-02)
- File: `main.tex`, L511–515 (Theorem 1 part (i), the line defining $\mathbf m^{b'}_j$).

**Change:** after "where $P^{b'}_{j,\btheta}$ is the tilt of Definition~\ref{def:tilt} with
$b'$ in place of $b''$" append:

> (well defined as a probability measure because $b'=\mu>0$ on the fitted range for the
> families in scope).

- Repair class: Notation-Fix (L1). Within-file, no ref change.

---

## Patch 3 — Prop 1(iii): disambiguate $\bar\bx_j$ (Issue I-03)
- File: `main.tex`, L431–439 (Prop 1 part iii) and L1464–1468 (its proof).

**Change (statement, L432–433):** replace "$e^{\bgamma'\bar\bx_j}$" context so $\bar\bx_j$ is
identified as the $\nu_j$-weighted mean:

> $\mathcal{L}_{\nu_j}(\bgamma) \ge \nu_j(\R^p)\, e^{\bgamma' \bar\bx_j^{\nu}}$, where
> $\bar\bx_j^{\nu} = \nu_j(\R^p)^{-1}\!\int \bx\,\nu_j(\dd\bx)$ is the $\nu_j$-weighted cell
> mean (which coincides with the plain cell mean of part (i) only when the offsets are
> constant within the cell),

and use $\bar\bx_j^{\nu}$ correspondingly at L1466 in the proof.

- Repair class: Notation-Fix (L1). Within-file.

---

## Patch 4 — Thm S4: name the Fahrmeir–Kaufmann conditions (Issue I-04)
- File: `supplement.tex`, L1357–1363 (proof of Theorem S4).

**Change:** replace "the normalizing regularity and divergence conditions of
\citet[Sec.~3]{fahrmeir1985} hold for the canonical-link exponential family with bounded
covariates" with an explicit statement of which conditions and why they hold:

> the divergence condition $\lambda_{\min}(\bF_N)\to\infty$ and the normalization
> $\bF_N^{-1/2}\bF_N(\btheta)\bF_N^{-1/2}\to\bI$ uniformly on shrinking $\bF_N$-neighborhoods
> of $\btheta_0$ — Fahrmeir--Kaufmann's (N1)–(N3) — follow from
> Assumption~\ref{ass:classical}(i)–(ii) (bounded covariates give the uniform third-derivative
> control; $N^{-1}\bF_N\to\bF_\infty\succ\bzero$ gives the normalization), so their Theorems~1–2
> apply.

- Repair class: Citation-Fix (L1). Within-file.

---

## Patch 5 — Thm S1: inner-loop / essential-cyclicity equivalence (Issue I-05)
- File: `supplement.tex`, L1167–1172 (proof of Theorem S1, part iii paragraph).

**Change:** after "all invariant to replacing the $\bbeta$-step by exact maximization and to
essentially-cyclic orderings (every block updated at least once every $B$ iterations)" append:

> In particular, running the inner $\balpha$-loop to tolerance before each outer $\bbeta$-step
> (as Algorithm~\ref{alg:absorb} does) is an instance of an essentially-cyclic ordering: each
> $\balpha$-block is updated at least once, and the $\bbeta$-block exactly once, per outer
> cycle, so the vanishing-per-block-increment argument applies verbatim with $B$ equal to the
> (finite, since the inner loop converges $R$-linearly by Theorem~\ref{thm:rate}(i)) inner-loop
> length.

- Repair class: Fill-Skipped-Steps (L1). Within-file, references Thm S2 (same file, `\ref` OK).

---

## Pre-patch validation
- [x] Every `\ref{}`/`\eqref{}` inside a patch resolves within the patch's target file
      (Patch 5's `\ref{thm:rate}` and `\ref{alg:absorb}` are both in `supplement.tex`).
- [x] No cross-file reference introduced by any patch.
- [x] No new `\cite{}` keys — `references.bib` unchanged; `repair_references.bib` is empty.
- [x] Supplement S-numbering untouched (no new supplement theorem-like environments added).
