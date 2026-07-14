# stpois

**Poisson event-history regression for `stset` data in Stata**

`stpois` fits Poisson event-history (exponential hazard) models directly on `stset` survival data, using the event indicator `_d` as outcome and time at risk `_t - _t0` as exposure. It reproduces `streg, distribution(exponential)` exactly on the default path, and adds two things that command does not have — with **no approximations anywhere**: every path returns the exact MLE with exact OIM/robust/cluster standard errors.

1. **High-dimensional fixed effects** — `absorb(terms)` absorbs one or more sets of fixed effects (including interactions such as `absorb(region#period)`) via iteratively reweighted alternating projections (the `ppmlhdfe` approach), implemented in pure Mata with no dependencies. Each demeaning pass is a single O(n) sweep using precomputed sort permutations and group boundaries, independent of the number of levels. Coefficients and standard errors are numerically exact relative to the equivalent dummy-variable model: at 1M observations with a 1,000-level fixed effect, `absorb()` takes ~3s where explicit dummies take ~33s.
2. **Fast exact estimation on collapsed cells** — `fast` fits the exact individual-level MLE by Newton iterations on exponentially tilted within-cell moments: one O(n·p²) microdata pass per iteration (p = number of individual-level columns), all remaining algebra at the cell level. When every predictor is categorical, iterations run entirely on the collapsed (events, exposure) cell sums. Group structures come from a single sort with running-sum segment aggregation — no data copies, no collapse round-trips.

Factor-variable interactions are supported on all paths: continuous#continuous, continuous#categorical, categorical#categorical, and interacted fixed effects in `absorb()`.

## Installation

From Stata (requires Stata 14 or later):

```stata
net install stpois, from("https://raw.githubusercontent.com/andreasljungstrom/stpois/main/")
```

Or clone this repository and add it to your adopath:

```stata
adopath ++ "path/to/stpois"
```

## Quick start

```stata
webuse stan3, clear
stset t1, failure(died) id(id)

* Standard exponential-hazard Poisson (matches streg, dist(exponential))
stpois age posttran, irr

* Absorb fixed effects (exact coefficients and SEs)
stpois age posttran, absorb(surgery) vce(robust)

* Fast exact estimation on collapsed cells
* bare names = continuous, i. terms = categorical (define the cells)
stpois age i.surgery, fast

* Interactions work on every path
stpois c.age##i.posttran i.surgery, fast
stpois age, absorb(posttran#surgery)
```

See `help stpois` after installation for full syntax and the mathematical basis of the methods.

## Syntax

```stata
stpois varlist [if] [in] [weight] [, irr nolog level(#) noconstant
       vce(vcetype) robust cluster(varname)
       absorb(terms) tolerance(#) maxiter(#)
       fast ]
```

- `fast` and `absorb()` may not be combined.
- With `fast`, terms built only from factor-variable pieces enter at the cell level (and define the cells); terms involving a continuous piece enter at the individual level. At least one categorical term is required.
- `absorb()` accepts variable names and `var1#var2` interactions.

## Files

| File | Purpose |
|---|---|
| `stpois.ado` | Main command: parsing, dispatch, ereturn, display |
| `stpois_p.ado` | `predict` after `stpois` (xb, n, hazard, survival) |
| `_stpois_hdfe.ado` | IRLS + alternating-projections FE absorption (Mata) |
| `_stpois_fast.ado` | Cell-accelerated exact MLE (Mata Newton) |
| `stpois.sthlp` | Stata help file |
| `stpois.pkg`, `stata.toc` | `net install` metadata |
| `tests/` | Validation suite (run from repo root) |
| `paper/` | Companion paper (Quarto/Typst source and PDF) plus benchmark scripts |

## Testing

From the repository root, in Stata:

```stata
do tests/test_stpois.do
do tests/test_fast_exact.do
```

The suite checks that (i) the standard path matches `streg, distribution(exponential)`; (ii) `absorb()` matches explicit dummy-variable Poisson on coefficients *and* standard errors, including interacted fixed effects; (iii) `fast` matches the full MLE on coefficients, standard errors (OIM, robust, clustered), and log likelihood to 1e-6, for all interaction types and for all-categorical models.

## Citation

If you use `stpois`, please cite the companion paper (see `paper/`):

> Ljungström, A. (2026). stpois: Fast Poisson event-history regression with
> high-dimensional fixed effects in Stata. SocArXiv preprint.

## Author

Andreas Ljungström, Swedish Institute for Social Research (SOFI), Stockholm University — andreas.ljungstrom@sofi.su.se

## License

MIT — see [LICENSE](LICENSE).
