# stpois

**Poisson event-history regression for `stset` data in Stata**

`stpois` fits Poisson event-history (exponential hazard) models directly on `stset` survival data, using the event indicator `_d` as outcome and time at risk `_t - _t0` as exposure. It reproduces `streg, distribution(exponential)` exactly on the default path, and adds two things that command does not have:

1. **High-dimensional fixed effects** — `absorb(varlist)` absorbs one or more sets of fixed effects via iteratively reweighted alternating projections (the `ppmlhdfe` approach), implemented in pure Mata with no dependencies. Coefficients and standard errors (OIM, robust HC1, clustered) are numerically exact relative to the equivalent dummy-variable model.
2. **Fast approximate estimation on collapsed data** — `fast(offset)` and `fast(moments)` collapse millions of individual spell records to a small number of categorical risk cells while retaining continuous covariates, either by folding them into a mutated exposure (two-stage offset) or by including within-cell moments as aggregation-bias corrections (CGF/Jensen expansion).

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

* Fast approximate estimation on collapsed cells
* bare names = continuous, i. terms = categorical (collapsed by)
stpois age i.surgery, fast(offset)
stpois age i.surgery, fast(moments)
```

See `help stpois` after installation for full syntax, the mathematical basis of the fast methods, and approximation warnings.

## Syntax

```stata
stpois varlist [if] [in] [weight] [, irr nolog level(#) noconstant
       vce(vcetype) robust cluster(varname)
       absorb(varlist) tolerance(#) maxiter(#)
       fast(offset|moments) skewness ]
```

- `fast()` and `absorb()` may not be combined.
- In the fast paths, bare variable names are treated as continuous and `i.`-prefixed variables as categorical grouping variables.
- Fast-path standard errors are conditional on first-stage estimates and do not propagate first-stage uncertainty.

## Files

| File | Purpose |
|---|---|
| `stpois.ado` | Main command: parsing, dispatch, ereturn, display |
| `stpois_p.ado` | `predict` after `stpois` (xb, n, hazard, survival) |
| `_stpois_hdfe.ado` | IRLS + alternating-projections FE absorption (Mata) |
| `_stpois_fast_offset.ado` | Two-stage multiplicative-offset fast path |
| `_stpois_fast_moments.ado` | CGF/Jensen moment-correction fast path |
| `stpois.sthlp` | Stata help file |
| `stpois.pkg`, `stata.toc` | `net install` metadata |
| `tests/test_stpois.do` | Validation suite (run from repo root) |
| `paper/` | Companion paper (Quarto/Typst source and PDF) |

## Testing

From the repository root, in Stata:

```stata
do tests/test_stpois.do
```

The suite checks (i) convergence of the standard path, (ii) that `absorb()` matches explicit dummy-variable Poisson on coefficients *and* standard errors to 1e-4, and (iii) that both fast paths recover the structural parameters on synthetic data.

## Citation

If you use `stpois`, please cite the companion paper (see `paper/`):

> Ljungström, A. (2026). stpois: Fast Poisson event-history regression with
> high-dimensional fixed effects in Stata. SocArXiv preprint.

## Author

Andreas Ljungström Swedish Institute for Social Research (SOFI), Stockholm University andreas.ljungstrom@sofi.su.se

## License

MIT — see [LICENSE](LICENSE).
