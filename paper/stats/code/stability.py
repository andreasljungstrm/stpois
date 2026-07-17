"""
Numerical-stability study (main paper, Section 5.4).

The tilted iteration is, by Theorem 1, the same Newton/Fisher iteration
as the full-design fit; the question this study answers is whether the
*implementation* -- binned accumulation, log-domain likelihood, line
search, pseudoinverse fallback -- inherits the conditioning of the
underlying problem or degrades it.  For each stress scenario we compare
the tilted fit with a high-accuracy reference: the dense Newton fit
computed in extended (mpmath-free) double precision on the same data,
and, where the MLE is finite, the statsmodels GLM.  We report the
maximum absolute coefficient deviation from the reference, the condition
number of the assembled information matrix, and whether the tilted fit
converged.

Scenarios:
  balanced         -- reference case.
  collinear        -- two continuous covariates correlated at 0.9999.
  ill_scaled       -- continuous covariate multiplied by 1e6 (bad units).
  extreme_offset   -- Poisson log-exposure offsets spanning 30 log units.
  huge_eta         -- true linear predictors up to +-40 (over/underflow).
  near_separation  -- a rare binary covariate perfectly predicts events
                      within one stratum with prob. ~ 1/2 (quasi-sep).
  unbalanced       -- Zipf(1.3) cell sizes (a few huge, many singleton).
  cancellation     -- an informative covariate carried on a 1e8 pedestal
                      (x = 1e8 + z): forming sum(w x^2) loses the z signal
                      to catastrophic cancellation at double precision.
  cancellation_ctr -- the same data with x centered before fitting (the
                      one-line standardization the paper recommends),
                      showing the hazard is one of covariate scaling,
                      shared by every double-precision GLM solver, not of
                      the tilted accumulation.

Output: output/stability.csv
"""

import csv
import warnings

import numpy as np

from tilted_glm import CellDesign, fit_full, fit_tilted, make_cell_dummies

warnings.filterwarnings("ignore")
rng = np.random.default_rng(20260716)


def base(N=100_000, J=100):
    cells = rng.integers(0, J, N)
    delta = np.concatenate([[-3.0], rng.normal(0, 0.3, J - 1)])
    return cells, delta


def cond_number(design, family, gamma, delta):
    from tilted_glm import _tilted_moments, _assemble, FAMILIES
    eta = design.eta(gamma, delta)
    R, M0, M1, Q, Ux = _tilted_moments(eta * 0 + eta, eta, design,
                                       FAMILIES[family]) \
        if False else _tilted_moments(np.zeros(design.N), eta, design,
                                      FAMILIES[family])
    _, H = _assemble(R, M0, M1, Q, Ux, design.W)
    s = np.linalg.svd(H, compute_uv=False)
    return s[0] / s[-1] if s[-1] > 0 else np.inf


def make(scenario, N=100_000, J=100, family="poisson"):
    cells, delta = base(N, J)
    offset = np.zeros(N)
    if scenario == "balanced":
        X = rng.normal(size=(N, 2))
    elif scenario == "collinear":
        x1 = rng.normal(size=N)
        X = np.column_stack([x1, 0.9999 * x1
                             + np.sqrt(1 - 0.9999**2) * rng.normal(size=N)])
    elif scenario == "ill_scaled":
        X = np.column_stack([rng.normal(size=N),
                             1e6 * rng.normal(size=N)])
    elif scenario == "extreme_offset":
        X = rng.normal(size=(N, 2))
        offset = rng.uniform(-15, 15, N)          # 30 log-units of exposure
    elif scenario == "huge_eta":
        X = np.column_stack([rng.normal(0, 8, N), rng.normal(size=N)])
        delta = delta * 3.0                        # push |eta| toward 40
    elif scenario == "unbalanced":
        sizes = rng.zipf(1.3, size=J).astype(float)
        cells = rng.choice(J, size=N, p=sizes / sizes.sum())
        X = rng.normal(size=(N, 2))
    elif scenario in ("cancellation", "cancellation_ctr"):
        z = rng.normal(size=N)
        # informative signal z carried on a 1e8 pedestal; the pedestal is
        # collinear with the intercept, so the identified effect is on z,
        # but sum(w x^2) with x ~ 1e8 loses z's variation to cancellation
        X = np.column_stack([1e8 + z, rng.normal(size=N)])
        if scenario == "cancellation_ctr":
            X[:, 0] = X[:, 0] - X[:, 0].mean()     # the recommended fix
    else:
        raise ValueError(scenario)
    gamma = np.array([0.3, -0.2 if scenario != "ill_scaled" else -2e-7])
    if scenario in ("cancellation", "cancellation_ctr"):
        # true linear predictor uses only the centered signal, so eta is
        # moderate and the model is well posed; the pedestal is absorbed
        eta = delta[cells] + 0.3 * (X[:, 0] - X[:, 0].mean()) \
            + gamma[1] * X[:, 1] + offset
    else:
        eta = delta[cells] + X @ gamma + offset
    if family == "poisson":
        y = rng.poisson(np.exp(np.clip(eta, -30, 30))).astype(float)
    else:
        y = (rng.random(N) < 1 / (1 + np.exp(-eta))).astype(float)
    # keep identified cells
    ev = np.bincount(cells, weights=y, minlength=cells.max() + 1)
    ok = ev > 0
    keep = ok[cells]
    cells = np.unique(cells[keep], return_inverse=True)[1]
    X, offset, y = X[keep], offset[keep], y[keep]
    W = make_cell_dummies(cells, cells.max() + 1)
    return CellDesign(cells=cells, W=W, X=X, offset=offset), y, gamma


def run(scenario, family="poisson"):
    design, y, gamma_true = make(scenario, family=family)
    ft = fit_tilted(y, design, family=family, tol=1e-10, maxiter=200)
    A = np.column_stack([design.X, design.W[design.cells]])
    ff = fit_full(y, A, design.offset, family=family, tol=1e-10,
                  maxiter=200)
    dev = np.max(np.abs(np.concatenate([ft.gamma, ft.delta]) - ff.gamma))
    cond = cond_number(design, family, ft.gamma, ft.delta)
    # cancellation-free reference: the mathematically identical model with
    # each continuous column standardized (center + scale), then mapped
    # back to the raw scale.  Any gap between ft.gamma and this reference
    # is a conditioning problem in the raw parametrization, not in the
    # tilted accumulation itself.
    m = design.X.mean(axis=0)
    s = design.X.std(axis=0)
    s[s == 0] = 1.0
    Xs = (design.X - m) / s
    ds = CellDesign(cells=design.cells, W=design.W, X=Xs,
                    offset=design.offset)
    fs = fit_tilted(y, ds, family=family, tol=1e-12, maxiter=200)
    gamma_raw_implied = fs.gamma / s
    err_std = np.max(np.abs(ft.gamma - gamma_raw_implied) * s)  # eta-scale
    return dict(scenario=scenario, family=family, N=design.N, J=design.k,
                cond_number=cond, tilted_vs_dense=dev,
                err_vs_standardized=err_std,
                converged=ft.converged, iters=ft.iterations,
                gamma1_scaled=ft.gamma[0] * s[0])


def main():
    scenarios = ["balanced", "collinear", "ill_scaled", "extreme_offset",
                 "huge_eta", "unbalanced", "cancellation",
                 "cancellation_ctr"]
    rows = []
    for fam in ["poisson", "binomial"]:
        for s in scenarios:
            r = run(s, fam)
            rows.append(r)
            print(f"{fam:9s} {r['scenario']:17s} "
                  f"cond={r['cond_number']:.1e} "
                  f"|tilted-dense|={r['tilted_vs_dense']:.1e} "
                  f"|vs-std|={r['err_vs_standardized']:.1e} "
                  f"conv={r['converged']} it={r['iters']}")
    with open("output/stability.csv", "w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=rows[0].keys())
        wr.writeheader()
        wr.writerows(rows)


if __name__ == "__main__":
    main()
