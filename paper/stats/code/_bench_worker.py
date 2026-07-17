"""
Subprocess worker for the benchmark harness (main paper, Section 5.2).

Fits ONE method on ONE generated data set in a fresh process, so that
peak resident memory (``ru_maxrss``) isolates that method's footprint
rather than the union of all competitors loaded in one interpreter.
Prints a single JSON line: time (median of `reps` timed fits), peak RSS
in MB, iteration count, and the two continuous-covariate coefficients
(for the exactness check the driver performs against the tilted fit).

Usage (invoked by bench_baselines.py / bench_hdfe.py):
    python3 _bench_worker.py <method> <family> <N> <J> <p> <seed> <reps>
                             [absorb_levels]

method in {tilted, dense, statsmodels, sparse, glum, pyfixest}.
"""

import json
import resource
import sys
import time

import numpy as np


def peak_rss_mb():
    # ru_maxrss is kilobytes on Linux, bytes on macOS
    r = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return r / 1024.0 if sys.platform != "darwin" else r / 1024.0 / 1024.0


def gen(N, J, p, family, seed, absorb_levels=0):
    rng = np.random.default_rng(seed)
    cells = rng.integers(0, J, size=N)
    X = rng.normal(size=(N, p))
    delta = np.concatenate([[-3.0 if family == "poisson" else 0.0],
                            rng.normal(0, 0.3, J - 1)])
    gamma = np.linspace(0.3, -0.2, p)
    eta = delta[cells] + X @ gamma
    offset = np.zeros(N)
    absorb = None
    if absorb_levels:
        absorb = rng.integers(0, absorb_levels, size=N)
        # raise the baseline so zero-event groups are rare and the MLE
        # exists cleanly for every competitor (dense dummies included)
        delta[0] = -1.0 if family == "poisson" else 0.0
        eta = delta[cells] + X @ gamma + rng.normal(0, 0.3,
                                                    absorb_levels)[absorb]
    if family == "poisson":
        offset = np.log(0.5 + 0.5 * rng.random(N))
        y = rng.poisson(np.exp(eta + offset)).astype(float)
    else:
        eta = eta + (3.0 if family == "binomial" else 0.0)
        y = (rng.random(N) < 1 / (1 + np.exp(-eta))).astype(float)

    # iteratively drop non-identified cells and absorbed groups (all
    # baselines need the MLE to exist / the design to be connected)
    def evmask(idx, K):
        e = np.bincount(idx, weights=y, minlength=K)
        if family == "binomial":
            nc = np.bincount(idx, minlength=K)
            return (e > 0) & (e < nc)
        return e > 0

    for _ in range(20):
        ok = evmask(cells, cells.max() + 1)
        keep = ok[cells]
        if absorb is not None:
            oka = evmask(absorb, absorb.max() + 1)
            keep = keep & oka[absorb]
        if keep.all():
            break
        cells = np.unique(cells[keep], return_inverse=True)[1]
        X, offset, y = X[keep], offset[keep], y[keep]
        if absorb is not None:
            absorb = np.unique(absorb[keep], return_inverse=True)[1]
    return cells, X, offset, y, int(cells.max()) + 1, absorb


def main():
    method, family = sys.argv[1], sys.argv[2]
    N, J, p, seed, reps = (int(sys.argv[i]) for i in range(3, 8))
    absorb_levels = int(sys.argv[8]) if len(sys.argv) > 8 else 0
    cells, X, offset, y, J2, absorb = gen(N, J, p, family, seed,
                                          absorb_levels)

    # peak RSS is a monotone high-water mark, so we read it right after a
    # single fit (before the timing repeats) to isolate one fit's
    # footprint rather than the union of repeated fits.
    single_peak = {}

    def timed(fit_once):
        coef, iters = fit_once()               # one fit for correctness
        single_peak["mb"] = peak_rss_mb()      # single-fit peak
        ts = []
        for _ in range(reps):
            t0 = time.perf_counter()
            fit_once()
            ts.append(time.perf_counter() - t0)
        return coef, iters, float(np.median(ts))

    import warnings
    warnings.filterwarnings("ignore")

    if method == "tilted":
        from tilted_glm import CellDesign, fit_tilted, make_cell_dummies
        W = make_cell_dummies(cells, J2)
        design = CellDesign(cells=cells, W=W, X=X, offset=offset)

        def f():
            r = fit_tilted(y, design, family=family,
                           absorb=[absorb] if absorb is not None else None,
                           tol=1e-9)
            return np.asarray(r.gamma[:p]), r.iterations
        coef, iters, t = timed(f)

    elif method == "dense":
        from tilted_glm import fit_full, make_cell_dummies
        W = make_cell_dummies(cells, J2)
        cols = [X, W[cells]]
        if absorb is not None:
            cols.append(np.eye(int(absorb.max()) + 1)[absorb][:, 1:])
        A = np.column_stack(cols)

        def f():
            r = fit_full(y, A, offset, family=family, tol=1e-9)
            return np.asarray(r.gamma[:p]), r.iterations
        coef, iters, t = timed(f)

    elif method == "statsmodels":
        import statsmodels.api as sm
        from tilted_glm import make_cell_dummies
        W = make_cell_dummies(cells, J2)
        A = np.column_stack([X, W[cells]])
        fam = (sm.families.Poisson() if family == "poisson"
               else sm.families.Binomial())

        def f():
            r = sm.GLM(y, A, family=fam, offset=offset).fit(tol=1e-9)
            return np.asarray(r.params[:p]), int(r.fit_history["iteration"])
        coef, iters, t = timed(f)

    elif method == "sparse":
        import scipy.sparse as sp
        import scipy.sparse.linalg as spla
        from tilted_glm import FAMILIES, make_cell_dummies
        W = make_cell_dummies(cells, J2)
        S = sp.hstack([sp.csr_matrix(X),
                       sp.csr_matrix((np.ones(len(y)),
                                      (np.arange(len(y)), cells)),
                                     shape=(len(y), J2)) @ sp.csr_matrix(W)
                       ]).tocsr()
        fam = FAMILIES[family]

        def f():
            theta = np.zeros(S.shape[1])
            if family == "poisson":
                theta[p] = np.log(max(y.mean(), 1e-8)
                                  / max(np.exp(offset).mean(), 1e-12))
            for it in range(100):
                eta = S @ theta + offset
                w = fam.weight(eta)
                U = S.T @ (y - fam.mean(eta))
                H = (S.T @ sp.diags(w) @ S).tocsc()
                step = spla.spsolve(H, U)
                theta = theta + step
                if np.max(np.abs(step) / (1 + np.abs(theta))) < 1e-9:
                    break
            return np.asarray(theta[:p]), it + 1
        coef, iters, t = timed(f)

    elif method == "glum":
        import pandas as pd
        from glum import GeneralizedLinearRegressor
        cols = {f"x{j}": X[:, j] for j in range(p)}
        cols["cell"] = pd.Categorical(cells)
        cols["y"] = y
        if offset.any():
            cols["off"] = offset
        df = pd.DataFrame(cols)
        rhs = " + ".join([f"x{j}" for j in range(p)] + ["C(cell)"])
        if absorb is not None:
            df["absorb"] = pd.Categorical(absorb)
            rhs += " + C(absorb)"
        formula = f"y ~ {rhs}"
        kw = dict(family=family, alpha=0.0, fit_intercept=True,
                  formula=formula, gradient_tol=1e-9)

        def f():
            m = GeneralizedLinearRegressor(**kw)
            m.fit(df, offset=offset if offset.any() else None)
            return np.asarray(m.coef_[:p]), int(np.max(m.n_iter_))
        coef, iters, t = timed(f)

    elif method == "pyfixest":
        import pandas as pd
        import pyfixest as pf
        d = {f"x{j}": X[:, j] for j in range(p)}
        d["y"] = y
        d["cell"] = cells
        if absorb is not None:
            d["absorb"] = absorb
        df = pd.DataFrame(d)
        xs = " + ".join(f"x{j}" for j in range(p))
        fe = "cell" + (" + absorb" if absorb is not None else "")
        formula = f"y ~ {xs} | {fe}"

        def f():
            if family == "poisson":
                r = pf.fepois(formula, data=df, iwls_tol=1e-9,
                              fixef_tol=1e-9, separation_check=[])
            else:
                r = pf.feglm(formula, data=df, family="logit")
            return np.asarray(r.coef().values[:p]), np.nan
        coef, iters, t = timed(f)

    else:
        raise SystemExit(f"unknown method {method}")

    print(json.dumps(dict(method=method, family=family, N=int(len(y)),
                          J=J2, p=p, absorb=absorb_levels,
                          time=t, peak_rss_mb=single_peak.get("mb",
                                                              peak_rss_mb()),
                          iters=(None if iters is None or
                                 (isinstance(iters, float)
                                  and np.isnan(iters)) else int(iters)),
                          coef=list(map(float, coef)))))


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:  # report a structured failure reason
        import json as _json
        etype = type(exc).__name__
        msg = str(exc)[:120]
        reason = ("singular" if "singular" in msg.lower()
                  else "oom" if isinstance(exc, MemoryError)
                  else etype)
        print(_json.dumps(dict(status="error", error=etype,
                               reason=reason, msg=msg)))
