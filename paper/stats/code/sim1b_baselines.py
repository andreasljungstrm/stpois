"""
Study: external and sparse baselines (main paper, Section 5.2), and
sensitivity to cell imbalance and rare events (Section 5.3).

Baselines for the same Poisson / logistic model (N = 200,000, J = 200
cells, p = 2 continuous covariates), all run to the same convergence
tolerance from comparable starting values:

  * tilted        -- tilted-moment Newton (this paper);
  * dense         -- textbook dense Newton on the full design (NumPy);
  * statsmodels   -- statsmodels GLM (IRLS on the dense design), the
                     standard general-purpose Python implementation;
  * sparse        -- IRLS with the design stored as a scipy.sparse CSR
                     matrix (continuous columns + cell dummies) and the
                     normal equations formed and solved sparsely: the
                     natural "just use sparse matrices" alternative.

We report total wall-clock time to convergence and the maximum absolute
coefficient deviation from the tilted fit.

The sensitivity panel re-runs the tilted-vs-dense comparison under (a)
heavily imbalanced cells (Zipf-distributed cell sizes) and (b) rare
events (Poisson mean ~ 0.002; logistic prevalence ~ 0.5%).

Output: code/output/sim1b_baselines.csv, sim1b_sensitivity.csv
"""

import csv
import time

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import statsmodels.api as sm

from tilted_glm import CellDesign, fit_full, fit_tilted, make_cell_dummies

rng = np.random.default_rng(20260716)


def gen(N, J, family, imbalance=False, rare=False):
    if imbalance:
        sizes = rng.zipf(1.6, size=J).astype(float)
        probs = sizes / sizes.sum()
        cells = rng.choice(J, size=N, p=probs)
    else:
        cells = rng.integers(0, J, size=N)
    X = rng.normal(size=(N, 2))
    base = -6.5 if rare else -3.0
    delta = np.concatenate([[base], rng.normal(0, 0.3, J - 1)])
    W = make_cell_dummies(cells, J)
    eta = (W @ delta)[cells] + X @ [0.3, -0.2]
    if family == "poisson":
        offset = np.log(0.5 + 0.5 * rng.random(N))
        y = rng.poisson(np.exp(eta + offset)).astype(float)
    else:
        offset = np.zeros(N)
        if not rare:
            eta = eta + 3.0
        y = (rng.random(N) < 1 / (1 + np.exp(-eta))).astype(float)
    # drop cells with no events (all baselines require the MLE to exist)
    ev = np.bincount(cells, weights=y, minlength=J)
    if family == "binomial":
        n_c = np.bincount(cells, minlength=J)
        ok = (ev > 0) & (ev < n_c)
    else:
        ok = ev > 0
    keep = ok[cells]
    relab = np.cumsum(ok) - 1
    cells, X, offset, y = relab[cells[keep]], X[keep], offset[keep], y[keep]
    J2 = int(ok.sum())
    return cells, X, offset, y, make_cell_dummies(cells, J2), J2


def run_baselines(N=200_000, J=200):
    rows = []
    for family in ["poisson", "binomial"]:
        cells, X, offset, y, W, J2 = gen(N, J, family)
        design = CellDesign(cells=cells, W=W, X=X, offset=offset)

        t0 = time.perf_counter()
        ft = fit_tilted(y, design, family=family, tol=1e-9)
        t_tilt = time.perf_counter() - t0
        ref = np.concatenate([ft.gamma, ft.delta])

        A = np.column_stack([X, W[cells]])
        t0 = time.perf_counter()
        ff = fit_full(y, A, offset, family=family, tol=1e-9)
        t_dense = time.perf_counter() - t0
        d_dense = np.max(np.abs(ff.gamma - ref))

        t0 = time.perf_counter()
        fam_sm = sm.families.Poisson() if family == "poisson" \
            else sm.families.Binomial()
        res = sm.GLM(y, A, family=fam_sm, offset=offset).fit(tol=1e-9)
        t_sm = time.perf_counter() - t0
        d_sm = np.max(np.abs(res.params - ref))

        # sparse IRLS: CSR design, sparse normal equations
        S = sp.hstack([sp.csr_matrix(X),
                       sp.csr_matrix((np.ones(len(y)),
                                      (np.arange(len(y)), cells)),
                                     shape=(len(y), J2)) @ sp.csr_matrix(W)
                       ]).tocsr()
        from tilted_glm import FAMILIES
        fam = FAMILIES[family]
        theta = np.zeros(S.shape[1])
        if family == "poisson":
            theta[2] = np.log(max(y.mean(), 1e-8) /
                              max(np.exp(offset).mean(), 1e-12))
        t0 = time.perf_counter()
        for it in range(100):
            eta = S @ theta + offset
            w = fam.weight(eta)
            resid = y - fam.mean(eta)
            U = S.T @ resid
            H = (S.T @ sp.diags(w) @ S).tocsc()
            step = spla.spsolve(H, U)
            theta = theta + step
            if np.max(np.abs(step) / (1 + np.abs(theta))) < 1e-9:
                break
        t_sparse = time.perf_counter() - t0
        d_sparse = np.max(np.abs(theta - ref))

        rows.append(dict(family=family, N=len(y), J=J2,
                         t_tilted=t_tilt, t_dense=t_dense,
                         t_statsmodels=t_sm, t_sparse=t_sparse,
                         dev_dense=d_dense, dev_statsmodels=d_sm,
                         dev_sparse=d_sparse))
        print(f"{family:9s} tilted={t_tilt:.2f}s dense={t_dense:.2f}s "
              f"statsmodels={t_sm:.2f}s sparse={t_sparse:.2f}s | "
              f"dev: dense={d_dense:.1e} sm={d_sm:.1e} "
              f"sparse={d_sparse:.1e}")
    with open("output/sim1b_baselines.csv", "w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=rows[0].keys())
        wr.writeheader()
        wr.writerows(rows)


def run_sensitivity(N=200_000, J=200):
    rows = []
    for family in ["poisson", "binomial"]:
        for label, imb, rare in [("balanced", False, False),
                                 ("imbalanced", True, False),
                                 ("rare", False, True),
                                 ("imbalanced+rare", True, True)]:
            cells, X, offset, y, W, J2 = gen(N, J, family,
                                             imbalance=imb, rare=rare)
            design = CellDesign(cells=cells, W=W, X=X, offset=offset)
            t0 = time.perf_counter()
            ft = fit_tilted(y, design, family=family, tol=1e-9)
            t_tilt = time.perf_counter() - t0
            A = np.column_stack([X, W[cells]])
            ffull = fit_full(y, A, offset, family=family, tol=1e-9)
            dev = np.max(np.abs(ffull.gamma
                                - np.concatenate([ft.gamma, ft.delta])))
            share = np.mean(y > 0) if family == "poisson" else y.mean()
            rows.append(dict(family=family, scenario=label, N=len(y),
                             J=J2, event_share=share,
                             iters=ft.iterations, converged=ft.converged,
                             dev=dev, t_tilted=t_tilt))
            print(f"{family:9s} {label:16s} J={J2:>4} "
                  f"events={share:.4f} iters={ft.iterations} "
                  f"conv={ft.converged} dev={dev:.1e} t={t_tilt:.2f}s")
    with open("output/sim1b_sensitivity.csv", "w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=rows[0].keys())
        wr.writeheader()
        wr.writerows(rows)


if __name__ == "__main__":
    run_baselines()
    run_sensitivity()
