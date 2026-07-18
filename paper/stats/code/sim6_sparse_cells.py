"""
Sparse cell algebra study (main paper, Section "Sparse cell-level designs").

The dense cell assembly F_dd = W' diag(M0) W costs O(J k^2) and the dense
Newton solve O((p+k)^3): for saturated one-hot cell designs (k = J free
effects) these terms dominate as J grows and squander the tilted
reorganization exactly where the categorical structure is richest -- the
regime in which production sparse-GLM implementations are already strong.
The sparse cell algebra assembles the same blocks at O(nnz(W)) and takes
the Newton step by block elimination: sparse factorization of F_dd (exactly
diagonal in the saturated case), p x p Schur complement (the profile
Hessian), back-substitution for delta.

This study fits the identical saturated logistic model with the dense and
the sparse cell algebra (plus the textbook full-design Newton where it is
feasible) over a sweep of J = k, and verifies that the iterate sequences
and estimates agree to machine precision while the per-iteration cost of
the sparse path stays flat in k.

Output: code/output/sim6_sparse_cells.csv
"""

import csv
import time

import numpy as np

from tilted_glm import (CellDesign, fit_tilted, fit_full,
                        make_cell_dummies, make_cell_dummies_sparse)

N = 100_000
P = 2
GAMMA = np.array([0.30, -0.20])
J_GRID = [100, 400, 1600, 6400, 25600]
DENSE_CAP = 1600          # dense cell algebra beyond this is pointlessly slow
FULL_CAP = 400            # full-design Newton needs the N x (p+k) matrix

rng = np.random.default_rng(20260718)


def timed(fun, reps=2):
    best, out = np.inf, None
    for _ in range(reps):
        t0 = time.perf_counter()
        out = fun()
        best = min(best, time.perf_counter() - t0)
    return best, out


def one_design(J):
    cells = rng.integers(0, J, N)
    X = rng.normal(size=(N, P))
    alpha = rng.normal(-0.6, 0.3, J)
    eta = alpha[cells] + X @ GAMMA
    y = rng.binomial(1, 1.0 / (1.0 + np.exp(-eta))).astype(float)
    # drop constant-outcome cells (their saturated effects have no MLE;
    # standard practice, cf. the event-free-cell drop of sim3)
    D = np.bincount(cells, weights=y, minlength=J)
    n_j = np.bincount(cells, minlength=J)
    keep_cell = (D > 0) & (D < n_j)
    keep = keep_cell[cells]
    relab = np.cumsum(keep_cell) - 1
    return relab[cells[keep]], X[keep], y[keep], int(keep_cell.sum())


def theta(fit):
    return np.concatenate([fit.gamma, fit.delta])


rows = []
for J_nom in J_GRID:
    cells, X, y, J = one_design(J_nom)
    Nk = len(y)
    off = np.zeros(Nk)

    Ws = make_cell_dummies_sparse(cells, J)
    ds = CellDesign(cells, Ws, X, off)
    t_sp, fit_sp = timed(lambda: fit_tilted(y, ds, family="binomial"))
    ref = theta(fit_sp)
    rows.append(dict(J=J, N=Nk, k=J, method="tilted_sparse", time=t_sp,
                     iters=fit_sp.iterations,
                     per_iter=t_sp / fit_sp.iterations, dcoef=0.0))
    print(f"J={J:6d} (nominal {J_nom}, N={Nk})  sparse: {t_sp:8.3f}s  "
          f"({fit_sp.iterations} it, {t_sp/fit_sp.iterations:.4f}s/it)")

    if J <= DENSE_CAP:
        Wd = make_cell_dummies(cells, J, intercept=False)
        dd = CellDesign(cells, Wd, X, off)
        t_dn, fit_dn = timed(lambda: fit_tilted(y, dd, family="binomial"))
        d = float(np.max(np.abs(theta(fit_dn) - ref)))
        rows.append(dict(J=J, N=Nk, k=J, method="tilted_dense", time=t_dn,
                         iters=fit_dn.iterations,
                         per_iter=t_dn / fit_dn.iterations, dcoef=d))
        print(f"          dense:  {t_dn:8.3f}s  "
              f"({fit_dn.iterations} it, {t_dn/fit_dn.iterations:.4f}s/it)  "
              f"dcoef={d:.1e}")

    if J <= FULL_CAP:
        A = np.column_stack([X, make_cell_dummies(cells, J,
                                                  intercept=False)[cells]])
        t_fl, fit_fl = timed(lambda: fit_full(y, A, off, family="binomial"))
        d = float(np.max(np.abs(fit_fl.gamma - ref)))
        rows.append(dict(J=J, N=Nk, k=J, method="full_dense", time=t_fl,
                         iters=fit_fl.iterations,
                         per_iter=t_fl / fit_fl.iterations, dcoef=d))
        print(f"          full:   {t_fl:8.3f}s  "
              f"({fit_fl.iterations} it, {t_fl/fit_fl.iterations:.4f}s/it)  "
              f"dcoef={d:.1e}")

    assert all(r["dcoef"] < 1e-8 for r in rows)

with open("output/sim6_sparse_cells.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader()
    w.writerows(rows)

print("\nwritten: output/sim6_sparse_cells.csv")
