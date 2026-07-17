"""
Simulation study 5 (Section 6.5): the logit failure mode in the
proportional regime, and its conditional-likelihood repair.

Theorem 7 shows that the Poisson profile MLE is exactly unbiased at
J/N -> c because profiling coincides with conditioning.  For the logit
the two operations differ: the profile (unconditional fixed-effects) MLE
carries incidental-parameter bias of order 1/m -- for m = 2 the classic
result of Andersen (1970) gives plim gamma-hat = 2 * gamma -- while the
conditional MLE (conditioning on cell totals, Chamberlain 1980) removes
the bias entirely at the cost of a combinatorial denominator.

Both estimators are computed with cell-collapsed machinery:

  * `profile_logit`: concentrates one logit fixed effect per cell by a
    vectorised inner Newton (all cells simultaneously, O(N) per inner
    step) and takes outer Newton steps on the profile score with the
    Schur-complement Hessian -- the same tilted-moment algebra as the
    Poisson case;
  * `conditional_logit`: exact conditional likelihood; within-cell
    subset sums are enumerated per (cell size, total) group, vectorised
    across cells, feasible for the bounded cell sizes of the
    proportional regime (here m <= 8).

We report mean bias over REPS replications for m in {2, 4, 8} with
gamma_0 = 1, together with the Poisson profile bias on matched designs
(from Study 3) for contrast.

Output: code/output/sim5_logit_bias.csv
"""

import csv
from itertools import combinations

import numpy as np


def _sigmoid(z):
    return 1.0 / (1.0 + np.exp(-np.clip(z, -35.0, 35.0)))


REPS = 500
GAMMA0 = 1.0
rng = np.random.default_rng(20260716)


# ----------------------------------------------------------------------
# profile (unconditional fixed-effects) logit, one effect per cell
# ----------------------------------------------------------------------


def profile_logit(y, x, cells, tol=1e-9, maxiter=60):
    """Profile MLE of gamma with one logit fixed effect per cell.

    Cells with all-0 or all-1 outcomes are separated (alpha-hat = -+inf)
    and must be dropped by the caller.  The inner concentration solves
    every cell's one-dimensional MLE simultaneously by vectorised Newton;
    the outer step uses the Schur-complement profile Hessian.
    """
    J = int(cells.max()) + 1
    D = np.bincount(cells, weights=y, minlength=J)
    gamma = 0.0
    alpha = np.zeros(J)
    for _ in range(maxiter):
        # inner: alpha_j solves sum_i sigma(alpha_j + gamma x_i) = D_j
        for _ in range(80):
            mu = _sigmoid(alpha[cells] + gamma * x)
            f = np.bincount(cells, weights=mu, minlength=J) - D
            fp = np.bincount(cells, weights=mu * (1 - mu), minlength=J)
            step_a = f / np.maximum(fp, 1e-12)
            alpha -= step_a
            if np.max(np.abs(step_a)) < 1e-12:
                break
        mu = _sigmoid(alpha[cells] + gamma * x)
        w = mu * (1 - mu)
        U = float(np.sum((y - mu) * x))
        Q = float(np.sum(w * x * x))
        m1 = np.bincount(cells, weights=w * x, minlength=J)
        M0 = np.bincount(cells, weights=w, minlength=J)
        Hp = Q - float(np.sum(m1 * m1 / np.maximum(M0, 1e-12)))
        step = U / Hp
        gamma += step
        if abs(step) < tol * (1 + abs(gamma)):
            break
    return gamma


# ----------------------------------------------------------------------
# exact conditional logit for bounded cell sizes (Chamberlain 1980)
# ----------------------------------------------------------------------


def conditional_logit(y, x, cells, m, tol=1e-10, maxiter=60):
    """Conditional MLE of gamma given cell totals; cells all of size m.

    Cells are grouped by their total D in (0, m); for each group the
    C(m, D) within-cell subsets are enumerated once and all cells of the
    group are processed as one vectorised batch.
    """
    J = int(cells.max()) + 1
    order = np.argsort(cells, kind="stable")
    xs = x[order].reshape(J, m)
    ys = y[order].reshape(J, m)
    D = ys.sum(axis=1).astype(int)
    sx_obs = (ys * xs).sum(axis=1)
    groups = {}
    for d in range(1, m):
        sel = D == d
        if sel.any():
            combos = np.array([np.isin(np.arange(m), c)
                               for c in combinations(range(m), d)],
                              dtype=float)          # (n_combos, m)
            groups[d] = (xs[sel] @ combos.T, sx_obs[sel])  # subset sums
    gamma = 0.0
    for _ in range(maxiter):
        U, H = 0.0, 0.0
        for d, (S, sobs) in groups.items():
            # log-denominator: logsumexp over subsets of size d
            Z = gamma * S                            # (n_cells, n_combos)
            Zm = Z.max(axis=1, keepdims=True)
            wgt = np.exp(Z - Zm)
            wgt /= wgt.sum(axis=1, keepdims=True)
            Es = (wgt * S).sum(axis=1)               # E[sum_S x | D]
            Es2 = (wgt * S * S).sum(axis=1)
            U += float(np.sum(sobs - Es))
            H += float(np.sum(Es2 - Es * Es))
        step = U / H
        gamma += step
        if abs(step) < tol * (1 + abs(gamma)):
            break
    return gamma


def one_rep(N, m):
    J = N // m
    cells = np.repeat(np.arange(J), m)
    alpha = rng.normal(0.0, 1.0, J)
    x = rng.normal(size=N)
    eta = alpha[cells] + GAMMA0 * x
    y = (rng.random(N) < 1.0 / (1.0 + np.exp(-eta))).astype(float)
    D = np.bincount(cells, weights=y, minlength=J)
    keep_c = (D > 0) & (D < m)                # informative cells only
    keep = keep_c[cells]
    relab = np.cumsum(keep_c) - 1
    c2, x2, y2 = relab[cells[keep]], x[keep], y[keep]
    g_prof = profile_logit(y2, x2, c2)
    g_cond = conditional_logit(y2, x2, c2, m)
    return g_prof, g_cond


def main():
    rows = []
    for (N, m) in [(8_000, 2), (8_000, 4), (8_000, 8)]:
        gp = np.empty(REPS)
        gc = np.empty(REPS)
        for r in range(REPS):
            gp[r], gc[r] = one_rep(N, m)
        rows.append(dict(N=N, m=m, c=1.0 / m,
                         mean_profile=gp.mean(), sd_profile=gp.std(),
                         bias_profile=gp.mean() - GAMMA0,
                         mean_cond=gc.mean(), sd_cond=gc.std(),
                         bias_cond=gc.mean() - GAMMA0))
        r = rows[-1]
        print(f"N={N} m={m}  profile: mean={r['mean_profile']:.4f} "
              f"bias={r['bias_profile']:+.4f}   conditional: "
              f"mean={r['mean_cond']:.4f} bias={r['bias_cond']:+.4f}")
    with open("output/sim5_logit_bias.csv", "w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=rows[0].keys())
        wr.writeheader()
        wr.writerows(rows)


if __name__ == "__main__":
    main()
