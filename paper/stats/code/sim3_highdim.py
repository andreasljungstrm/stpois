"""
Proportional-regime study (Supplement, Section S4): J/N -> c.

Poisson model with one fixed effect per cell (absorbed / profiled out) and
p = 2 individual-level covariates.  Cells have bounded size m, so the
number of incidental parameters J = N/m grows proportionally with N:
J/N = 1/m in (0, 1).  Supplement Theorem S6 states that the profile MLE gamma-hat
remains consistent and asymptotically normal (no incidental-parameter bias,
by the coincidence of the Poisson profile likelihood with the conditional
multinomial likelihood), and that both the profile-information and the
cell-clustered sandwich variance estimators are consistent.

We report, over `REPS` Monte Carlo replications:
  * mean bias of gamma-hat_1,
  * sd of gamma-hat_1 across replications (the truth to be estimated),
  * mean of the profile-information SE and of the cluster SE,
  * empirical coverage of nominal 95% intervals from each,
  * for contrast, the coverage of the *naive* SE that ignores the fixed
    effects (fits gamma with a single intercept only) -- demonstrating that
    ignoring the incidental structure fails while profiling does not.

Output: code/output/sim3_highdim.csv
"""

import csv
import numpy as np

from tilted_glm import profile_poisson

REPS = 1000
GAMMA = np.array([0.30, -0.20])
Z975 = 1.959963984540054

rng = np.random.default_rng(20260716)


def one_rep(N, m):
    J = N // m
    cells = np.repeat(np.arange(J), m)
    alpha = rng.normal(-1.2, 0.8, J)            # heterogeneous cell effects
    X = np.column_stack([rng.normal(size=N),
                         rng.normal(size=N) + 0.5 * alpha[cells]])
    offset = np.log(0.5 + rng.random(N))
    eta = alpha[cells] + X @ GAMMA + offset
    y = rng.poisson(np.exp(eta)).astype(float)
    if np.bincount(cells, weights=y, minlength=J).min() == 0:
        # drop event-free cells (their alpha-hat = -inf; standard practice)
        keep_cell = np.bincount(cells, weights=y, minlength=J) > 0
        keep = keep_cell[cells]
        relab = np.cumsum(keep_cell) - 1
        cells2, X2, off2, y2 = relab[cells[keep]], X[keep], offset[keep], \
            y[keep]
    else:
        cells2, X2, off2, y2 = cells, X, offset, y
    g, Hp, Vp, Vcl = profile_poisson(y2, X2, cells2, off2)
    # naive contrast: ignore the cell effects entirely (single intercept)
    gn, sn = naive_poisson(y, X, offset)
    return g, np.sqrt(np.diag(Vp)), np.sqrt(np.diag(Vcl)), gn, sn


def naive_poisson(y, X, offset):
    """Poisson MLE with intercept + X, no cell effects (misspecified)."""
    A = np.column_stack([np.ones(len(y)), X])
    theta = np.zeros(3)
    theta[0] = np.log(y.mean()) - np.log(np.exp(offset).mean())
    for _ in range(50):
        mu = np.exp(A @ theta + offset)
        U = A.T @ (y - mu)
        H = (A * mu[:, None]).T @ A
        step = np.linalg.solve(H, U)
        theta += step
        if np.max(np.abs(step)) < 1e-10:
            break
    se = np.sqrt(np.diag(np.linalg.inv(H)))
    return theta[1:], se[1:]


def main():
    rows = []
    for (N, m) in [(20_000, 2), (20_000, 5), (20_000, 20),
                   (80_000, 2), (80_000, 5), (80_000, 20)]:
        est = np.empty((REPS, 2))
        se_p = np.empty((REPS, 2))
        se_c = np.empty((REPS, 2))
        est_n = np.empty((REPS, 2))
        se_n = np.empty((REPS, 2))
        for r in range(REPS):
            g, sp, sc, gn, sn = one_rep(N, m)
            est[r], se_p[r], se_c[r], est_n[r], se_n[r] = g, sp, sc, gn, sn
        bias = est.mean(axis=0) - GAMMA
        sd = est.std(axis=0)
        cov_p = np.mean(np.abs(est - GAMMA) <= Z975 * se_p, axis=0)
        cov_c = np.mean(np.abs(est - GAMMA) <= Z975 * se_c, axis=0)
        cov_n = np.mean(np.abs(est_n - GAMMA) <= Z975 * se_n, axis=0)
        bias_n = est_n.mean(axis=0) - GAMMA
        rows.append(dict(N=N, m=m, c=1.0 / m,
                         bias1=bias[0], sd1=sd[0],
                         mean_se_profile1=se_p[:, 0].mean(),
                         mean_se_cluster1=se_c[:, 0].mean(),
                         cover_profile1=cov_p[0], cover_cluster1=cov_c[0],
                         bias2=bias[1], sd2=sd[1],
                         cover_profile2=cov_p[1], cover_cluster2=cov_c[1],
                         naive_bias2=bias_n[1], cover_naive2=cov_n[1]))
        print(f"N={N:>6} m={m:>2} c={1/m:5.2f}  bias1={bias[0]:+.5f} "
              f"sd1={sd[0]:.5f} se_p={se_p[:,0].mean():.5f} "
              f"cov_p={cov_p[0]:.3f} cov_c={cov_c[0]:.3f}  "
              f"bias2={bias[1]:+.5f} cov_p2={cov_p[1]:.3f}  "
              f"naive: bias2={bias_n[1]:+.4f} cov2={cov_n[1]:.3f}")
    with open("output/sim3_highdim.csv", "w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=rows[0].keys())
        wr.writeheader()
        wr.writerows(rows)


if __name__ == "__main__":
    main()
