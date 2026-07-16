"""
Simulation study 2 (Section 6.2): linear convergence rate of the weighted
alternating projections and its spectral prediction (Theorem 3).

Two absorbed factors (J1 "regions" and J2 "cohorts") with an assignment
mechanism whose mixing parameter rho controls how strongly the factors are
associated: with probability rho a unit's cohort is deterministically tied
to its region block, otherwise uniform.  As rho -> 1 the weighted bipartite
graph of the two factors becomes nearly disconnected, the Friedrichs cosine
approaches 1, and the demeaning slows down exactly as predicted.

For each rho we
  * compute the predicted asymptotic rate sigma_2^2 from the weighted
    bipartite transition kernel (Corollary 2), and
  * measure the empirical per-sweep contraction factor of the weighted
    alternating projections on the IRLS working variable at the Poisson
    optimum weights.

Output: code/output/sim2_convergence.csv
"""

import csv
import numpy as np

from tilted_glm import demean_ap, friedrichs_rate

rng = np.random.default_rng(20260716)


def gen_factors(N, J1, J2, rho):
    f1 = rng.integers(0, J1, size=N)
    # block-aligned cohort: cohort index tied to region block when aligned
    aligned = rng.random(N) < rho
    f2 = np.where(aligned,
                  (f1 * J2) // J1,        # deterministic tie
                  rng.integers(0, J2, size=N))
    return f1, f2


def empirical_rate(history, tail=10):
    h = np.asarray(history)
    h = h[h > 1e-13]
    if len(h) < tail + 2:
        tail = max(2, len(h) - 2)
    ratios = h[-tail:] / h[-tail - 1:-1]
    return float(np.mean(ratios))


def main():
    N, J1, J2 = 200_000, 100, 80
    rows = []
    for rho in [0.0, 0.3, 0.6, 0.8, 0.9, 0.95, 0.99]:
        f1, f2 = gen_factors(N, J1, J2, rho)
        # Poisson-style weights at a "converged" fit: mu_i from a synthetic
        # linear predictor with both effects present
        a1 = rng.normal(0, 0.4, J1)
        a2 = rng.normal(0, 0.4, J2)
        eta = -3.0 + a1[f1] + a2[f2] + 0.3 * rng.normal(size=N)
        w = np.exp(eta)
        pred = friedrichs_rate(f1, f2, w)
        # demean a generic covariate; track weighted residual norms
        x = rng.normal(size=N) + 0.5 * a1[f1] - 0.3 * a2[f2]
        _, hist, iters = demean_ap(x, [f1, f2], w, tol=1e-13,
                                   maxiter=20000, return_history=True)
        emp = empirical_rate(hist)
        rows.append(dict(rho=rho, predicted=pred, empirical=emp,
                         sweeps=iters))
        print(f"rho={rho:4.2f} predicted={pred:.6f} empirical={emp:.6f} "
              f"sweeps={iters}")
    with open("output/sim2_convergence.csv", "w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=rows[0].keys())
        wr.writeheader()
        wr.writerows(rows)


if __name__ == "__main__":
    main()
