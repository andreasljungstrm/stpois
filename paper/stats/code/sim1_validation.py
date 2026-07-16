"""
Simulation study 1 (Section 6.1): exactness and cost of tilted-moment Newton.

For Poisson, logistic, and gamma models with a cell structure of J cells and
p continuous covariates, we verify that

  (a) the tilted-moment Newton/Fisher iterates coincide with the full dense
      Newton iterates started from the same point, iteration by iteration,
      to machine precision (Theorem 1: the decomposition is exact, not an
      approximation); and
  (b) the measured per-iteration wall-clock cost tracks the leading-order
      FLOP ratio (p+k)^2 / p^2 of Corollary 1.

Output: code/output/sim1_validation.csv and a console summary.
"""

import time
import numpy as np

from tilted_glm import (CellDesign, fit_full, fit_tilted, make_cell_dummies,
                        flops_per_iteration)

rng = np.random.default_rng(20260716)


def gen_data(N, J, p, family):
    cells = rng.integers(0, J, size=N)
    X = rng.normal(size=(N, p))
    X[:, 0] = 30 + 10 * rng.random(N)          # an age-like covariate
    X[:, 0] = (X[:, 0] - X[:, 0].mean()) / X[:, 0].std()
    delta_true = np.concatenate([[-3.0], rng.normal(0, 0.3, J - 1)])
    gamma_true = np.linspace(0.2, -0.2, p)
    W = make_cell_dummies(cells, J)
    offset = np.log(0.5 + 0.5 * rng.random(N)) if family == "poisson" \
        else np.zeros(N)
    eta = (W @ delta_true)[cells] + X @ gamma_true + offset
    if family == "poisson":
        y = rng.poisson(np.exp(eta)).astype(float)
    elif family == "binomial":
        eta += 3.0                              # keep event rate moderate
        delta_true[0] += 3.0
        y = (rng.random(N) < 1 / (1 + np.exp(-eta))).astype(float)
    else:                                       # gamma, shape nu = 2
        nu = 2.0
        y = rng.gamma(nu, np.exp(eta + 3.0) / nu)
        delta_true[0] += 3.0
    return CellDesign(cells=cells, W=W, X=X, offset=offset), y


def run_one(family, N, J, p, reps_timing=3):
    design, y = gen_data(N, J, p, family)
    # (a) exactness: run both with iterate tracing from the same start
    fit_t = fit_tilted(y, design, family=family, tol=1e-10,
                       trace_iterates=True)
    A = np.column_stack([design.X, design.W[design.cells]])
    start = np.zeros(p + design.k)
    if family == "poisson":
        mu0 = max(np.mean(y), 1e-8) / max(np.mean(np.exp(design.offset)),
                                          1e-12)
        start[p] = np.log(mu0)
    fit_f = fit_full(y, A, design.offset, family=family, tol=1e-10,
                     start=start, trace_iterates=True)
    m = min(len(fit_t.trace), len(fit_f.trace))
    # both traces store the stacked (gamma, delta) vector in the same order
    max_dev = 0.0
    for a in range(m):
        max_dev = max(max_dev, np.max(np.abs(fit_t.trace[a]
                                             - fit_f.trace[a])))
    final_dev = np.max(np.abs(np.concatenate([fit_t.gamma, fit_t.delta])
                              - fit_f.gamma))

    # (b) timing of one Newton pass, averaged
    from tilted_glm import _tilted_moments, _assemble, FAMILIES
    fam = FAMILIES[family]
    eta = design.eta(fit_t.gamma, fit_t.delta)
    t_tilt = []
    for _ in range(reps_timing):
        t0 = time.perf_counter()
        R, M0, M1, Q, Ux = _tilted_moments(y, eta, design, fam)
        U, H = _assemble(R, M0, M1, Q, Ux, design.W)
        np.linalg.solve(H, U)
        t_tilt.append(time.perf_counter() - t0)
    etaf = A @ fit_f.gamma + design.offset
    t_full = []
    for _ in range(reps_timing):
        t0 = time.perf_counter()
        omega = fam.weight(etaf)
        resid = (y - fam.mean(etaf)) if fam.canonical else \
            omega * fam.working_residual(y, etaf)
        U = A.T @ resid
        H = (A * omega[:, None]).T @ A
        np.linalg.solve(H, U)
        t_full.append(time.perf_counter() - t0)
    fl_full, fl_tilt = flops_per_iteration(N, p, design.k)
    return dict(family=family, N=N, J=J, p=p, iters=fit_t.iterations,
                max_iterate_dev=max_dev, final_dev=final_dev,
                t_tilted=np.median(t_tilt), t_full=np.median(t_full),
                speedup=np.median(t_full) / np.median(t_tilt),
                flop_ratio=fl_full / fl_tilt)


def main():
    rows = []
    for family in ["poisson", "binomial", "gamma"]:
        for (N, J, p) in [(200_000, 50, 2), (200_000, 200, 2),
                          (200_000, 500, 2), (1_000_000, 200, 2),
                          (200_000, 200, 5)]:
            r = run_one(family, N, J, p)
            rows.append(r)
            print(f"{family:9s} N={N:>9,} J={J:>4} p={p} "
                  f"iter-dev={r['max_iterate_dev']:.2e} "
                  f"final-dev={r['final_dev']:.2e} "
                  f"t_full={r['t_full']:.3f}s t_tilt={r['t_tilted']:.3f}s "
                  f"speedup={r['speedup']:.1f}x flop-ratio="
                  f"{r['flop_ratio']:.1f}x")
    import csv
    with open("output/sim1_validation.csv", "w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=rows[0].keys())
        wr.writeheader()
        wr.writerows(rows)


if __name__ == "__main__":
    main()
