"""
Simulation study 4 (Section 6.4): separation and Firth-adjusted tilted
moments.

Logistic model with a cell structure and a rare binary covariate z that is
quasi-separating with non-negligible probability in small samples: in cells
of one stratum, z = 1 occurs only among events.  The unadjusted MLE then
diverges (|gamma_z| -> infinity along the iterations, likelihood approaches
a finite supremum that is not attained); Firth's bias-reducing penalty,
computed inside the tilted-moment accumulation (Proposition 4), keeps the
maximiser finite and second-order unbiased.

We report, over REPS replications at each N:
  * fraction of replications flagged as separated for the MLE
    (|gamma_z-hat| > 10 or non-convergence),
  * median gamma_z-hat for MLE (on non-separated reps) and Firth (all reps),
  * median absolute error of the two estimators around the truth.

Output: code/output/sim4_separation.csv
"""

import csv
import numpy as np

from tilted_glm import CellDesign, fit_tilted, make_cell_dummies

REPS = 500
GAMMA_Z = 1.5      # true log-odds effect of the rare exposure
rng = np.random.default_rng(20260716)


def one_rep(N, J=10):
    cells = rng.integers(0, J, size=N)
    x = rng.normal(size=N)
    z = (rng.random(N) < 0.03).astype(float)     # rare exposure
    delta = np.concatenate([[-2.2], rng.normal(0, 0.3, J - 1)])
    W = make_cell_dummies(cells, J)
    eta = (W @ delta)[cells] + 0.3 * x + GAMMA_Z * z
    y = (rng.random(N) < 1 / (1 + np.exp(-eta))).astype(float)
    X = np.column_stack([x, z])
    design = CellDesign(cells=cells, W=W, X=X, offset=np.zeros(N))
    mle = fit_tilted(y, design, family="binomial", maxiter=100)
    fir = fit_tilted(y, design, family="binomial", maxiter=200,
                     tol=1e-6, firth=True)
    sep = (not mle.converged) or (abs(mle.gamma[1]) > 10)
    return mle.gamma[1], fir.gamma[1], sep, fir.converged


def main():
    rows = []
    for N in [300, 600, 1200, 2400]:
        g_mle, g_fir, seps, fconv = [], [], 0, 0
        for r in range(REPS):
            gm, gf, sep, fc = one_rep(N)
            seps += sep
            fconv += fc
            if not sep:
                g_mle.append(gm)
            g_fir.append(gf)
        g_mle, g_fir = np.array(g_mle), np.array(g_fir)
        rows.append(dict(
            N=N, sep_frac=seps / REPS, firth_conv_frac=fconv / REPS,
            med_mle=np.median(g_mle) if len(g_mle) else np.nan,
            med_firth=np.median(g_fir),
            mae_mle=np.median(np.abs(g_mle - GAMMA_Z)) if len(g_mle)
            else np.nan,
            mae_firth=np.median(np.abs(g_fir - GAMMA_Z))))
        r = rows[-1]
        print(f"N={N:>5} separated={r['sep_frac']:.3f} "
              f"firth-conv={r['firth_conv_frac']:.3f} "
              f"med(MLE|finite)={r['med_mle']:.3f} "
              f"med(Firth)={r['med_firth']:.3f} "
              f"MAE mle={r['mae_mle']:.3f} firth={r['mae_firth']:.3f}")
    with open("output/sim4_separation.csv", "w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=rows[0].keys())
        wr.writeheader()
        wr.writerows(rows)


if __name__ == "__main__":
    main()
