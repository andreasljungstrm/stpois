"""
Case study (Section 7): register-scale spatial epidemiology.

A fully reproducible, register-scale demonstration: N = 5 million
person-year episodes, mortality modelled with

  * 10,000 municipality x year fixed effects (500 municipalities x 20
    years), absorbed by the nonlinear Gauss--Seidel of Algorithm 1;
  * a cell-level design from education (4) x sex (2) x age band (12)
    = 96 categorical cells estimated via tilted moments;
  * two individual-level continuous covariates (standardised log income,
    and a pollution exposure index that varies within municipalities).

Rates, covariate effects, and the age gradient are calibrated to published
Nordic register-based mortality studies; the data generator is part of the
script, so every number in Section 7 of the paper can be reproduced
end-to-end on a laptop.  A dense dummy-variable fit would require a design
matrix with ~9,700 columns (about 0.4 TB of doubles at N = 5 x 10^6), so
the comparison fit is infeasible at full scale; exactness is instead
established on the subpopulation of the first 10 municipalities against
the dense fit.

Output: code/output/app_register.txt (timings, estimates, SEs).
"""

import time
import numpy as np

from tilted_glm import CellDesign, fit_tilted, fit_full

rng = np.random.default_rng(20260716)

N = 5_000_000
N_MUNI, N_YEAR = 500, 20
N_EDU, N_SEX, N_AGE = 4, 2, 12

# --- true parameters, loosely calibrated to Nordic all-cause mortality ---
BETA_INCOME = -0.18          # per SD of log income
BETA_POLLUT = 0.06           # per SD of exposure index
EDU_EFF = np.array([0.0, -0.12, -0.25, -0.40])
SEX_EFF = np.array([0.0, -0.45])                    # female advantage
AGE_EFF = np.linspace(0.0, 4.4, N_AGE)              # 35-39 ... 90+, Gompertz
BASE = -6.74                 # ~12 deaths / 1000 person-years overall


def generate():
    muni = rng.integers(0, N_MUNI, size=N)
    year = rng.integers(0, N_YEAR, size=N)
    fe_muni = rng.normal(0, 0.15, N_MUNI)
    trend = -0.012 * np.arange(N_YEAR)              # secular decline
    fe_year = trend + rng.normal(0, 0.03, N_YEAR)
    edu = rng.choice(N_EDU, size=N, p=[0.25, 0.35, 0.25, 0.15])
    sex = rng.integers(0, N_SEX, size=N)
    age = rng.choice(N_AGE, size=N,
                     p=np.array([12, 12, 11, 11, 10, 9, 8, 8, 7, 5, 4, 3])
                     / 100)
    income = rng.normal(-0.2 * (edu == 0) + 0.3 * (edu == 3), 1.0, N)
    pollut = 0.6 * fe_muni[muni] / 0.15 * 0.3 + rng.normal(0, 1, N)
    pollut = (pollut - pollut.mean()) / pollut.std()
    ptime = 0.5 + 0.5 * rng.random(N)
    eta = (BASE + fe_muni[muni] + fe_year[year] + EDU_EFF[edu]
           + SEX_EFF[sex] + AGE_EFF[age] + BETA_INCOME * income
           + BETA_POLLUT * pollut + np.log(ptime))
    y = rng.poisson(np.exp(eta)).astype(float)
    return muni, year, edu, sex, age, income, pollut, ptime, y


def build_design(edu, sex, age, income, pollut, ptime):
    cells = (edu * N_SEX + sex) * N_AGE + age
    Jc = N_EDU * N_SEX * N_AGE
    # cell-level design: intercept + edu, sex, age main-effect dummies
    W = np.zeros((Jc, 1 + (N_EDU - 1) + (N_SEX - 1) + (N_AGE - 1)))
    for j in range(Jc):
        e, rem = divmod(j, N_SEX * N_AGE)
        s, a = divmod(rem, N_AGE)
        W[j, 0] = 1.0
        if e > 0:
            W[j, e] = 1.0
        if s > 0:
            W[j, N_EDU - 1 + s] = 1.0
        if a > 0:
            W[j, N_EDU - 1 + N_SEX - 1 + a] = 1.0
    X = np.column_stack([income, pollut])
    return CellDesign(cells=cells, W=W, X=X,
                      offset=np.log(ptime)), cells


def main():
    out = []

    def log(msg):
        print(msg)
        out.append(msg)

    t0 = time.perf_counter()
    muni, year, edu, sex, age, income, pollut, ptime, y = generate()
    log(f"generated N={N:,} episodes, {int(y.sum()):,} deaths "
        f"(crude rate {y.sum()/ptime.sum()*1000:.2f}/1000 py) "
        f"in {time.perf_counter()-t0:.1f}s")

    # drop observations in zero-event absorbed groups (their fixed effects
    # diverge to -infinity: separation in the MLE sense) and relabel
    grp = muni * N_YEAR + year                      # 10,000 groups
    G0 = int(grp.max()) + 1
    ev = np.bincount(grp, weights=y, minlength=G0) > 0
    keep0 = ev[grp]
    n_drop = int(N - keep0.sum())
    relab0 = np.cumsum(ev) - 1
    muni, year, edu, sex, age = (v[keep0] for v in
                                 (muni, year, edu, sex, age))
    income, pollut, ptime, y = (v[keep0] for v in
                                (income, pollut, ptime, y))
    grp = relab0[grp[keep0]]
    log(f"dropped {n_drop:,} episodes in {int(G0 - ev.sum()):,} "
        f"zero-event municipality-year groups (separation)")

    design, cells = build_design(edu, sex, age, income, pollut, ptime)
    absorb = [grp]
    log(f"absorbed groups: {int(absorb[0].max())+1:,}; "
        f"cell-level cells: {design.J}; k={design.k}, p={design.p}")

    t0 = time.perf_counter()
    fit = fit_tilted(y, design, family="poisson", absorb=absorb,
                     tol=1e-8, vce="cluster", cluster=muni)
    t_fit = time.perf_counter() - t0
    log(f"tilted+absorbed fit: {t_fit:.1f}s, {fit.iterations} Newton "
        f"iterations, converged={fit.converged}, ll={fit.loglik:.1f}")
    names = (["income", "pollution"]
             + ["_cons"] + [f"edu{e}" for e in range(1, N_EDU)]
             + ["female"] + [f"age{a}" for a in range(1, N_AGE)])
    est = np.concatenate([fit.gamma, fit.delta])
    truth = np.concatenate([[BETA_INCOME, BETA_POLLUT],
                            [np.nan], EDU_EFF[1:], SEX_EFF[1:], AGE_EFF[1:]])
    for nm, b, se, tr in zip(names, est, fit.se, truth):
        log(f"  {nm:10s} {b:+.4f}  (se {se:.4f})   truth "
            f"{'' if np.isnan(tr) else f'{tr:+.4f}'}")

    # exactness cross-check against a dense dummy-variable fit, on the
    # subpopulation of the first 10 municipalities (the dense design for
    # the full data would need ~10,000 fixed-effect columns / ~0.4 TB)
    sel = muni < 10
    sub_grp = grp[sel]
    uniq, inv = np.unique(sub_grp, return_inverse=True)
    ysub = y[sel]
    dsub = CellDesign(cells=cells[sel], W=design.W,
                      X=design.X[sel], offset=design.offset[sel])
    Gsub = len(uniq)
    t0 = time.perf_counter()
    fsub = fit_tilted(ysub, dsub, family="poisson", absorb=[inv],
                      tol=1e-9)
    t_sub = time.perf_counter() - t0
    A = np.column_stack([dsub.X, dsub.W[dsub.cells][:, 1:],
                         np.eye(Gsub)[inv]])
    t0 = time.perf_counter()
    ffull = fit_full(ysub, A, dsub.offset, family="poisson", tol=1e-9)
    t_dense = time.perf_counter() - t0
    dev = np.max(np.abs(np.concatenate([fsub.gamma, fsub.delta[1:]])
                        - ffull.gamma[:dsub.p + dsub.k - 1]))
    log(f"subpopulation cross-check (n={sel.sum():,}, {Gsub:,} absorbed "
        f"groups): max |tilted - dense| = {dev:.2e}; "
        f"tilted {t_sub:.1f}s vs dense {t_dense:.1f}s "
        f"({t_dense/t_sub:.0f}x)")

    with open("output/app_register.txt", "w") as fh:
        fh.write("\n".join(out) + "\n")


if __name__ == "__main__":
    main()
