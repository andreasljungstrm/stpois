"""
Application (main paper, Section 6.1): on-time performance of NYC
flights, 2013 -- a public-data validation of the tilted-moment framework
against a standard general-purpose implementation.

Data: the nycflights13 data set (Bureau of Transportation Statistics via
the R package nycflights13; 336,776 departures from EWR/JFK/LGA in
2013).  Model: logistic regression for a late arrival (arrival delay
>= 15 minutes, the FAA on-time criterion) with

  * cell structure: carrier x origin x destination x month (all observed
    combinations), entering through a main-effects cell-level design
    (intercept + carrier + origin + destination + month dummies);
  * individual-level continuous covariates: scheduled departure hour
    (centred) and log distance (standardised).

The tilted logistic fit is compared with statsmodels GLM on the
identical dense design: identical coefficients to ~1e-9 and the wall
clock ratio.  Expect data/flights.csv (see README for the download URL).

Output: code/output/app_flights.txt
"""

import csv
import time

import numpy as np

from tilted_glm import CellDesign, fit_tilted

DATA = "data/flights.csv"
URL = ("https://raw.githubusercontent.com/byuidatascience/"
       "data4python4ds/master/data-raw/flights/flights.csv")


def load():
    import os
    if not os.path.exists(DATA):
        os.makedirs("data", exist_ok=True)
        print(f"downloading {URL} ...")
        import urllib.request
        urllib.request.urlretrieve(URL, DATA)
    carrier, origin, dest, month = [], [], [], []
    hour, dist, late = [], [], []
    with open(DATA) as fh:
        rd = csv.DictReader(fh)
        for row in rd:
            if row["arr_delay"] in ("", "NA") or row["distance"] in ("", "NA"):
                continue
            carrier.append(row["carrier"])
            origin.append(row["origin"])
            dest.append(row["dest"])
            month.append(int(row["month"]))
            hour.append(float(row["hour"]) + float(row["minute"]) / 60.0)
            dist.append(float(row["distance"]))
            late.append(1.0 if float(row["arr_delay"]) >= 15.0 else 0.0)
    return (np.array(carrier), np.array(origin), np.array(dest),
            np.array(month), np.array(hour), np.array(dist),
            np.array(late))


def codes(a):
    uniq, inv = np.unique(a, return_inverse=True)
    return uniq, inv


def main():
    out = []

    def log(msg):
        print(msg)
        out.append(msg)

    t0 = time.perf_counter()
    carrier, origin, dest, month, hour, dist, late = load()
    N = len(late)
    log(f"nycflights13: N={N:,} flights with observed arrival delay; "
        f"{int(late.sum()):,} late arrivals "
        f"({100*late.mean():.1f}%); loaded in "
        f"{time.perf_counter()-t0:.1f}s")

    # drop separated factor levels: a dummy level whose flights are all
    # late or all on time sends its coefficient to +-infinity (the MLE
    # does not exist; Section 5.3 of the supplement).  In these data one
    # destination (LEX, a single on-time flight) is separated.
    changed = True
    while changed:
        changed = False
        for name in ("carrier", "origin", "dest", "month"):
            fac = {"carrier": carrier, "origin": origin,
                   "dest": dest, "month": month}[name]
            lv, iv = codes(fac)
            tot = np.bincount(iv)
            lat = np.bincount(iv, weights=late)
            bad = (lat == 0) | (lat == tot)
            if bad.any():
                keep = ~bad[iv]
                carrier, origin, dest, month = (v[keep] for v in
                                                (carrier, origin,
                                                 dest, month))
                hour, dist, late = hour[keep], dist[keep], late[keep]
                log(f"dropped {int((~keep).sum())} flights in separated "
                    f"{name} levels: {[str(l) for l in lv[bad]]}")
                changed = True

    uc, ic = codes(carrier)
    uo, io = codes(origin)
    ud, id_ = codes(dest)
    um, im = codes(month)
    nC, nO, nD, nM = len(uc), len(uo), len(ud), len(um)
    N = len(late)

    # cells: observed carrier x origin x dest x month combinations
    key = ((ic * nO + io) * nD + id_) * nM + im
    ukey, cells = np.unique(key, return_inverse=True)
    J = len(ukey)

    # cell-level main-effects design W (J x k)
    kc = ukey // (nO * nD * nM)
    rem = ukey % (nO * nD * nM)
    ko = rem // (nD * nM)
    rem = rem % (nD * nM)
    kd = rem // nM
    km = rem % nM
    k = 1 + (nC - 1) + (nO - 1) + (nD - 1) + (nM - 1)
    W = np.zeros((J, k))
    W[:, 0] = 1.0
    col = 1
    for lev in range(1, nC):
        W[kc == lev, col] = 1.0
        col += 1
    for lev in range(1, nO):
        W[ko == lev, col] = 1.0
        col += 1
    for lev in range(1, nD):
        W[kd == lev, col] = 1.0
        col += 1
    for lev in range(1, nM):
        W[km == lev, col] = 1.0
        col += 1

    X = np.column_stack([hour - hour.mean(),
                         (np.log(dist) - np.log(dist).mean())
                         / np.log(dist).std()])
    log(f"cells J={J:,} (carrier x origin x dest x month); "
        f"k={k} cell-level parameters "
        f"({nC} carriers, {nO} origins, {nD} destinations, {nM} months); "
        f"p=2 continuous (dep hour, log distance)")

    design = CellDesign(cells=cells, W=W, X=X, offset=np.zeros(N))
    t0 = time.perf_counter()
    fit = fit_tilted(late, design, family="binomial", tol=1e-9,
                     vce="robust")
    t_tilt = time.perf_counter() - t0
    log(f"tilted logistic fit: {t_tilt:.2f}s, {fit.iterations} Newton "
        f"iterations, converged={fit.converged}, ll={fit.loglik:.1f}")
    log(f"  dep hour (per hour) : {fit.gamma[0]:+.5f} "
        f"(robust se {fit.se[0]:.5f}) -> OR/hour "
        f"{np.exp(fit.gamma[0]):.4f}")
    log(f"  log distance (per sd): {fit.gamma[1]:+.5f} "
        f"(robust se {fit.se[1]:.5f})")

    # external check: statsmodels GLM on the identical dense design
    import statsmodels.api as sm
    A = np.column_stack([X, W[cells]])
    log(f"dense design for statsmodels: {A.shape[0]:,} x {A.shape[1]} "
        f"({A.nbytes/1e9:.2f} GB)")
    t0 = time.perf_counter()
    res = sm.GLM(late, A, family=sm.families.Binomial()).fit(tol=1e-9)
    t_sm = time.perf_counter() - t0
    ref = np.concatenate([fit.gamma, fit.delta])
    dev = np.max(np.abs(res.params - ref))
    log(f"statsmodels GLM: {t_sm:.2f}s; max |tilted - statsmodels| "
        f"= {dev:.2e}; speed ratio {t_sm/t_tilt:.0f}x")

    with open("output/app_flights.txt", "w") as fh:
        fh.write("\n".join(out) + "\n")


if __name__ == "__main__":
    main()
