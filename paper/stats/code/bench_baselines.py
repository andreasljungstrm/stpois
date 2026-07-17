"""
Benchmark study (main paper, Section 5.2): serious external competitors.

Each (method, configuration) is fitted in a fresh subprocess
(``_bench_worker.py``) so that peak resident memory isolates that
method's footprint.  We report, for the identical model and tolerance
($10^{-9}$): median wall-clock over repeated fits, peak RSS, iteration
count, and the maximum absolute deviation of the continuous-covariate
coefficients from the tilted fit (the exactness check).

Panel A (no absorbed factor): the cell structure enters as ordinary
covariates; competitors are dense Newton, ``statsmodels`` GLM (the
standard general-purpose Python IRLS), a sparse-design IRLS
(``scipy.sparse``), and ``glum`` (QuantCo's optimized GLM with the
``tabmat`` categorical-matrix backend -- among the fastest Python GLMs).
Panel A also sweeps the continuous dimension p in {2, 20, 50}.

Panel B (one high-dimensional absorbed factor): competitors are the
Poisson HDFE estimator ``pyfixest.fepois`` (the Python implementation of
the ``fixest`` alternating-projections algorithm), dense-dummy Newton,
and ``glum`` with dense dummies.

Output: output/bench_baselines_A.csv, bench_baselines_B.csv
"""

import csv
import json
import subprocess
import sys

WORKER = "_bench_worker.py"


def run(method, family, N, J, p, seed=20260716, reps=5, absorb=0):
    cmd = [sys.executable, WORKER, method, family, str(N), str(J),
           str(p), str(seed), str(reps), str(absorb)]
    try:
        out = subprocess.run(cmd, capture_output=True, text=True,
                             timeout=1200)
    except subprocess.TimeoutExpired:
        return None
    for line in out.stdout.splitlines()[::-1]:
        line = line.strip()
        if line.startswith("{"):
            try:
                obj = json.loads(line)
            except json.JSONDecodeError:
                continue
            if obj.get("status") == "error":
                obj["_failed"] = True
                return obj
            return obj
    # no JSON at all: subprocess died (timeout / OOM kill / segfault)
    return {"_failed": True, "reason": "timeout/oom",
            "msg": (out.stderr or "")[-120:]}


def _emit(rows, panel, family, N, J, p, m, r, ref):
    if r.get("_failed"):
        rows.append(dict(panel=panel, family=family, N=N, J=J, p=p,
                         method=m, time=None, peak_rss_mb=None, iters=None,
                         coef_dev=None, status="fail",
                         reason=r.get("reason", "?")))
        print(f"  {family:8s} p={p:<3} {m:12s} FAILED ({r.get('reason')})")
        return
    dev = (max(abs(a - b) for a, b in zip(r["coef"], ref))
           if ref else None)
    rows.append(dict(panel=panel, family=family, N=r["N"], J=r["J"], p=p,
                     method=m, time=r["time"], peak_rss_mb=r["peak_rss_mb"],
                     iters=r["iters"], coef_dev=dev, status="ok",
                     reason=""))
    print(f"  {family:8s} p={p:<3} {m:12s} t={r['time']:.3f}s "
          f"mem={r['peak_rss_mb']:.0f}MB iters={r['iters']} "
          f"dev={dev:.1e}")


def panel_A():
    rows = []
    # A1: full five-way comparison at p = 2 (Poisson and logistic)
    for family in ["poisson", "binomial"]:
        methods = ["tilted", "dense", "statsmodels", "sparse", "glum"]
        res = {m: run(m, family, 200_000, 200, 2) for m in methods}
        ref = res["tilted"]["coef"] if not res["tilted"].get("_failed") \
            else None
        for m in methods:
            _emit(rows, "A1", family, 200_000, 200, 2, m, res[m], ref)
    # A2: continuous-dimension sweep for the robust exact-MLE methods
    # (statsmodels' dense IRLS and glum's unpenalized solver are excluded:
    #  the former needs tens of seconds and gigabytes per fit, the latter
    #  errors with a singular matrix under alpha=0 for p >= 20)
    for p in [2, 20, 50]:
        methods = ["tilted", "dense", "sparse"]
        res = {m: run(m, "poisson", 200_000, 200, p) for m in methods}
        ref = res["tilted"]["coef"] if not res["tilted"].get("_failed") \
            else None
        for m in methods:
            _emit(rows, "A2", "poisson", 200_000, 200, p, m, res[m], ref)
    return rows


def panel_B():
    rows = []
    methods = ["tilted", "pyfixest", "dense", "glum"]
    for (N, absorb) in [(200_000, 500), (200_000, 5_000),
                        (500_000, 20_000)]:
        res = {m: run(m, "poisson", N, 100, 2, absorb=absorb,
                      reps=3) for m in methods}
        ref = res["tilted"]["coef"] if not res["tilted"].get("_failed") \
            else None
        for m in methods:
            r = res[m]
            if r.get("_failed"):
                rows.append(dict(panel="B", family="poisson", N=N, J=100,
                                 p=2, absorb=absorb, method=m, time=None,
                                 peak_rss_mb=None, iters=None,
                                 coef_dev=None, status="fail",
                                 reason=r.get("reason", "?")))
                print(f"  absorb={absorb:<6} {m:12s} FAILED "
                      f"({r.get('reason')})")
                continue
            dev = (max(abs(a - b) for a, b in zip(r["coef"], ref))
                   if ref else None)
            rows.append(dict(panel="B", family="poisson", N=r["N"], J=100,
                             p=2, absorb=absorb, method=m, time=r["time"],
                             peak_rss_mb=r["peak_rss_mb"], iters=r["iters"],
                             coef_dev=dev, status="ok", reason=""))
            print(f"  absorb={absorb:<6} {m:12s} "
                  f"t={r['time']:.3f}s mem={r['peak_rss_mb']:.0f}MB "
                  f"dev={dev:.1e}")
    return rows


def main():
    print("Panel A (no absorbed factor):")
    a = panel_A()
    print("Panel B (one high-dimensional absorbed factor):")
    b = panel_B()
    with open("output/bench_baselines_A.csv", "w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=a[0].keys())
        wr.writeheader()
        wr.writerows(a)
    with open("output/bench_baselines_B.csv", "w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=b[0].keys())
        wr.writeheader()
        wr.writerows(b)


if __name__ == "__main__":
    main()
