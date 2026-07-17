"""
Scaling study (main paper, Section 5.2): O(N) time and memory in the
sample size, and the memory gap to a dense design.

The tilted fit is run at N from 1e6 to 5e7 (Poisson, J = 500 cells,
p = 2), each in a fresh subprocess so peak resident memory is isolated.
We report wall-clock, peak RSS, and -- for contrast -- the size a dense
N x (p + J) double-precision design matrix would occupy (which a dense or
IRLS solver must form and, in the sparse case, whose Gram matrix it must
refactor each iteration).  The tilted pass forms no design matrix: its
footprint is the loaded columns O(N (p + 1)) plus O(J (k^2 + kp)) cell
accumulators.

Output: output/scaling.csv
"""

import csv
import json
import subprocess
import sys


def run_tilted(N, J=500, p=2):
    cmd = [sys.executable, "_bench_worker.py", "tilted", "poisson",
           str(N), str(J), str(p), "20260716", "1", "0"]
    out = subprocess.run(cmd, capture_output=True, text=True, timeout=2400)
    for line in out.stdout.splitlines()[::-1]:
        if line.strip().startswith("{"):
            return json.loads(line)
    raise RuntimeError(out.stderr[-500:])


def main():
    rows = []
    J, p = 500, 2
    for N in [1_000_000, 5_000_000, 20_000_000, 50_000_000]:
        r = run_tilted(N, J, p)
        dense_gb = N * (p + J) * 8 / 1e9
        rows.append(dict(N=N, J=r["J"], p=p, time=r["time"],
                         peak_rss_mb=r["peak_rss_mb"],
                         peak_rss_gb=r["peak_rss_mb"] / 1024,
                         dense_design_gb=dense_gb,
                         iters=r["iters"],
                         mem_ratio=dense_gb * 1024 / r["peak_rss_mb"]))
        rr = rows[-1]
        print(f"N={N:>11,} t={rr['time']:.2f}s peak={rr['peak_rss_gb']:.2f}GB "
              f"iters={rr['iters']} | dense design would be "
              f"{dense_gb:.1f}GB ({rr['mem_ratio']:.0f}x)")
    # linearity check: time and (memory - baseline) per million rows
    print("\ntime per 1e6 rows:",
          [round(r["time"] / (r["N"] / 1e6), 3) for r in rows])
    with open("output/scaling.csv", "w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=rows[0].keys())
        wr.writeheader()
        wr.writerows(rows)


if __name__ == "__main__":
    main()
