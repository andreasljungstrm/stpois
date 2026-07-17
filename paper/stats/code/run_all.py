"""
Run every numerical study in the manuscript, in order.

    python3 run_all.py

Total runtime is roughly 15-30 minutes on a single laptop core, dominated
by the Monte Carlo loops.  Results land in ./output/; the tables in the
main paper (Sections 5-6) and the supplement (Section S4) are
transcriptions of these files.  Requires numpy and scipy; statsmodels
for the baseline study and the flights application (which also
downloads data/flights.csv on first run).
"""

import os
import runpy
import time

os.makedirs("output", exist_ok=True)

for script in ["sim1_validation.py", "sim1b_baselines.py",
               "sim2_convergence.py",
               "sim3_highdim.py", "sim4_separation.py", "sim5_logit_bias.py",
               "app_register.py", "app_flights.py"]:
    print(f"\n===== {script} =====")
    t0 = time.perf_counter()
    runpy.run_path(script, run_name="__main__")
    print(f"----- {script}: {time.perf_counter() - t0:.1f}s -----")
