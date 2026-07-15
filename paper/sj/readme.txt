NOTE:  readme.txt template -- do not remove empty entries, but you may
                              add entries for additional authors
------------------------------------------------------------------------------

Package name:   <leave blank>

Title:  stpois: Fast Poisson event-history regression with
        high-dimensional fixed effects

Author 1 name:  Andreas Ljungström
Author 1 from:  Swedish Institute for Social Research (SOFI),
                Stockholm University, Stockholm, Sweden
Author 1 email: andreas.ljungstrom@sofi.su.se

Author 2 name:
Author 2 from:
Author 2 email:

Author 3 name:
Author 3 from:
Author 3 email:

Author 4 name:
Author 4 from:
Author 4 email:

Author 5 name:
Author 5 from:
Author 5 email:

Help keywords:  stpois, event-history analysis, survival analysis,
                Poisson regression, piecewise exponential model,
                high-dimensional fixed effects, register data

File list:      stpois.ado
                stpois_p.ado
                _stpois_hdfe.ado
                _stpois_fast.ado
                stpois.sthlp
                tests/test_stpois.do
                tests/test_fast_exact.do
                sj_examples.do
                benchmark.do
                benchmark_hdfe.do
                benchmark_results.csv
                benchmark_hdfe_results.csv

Notes:  Requires Stata 14 or later.  No dependencies beyond official
        Stata (both estimation engines are implemented in Mata).  The
        test suite (tests/) verifies that all three estimation paths
        reproduce streg/poisson coefficients and OIM, robust, and
        cluster-robust standard errors to 1e-6, including weighted
        estimation (fw/iw/pw), the combined fast + absorb() engine,
        and predictions after absorb().  All examples in the article
        use webuse stan3 or simulated data with a fixed seed
        (20260713), so no proprietary data are needed; sj_examples.do
        regenerates every log shown in the article via sjlog.  The
        benchmark scripts write the timing results reproduced in
        tables 1-4 to the two CSV files listed above.
