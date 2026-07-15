// make_all.do -- full pipeline for the stpois package and paper
// Run from the repository root:  do paper/make_all.do
//
// Order: validation suite -> SJ log output -> timing benchmarks.
// After this completes, render both paper versions with:
//   bash paper/build.sh
//
// To regenerate only the SJ Stata log output (ex_*.log.tex):
//   do paper/sj/sj_examples.do

// sj_examples.do changes the working directory; save and restore
local repo_root = "`c(pwd)'"

do tests/test_stpois.do
do tests/test_fast_exact.do
do paper/sj/sj_examples.do    // does cd paper/sj internally
cd "`repo_root'"
do paper/benchmark.do
do paper/benchmark_hdfe.do
do paper/benchmark_ppmlhdfe.do

di _n "===== make_all: pipeline complete ====="
di "  - tests passed"
di "  - paper/sj/ex_*.log.tex regenerated"
di "  - paper/benchmark_results.csv written"
di "  - paper/benchmark_hdfe_results.csv written"
di "  - paper/benchmark_ppmlhdfe_results.csv written"
di ""
di "  Render both paper versions with:"
di "    bash paper/build.sh"
