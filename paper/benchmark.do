// benchmark.do -- verification and timing benchmarks for the stpois paper
// Run from repo root. Writes paper/benchmark_results.csv (versioned), which
// backs Tables 1-3 of the paper. Timings are the MEDIAN of `reps' runs; the
// CSV also records min and max so run-to-run dispersion can be reported.
discard
adopath ++ "`c(pwd)'"
set more off

local reps 10

// median/min/max timer for a command string: returns r(med), r(lo), r(hi)
capture program drop _medtime
program _medtime, rclass
    args reps cmd
    tempname T
    matrix `T' = J(`reps', 1, .)
    forvalues r = 1/`reps' {
        timer clear 1
        timer on 1
        quietly `cmd'
        timer off 1
        quietly timer list 1
        matrix `T'[`r', 1] = r(t1)
    }
    preserve
        clear
        quietly svmat double `T'
        quietly summarize `T'1, detail
        return scalar med = r(p50)
        return scalar lo  = r(min)
        return scalar hi  = r(max)
    restore
end

tempname bres
file open `bres' using paper/benchmark_results.csv, write replace text
file write `bres' "section,n_records,flavor,reps,streg_med,streg_lo,streg_hi,stpois_med,fast_med,fast_lo,fast_hi" _n
local flav = c(flavor)

// ── A. streg verification (stan3) ──────────────────────────────────────────
di _n "===== A. streg verification ====="
webuse stan3, clear
stset t1, failure(died) id(id)
qui streg age posttran, distribution(exponential) nohr nolog
local b_age_st  = _b[age]
local se_age_st = _se[age]
stpois age posttran, nolog
di "age: streg b=" %10.7f `b_age_st' " se=" %10.7f `se_age_st' ///
   "   stpois b=" %10.7f _b[age] " se=" %10.7f _se[age]
assert abs(_b[age]  - `b_age_st')  < 1e-6
assert abs(_se[age] - `se_age_st') < 1e-6
di "PASS: stpois == streg, dist(exponential)"

// ── B. Mixed model: continuous + categorical ────────────────────────────────
di _n "===== B. Mixed continuous + categorical ====="
capture program drop bench_dgp
program bench_dgp
    args n
    clear
    qui set obs `n'
    set seed 20260713
    qui gen edu    = ceil(4 * runiform())
    qui gen region = ceil(5 * runiform())
    qui gen period = ceil(10 * runiform())
    qui gen age    = 30 + 20 * runiform()
    qui gen income = rnormal(10, 2)
    qui gen id2    = _n
    qui gen lambda = exp(-6 + 0.03*age + 0.05*income ///
        + 0.3*(edu==2) + 0.5*(edu==3) - 0.2*(edu==4) ///
        + 0.1*(region==2) - 0.2*(region==3) + 0.2*(period >= 6))
    qui gen ptime  = 0.5 + 0.5 * runiform()
    qui gen failed = (runiform() < 1 - exp(-lambda * ptime))
    qui stset ptime, failure(failed) id(id2)
end

foreach n in 100000 1000000 5000000 {
    di _n "--- n = `n' ---"
    bench_dgp `n'
    qui streg age income i.edu i.region i.period, distribution(exponential) nolog
    local b_check = _b[age]
    qui stpois age income i.edu i.region i.period, nolog
    assert abs(_b[age] - `b_check') < 1e-6
    qui stpois age income i.edu i.region i.period, fast nolog
    assert abs(_b[age] - `b_check') < 1e-6

    _medtime `reps' "streg age income i.edu i.region i.period, distribution(exponential) nolog"
    local sr_med = r(med)
    local sr_lo  = r(lo)
    local sr_hi  = r(hi)
    _medtime `reps' "stpois age income i.edu i.region i.period, nolog"
    local sp_med = r(med)
    _medtime `reps' "stpois age income i.edu i.region i.period, fast nolog"
    local ft_med = r(med)
    local ft_lo  = r(lo)
    local ft_hi  = r(hi)
    di "n=`n'  streg=" `sr_med' "  stpois=" `sp_med' "  fast=" `ft_med'
    file write `bres' "B_mixed,`n',`flav',`reps',`=string(`sr_med',"%9.4f")',`=string(`sr_lo',"%9.4f")',`=string(`sr_hi',"%9.4f")',`=string(`sp_med',"%9.4f")',`=string(`ft_med',"%9.4f")',`=string(`ft_lo',"%9.4f")',`=string(`ft_hi',"%9.4f")'" _n
}

// ── C. Episode-split workflow ────────────────────────────────────────────────
di _n "===== C. Episode-split workflow (stsplit) ====="
clear
qui set obs 100000
set seed 20260713
qui gen edu    = ceil(4 * runiform())
qui gen region = ceil(5 * runiform())
qui gen age0   = 20 + 30 * runiform()
qui gen income = rnormal(10, 2)
qui gen id2    = _n
qui gen lambda = exp(-4.5 + 0.02*age0 + 0.03*income ///
    + 0.3*(edu==2) + 0.5*(edu==3) - 0.2*(edu==4)    ///
    + 0.1*(region==2) - 0.2*(region==3))
qui gen tfail  = -ln(runiform()) / lambda
qui gen failed = (tfail <= 30)
qui gen t_end  = min(tfail, 30)
qui stset t_end, failure(failed) id(id2)
qui stsplit fu, at(1(1)30)
qui count
local n_episodes = r(N)
di "After stsplit: `n_episodes' episode records"

qui streg age0 income i.edu i.region i.fu, distribution(exponential) nolog
local b_check = _b[age0]
qui stpois age0 income i.edu i.region i.fu, nolog
assert abs(_b[age0] - `b_check') < 1e-6
qui stpois age0 income i.edu i.region i.fu, fast nolog
assert abs(_b[age0] - `b_check') < 1e-6

_medtime `reps' "streg age0 income i.edu i.region i.fu, distribution(exponential) nolog"
local sr_med = r(med)
local sr_lo  = r(lo)
local sr_hi  = r(hi)
_medtime `reps' "stpois age0 income i.edu i.region i.fu, nolog"
local sp_med = r(med)
_medtime `reps' "stpois age0 income i.edu i.region i.fu, fast nolog"
local ft_med = r(med)
local ft_lo  = r(lo)
local ft_hi  = r(hi)
di "episodes: streg=" `sr_med' "  stpois=" `sp_med' "  fast=" `ft_med'
file write `bres' "C_episode_split,`n_episodes',`flav',`reps',`=string(`sr_med',"%9.4f")',`=string(`sr_lo',"%9.4f")',`=string(`sr_hi',"%9.4f")',`=string(`sp_med',"%9.4f")',`=string(`ft_med',"%9.4f")',`=string(`ft_lo',"%9.4f")',`=string(`ft_hi',"%9.4f")'" _n

// ── D. All-categorical predictors ────────────────────────────────────────────
di _n "===== D. All-categorical predictors ====="
foreach n in 1000000 5000000 {
    di _n "--- n = `n' ---"
    bench_dgp `n'
    qui streg i.edu i.region i.period, distribution(exponential) nolog
    local b_check = _b[3.edu]
    qui stpois i.edu i.region i.period, nolog
    assert abs(_b[3.edu] - `b_check') < 1e-6
    qui stpois i.edu i.region i.period, fast nolog
    assert abs(_b[3.edu] - `b_check') < 1e-6

    _medtime `reps' "streg i.edu i.region i.period, distribution(exponential) nolog"
    local sr_med = r(med)
    local sr_lo  = r(lo)
    local sr_hi  = r(hi)
    _medtime `reps' "stpois i.edu i.region i.period, nolog"
    local sp_med = r(med)
    _medtime `reps' "stpois i.edu i.region i.period, fast nolog"
    local ft_med = r(med)
    local ft_lo  = r(lo)
    local ft_hi  = r(hi)
    di "n=`n'  streg=" `sr_med' "  stpois=" `sp_med' "  fast=" `ft_med'
    file write `bres' "D_all_categorical,`n',`flav',`reps',`=string(`sr_med',"%9.4f")',`=string(`sr_lo',"%9.4f")',`=string(`sr_hi',"%9.4f")',`=string(`sp_med',"%9.4f")',`=string(`ft_med',"%9.4f")',`=string(`ft_lo',"%9.4f")',`=string(`ft_hi',"%9.4f")'" _n
}

file close `bres'
di _n "===== BENCHMARK DONE (results in paper/benchmark_results.csv) ====="
