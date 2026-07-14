// benchmark.do -- verification and timing benchmarks for the stpois paper
// Run from repo root
discard
adopath ++ "`c(pwd)'"

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
    timer clear
    timer on 1
    qui streg age income i.edu i.region i.period, distribution(exponential) nolog
    timer off 1
    local b_check = _b[age]
    timer on 2
    qui stpois age income i.edu i.region i.period, nolog
    timer off 2
    assert abs(_b[age] - `b_check') < 1e-6
    timer on 3
    qui stpois age income i.edu i.region i.period, fast nolog
    timer off 3
    assert abs(_b[age] - `b_check') < 1e-6
    timer list
    di "n=`n'  streg=" r(t1) "  stpois=" r(t2) "  fast=" r(t3)
}

// ── C. Episode-split workflow ────────────────────────────────────────────────
// 100,000 subjects observed up to 30 years, split into 1-year episodes with
// a time-varying period effect: the canonical piecewise-exponential setup.
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

timer clear
timer on 4
qui stsplit fu, at(1(1)30)
timer off 4
qui count
di "After stsplit: " r(N) " episode records (stsplit took " r(t4) "s)"

timer on 1
qui streg age0 income i.edu i.region i.fu, distribution(exponential) nolog
timer off 1
local b_check = _b[age0]
timer on 2
qui stpois age0 income i.edu i.region i.fu, nolog
timer off 2
assert abs(_b[age0] - `b_check') < 1e-6
timer on 3
qui stpois age0 income i.edu i.region i.fu, fast nolog
timer off 3
assert abs(_b[age0] - `b_check') < 1e-6
timer list
di "episodes: streg=" r(t1) "  stpois=" r(t2) "  fast=" r(t3)

// ── D. All-categorical predictors ────────────────────────────────────────────
// When every predictor is categorical, fast iterates on the collapsed
// (events, exposure) cell sums only.
di _n "===== D. All-categorical predictors ====="
foreach n in 1000000 5000000 {
    di _n "--- n = `n' ---"
    bench_dgp `n'
    timer clear
    timer on 1
    qui streg i.edu i.region i.period, distribution(exponential) nolog
    timer off 1
    local b_check = _b[3.edu]
    timer on 2
    qui stpois i.edu i.region i.period, nolog
    timer off 2
    assert abs(_b[3.edu] - `b_check') < 1e-6
    timer on 3
    qui stpois i.edu i.region i.period, fast nolog
    timer off 3
    assert abs(_b[3.edu] - `b_check') < 1e-6
    timer list
    di "n=`n'  streg=" r(t1) "  stpois=" r(t2) "  fast=" r(t3)
}

di _n "===== BENCHMARK DONE ====="
