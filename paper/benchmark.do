// benchmark.do -- streg verification, mtopel acid test, and timing benchmarks
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
local b_pt_st   = _b[posttran]
local se_pt_st  = _se[posttran]
local b_c_st    = _b[_cons]
local se_c_st   = _se[_cons]

stpois age posttran, nolog
di "age:      streg b=" %10.7f `b_age_st' " se=" %10.7f `se_age_st' ///
   "   stpois b=" %10.7f _b[age] " se=" %10.7f _se[age]
di "posttran: streg b=" %10.7f `b_pt_st'  " se=" %10.7f `se_pt_st'  ///
   "   stpois b=" %10.7f _b[posttran] " se=" %10.7f _se[posttran]
di "_cons:    streg b=" %10.7f `b_c_st'   " se=" %10.7f `se_c_st'   ///
   "   stpois b=" %10.7f _b[_cons] " se=" %10.7f _se[_cons]
assert abs(_b[age]      - `b_age_st') < 1e-6
assert abs(_se[age]     - `se_age_st') < 1e-6
assert abs(_b[posttran] - `b_pt_st')  < 1e-6
assert abs(_se[posttran]- `se_pt_st') < 1e-6
assert abs(_b[_cons]    - `b_c_st')   < 1e-5
assert abs(_se[_cons]   - `se_c_st')  < 1e-5
di "PASS: stpois == streg, dist(exponential) on coefs and SEs"

// ── B. Murphy–Topel acid test ──────────────────────────────────────────────
// With mtopel, fast(offset) SEs must equal the full joint model's SEs
di _n "===== B. Murphy–Topel acid test ====="
clear
set obs 50000
set seed 20260713
gen edu    = ceil(4 * runiform())
gen region = ceil(5 * runiform())
gen age    = 30 + 20 * runiform() + 2*edu
gen income = rnormal(10, 2) - 0.5*region
gen id2    = _n
gen lambda = exp(-6 + 0.03*age + 0.05*income ///
    + 0.3*(edu==2) + 0.5*(edu==3) - 0.2*(edu==4) ///
    + 0.1*(region==2) - 0.2*(region==3))
gen ptime  = 0.5 + 0.5 * runiform()
gen failed = (runiform() < 1 - exp(-lambda * ptime))
stset ptime, failure(failed) id(id2)

tempvar expv
qui gen double `expv' = _t - _t0
qui poisson _d age income i.edu i.region if _st, exposure(`expv') nolog
local se_edu3_full = _se[3.edu]
local se_cons_full = _se[_cons]

// Conditional (default): constant SE understated
stpois age income i.edu i.region, fast(offset) nolog
local se_cons_cond = _se[_cons]
di "se(_cons): full=" %9.6f `se_cons_full' "  conditional=" %9.6f `se_cons_cond'
assert `se_cons_cond' < `se_cons_full'

// Murphy–Topel: matches full joint SEs
stpois age income i.edu i.region, fast(offset) mtopel nolog
di "se(3.edu): full=" %9.6f `se_edu3_full' "  mtopel=" %9.6f _se[3.edu]
di "se(_cons): full=" %9.6f `se_cons_full' "  mtopel=" %9.6f _se[_cons]
assert abs(_se[3.edu] - `se_edu3_full') < 1e-6
assert abs(_se[_cons] - `se_cons_full') < 1e-6
assert e(mtopel) == 1
di "PASS: mtopel SEs equal the full joint-model SEs"

// ── C. Timing benchmarks ───────────────────────────────────────────────────
di _n "===== C. Timing benchmarks ====="
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
    qui stpois age income i.edu i.region i.period, fast(moments) nolog
    timer off 3
    assert abs(_b[age] - `b_check') < 1e-6

    timer on 4
    qui stpois age income, absorb(edu region period) nolog
    timer off 4

    timer list
    di "n=`n'  streg=" r(t1) "  stpois=" r(t2) "  fast(moments)=" r(t3) "  absorb=" r(t4)
}
di _n "===== BENCHMARK DONE ====="
