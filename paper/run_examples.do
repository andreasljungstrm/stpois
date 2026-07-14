// run_examples.do -- generate real output for the stpois paper
// Run from repo root
discard
adopath ++ "`c(pwd)'"

capture log close
log using paper/examples_output.log, replace text

di _n "===== EXAMPLE 1: standard path vs streg ====="
webuse stan3, clear
stset t1, failure(died) id(id)
stpois age posttran, nolog
streg age posttran, distribution(exponential) nolog nohr

di _n "===== EXAMPLE 2: HDFE ====="
tempvar expv
qui gen double `expv' = _t - _t0
poisson _d age posttran i.surgery if _st, exposure(`expv') nolog
stpois age posttran, absorb(surgery) nolog tol(1e-10)
stpois age posttran, absorb(surgery) vce(robust) nolog

di _n "===== EXAMPLE 3: fast on 500,000 observations ====="
clear
set obs 500000
set seed 20260713
gen edu    = ceil(4 * runiform())
gen region = ceil(5 * runiform())
gen period = ceil(10 * runiform())
gen age    = 30 + 20 * runiform()
gen income = rnormal(10, 2)
gen id2    = _n
gen lambda = exp(-6 + ///
    0.03 * age + 0.05 * income                       + ///
    0.3 * (edu==2) + 0.5 * (edu==3) - 0.2 * (edu==4) + ///
    0.1 * (region==2) - 0.2 * (region==3)            + ///
    0.2 * (period >= 6))
gen ptime  = 0.5 + 0.5 * runiform()
gen failed = (runiform() < 1 - exp(-lambda * ptime))
stset ptime, failure(failed) id(id2)

stpois age income i.edu i.region i.period, fast nolog
local b_age_f = _b[age]
local ll_f    = e(ll)
qui stpois age income i.edu i.region i.period, nolog
di "standard: age=" %12.9f _b[age] "  fast: age=" %12.9f `b_age_f'
di "standard: ll=" %14.5f e(ll)   "  fast: ll=" %14.5f `ll_f'
assert abs(_b[age] - `b_age_f') < 1e-8
assert abs(e(ll) - `ll_f') < 1e-4

di _n "===== EXAMPLE 4: interactions with fast ====="
stpois c.age##i.edu income i.region i.period, fast nolog

log close
