// run_examples.do -- generate real output for the stpois paper
// Run from repo root
discard
adopath ++ "`c(pwd)'"

log using paper/examples_output.log, replace text

di _n "===== EXAMPLE 1: standard path vs streg ====="
webuse stan3, clear
stset t1, failure(died) id(id)
stpois age posttran, nolog
streg age posttran, distribution(exponential) nolog nohr

di _n "===== EXAMPLE 2: HDFE acid test ====="
tempvar expv
qui gen double `expv' = _t - _t0
poisson _d age posttran i.surgery if _st, exposure(`expv') nolog
stpois age posttran, absorb(surgery) nolog tol(1e-10)
stpois age posttran, absorb(surgery) vce(robust) nolog

di _n "===== EXAMPLE 3: synthetic data, fast paths ====="
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
gen prob_e = 1 - exp(-lambda * ptime)
gen failed = (runiform() < prob_e)
gen t_end  = ptime
stset t_end, failure(failed) id(id2)
qui count if _d
di "Events: " r(N) " of " _N

timer clear
timer on 1
stpois age income i.edu i.region i.period, nolog
timer off 1
estimates store full
local b_age_full = _b[age]
local b_inc_full = _b[income]

timer on 2
stpois age income i.edu i.region i.period, fast(offset) nolog
timer off 2
estimates store fast_off

timer on 3
stpois age income i.edu i.region i.period, fast(moments) nolog
timer off 3
estimates store fast_mom

timer list

di _n "Full MLE:      age=" %9.6f `b_age_full' "  income=" %9.6f `b_inc_full'

estimates restore fast_mom
di "fast(moments): age=" %9.6f _b[age] "  income=" %9.6f _b[income]

estimates table full fast_off fast_mom, b(%9.6f) se(%9.6f) stats(N)

log close
