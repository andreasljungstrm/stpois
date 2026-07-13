// benchmark_hdfe.do -- absorb() vs explicit dummies, high-dimensional FE
discard
adopath ++ "`c(pwd)'"

clear
set obs 1000000
set seed 20260713
gen g      = ceil(1000 * runiform())    // 1,000-level fixed effect
gen age    = 30 + 20 * runiform()
gen income = rnormal(10, 2)
gen id2    = _n
gen geff   = (mod(g, 7) - 3) * 0.1
gen lambda = exp(-6 + 0.03*age + 0.05*income + geff)
gen ptime  = 0.5 + 0.5 * runiform()
gen failed = (runiform() < 1 - exp(-lambda * ptime))
qui stset ptime, failure(failed) id(id2)

tempvar expv
qui gen double `expv' = _t - _t0

timer clear
timer on 1
qui poisson _d age income i.g if _st, exposure(`expv') nolog
timer off 1
local b_ref  = _b[age]
local se_ref = _se[age]

timer on 2
qui stpois age income, absorb(g) nolog
timer off 2

timer list
di "poisson i.g (1000 dummies): " r(t1) "s   stpois absorb(g): " r(t2) "s"
di "age: ref=" %10.7f `b_ref' " (se " %10.7f `se_ref' ")   absorb=" %10.7f _b[age] " (se " %10.7f _se[age] ")"
assert abs(_b[age]  - `b_ref')  < 1e-5
assert abs(_se[age] - `se_ref') < 1e-6
di "PASS: absorb matches 1000-dummy poisson"
