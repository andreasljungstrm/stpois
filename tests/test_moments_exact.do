// test_moments_exact.do -- fast(moments) must match full MLE on coefs AND SEs
discard
adopath ++ "`c(pwd)'"

clear
set obs 50000
set seed 20260713
gen edu    = ceil(4 * runiform())
gen region = ceil(5 * runiform())
gen age    = 30 + 20 * runiform() + 2*edu   // correlate with cells
gen income = rnormal(10, 2) - 0.5*region
gen id2    = _n
gen lambda = exp(-6 + 0.03*age + 0.05*income ///
    + 0.3*(edu==2) + 0.5*(edu==3) - 0.2*(edu==4) ///
    + 0.1*(region==2) - 0.2*(region==3))
gen ptime  = 0.5 + 0.5 * runiform()
gen failed = (runiform() < 1 - exp(-lambda * ptime))
stset ptime, failure(failed) id(id2)

// Reference full MLE
tempvar expv
qui gen double `expv' = _t - _t0
qui poisson _d age income i.edu i.region if _st, exposure(`expv') nolog
local b_age = _b[age]
local se_age = _se[age]
local b_inc = _b[income]
local se_inc = _se[income]
local b_edu3 = _b[3.edu]
local se_edu3 = _se[3.edu]
local b_cons = _b[_cons]
local se_cons = _se[_cons]
local ll_ref = e(ll)

stpois age income i.edu i.region, fast(moments) nolog
di _n "age:    ref b=" `b_age'  " se=" `se_age'  "  new b=" _b[age]    " se=" _se[age]
di    "income: ref b=" `b_inc'  " se=" `se_inc'  "  new b=" _b[income] " se=" _se[income]
di    "3.edu:  ref b=" `b_edu3' " se=" `se_edu3' "  new b=" _b[3.edu]  " se=" _se[3.edu]
di    "_cons:  ref b=" `b_cons' " se=" `se_cons' "  new b=" _b[_cons]  " se=" _se[_cons]
di    "ll:     ref=" `ll_ref' "  new=" e(ll)
assert abs(_b[age]      - `b_age')   < 1e-6
assert abs(_se[age]     - `se_age')  < 1e-6
assert abs(_b[income]   - `b_inc')   < 1e-6
assert abs(_se[income]  - `se_inc')  < 1e-6
assert abs(_b[3.edu]    - `b_edu3')  < 1e-6
assert abs(_se[3.edu]   - `se_edu3') < 1e-6
assert abs(_b[_cons]    - `b_cons')  < 1e-5
assert abs(_se[_cons]   - `se_cons') < 1e-6
assert abs(e(ll) - `ll_ref') < 1e-4
di "PASS: exact match on coefficients, SEs, and ll"

// robust VCE
qui poisson _d age income i.edu i.region if _st, exposure(`expv') vce(robust) nolog
local se_age_r = _se[age]
local se_edu3_r = _se[3.edu]
stpois age income i.edu i.region, fast(moments) vce(robust) nolog
di _n "robust se(age):   ref=" `se_age_r' " new=" _se[age]
di    "robust se(3.edu): ref=" `se_edu3_r' " new=" _se[3.edu]
assert abs(_se[age]   - `se_age_r')  < 1e-6
assert abs(_se[3.edu] - `se_edu3_r') < 1e-6
di "PASS: robust VCE exact"

// cluster VCE
qui poisson _d age income i.edu i.region if _st, exposure(`expv') vce(cluster region) nolog
local se_age_c = _se[age]
stpois age income i.edu i.region, fast(moments) vce(cluster region) nolog
di _n "cluster se(age): ref=" `se_age_c' " new=" _se[age]
assert abs(_se[age] - `se_age_c') < 1e-6
di "PASS: cluster VCE exact"

// margins/predict smoke test
qui margins, dydx(age)
di "PASS: margins runs"

di _n "=== ALL MOMENTS ACID TESTS PASSED ==="
