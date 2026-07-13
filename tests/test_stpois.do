// test_stpois.do  -- validation for stpois
// Run from the repository root:  do tests/test_stpois.do
discard
adopath ++ "`c(pwd)'"

// ── 1. Standard path (regression test) ─────────────────────────────────────
di _n "=== 1. Standard path ==="
webuse stan3, clear
stset t1, failure(died) id(id)
stpois age posttran, nolog
assert e(converged) == 1
di "PASS: N=" e(N) " ll=" e(ll)

// ── 2. HDFE acid test ──────────────────────────────────────────────────────
// stpois age posttran, absorb(surgery) must match
// poisson _d age posttran i.surgery, exposure(ptime) on coefs AND SEs
di _n "=== 2. HDFE acid test ==="
tempvar expv
qui gen double `expv' = _t - _t0

poisson _d age posttran i.surgery if _st, exposure(`expv') nolog
local b_age_ref  = _b[age]
local b_post_ref = _b[posttran]
local se_age_ref = _se[age]

stpois age posttran, absorb(surgery) nolog tol(1e-10)
di "age:      ref=" `b_age_ref'  " hdfe=" _b[age]
di "posttran: ref=" `b_post_ref' " hdfe=" _b[posttran]
di "se(age):  ref=" `se_age_ref' " hdfe=" _se[age]
assert abs(_b[age]      - `b_age_ref')  < 1e-4
assert abs(_b[posttran] - `b_post_ref') < 1e-4
assert abs(_se[age]     - `se_age_ref') < 1e-4
di "PASS: coefs and SEs match to 1e-4"

// HDFE robust
stpois age posttran, absorb(surgery) vce(robust) nolog
assert _se[age] > 0
di "PASS: HDFE + robust VCE runs"

// ── 3. Synthetic data for fast methods ─────────────────────────────────────
// 4 edu × 5 region × 10 period = 200 cells
// Binary outcome: person either fails or not during their period
di _n "=== 3. Synthetic data ==="
clear
set obs 5000
set seed 20260713

gen edu    = ceil(4 * runiform())
gen region = ceil(5 * runiform())
gen period = ceil(10 * runiform())
gen age    = 30 + 20 * runiform()
gen income = rnormal(10, 2)
gen id2    = _n

// True hazard (rate per unit time)
gen lambda = exp( ///
    0.03 * age + 0.05 * income                         + ///
    0.3 * (edu==2) + 0.5 * (edu==3) - 0.2 * (edu==4) + ///
    0.1 * (region==2) - 0.2 * (region==3)             + ///
    0.2 * (period >= 6))

// Each person observed for ptime; failure = Bernoulli with P = 1-exp(-lambda*ptime)
gen ptime  = 0.5 + 0.5 * runiform()
gen prob_e = 1 - exp(-lambda * ptime)
gen failed = (runiform() < prob_e)
gen t_end  = ptime

stset t_end, failure(failed) id(id2)

qui count if _d
di "Events: " r(N) " out of " _N

// Full MLE reference
tempvar expv2
qui gen double `expv2' = _t - _t0
poisson _d age income i.edu i.region i.period if _st, exposure(`expv2') nolog
local b_age_full  = _b[age]
local b_inc_full  = _b[income]
di "Full MLE: age=" `b_age_full' " income=" `b_inc_full'

// ── 4. fast(offset) ─────────────────────────────────────────────────────────
di _n "=== 4. fast(offset) ==="
stpois age income i.edu i.region i.period, fast(offset) nolog
di "PASS: N=" e(N) " N_cells=" e(N_cells) " ll=" e(ll)
assert e(N_cells) < e(N)
assert e(N_cells) >= 100   // most of 200 possible cells should have data
assert e(converged) == 1

// ── 5. fast(moments): exact MLE acid test ──────────────────────────────────
di _n "=== 5. fast(moments) ==="
// Reference SEs from the full model
qui poisson _d age income i.edu i.region i.period if _st, exposure(`expv2') nolog
local se_age_full = _se[age]
local b_edu3_full = _b[3.edu]
local se_edu3_full = _se[3.edu]

stpois age income i.edu i.region i.period, fast(moments) nolog
di "PASS: N=" e(N) " N_cells=" e(N_cells)
assert e(N_cells) >= 100
assert e(converged) == 1
di "age:   full=" `b_age_full' " moments=" _b[age]
di "3.edu: full=" `b_edu3_full' " moments=" _b[3.edu]
// Exact: coefficients AND SEs must match the full MLE
assert abs(_b[age]     - `b_age_full')   < 1e-6
assert abs(_se[age]    - `se_age_full')  < 1e-6
assert abs(_b[income]  - `b_inc_full')   < 1e-6
assert abs(_b[3.edu]   - `b_edu3_full')  < 1e-6
assert abs(_se[3.edu]  - `se_edu3_full') < 1e-6
di "PASS: fast(moments) matches full MLE on coefficients and SEs (1e-6)"

// ── 6. fast methods on stan3 ────────────────────────────────────────────────
di _n "=== 6. fast methods on stan3 ==="
webuse stan3, clear
stset t1, failure(died) id(id)
stpois age i.surgery, fast(offset) nolog
assert e(N_cells) == 2
di "PASS: fast(offset) stan3 → 2 cells"

// fast(moments) is exact even with 2 cells: age identified from
// within-cell variation, matching poisson _d age i.surgery exactly
tempvar expv3
qui gen double `expv3' = _t - _t0
qui poisson _d age i.surgery if _st, exposure(`expv3') nolog
local b_age_s3  = _b[age]
local se_age_s3 = _se[age]
stpois age i.surgery, fast(moments) nolog
assert e(N_cells) == 2
assert abs(_b[age]  - `b_age_s3')  < 1e-6
assert abs(_se[age] - `se_age_s3') < 1e-6
di "PASS: fast(moments) stan3 → 2 cells, exact match on age coef and SE"

// ── 7. Replay and e() locals ─────────────────────────────────────────────────
di _n "=== 7. e() locals and replay ==="
webuse stan3, clear
stset t1, failure(died) id(id)
stpois age posttran, nolog
assert "`e(cmd)'" == "stpois"
assert "`e(predict)'" == "stpois_p"
stpois   // replay
di "PASS: cmd locals and replay"

di _n "=== ALL TESTS PASSED ==="
