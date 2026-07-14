// test_stpois.do  -- validation for stpois
// Run from the repository root:  do tests/test_stpois.do
discard
adopath ++ "`c(pwd)'"

// ── 1. Standard path (regression test + streg equivalence) ─────────────────
di _n "=== 1. Standard path ==="
webuse stan3, clear
stset t1, failure(died) id(id)
stpois age posttran, nolog
assert e(converged) == 1
di "PASS: N=" e(N) " ll=" e(ll)

local b_age_p  = _b[age]
local se_age_p = _se[age]
qui streg age posttran, distribution(exponential) nohr nolog
assert abs(`b_age_p'  - _b[age])  < 1e-6
assert abs(`se_age_p' - _se[age]) < 1e-6
di "PASS: matches streg, dist(exponential)"

// ── 2. HDFE acid test ──────────────────────────────────────────────────────
di _n "=== 2. HDFE acid test ==="
webuse stan3, clear
stset t1, failure(died) id(id)
tempvar expv
qui gen double `expv' = _t - _t0

poisson _d age posttran i.surgery if _st, exposure(`expv') nolog
local b_age_ref  = _b[age]
local b_post_ref = _b[posttran]
local se_age_ref = _se[age]

stpois age posttran, absorb(surgery) nolog tol(1e-10)
di "age:      ref=" `b_age_ref'  " hdfe=" _b[age]
di "se(age):  ref=" `se_age_ref' " hdfe=" _se[age]
assert abs(_b[age]      - `b_age_ref')  < 1e-4
assert abs(_b[posttran] - `b_post_ref') < 1e-4
assert abs(_se[age]     - `se_age_ref') < 1e-4
di "PASS: coefs and SEs match to 1e-4"

// HDFE robust
stpois age posttran, absorb(surgery) vce(robust) nolog
assert _se[age] > 0
di "PASS: HDFE + robust VCE runs"

// HDFE with interaction in the varlist
qui poisson _d c.age##i.posttran i.surgery if _st, exposure(`expv') nolog
local b_int_ref  = _b[c.age#1.posttran]
local se_int_ref = _se[c.age#1.posttran]
stpois c.age##i.posttran, absorb(surgery) nolog tol(1e-10)
di "c.age#1.posttran: ref=" `b_int_ref' " hdfe=" _b[c.age#1.posttran]
assert abs(_b[c.age#1.posttran]  - `b_int_ref')  < 1e-4
assert abs(_se[c.age#1.posttran] - `se_int_ref') < 1e-4
di "PASS: HDFE with varlist interaction matches"

// ── 3. Synthetic data ───────────────────────────────────────────────────────
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

gen lambda = exp( ///
    0.03 * age + 0.05 * income                       + ///
    0.3 * (edu==2) + 0.5 * (edu==3) - 0.2 * (edu==4) + ///
    0.1 * (region==2) - 0.2 * (region==3)            + ///
    0.2 * (period >= 6))

gen ptime  = 0.5 + 0.5 * runiform()
gen prob_e = 1 - exp(-lambda * ptime)
gen failed = (runiform() < prob_e)
gen t_end  = ptime

stset t_end, failure(failed) id(id2)
tempvar expv2
qui gen double `expv2' = _t - _t0

// ── 4. fast: exact MLE acid test ────────────────────────────────────────────
di _n "=== 4. fast exact MLE ==="
qui poisson _d age income i.edu i.region i.period if _st, exposure(`expv2') nolog
local b_age_full  = _b[age]
local se_age_full = _se[age]
local b_inc_full  = _b[income]
local b_edu3_full = _b[3.edu]
local se_edu3_full = _se[3.edu]

stpois age income i.edu i.region i.period, fast nolog
di "PASS: N=" e(N) " N_cells=" e(N_cells)
assert e(N_cells) >= 100
assert e(converged) == 1
assert abs(_b[age]     - `b_age_full')   < 1e-6
assert abs(_se[age]    - `se_age_full')  < 1e-6
assert abs(_b[income]  - `b_inc_full')   < 1e-6
assert abs(_b[3.edu]   - `b_edu3_full')  < 1e-6
assert abs(_se[3.edu]  - `se_edu3_full') < 1e-6
di "PASS: fast matches full MLE on coefficients and SEs (1e-6)"

// ── 5. fast with interactions ───────────────────────────────────────────────
di _n "=== 5. fast with interactions ==="

// cont#cont
qui poisson _d c.age#c.income age income i.edu if _st, exposure(`expv2') nolog
local b_cc  = _b[c.age#c.income]
local se_cc = _se[c.age#c.income]
stpois c.age#c.income age income i.edu, fast nolog
di "c.age#c.income: ref=" `b_cc' " fast=" _b[c.age#c.income]
assert abs(_b[c.age#c.income]  - `b_cc')  < 1e-6
assert abs(_se[c.age#c.income] - `se_cc') < 1e-6
di "PASS: cont#cont interaction exact"

// cont#cat
qui poisson _d c.age##i.edu i.region if _st, exposure(`expv2') nolog
local b_ce  = _b[c.age#3.edu]
local se_ce = _se[c.age#3.edu]
stpois c.age##i.edu i.region, fast nolog
di "c.age#3.edu: ref=" `b_ce' " fast=" _b[c.age#3.edu]
assert abs(_b[c.age#3.edu]  - `b_ce')  < 1e-6
assert abs(_se[c.age#3.edu] - `se_ce') < 1e-6
di "PASS: cont#cat interaction exact"

// cat#cat
qui poisson _d age i.edu##i.region if _st, exposure(`expv2') nolog
local b_ee  = _b[3.edu#2.region]
local se_ee = _se[3.edu#2.region]
stpois age i.edu##i.region, fast nolog
di "3.edu#2.region: ref=" `b_ee' " fast=" _b[3.edu#2.region]
assert abs(_b[3.edu#2.region]  - `b_ee')  < 1e-6
assert abs(_se[3.edu#2.region] - `se_ee') < 1e-6
di "PASS: cat#cat interaction exact"

// ── 6. fast with all-categorical predictors ─────────────────────────────────
di _n "=== 6. fast all-categorical ==="
qui poisson _d i.edu i.region i.period if _st, exposure(`expv2') nolog
local b_edu3_c  = _b[3.edu]
local se_edu3_c = _se[3.edu]
local ll_c      = e(ll)
stpois i.edu i.region i.period, fast nolog
assert abs(_b[3.edu]  - `b_edu3_c')  < 1e-6
assert abs(_se[3.edu] - `se_edu3_c') < 1e-6
assert abs(e(ll) - `ll_c') < 1e-4
di "PASS: all-categorical fast exact (pure cell-level Newton)"

// ── 7. absorb() with interaction FE ─────────────────────────────────────────
di _n "=== 7. absorb(edu#region) ==="
qui poisson _d age income i.edu#i.region if _st, exposure(`expv2') nolog
local b_age_fe  = _b[age]
local se_age_fe = _se[age]
stpois age income, absorb(edu#region) nolog tol(1e-10)
di "age: ref=" `b_age_fe' " hdfe=" _b[age]
assert abs(_b[age]  - `b_age_fe')  < 1e-4
assert abs(_se[age] - `se_age_fe') < 1e-4
di "PASS: interaction FE matches explicit i.edu#i.region dummies"

// ── 8. fast on stan3 (2 cells) ──────────────────────────────────────────────
di _n "=== 8. fast on stan3 ==="
webuse stan3, clear
stset t1, failure(died) id(id)
tempvar expv3
qui gen double `expv3' = _t - _t0
qui poisson _d age i.surgery if _st, exposure(`expv3') nolog
local b_age_s3  = _b[age]
local se_age_s3 = _se[age]
stpois age i.surgery, fast nolog
assert e(N_cells) == 2
assert abs(_b[age]  - `b_age_s3')  < 1e-6
assert abs(_se[age] - `se_age_s3') < 1e-6
di "PASS: fast stan3 -> 2 cells, exact on age coef and SE"

// ── 9. Replay and e() locals ────────────────────────────────────────────────
di _n "=== 9. e() locals and replay ==="
webuse stan3, clear
stset t1, failure(died) id(id)
stpois age posttran, nolog
assert "`e(cmd)'" == "stpois"
assert "`e(predict)'" == "stpois_p"
stpois   // replay
di "PASS: cmd locals and replay"

di _n "=== ALL TESTS PASSED ==="
