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
assert abs(_b[age]      - `b_age_ref')  < 1e-6
assert abs(_b[posttran] - `b_post_ref') < 1e-6
assert abs(_se[age]     - `se_age_ref') < 1e-6
assert "`e(chi2type)'" == "Wald"
di "PASS: coefs and SEs match to 1e-6; Wald chi2 reported"

// Wald chi2 matches the joint test on the dummy-variable model
qui poisson _d age posttran i.surgery if _st, exposure(`expv') nolog
qui test age posttran
local wald_ref = r(chi2)
qui stpois age posttran, absorb(surgery) nolog tol(1e-10)
di "Wald chi2: ref=" `wald_ref' " hdfe=" e(chi2)
assert abs(e(chi2) - `wald_ref') < 1e-4
assert e(df_m) == 2
di "PASS: HDFE Wald chi2 matches joint test on dummy model"

// HDFE robust: exact vs dummy-variable poisson with vce(robust)
qui poisson _d age posttran i.surgery if _st, exposure(`expv') vce(robust) nolog
local se_age_rr  = _se[age]
local se_post_rr = _se[posttran]
stpois age posttran, absorb(surgery) vce(robust) nolog tol(1e-10)
di "robust se(age): ref=" `se_age_rr' " hdfe=" _se[age]
assert abs(_se[age]      - `se_age_rr')  < 1e-6
assert abs(_se[posttran] - `se_post_rr') < 1e-6
di "PASS: HDFE + robust VCE exact"

// HDFE with interaction in the varlist
qui poisson _d c.age##i.posttran i.surgery if _st, exposure(`expv') nolog
local b_int_ref  = _b[c.age#1.posttran]
local se_int_ref = _se[c.age#1.posttran]
stpois c.age##i.posttran, absorb(surgery) nolog tol(1e-10)
di "c.age#1.posttran: ref=" `b_int_ref' " hdfe=" _b[c.age#1.posttran]
assert abs(_b[c.age#1.posttran]  - `b_int_ref')  < 1e-6
assert abs(_se[c.age#1.posttran] - `se_int_ref') < 1e-6
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
assert e(N_cells) < .
di "PASS: all-categorical fast exact (collapse + poisson route)"

// 6b. all-categorical route under iweight (weighted ll reconstruction) and
//     under robust (which must fall through to the Mata cell-sum Newton)
tempvar iw6
qui gen double `iw6' = 0.5 + runiform() if _st
qui poisson _d i.edu i.region i.period [iw=`iw6'] if _st, exposure(`expv2') nolog
local b_iw6  = _b[3.edu]
local se_iw6 = _se[3.edu]
local ll_iw6 = e(ll)
stpois i.edu i.region i.period [iw=`iw6'], fast nolog
assert abs(_b[3.edu]  - `b_iw6')  < 1e-6
assert abs(_se[3.edu] - `se_iw6') < 1e-6
assert abs(e(ll) - `ll_iw6') < 1e-3
di "PASS: all-categorical fast + iweight exact (b, se, ll)"

qui poisson _d i.edu i.region i.period if _st, exposure(`expv2') vce(robust) nolog
local se_edu3_r = _se[3.edu]
stpois i.edu i.region i.period, fast vce(robust) nolog
assert abs(_se[3.edu] - `se_edu3_r') < 1e-6
di "PASS: all-categorical fast + robust exact (Mata fallthrough)"

// ── 7. absorb() with interaction FE ─────────────────────────────────────────
di _n "=== 7. absorb(edu#region) ==="
qui poisson _d age income i.edu#i.region if _st, exposure(`expv2') nolog
local b_age_fe  = _b[age]
local se_age_fe = _se[age]
stpois age income, absorb(edu#region) nolog tol(1e-10)
di "age: ref=" `b_age_fe' " hdfe=" _b[age]
assert abs(_b[age]  - `b_age_fe')  < 1e-6
assert abs(_se[age] - `se_age_fe') < 1e-6
di "PASS: interaction FE matches explicit i.edu#i.region dummies"

// ── 7b. Multi-term absorb ───────────────────────────────────────────────────
di _n "=== 7b. absorb(edu region period) ==="
qui poisson _d age income i.edu i.region i.period if _st, exposure(`expv2') nolog
local b_age_m  = _b[age]
local se_age_m = _se[age]
stpois age income, absorb(edu region period) nolog tol(1e-10)
di "age: ref=" `b_age_m' " hdfe=" _b[age]
assert abs(_b[age]  - `b_age_m')  < 1e-6
assert abs(_se[age] - `se_age_m') < 1e-6
di "PASS: multi-term absorb matches explicit dummies"

// ── 7c. absorb() with cluster VCE ───────────────────────────────────────────
di _n "=== 7c. absorb() + vce(cluster) ==="
qui poisson _d age income i.edu if _st, exposure(`expv2') vce(cluster region) nolog
local b_age_cl  = _b[age]
local se_age_cl = _se[age]
stpois age income, absorb(edu) vce(cluster region) nolog tol(1e-10)
di "cluster se(age): ref=" `se_age_cl' " hdfe=" _se[age]
assert abs(_b[age]  - `b_age_cl')  < 1e-6
assert abs(_se[age] - `se_age_cl') < 1e-6
assert e(N_clust) == 5
di "PASS: clustered absorb SEs match dummy poisson vce(cluster)"

// ── 7d. Separation drop: header scalars on the final sample ────────────────
di _n "=== 7d. separation drop consistency ==="
clear
set obs 1000
set seed 20260713
gen g      = 1 + (_n > 200)                    // group 1: first 200 obs
gen x      = rnormal()
gen ptime  = 1 + runiform()
gen failed = (runiform() < 0.3) * (g == 2)     // group 1 has no events
gen idz    = _n
qui stset ptime, failure(failed) id(idz)
stpois x, absorb(g) nolog
assert e(N)     == 800
assert e(N_sub) == 800
qui count if failed == 1
assert e(N_fail) == r(N)
di "PASS: zero-event group dropped; header scalars computed on final sample"

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

// ── 8b. Continuous coefficient identified from within-cell variation only ──
// agew has, by construction, mean zero within every surgery cell, so it has
// zero between-cell variation; fast must still recover the exact MLE.
di _n "=== 8b. fast, zero between-cell variation ==="
qui bysort surgery: egen double agebar = mean(age)
qui gen double agew = age - agebar
qui summ agew if surgery == 0, meanonly
assert abs(r(mean)) < 1e-10
qui summ agew if surgery == 1, meanonly
assert abs(r(mean)) < 1e-10
qui poisson _d agew i.surgery if _st, exposure(`expv3') nolog
local b_agew  = _b[agew]
local se_agew = _se[agew]
stpois agew i.surgery, fast nolog
di "agew: ref=" `b_agew' " fast=" _b[agew]
assert abs(_b[agew]  - `b_agew')  < 1e-6
assert abs(_se[agew] - `se_agew') < 1e-6
di "PASS: continuous coef identified purely within cells, exact"

// ── 8c. Weights on fast: exact vs weighted poisson ─────────────────────────
di _n "=== 8c. weights on fast ==="
gen int fwv = 1 + mod(_n, 3)
gen double pwv = 0.5 + runiform()

// fweight
qui poisson _d age i.surgery [fw=fwv] if _st, exposure(`expv3') nolog
local b_fw  = _b[age]
local se_fw = _se[age]
local N_fw  = e(N)
local ll_fw = e(ll)
stpois age i.surgery [fw=fwv], fast nolog
assert abs(_b[age]  - `b_fw')  < 1e-6
assert abs(_se[age] - `se_fw') < 1e-6
assert e(N) == `N_fw'
di "PASS: fast + fweight exact (coef, SE, N)"

// pweight (implies robust)
qui poisson _d age i.surgery [pw=pwv] if _st, exposure(`expv3') nolog
local b_pw  = _b[age]
local se_pw = _se[age]
stpois age i.surgery [pw=pwv], fast nolog
assert abs(_b[age]  - `b_pw')  < 1e-6
assert abs(_se[age] - `se_pw') < 1e-6
di "PASS: fast + pweight exact (robust SEs implied)"

// iweight
qui poisson _d age i.surgery [iw=pwv] if _st, exposure(`expv3') nolog
local b_iw  = _b[age]
local se_iw = _se[age]
stpois age i.surgery [iw=pwv], fast nolog
assert abs(_b[age]  - `b_iw')  < 1e-6
assert abs(_se[age] - `se_iw') < 1e-6
di "PASS: fast + iweight exact"

// ── 8d. Weights on absorb(): exact vs weighted dummy poisson ───────────────
di _n "=== 8d. weights on absorb() ==="
qui poisson _d age posttran i.surgery [fw=fwv] if _st, exposure(`expv3') nolog
local b_afw  = _b[age]
local se_afw = _se[age]
stpois age posttran [fw=fwv], absorb(surgery) nolog tol(1e-10)
assert abs(_b[age]  - `b_afw')  < 1e-6
assert abs(_se[age] - `se_afw') < 1e-6
di "PASS: absorb + fweight exact"

qui poisson _d age posttran i.surgery [pw=pwv] if _st, exposure(`expv3') nolog
local b_apw  = _b[age]
local se_apw = _se[age]
stpois age posttran [pw=pwv], absorb(surgery) nolog tol(1e-10)
assert abs(_b[age]  - `b_apw')  < 1e-6
assert abs(_se[age] - `se_apw') < 1e-6
di "PASS: absorb + pweight exact (robust SEs implied)"

// ── 8e. fast + absorb() combined ────────────────────────────────────────────
di _n "=== 8e. fast + absorb() ==="
qui poisson _d age i.posttran i.surgery if _st, exposure(`expv3') nolog
local b_c  = _b[age]
local se_c = _se[age]
local b_cp  = _b[1.posttran]
local se_cp = _se[1.posttran]
stpois age i.posttran, fast absorb(surgery) nolog tol(1e-10)
di "age: ref=" `b_c' " fast+absorb=" _b[age]
assert abs(_b[age]         - `b_c')   < 1e-6
assert abs(_se[age]        - `se_c')  < 1e-6
assert abs(_b[1.posttran]  - `b_cp')  < 1e-6
assert abs(_se[1.posttran] - `se_cp') < 1e-6
assert "`e(chi2type)'" == "Wald"
di "PASS: fast+absorb matches dummy poisson (coefs and SEs)"

// robust
qui poisson _d age i.posttran i.surgery if _st, exposure(`expv3') vce(robust) nolog
local se_cr = _se[age]
stpois age i.posttran, fast absorb(surgery) vce(robust) nolog tol(1e-10)
assert abs(_se[age] - `se_cr') < 1e-6
di "PASS: fast+absorb robust VCE exact"

// with fweights
qui poisson _d age i.posttran i.surgery [fw=fwv] if _st, exposure(`expv3') nolog
local b_cw  = _b[age]
local se_cw = _se[age]
stpois age i.posttran [fw=fwv], fast absorb(surgery) nolog tol(1e-10)
assert abs(_b[age]  - `b_cw')  < 1e-6
assert abs(_se[age] - `se_cw') < 1e-6
di "PASS: fast+absorb + fweight exact"

// ── 8f. fast + multi-term absorb + cluster VCE ─────────────────────────────
di _n "=== 8f. fast + absorb(region period) + cluster ==="
clear
set obs 20000
set seed 20260713
gen edu    = ceil(4 * runiform())
gen region = ceil(5 * runiform())
gen period = ceil(10 * runiform())
gen age    = 30 + 20 * runiform() + edu          // correlate with FE space
gen lambda = exp(-5 + 0.03*age + 0.3*(edu==2) + 0.5*(edu==3) ///
    + 0.1*(region==2) + 0.2*(period >= 6))
gen ptime  = 0.5 + 0.5 * runiform()
gen failed = (runiform() < 1 - exp(-lambda * ptime))
gen id3    = _n
qui stset ptime, failure(failed) id(id3)
tempvar expv8
qui gen double `expv8' = _t - _t0
qui poisson _d age i.edu i.region i.period if _st, exposure(`expv8') vce(cluster region) nolog
local b_m  = _b[age]
local se_m = _se[age]
local b_m3  = _b[3.edu]
local se_m3 = _se[3.edu]
stpois age i.edu, fast absorb(region period) vce(cluster region) nolog tol(1e-10)
di "age:   ref=" `b_m'  " combo=" _b[age]
di "3.edu: ref=" `b_m3' " combo=" _b[3.edu]
assert abs(_b[age]    - `b_m')   < 1e-6
assert abs(_se[age]   - `se_m')  < 1e-6
assert abs(_b[3.edu]  - `b_m3')  < 1e-6
assert abs(_se[3.edu] - `se_m3') < 1e-6
assert e(N_clust) == 5
di "PASS: fast + multi-term absorb + clustered SEs exact"

// ── 9. Replay and e() locals ────────────────────────────────────────────────
di _n "=== 9. e() locals and replay ==="
webuse stan3, clear
stset t1, failure(died) id(id)
stpois age posttran, nolog
assert "`e(cmd)'" == "stpois"
assert "`e(predict)'" == "stpois_p"
stpois   // replay
di "PASS: cmd locals and replay"

// predict after absorb(): full linear predictor and level statistics,
// exact against the dummy-variable poisson predictions
tempvar expv9
qui gen double `expv9' = _t - _t0
qui poisson _d age posttran i.surgery if _st, exposure(`expv9') nolog
predict double n_ref, n
qui stpois age posttran, absorb(surgery) nolog tol(1e-10)
predict double n_hdfe, n
qui gen double n_diff = abs(n_hdfe - n_ref)
qui summ n_diff
assert r(max) < 1e-5
predict double hz_hdfe, hazard
predict double xb_hdfe, xb
predict double xb0_hdfe, xb nooffset
qui gen double hz_chk = abs(hz_hdfe - exp(xb0_hdfe))
qui summ hz_chk
assert r(max) < 1e-10
di "PASS: predict after absorb() exact (n vs dummy poisson; hazard = exp(xb))"

// e(fes_var) points at the stored FE variable
assert "`e(fes_var)'" == "_stpois_fe"
confirm numeric variable _stpois_fe
di "PASS: FE contribution stored in e(fes_var)"

// predict after fast+absorb
qui stpois age i.posttran, fast absorb(surgery) nolog tol(1e-10)
predict double n_fa, n
qui poisson _d age i.posttran i.surgery if _st, exposure(`expv9') nolog
predict double n_far, n
qui gen double n_fad = abs(n_fa - n_far)
qui summ n_fad
assert r(max) < 1e-5
di "PASS: predict after fast+absorb exact"

di _n "=== ALL TESTS PASSED ==="
