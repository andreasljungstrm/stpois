// benchmark_hdfe.do -- absorb() vs explicit dummies, high-dimensional FE
// Run from repo root. Writes paper/benchmark_hdfe_results.csv (versioned),
// which backs Table 4 of the paper. Timings are the MEDIAN of `reps' runs;
// the CSV also records the min and max so dispersion can be reported.
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

// agreement check (single fit)
qui poisson _d age income i.g if _st, exposure(`expv') nolog
local b_ref  = _b[age]
local se_ref = _se[age]
qui stpois age income, absorb(g) nolog
local b_st  = _b[age]
local se_st = _se[age]
di "age: ref=" %10.7f `b_ref' " (se " %10.7f `se_ref' ")   absorb=" %10.7f `b_st' " (se " %10.7f `se_st' ")"
assert abs(`b_st'  - `b_ref')  < 1e-5
assert abs(`se_st' - `se_ref') < 1e-6
di "PASS: absorb matches 1000-dummy poisson"

// median timings
capture program drop _poisdum
program _poisdum
    tempvar e
    qui gen double `e' = _t - _t0
    qui poisson _d age income i.g if _st, exposure(`e') nolog
end
_medtime `reps' "_poisdum"
local t1_med = r(med)
local t1_lo  = r(lo)
local t1_hi  = r(hi)
_medtime `reps' "stpois age income, absorb(g) nolog"
local t2_med = r(med)
local t2_lo  = r(lo)
local t2_hi  = r(hi)

// ppmlhdfe on the same design (if installed), for a consistent head-to-head
local have_ppml = 0
capture which ppmlhdfe
if _rc == 0 {
    local have_ppml = 1
    tempvar lnexp
    qui gen double `lnexp' = ln(_t - _t0)
    qui ppmlhdfe _d age income if _st, absorb(g) offset(`lnexp')
    local b_pp  = _b[age]
    di "ppmlhdfe age=" %10.7f `b_pp' "   (default vce = robust)"
    assert abs(`b_pp' - `b_ref') < 1e-5
    _medtime `reps' "ppmlhdfe _d age income if _st, absorb(g) offset(`lnexp')"
    local t3_med = r(med)
    local t3_lo  = r(lo)
    local t3_hi  = r(hi)
    di "ppmlhdfe absorb(g):         median " `t3_med' "s [" `t3_lo' ", " `t3_hi' "]"
}

di "poisson i.g (1000 dummies): median " `t1_med' "s [" `t1_lo' ", " `t1_hi' "]"
di "stpois absorb(g):           median " `t2_med' "s [" `t2_lo' ", " `t2_hi' "]"

tempname bres
file open `bres' using paper/benchmark_hdfe_results.csv, write replace text
file write `bres' "method,n,flavor,reps,time_med_s,time_lo_s,time_hi_s,b_age,se_age" _n
file write `bres' "poisson_1000_dummies,1000000,`=c(flavor)',`reps',`=string(`t1_med',"%9.4f")',`=string(`t1_lo',"%9.4f")',`=string(`t1_hi',"%9.4f")',`=string(`b_ref',"%12.7f")',`=string(`se_ref',"%12.7f")'" _n
file write `bres' "stpois_absorb_g,1000000,`=c(flavor)',`reps',`=string(`t2_med',"%9.4f")',`=string(`t2_lo',"%9.4f")',`=string(`t2_hi',"%9.4f")',`=string(`b_st',"%12.7f")',`=string(`se_st',"%12.7f")'" _n
if `have_ppml' {
    file write `bres' "ppmlhdfe_absorb_g,1000000,`=c(flavor)',`reps',`=string(`t3_med',"%9.4f")',`=string(`t3_lo',"%9.4f")',`=string(`t3_hi',"%9.4f")',`=string(`b_pp',"%12.7f")',." _n
}
file close `bres'
di "(results in paper/benchmark_hdfe_results.csv)"
