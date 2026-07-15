// benchmark_ppmlhdfe.do -- head-to-head of stpois against ppmlhdfe and the
// collapse+poisson folk baseline. Verifies numerical AGREEMENT (flavor-
// independent) and records timings (median of `reps' runs). Run from repo root.
// Writes paper/benchmark_ppmlhdfe_results.csv (versioned).
//
// Requires: ppmlhdfe, reghdfe, ftools (ssc install). The agreement asserts are
// the substantive result; timings should be produced on the same Stata flavor
// as the other paper tables (Stata MP).
discard
adopath ++ "`c(pwd)'"
set more off

local reps 10          // repetitions for median timings

// median-of-reps timer helper: returns r(med) for a command string
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
    // median
    preserve
        clear
        quietly svmat double `T'
        quietly summarize `T'1, detail
        return scalar med = r(p50)
    restore
end

tempname R
file open `R' using paper/benchmark_ppmlhdfe_results.csv, write replace text
file write `R' "comparison,flavor,n,metric,stpois,other,absdiff,stpois_med_s,other_med_s" _n
local flav = c(flavor)

// ============================================================================
// A. absorb() vs ppmlhdfe -- high-dimensional FE (Table 4 design)
// ============================================================================
di _n "===== A. absorb() vs ppmlhdfe (1,000-level FE) ====="
clear
set obs 1000000
set seed 20260713
gen g      = ceil(1000 * runiform())
gen age    = 30 + 20 * runiform()
gen income = rnormal(10, 2)
gen id2    = _n
gen geff   = (mod(g, 7) - 3) * 0.1
gen lambda = exp(-6 + 0.03*age + 0.05*income + geff)
gen ptime  = 0.5 + 0.5 * runiform()
gen failed = (runiform() < 1 - exp(-lambda * ptime))
qui stset ptime, failure(failed) id(id2)
qui gen double expo = _t - _t0
qui gen double lnexp = ln(expo)

// -- coefficient agreement (ppmlhdfe's DEFAULT vce is robust, not oim; it has
//    no vce(oim), so SE agreement is checked robust-to-robust and cluster-to-
//    cluster below. stpois vce(oim) reproduces poisson's classical OIM.)
qui stpois age income, absorb(g) nolog
local sp_b  = _b[age]
local sp_se_oim = _se[age]        // stpois OIM (no ppmlhdfe counterpart)
qui ppmlhdfe _d age income if _st, absorb(g) offset(lnexp)
local pp_b  = _b[age]
di "COEF age:   stpois b=" %12.9f `sp_b' "   ppmlhdfe b=" %12.9f `pp_b'
di "  b diff=" %14.3e abs(`sp_b'-`pp_b')
di "  (stpois vce(oim) se=" %12.9f `sp_se_oim' "; ppmlhdfe has no oim)"
assert abs(`sp_b'  - `pp_b')  < 1e-5
file write `R' "absorb_vs_ppmlhdfe,`flav',1000000,b_age,`=string(`sp_b',"%14.10f")',`=string(`pp_b',"%14.10f")',`=string(abs(`sp_b'-`pp_b'),"%14.3e")',,"  _n

// -- robust
qui stpois age income, absorb(g) vce(robust) nolog
local sp_se_r = _se[age]
qui ppmlhdfe _d age income if _st, absorb(g) offset(lnexp) vce(robust)
local pp_se_r = _se[age]
di "ROBUST se: stpois=" %12.9f `sp_se_r' "  ppmlhdfe=" %12.9f `pp_se_r' ///
   "  diff=" %14.3e abs(`sp_se_r'-`pp_se_r')
file write `R' "absorb_vs_ppmlhdfe,`flav',1000000,se_age_robust,`=string(`sp_se_r',"%14.10f")',`=string(`pp_se_r',"%14.10f")',`=string(abs(`sp_se_r'-`pp_se_r'),"%14.3e")',,"  _n

// -- cluster
qui gen clu = mod(id2, 500)
qui stpois age income, absorb(g) vce(cluster clu) nolog
local sp_se_c = _se[age]
qui ppmlhdfe _d age income if _st, absorb(g) offset(lnexp) vce(cluster clu)
local pp_se_c = _se[age]
di "CLUSTER se: stpois=" %12.9f `sp_se_c' "  ppmlhdfe=" %12.9f `pp_se_c' ///
   "  diff=" %14.3e abs(`sp_se_c'-`pp_se_c')
file write `R' "absorb_vs_ppmlhdfe,`flav',1000000,se_age_cluster,`=string(`sp_se_c',"%14.10f")',`=string(`pp_se_c',"%14.10f")',`=string(abs(`sp_se_c'-`pp_se_c'),"%14.3e")',,"  _n
di "A: stpois absorb() coefficients == ppmlhdfe; SE differences recorded"

// -- timings (median)
_medtime `reps' "stpois age income, absorb(g) nolog"
local sp_t = r(med)
_medtime `reps' "ppmlhdfe _d age income if _st, absorb(g) offset(lnexp)"
local pp_t = r(med)
di "TIME absorb: stpois=" `sp_t' "s  ppmlhdfe=" `pp_t' "s  (median of `reps', flavor `flav')"
file write `R' "absorb_vs_ppmlhdfe,`flav',1000000,time,,,,`=string(`sp_t',"%9.4f")',`=string(`pp_t',"%9.4f")'" _n

// ============================================================================
// B. fast + absorb() vs ppmlhdfe -- combined mechanism
// ============================================================================
di _n "===== B. fast+absorb() vs ppmlhdfe ====="
clear
set obs 1000000
set seed 20260713
gen edu    = ceil(4 * runiform())
gen region = ceil(5 * runiform())
gen g      = ceil(200 * runiform())
gen age    = 30 + 20 * runiform()
gen income = rnormal(10, 2)
gen id2    = _n
gen geff   = (mod(g, 7) - 3) * 0.08
gen lambda = exp(-6 + 0.03*age + 0.05*income + 0.3*(edu==2) + 0.5*(edu==3) ///
    - 0.2*(edu==4) + 0.1*(region==2) - 0.2*(region==3) + geff)
gen ptime  = 0.5 + 0.5 * runiform()
gen failed = (runiform() < 1 - exp(-lambda * ptime))
qui stset ptime, failure(failed) id(id2)
qui gen double lnexp = ln(_t - _t0)

qui stpois age income i.edu i.region, fast absorb(g) nolog
local sp_b  = _b[age]
qui ppmlhdfe _d age income i.edu i.region if _st, absorb(g) offset(lnexp)
local pp_b  = _b[age]
di "COEF age: stpois b=" %12.9f `sp_b' "   ppmlhdfe b=" %12.9f `pp_b' ///
   "   diff=" %14.3e abs(`sp_b'-`pp_b')
assert abs(`sp_b'  - `pp_b')  < 1e-5
file write `R' "fastabsorb_vs_ppmlhdfe,`flav',1000000,b_age,`=string(`sp_b',"%14.10f")',`=string(`pp_b',"%14.10f")',`=string(abs(`sp_b'-`pp_b'),"%14.3e")',,"  _n
// robust-to-robust SE agreement
qui stpois age income i.edu i.region, fast absorb(g) vce(robust) nolog
local sp_se_r = _se[age]
qui ppmlhdfe _d age income i.edu i.region if _st, absorb(g) offset(lnexp) vce(robust)
local pp_se_r = _se[age]
di "ROBUST se: stpois=" %12.9f `sp_se_r' "  ppmlhdfe=" %12.9f `pp_se_r' ///
   "  diff=" %14.3e abs(`sp_se_r'-`pp_se_r')
file write `R' "fastabsorb_vs_ppmlhdfe,`flav',1000000,se_age_robust,`=string(`sp_se_r',"%14.10f")',`=string(`pp_se_r',"%14.10f")',`=string(abs(`sp_se_r'-`pp_se_r'),"%14.3e")',,"  _n
di "B: stpois fast+absorb coefficients == ppmlhdfe"

_medtime `reps' "stpois age income i.edu i.region, fast absorb(g) nolog"
local sp_t = r(med)
_medtime `reps' "ppmlhdfe _d age income i.edu i.region if _st, absorb(g) offset(lnexp)"
local pp_t = r(med)
di "TIME fast+absorb: stpois=" `sp_t' "s  ppmlhdfe=" `pp_t' "s"
file write `R' "fastabsorb_vs_ppmlhdfe,`flav',1000000,time,,,,`=string(`sp_t',"%9.4f")',`=string(`pp_t',"%9.4f")'" _n

// ============================================================================
// C. fast vs collapse+poisson -- the all-categorical folk baseline
// ============================================================================
di _n "===== C. fast vs collapse+poisson (all-categorical) ====="
clear
set obs 1000000
set seed 20260713
gen edu    = ceil(4 * runiform())
gen region = ceil(5 * runiform())
gen period = ceil(10 * runiform())
gen id2    = _n
gen lambda = exp(-6 + 0.3*(edu==2) + 0.5*(edu==3) - 0.2*(edu==4) ///
    + 0.1*(region==2) - 0.2*(region==3) + 0.2*(period >= 6))
gen ptime  = 0.5 + 0.5 * runiform()
gen failed = (runiform() < 1 - exp(-lambda * ptime))
qui stset ptime, failure(failed) id(id2)

qui stpois i.edu i.region i.period, fast nolog
local sp_b  = _b[3.edu]
local sp_se = _se[3.edu]

// folk baseline: collapse to cells, poisson with summed events/exposure
preserve
    qui gen double expo = _t - _t0
    qui collapse (sum) D=_d (sum) T=expo, by(edu region period)
    qui poisson D i.edu i.region i.period, exposure(T) nolog
    local cp_b  = _b[3.edu]
    local cp_se = _se[3.edu]
restore
di "3.edu: stpois-fast b=" %12.9f `sp_b' " se=" %12.9f `sp_se' ///
   "   collapse+poisson b=" %12.9f `cp_b' " se=" %12.9f `cp_se'
assert abs(`sp_b'  - `cp_b')  < 1e-6
assert abs(`sp_se' - `cp_se') < 1e-6
file write `R' "fast_vs_collapse,`flav',1000000,b_edu3,`=string(`sp_b',"%14.10f")',`=string(`cp_b',"%14.10f")',`=string(abs(`sp_b'-`cp_b'),"%14.3e")',,"  _n
file write `R' "fast_vs_collapse,`flav',1000000,se_edu3,`=string(`sp_se',"%14.10f")',`=string(`cp_se',"%14.10f")',`=string(abs(`sp_se'-`cp_se'),"%14.3e")',,"  _n
di "PASS C: stpois fast == collapse+poisson"

// timings: fast vs full collapse+poisson workflow (collapse cost included)
_medtime `reps' "stpois i.edu i.region i.period, fast nolog"
local sp_t = r(med)
capture program drop _collapsefit
program _collapsefit
    preserve
        qui gen double expo = _t - _t0
        qui collapse (sum) D=_d (sum) T=expo, by(edu region period)
        qui poisson D i.edu i.region i.period, exposure(T) nolog
    restore
end
_medtime `reps' "_collapsefit"
local cp_t = r(med)
di "TIME allcat: fast=" `sp_t' "s  collapse+poisson=" `cp_t' "s"
file write `R' "fast_vs_collapse,`flav',1000000,time,,,,`=string(`sp_t',"%9.4f")',`=string(`cp_t',"%9.4f")'" _n

file close `R'
di _n "===== DONE (flavor `flav'); results in paper/benchmark_ppmlhdfe_results.csv ====="
