// sj_examples.do -- generate sjlog output files for the Stata Journal
// version of the stpois paper. Run from the repository root:
//   do paper/sj/sj_examples.do
// Requires the sjlatex package (sjlog): net install sjlatex,
//   from(http://www.stata-journal.com/production)
discard
adopath ++ "`c(pwd)'"
cd paper/sj

set more off
capture log close

// ── Example 1: standard path ───────────────────────────────────────────────
webuse stan3, clear
qui stset t1, failure(died) id(id)
sjlog using ex_standard, replace
stpois age posttran, nolog
sjlog close, replace

// ── Example 2: absorb() ────────────────────────────────────────────────────
sjlog using ex_absorb, replace
stpois age posttran, absorb(surgery) nolog
sjlog close, replace

// ── Example 3: fast on 500,000 observations ───────────────────────────────
clear
qui set obs 500000
set seed 20260713
qui gen edu    = ceil(4 * runiform())
qui gen region = ceil(5 * runiform())
qui gen period = ceil(10 * runiform())
qui gen age    = 30 + 20 * runiform()
qui gen income = rnormal(10, 2)
qui gen lambda = exp(-6 + ///
    0.03 * age + 0.05 * income                       + ///
    0.3 * (edu==2) + 0.5 * (edu==3) - 0.2 * (edu==4) + ///
    0.1 * (region==2) - 0.2 * (region==3)            + ///
    0.2 * (period >= 6))
qui gen ptime  = 0.5 + 0.5 * runiform()
qui gen failed = (runiform() < 1 - exp(-lambda * ptime))
qui gen id2    = _n
qui stset ptime, failure(failed) id(id2)
sjlog using ex_fast, replace
stpois age income i.edu i.region i.period, fast nolog
sjlog close, replace

// ── Example 4: interactions with fast ──────────────────────────────────────
sjlog using ex_inter, replace
stpois c.age##i.edu income i.region i.period, fast nolog
sjlog close, replace

// ── Example 5: fast + absorb and predictions ───────────────────────────────
webuse stan3, clear
qui stset t1, failure(died) id(id)
sjlog using ex_combo, replace
stpois age i.posttran, fast absorb(surgery) nolog
sjlog close, replace

sjlog using ex_predict, replace
stpois age posttran, absorb(surgery) nolog
predict double haz, hazard
summarize haz
sjlog close, replace

cd ../..
