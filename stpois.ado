*! version 0.5.0  13jul2026  Andreas Ljungström, SOFI Stockholm University
*! Poisson event-history regression for stset data
program stpois, eclass properties(st)
    version 14

    if replay() {
        if "`e(cmd)'" != "stpois" error 301
        Display `0'
        exit
    }

    syntax varlist(fv ts) [if] [in]          ///
        [iweight pweight fweight]             ///
        [,                                    ///
        IRr                                   ///
        noLOg                                 ///
        Level(cilevel)                        ///
        noCONstant                            ///
        VCE(passthru)                         ///
        ROBust                                ///
        CLuster(varname)                      ///
        ABSorb(varlist)                       ///
        FAST(string)                          ///
        MTopel                                ///
        TOLerance(real 1e-8)                  ///
        MAXIter(integer 100)                  ///
        *                                     ///
        ]

    st_is 2 analysis

    // Validate
    if `"`fast'"' != "" & `"`absorb'"' != "" {
        di as error "fast() and absorb() may not be combined"
        exit 198
    }
    if `"`fast'"' != "" {
        local fast = lower("`fast'")
        if !inlist("`fast'", "offset", "moments") {
            di as error "fast() requires {bf:offset} or {bf:moments}"
            exit 198
        }
    }
    if "`mtopel'" != "" & `"`fast'"' != "offset" {
        di as error "mtopel may only be specified with fast(offset)"
        exit 198
    }

    // Build sample
    marksample touse
    if `"`absorb'"' != "" markout `touse' `absorb'
    qui replace `touse' = 0 if _st == 0

    // Exposure
    tempvar exposure
    qui gen double `exposure' = _t - _t0 if `touse'
    qui replace `touse' = 0 if (`exposure' <= 0 | `exposure' >= .) & `touse'

    // Survival summaries (read-only; computed before any collapse)
    qui count if _d == 1 & `touse'
    local N_fail = r(N)
    local idvar `_dta[st_id]'
    if `"`idvar'"' != "" {
        // count distinct ids in Mata; levelsof overflows the macro
        // buffer when the id variable has hundreds of thousands of levels
        capture confirm string variable `idvar'
        if _rc {
            mata: st_local("N_sub", ///
                strofreal(rows(uniqrows(st_data(., "`idvar'", "`touse'")))))
        }
        else {
            mata: st_local("N_sub", ///
                strofreal(rows(uniqrows(st_sdata(., "`idvar'", "`touse'")))))
        }
    }
    else {
        qui count if `touse'
        local N_sub = r(N)
    }
    qui summ `exposure' if `touse', meanonly
    local risk = r(sum)

    // Weight expression
    local wgtexpr
    if "`weight'" != "" local wgtexpr "[`weight'`exp']"

    // Common pass-through opts
    local poisopts `log' `constant'
    if `"`vce'"'    != "" local poisopts `poisopts' `vce'
    if "`robust'"   != "" local poisopts `poisopts' robust
    if "`cluster'"  != "" local poisopts `poisopts' cluster(`cluster')

    // ── Dispatch ──────────────────────────────────────────────────────────
    if `"`absorb'"' != "" {
        // Pass cluster() explicitly alongside poisopts
        local hdfe_cluster
        if "`cluster'" != "" local hdfe_cluster cluster(`cluster')
        _stpois_hdfe `varlist' `wgtexpr', ///
            touse(`touse')                 ///
            exposure(`exposure')           ///
            absorb(`absorb')               ///
            tol(`tolerance')               ///
            maxiter(`maxiter')             ///
            `hdfe_cluster'                 ///
            `poisopts' `options'
        local etitle "Poisson EHA — HDFE (absorb: `absorb')"
        ereturn local absorb "`absorb'"
    }
    else if `"`fast'"' == "offset" {
        _stpois_fast_offset `varlist' `wgtexpr', ///
            touse(`touse')                         ///
            exposure(`exposure')                   ///
            `mtopel'                               ///
            `poisopts' `options'
        if "`mtopel'" != "" {
            local etitle "Poisson EHA — fast: two-stage offset (Murphy–Topel SEs)"
        }
        else {
            local etitle "Poisson EHA — fast: two-stage offset (conditional SEs)"
        }
        ereturn local fast "offset"
    }
    else if `"`fast'"' == "moments" {
        _stpois_fast_moments `varlist' `wgtexpr', ///
            touse(`touse')                          ///
            exposure(`exposure')                    ///
            tol(`tolerance')                        ///
            maxiter(`maxiter')                      ///
            `poisopts' `options'
        local etitle "Poisson EHA — fast: iterated cell moments (exact MLE)"
        ereturn local fast "moments"
    }
    else {
        // ── Standard path ─────────────────────────────────────────────────
        qui poisson _d `varlist' `wgtexpr' if `touse', ///
            exposure(`exposure') `poisopts' `options'

        local N     = e(N)
        local ll    = e(ll)
        local ll_0  = e(ll_0)
        local chi2  = e(chi2)
        local df_m  = e(df_m)
        local rank  = e(rank)
        local k     = e(k)
        local vce_r = e(vce)
        local ic    = e(ic)
        local conv  = e(converged)

        tempname b V
        matrix `b' = e(b)
        matrix `V' = e(V)

        local ncols = colsof(`b')
        local blankeq
        forvalues i = 1/`ncols' {
            local blankeq `blankeq' ""
        }
        matrix coleq `b' = `blankeq'
        matrix coleq `V' = `blankeq'
        matrix roweq `V' = `blankeq'

        ereturn clear
        ereturn post `b' `V', esample(`touse') obs(`N')
        ereturn scalar ll        = `ll'
        ereturn scalar ll_0      = `ll_0'
        ereturn scalar chi2      = `chi2'
        ereturn scalar df_m      = `df_m'
        ereturn scalar rank      = `rank'
        ereturn scalar k         = `k'
        ereturn scalar ic        = `ic'
        ereturn scalar converged = `conv'
        ereturn local  vce       "`vce_r'"
        local etitle "Poisson event-history regression"
    }

    // ── Survival scalars (common) ──────────────────────────────────────────
    ereturn scalar N_fail = `N_fail'
    ereturn scalar N_sub  = `N_sub'
    ereturn scalar risk   = `risk'

    // ── Standard cmd locals ────────────────────────────────────────────────
    ereturn local cmd       "stpois"
    ereturn local cmd2      "stpois"
    ereturn local cmdline   "stpois `0'"
    ereturn local predict   "stpois_p"
    ereturn local depvar    "_d"
    ereturn local dead      "_d"
    ereturn local t0        "_t0"
    ereturn local chi2type  "LR"
    ereturn local title     "`etitle'"
    ereturn local properties "b V"

    Display, `irr' level(`level')
end


program Display
    syntax [, IRr Level(cilevel) *]

    di _n as txt `"`e(title)'"'

    if `"`e(fast)'"' == "offset" {
        if e(mtopel) == 1 {
            di as txt "  (Note: SEs are Murphy–Topel corrected; " ///
                      "first-stage uncertainty propagated)"
        }
        else {
            di as txt "  (Note: SEs are conditional on first-stage estimates; " ///
                      "specify mtopel to propagate first-stage uncertainty)"
        }
    }

    di _n                                                                    ///
        as txt "No. of subjects = "                                          ///
        as res %12.0fc e(N_sub)                                              ///
        _col(49) as txt "Number of obs"                                      ///
        _col(65) "=" as res %9.0fc e(N)
    di as txt "No. of failures = " as res %12.0fc e(N_fail)

    if e(chi2) < . {
        di as txt "Time at risk    = " as res %12.4g e(risk)                ///
            _col(49) as txt "LR chi2(" as res e(df_m) as txt ")"            ///
            _col(65) "=" as res %9.2f e(chi2)
        di as txt "Log likelihood  = " as res %12.5f e(ll)                  ///
            _col(49) as txt "Prob > chi2"                                    ///
            _col(65) "=" as res %9.4f chi2tail(e(df_m), e(chi2))
    }
    else {
        di as txt "Time at risk    = " as res %12.4g e(risk)
        di as txt "Log likelihood  = " as res %12.5f e(ll)
    }
    di

    if "`irr'" != "" {
        ereturn display, eform(IRR) level(`level')
    }
    else {
        ereturn display, level(`level')
    }
end
