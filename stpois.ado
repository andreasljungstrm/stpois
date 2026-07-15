*! version 0.8.1  15jul2026  Andreas Ljungström, SOFI Stockholm University
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
        ABSorb(string)                        ///
        FAST                                  ///
        D(name)                               ///
        TOLerance(real 1e-8)                  ///
        MAXIter(integer 100)                  ///
        *                                     ///
        ]

    st_is 2 analysis

    // absorb() accepts variable names and varname#varname interactions
    local absorb_raw ""
    if `"`absorb'"' != "" {
        foreach aterm of local absorb {
            local araw = subinstr("`aterm'", "#", " ", .)
            local araw = subinstr("`araw'", "i.", "", .)
            foreach av of local araw {
                capture confirm numeric variable `av'
                if _rc {
                    di as error "absorb(): `av' is not a numeric variable"
                    exit 198
                }
            }
            local absorb_raw `absorb_raw' `araw'
        }
    }

    // Build sample
    marksample touse
    if `"`absorb'"' != "" markout `touse' `absorb_raw'
    qui replace `touse' = 0 if _st == 0

    // Exposure
    tempvar exposure
    qui gen double `exposure' = _t - _t0 if `touse'
    qui replace `touse' = 0 if (`exposure' <= 0 | `exposure' >= .) & `touse'

    // Weight expression (standard path) and weight variable (Mata engines)
    local wgtexpr
    tempvar wv
    if "`weight'" != "" {
        local wgtexpr "[`weight'`exp']"
        qui gen double `wv' `exp' if `touse'
        markout `touse' `wv'
        qui replace `touse' = 0 if `wv' <= 0 & `touse'
        if "`weight'" == "fweight" {
            capture assert `wv' == int(`wv') if `touse'
            if _rc {
                di as error "may not use noninteger frequency weights"
                exit 401
            }
        }
    }

    // Cluster variable used by the Mata engines: also reachable via
    // vce(cluster clustvar); mark out its missings like poisson does
    local eng_clustv "`cluster'"
    if `"`vce'"' != "" {
        local vword1 = lower(word(`"`vce'"', 1))
        if substr("`vword1'", 1, 2) == "cl" local eng_clustv = word(`"`vce'"', 2)
    }
    if ("`fast'" != "" | `"`absorb'"' != "") & "`eng_clustv'" != "" {
        markout `touse' `eng_clustv'
    }

    // pweights imply robust standard errors on the Mata engines
    // (poisson handles this itself on the standard path)
    if "`weight'" == "pweight" & ("`fast'" != "" | `"`absorb'"' != "") ///
        & `"`vce'"' == "" & "`robust'" == "" & "`cluster'" == "" {
        local robust robust
    }

    // Separation: drop observations in absorbed-term groups with no
    // (weighted) events; their fixed effects diverge to -infinity
    if `"`absorb'"' != "" {
        foreach aterm of local absorb {
            local araw = subinstr("`aterm'", "#", " ", .)
            local araw = subinstr("`araw'", "i.", "", .)
            tempvar n_ev
            if "`weight'" != "" {
                qui bysort `araw': egen double `n_ev' = total(_d * `wv' * `touse')
            }
            else {
                qui bysort `araw': egen double `n_ev' = total(_d * `touse')
            }
            qui count if `touse' & `n_ev' == 0
            if r(N) > 0 {
                di as txt "  (absorb(`aterm'): dropped " r(N) " obs in zero-event groups)"
                qui replace `touse' = 0 if `n_ev' == 0
            }
            drop `n_ev'
        }
    }

    // FE-contribution variable: filled by the engines at convergence so
    // that predict can form full linear predictors after absorb()
    if `"`absorb'"' != "" {
        if "`d'" != "" {
            confirm new variable `d'
        }
        else {
            local d "_stpois_fe"
            capture drop _stpois_fe
        }
        qui gen double `d' = .
        label variable `d' "stpois absorb() fixed-effect contribution"
    }
    else if "`d'" != "" {
        di as error "d() requires absorb()"
        exit 198
    }

    // Options for the Mata engines
    local engopts
    if "`weight'" != "" local engopts wvar(`wv') wtype(`weight')
    if `"`absorb'"' != "" {
        local fespec ""
        foreach aterm of local absorb {
            local araw = subinstr("`aterm'", "#", " ", .)
            local araw = subinstr("`araw'", "i.", "", .)
            local fespec "`fespec';`araw'"
        }
        local fespec = substr("`fespec'", 2, .)
        local engopts `engopts' fespec(`fespec') dvar(`d')
    }

    // Common pass-through opts
    local poisopts `log' `constant'
    if `"`vce'"'    != "" local poisopts `poisopts' `vce'
    if "`robust'"   != "" local poisopts `poisopts' robust
    if "`cluster'"  != "" local poisopts `poisopts' cluster(`cluster')

    // ── Dispatch ──────────────────────────────────────────────────────────
    local eng_cluster
    if "`cluster'" != "" local eng_cluster cluster(`cluster')
    if "`fast'" != "" {
        // fast engine, with or without absorbed fixed effects
        _stpois_fast `varlist', ///
            touse(`touse')                  ///
            exposure(`exposure')            ///
            tol(`tolerance')                ///
            maxiter(`maxiter')              ///
            `engopts'                       ///
            `eng_cluster'                   ///
            `poisopts' `options'
        if `"`absorb'"' != "" {
            local etitle "Poisson EHA — fast + HDFE (absorb: `absorb')"
            ereturn local absorb "`absorb'"
        }
        else {
            local etitle "Poisson EHA — fast (cell-accelerated exact MLE)"
        }
        ereturn local fast "fast"
    }
    else if `"`absorb'"' != "" {
        _stpois_hdfe `varlist', ///
            touse(`touse')                 ///
            exposure(`exposure')           ///
            tol(`tolerance')               ///
            maxiter(`maxiter')             ///
            `engopts'                      ///
            `eng_cluster'                  ///
            `poisopts' `options'
        local etitle "Poisson EHA — HDFE (absorb: `absorb')"
        ereturn local absorb "`absorb'"
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
        // "_" is the null equation; assigning "" leaves eq names unchanged
        local blankeq
        forvalues i = 1/`ncols' {
            local blankeq `blankeq' _
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
    // Computed on the final estimation sample, so the header agrees with
    // e(N)/e(sample) even after absorb() drops zero-event groups
    tempvar esamp
    qui gen byte `esamp' = e(sample)
    if "`weight'" == "fweight" {
        // frequency weights scale failures and time at risk (as in streg)
        tempvar wfail
        qui gen double `wfail' = `wv' * _d if `esamp'
        qui summ `wfail' if `esamp', meanonly
        local N_fail = r(sum)
    }
    else {
        qui count if _d == 1 & `esamp'
        local N_fail = r(N)
    }
    local idvar `_dta[st_id]'
    if `"`idvar'"' != "" {
        // count distinct ids in Mata; levelsof overflows the macro
        // buffer when the id variable has hundreds of thousands of levels
        capture confirm string variable `idvar'
        if _rc {
            mata: st_local("N_sub", ///
                strofreal(rows(uniqrows(st_data(., "`idvar'", "`esamp'")))))
        }
        else {
            mata: st_local("N_sub", ///
                strofreal(rows(uniqrows(st_sdata(., "`idvar'", "`esamp'")))))
        }
    }
    else {
        qui count if `esamp'
        local N_sub = r(N)
    }
    if "`weight'" == "fweight" {
        tempvar wrisk
        qui gen double `wrisk' = `wv' * `exposure' if `esamp'
        qui summ `wrisk' if `esamp', meanonly
        local risk = r(sum)
    }
    else {
        qui summ `exposure' if `esamp', meanonly
        local risk = r(sum)
    }
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
    // absorb() reports a Wald test of the non-absorbed coefficients; an LR
    // test against a null without the fixed effects would not be valid
    if `"`absorb'"' != "" {
        ereturn local chi2type "Wald"
        ereturn local fes_var  "`d'"
    }
    else {
        ereturn local chi2type "LR"
    }
    ereturn local title     "`etitle'"
    ereturn local properties "b V"

    Display, `irr' level(`level')
end


program Display
    syntax [, IRr Level(cilevel) *]

    di _n as txt `"`e(title)'"'

    di _n                                                                    ///
        as txt "No. of subjects = "                                          ///
        as res %12.0fc e(N_sub)                                              ///
        _col(49) as txt "Number of obs"                                      ///
        _col(65) "=" as res %9.0fc e(N)
    di as txt "No. of failures = " as res %12.0fc e(N_fail)

    if e(chi2) < . {
        local c2t = cond(`"`e(chi2type)'"' == "", "LR", `"`e(chi2type)'"')
        di as txt "Time at risk    = " as res %12.4g e(risk)                ///
            _col(49) as txt "`c2t' chi2(" as res e(df_m) as txt ")"         ///
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
