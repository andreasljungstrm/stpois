*! version 0.9.0  16jul2026  Andreas Ljungström, SOFI Stockholm University
program stpois_p
    version 14

    // Quick first-pass to detect scores option.
    // svy calls: stpois_p double stub* if e(sample), scores
    // `anything' absorbs the type+stub without triggering newvarname rules.
    local full0 `"`0'"'
    syntax [anything] [if] [in] [, SCores *]

    if "`scores'" != "" {
        GenScores `full0'
        exit
    }

    // ── single-variable predictions ────────────────────────────────────────────
    local 0 `"`full0'"'
    syntax newvarname [if] [in] [,   ///
        XB                           ///
        N                            ///
        Hazard                       ///
        HR                           ///
        SURVival                     ///
        noOFFset                     ///
        ]

    // Exactly one statistic
    local nopt : word count `xb' `n' `hazard' `hr' `survival'
    if `nopt' > 1 {
        di as error "only one statistic may be specified"
        exit 198
    }
    if `nopt' == 0 local xb xb

    local newv `varlist'

    // After absorb(), e(b) does not contain the absorbed fixed effects.
    // Their total contribution was stored at estimation in the variable
    // named by e(fes_var) (default _stpois_fe, or the d() option), so the
    // full linear predictor is Xb + that variable.
    local fevar ""
    if `"`e(absorb)'"' != "" {
        local fevar `"`e(fes_var)'"'
        capture confirm numeric variable `fevar'
        if _rc {
            di as error "the fixed-effect variable `fevar' stored by absorb() was not found;"
            di as error "it is required for predictions after absorb() -- refit the model"
            exit 111
        }
    }

    // Copy e(b) to a named tempmatrix before scoring (avoids r(111) inside margins)
    tempname bmat
    matrix `bmat' = e(b)

    // Linear predictor Xβ (no offset), plus the absorbed FE contribution
    tempvar xbhat
    matrix score double `xbhat' = `bmat' `if' `in'
    if "`fevar'" != "" {
        qui replace `xbhat' = `xbhat' + `fevar' `if' `in'
    }

    if "`xb'" != "" {
        // `offset' = "nooffset" when user specifies nooffset; "" otherwise
        if "`offset'" != "" {
            gen double `newv' = `xbhat' `if' `in'
        }
        else {
            // xb + ln(exposure)
            gen double `newv' = `xbhat' + ln(_t - _t0) `if' `in'
        }
    }
    else if "`n'" != "" {
        // Expected events = exp(Xβ) * exposure
        gen double `newv' = exp(`xbhat') * (_t - _t0) `if' `in'
    }
    else if "`hazard'" != "" | "`hr'" != "" {
        // Hazard rate = exp(Xβ)
        gen double `newv' = exp(`xbhat') `if' `in'
    }
    else if "`survival'" != "" {
        // S(t) = exp(-exp(Xβ) * (t - t0))
        gen double `newv' = exp(-exp(`xbhat') * (_t - _t0)) `if' `in'
    }
end


// Score for svy linearization: dlnL_i/d(eta_i) = y_i - mu_i
// Called as: GenScores double stub* if e(sample), scores
program GenScores, sclass
    version 9, missing

    if `"`e(absorb)'"' != "" {
        di as error "predict, scores not supported after absorb()"
        di as error "for svy variance with absorb(), use: svy bootstrap: stpois ..., absorb(...)"
        exit 198
    }
    if `"`e(fast)'"' != "" {
        di as error "predict, scores not supported with fast"
        di as error "for svy variance with fast, use: svy bootstrap: stpois ..., fast"
        exit 198
    }

    syntax [anything] [if] [in] [, *]
    marksample touse

    // _score_spec creates the new variable(s) from the type+stub spec
    _score_spec `anything', `options'
    local scvar  `s(varlist)'
    local sctype `s(typlist)'

    // mu = exp(Xb + ln(exposure));  score = d - mu
    tempname bmat
    matrix `bmat' = e(b)
    tempvar xbhat
    qui matrix score double `xbhat' = `bmat' if `touse'
    qui gen `sctype' `scvar' = _d - exp(`xbhat' + ln(_t - _t0)) if `touse'

    sreturn local scorevars `scvar'
end
